#include "Transport.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TNamed.h"
#include "TTree.h"

#include <limits>
#include <map>

namespace supera::io::larcv2 {
namespace {
unsigned int track(std::optional<std::int32_t> value) {
  if (!value)
    return larcv::kINVALID_UINT;
  if (*value < 0)
    throw std::invalid_argument("Signed transport ID cannot be represented in LArCV2");
  return static_cast<unsigned int>(*value);
}

double coordinate(double value) {
  require_measurement(value);
  if (std::isnan(value))
    return unavailable_transport_coordinate;
  // Leave a margin below the signed integer limit used internally by Vertex.
  if (std::abs(value) >= 9.e12)
    return unavailable_transport_coordinate; // Exact value and loss count remain in the archive.
  return value;
}

void position(larcv::Particle &output, const GeneratorStep *step, int field) {
  std::array<double, 4> x{missing_value, missing_value, missing_value, missing_value};
  if (step)
    x = {step->position_cm[0], step->position_cm[1], step->position_cm[2], step->time_ns};
  for (auto &value : x)
    value = coordinate(value);
  if (field == 0)
    output.position(x[0], x[1], x[2], x[3]);
  else if (field == 1)
    output.end_position(x[0], x[1], x[2], x[3]);
  else if (field == 2)
    output.parent_position(x[0], x[1], x[2], x[3]);
  else
    output.ancestor_position(x[0], x[1], x[2], x[3]);
}

const GeneratorStep *initial(const TransportParticle &particle) {
  return particle.trajectory.empty() ? nullptr : &particle.trajectory.front();
}
} // namespace

std::vector<larcv::Particle> project_transport(const TransportCollection &input) {
  std::vector<larcv::Particle> result;
  result.reserve(input.particles.size());
  for (std::size_t i = 0; i < input.particles.size(); ++i) {
    const auto &source = input.particles[i];
    auto &particle = result.emplace_back();
    particle.id(i);
    particle.track_id(track(source.track_id));
    particle.pdg_code(source.pdg);
    particle.creation_process(source.process);

    // Inventory kinematics use MeV and MeV/c, as legacy transported particles do.
    const auto *start = initial(source);
    position(particle, start, 0);
    position(particle, source.trajectory.empty() ? nullptr : &source.trajectory.back(), 1);
    if (start) {
      const auto &p = start->momentum_gev_c;
      for (double value : {p[0] * 1000., p[1] * 1000., p[2] * 1000., start->energy_gev * 1000.})
        require_measurement(value);
      particle.momentum(p[0] * 1000., p[1] * 1000., p[2] * 1000.);
      particle.energy_init(start->energy_gev * 1000.);
    } else {
      particle.momentum(missing_value, missing_value, missing_value);
      particle.energy_init(missing_value);
    }

    // Direct ancestry and recovered primary affiliation occupy distinct fields.
    particle.parent_track_id(track(source.mother_id));
    position(particle, nullptr, 2);
    position(particle, nullptr, 3);
    const auto parent = source.primary ? record_index(i) : source.direct_parent;
    if (parent != missing_index) {
      require_index(parent, input.particles.size());
      const auto &record = input.particles[parent];
      particle.parent_id(parent);
      particle.parent_track_id(track(record.track_id));
      particle.parent_pdg_code(record.pdg);
      particle.parent_creation_process(record.process);
      position(particle, initial(record), 2);
    }
    if (source.primary_ancestor != missing_index) {
      require_index(source.primary_ancestor, input.particles.size());
      const auto &record = input.particles[source.primary_ancestor];
      particle.ancestor_track_id(track(record.track_id));
      particle.ancestor_pdg_code(record.pdg);
      particle.ancestor_creation_process(record.process);
      position(particle, initial(record), 3);
    }

    // No charge, detector entry, semantic grouping, or interaction is fabricated.
    particle.energy_deposit(missing_value);
    particle.first_step(unavailable_transport_coordinate, unavailable_transport_coordinate,
                        unavailable_transport_coordinate, unavailable_transport_coordinate);
    particle.last_step(unavailable_transport_coordinate, unavailable_transport_coordinate,
                       unavailable_transport_coordinate, unavailable_transport_coordinate);
  }
  std::vector<std::vector<larcv::InstanceID_t>> children(result.size());
  for (std::size_t i = 0; i < input.particles.size(); ++i) {
    const auto parent = input.particles[i].direct_parent;
    if (parent != missing_index && parent != record_index(i)) {
      require_index(parent, result.size());
      children[parent].push_back(i);
    }
  }
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i].children_id(children[i]);
  return result;
}

class TransportArchive::Impl {
public:
  std::unique_ptr<TFile> file;
  TTree *particles{};
  TTree *events{};
  unsigned long long run{}, subrun{}, event{};
  std::string collection, product, process, end_process, truth_product;
  long long index{}, source_index{}, direct_parent{}, available_parent{}, primary_ancestor{};
  long long truth_index{}, generated_index{};
  bool has_track{}, has_mother{}, has_status{}, has_truth{}, primary{};
  int track_id{}, mother_id{}, pdg{}, status{}, storage{}, ancestry{}, origin{};
  double mass{}, weight{};
  bool has_rescatter{};
  int rescatter{};
  unsigned int unrepresentable_trajectory_coordinates{};
  std::vector<double> polarization, generator_vertex;
  std::vector<long long> process_steps;
  std::vector<unsigned char> process_codes;
  std::vector<std::string> process_names;
  std::vector<int> daughters, map_ancestor, map_descendant, root_mothers;
  std::vector<double> trajectory;
};

TransportArchive::TransportArchive(const std::string &filename) : impl_(std::make_unique<Impl>()) {
  TDirectory::TContext directory;
  auto &s = *impl_;
  s.file.reset(TFile::Open(filename.c_str(), "RECREATE"));
  if (!s.file || s.file->IsZombie())
    throw std::runtime_error("Cannot open transport archive");
  s.particles = new TTree("transport", "Persisted transport inventory, schema version 1");
  s.events = new TTree("events", "Collection presence and complete suppression maps");
#define FIELD(name) s.particles->Branch(#name, &s.name)
  FIELD(run);
  FIELD(subrun);
  FIELD(event);
  FIELD(collection);
  FIELD(product);
  FIELD(index);
  FIELD(source_index);
  FIELD(has_track);
  FIELD(track_id);
  FIELD(has_mother);
  FIELD(mother_id);
  FIELD(pdg);
  FIELD(has_status);
  FIELD(status);
  FIELD(process);
  FIELD(end_process);
  FIELD(storage);
  FIELD(weight);
  FIELD(polarization);
  FIELD(generator_vertex);
  FIELD(has_rescatter);
  FIELD(rescatter);
  FIELD(process_steps);
  FIELD(process_codes);
  FIELD(process_names);
  FIELD(unrepresentable_trajectory_coordinates);
  FIELD(mass);
  FIELD(daughters);
  FIELD(trajectory);
  FIELD(has_truth);
  FIELD(truth_product);
  FIELD(truth_index);
  FIELD(generated_index);
  FIELD(origin);
  FIELD(primary);
  FIELD(direct_parent);
  FIELD(available_parent);
  FIELD(primary_ancestor);
  FIELD(ancestry);
#undef FIELD
#define FIELD(name) s.events->Branch(#name, &s.name)
  FIELD(run);
  FIELD(subrun);
  FIELD(event);
  FIELD(collection);
  FIELD(map_ancestor);
  FIELD(map_descendant);
  FIELD(root_mothers);
#undef FIELD
  TNamed("schema", "transport-v1; trajectory stride 8: x,y,z[cm],t[ns],px,py,pz[GeV/c],E[GeV]; "
                   "LArCV-only missing coordinate=-9e12; shared/archive measurements use NaN")
      .Write();
}

TransportArchive::~TransportArchive() = default;

void TransportArchive::write(const EventRecord &event) {
  auto &s = *impl_;
  if (!s.file->IsOpen())
    throw std::logic_error("Cannot write a closed transport archive");
  s.run = event.id.run;
  s.subrun = event.id.subrun;
  s.event = event.id.event;
  for (const auto &collection : event.transport) {
    s.collection = collection.name;
    s.root_mothers.assign(collection.root_mother_ids.begin(), collection.root_mother_ids.end());
    s.map_ancestor.clear();
    s.map_descendant.clear();
    for (const auto &[ancestor, descendants] : collection.suppression)
      for (auto descendant : descendants) {
        s.map_ancestor.push_back(ancestor);
        s.map_descendant.push_back(descendant);
      }
    if (s.events->Fill() < 0)
      throw std::runtime_error("Cannot save transport event metadata");
    for (std::size_t i = 0; i < collection.particles.size(); ++i) {
      const auto &p = collection.particles[i];
      s.index = record_index(i);
      s.product = p.product;
      s.source_index = p.source_index;
      s.has_track = p.track_id.has_value();
      s.track_id = p.track_id.value_or(0);
      s.has_mother = p.mother_id.has_value();
      s.mother_id = p.mother_id.value_or(0);
      s.pdg = p.pdg;
      s.has_status = p.status.has_value();
      s.status = p.status.value_or(0);
      s.process = p.process;
      s.end_process = p.end_process;
      s.storage = static_cast<int>(p.storage);
      s.mass = p.mass_gev;
      s.weight = p.weight;
      s.polarization.assign(p.polarization.begin(), p.polarization.end());
      s.generator_vertex.assign(p.generator_vertex.begin(), p.generator_vertex.end());
      s.has_rescatter = p.rescatter.has_value();
      s.rescatter = p.rescatter.value_or(0);
      s.process_steps.clear();
      s.process_codes.clear();
      s.process_names.clear();
      for (const auto &process : p.trajectory_processes) {
        s.process_steps.push_back(process.step);
        s.process_codes.push_back(process.source_code);
        s.process_names.push_back(process.name);
      }
      s.unrepresentable_trajectory_coordinates = 0;
      s.daughters.assign(p.daughter_ids.begin(), p.daughter_ids.end());
      s.primary = p.primary;
      s.direct_parent = p.direct_parent;
      s.available_parent = p.available_parent;
      s.primary_ancestor = p.primary_ancestor;
      s.ancestry = static_cast<int>(p.ancestry);
      s.has_truth = p.truth.has_value();
      s.truth_product = p.truth ? p.truth->product : "";
      s.truth_index = p.truth ? p.truth->truth_index : missing_index;
      s.generated_index = p.truth ? p.truth->generated_index : missing_index;
      s.origin = p.truth ? p.truth->origin : 0;
      s.trajectory.clear();
      for (const auto &step : p.trajectory) {
        for (double value :
             {step.position_cm[0], step.position_cm[1], step.position_cm[2], step.time_ns})
          s.unrepresentable_trajectory_coordinates +=
              std::isfinite(value) && std::abs(value) >= 9.e12;
        s.trajectory.insert(s.trajectory.end(), step.position_cm.begin(), step.position_cm.end());
        s.trajectory.push_back(step.time_ns);
        s.trajectory.insert(s.trajectory.end(), step.momentum_gev_c.begin(),
                            step.momentum_gev_c.end());
        s.trajectory.push_back(step.energy_gev);
      }
      if (s.particles->Fill() < 0)
        throw std::runtime_error("Cannot save transport particle");
    }
  }
}

void TransportArchive::close() {
  auto &s = *impl_;
  if (!s.file->IsOpen())
    return;
  s.file->cd();
  if (s.file->Write() <= 0)
    throw std::runtime_error("Cannot finalize transport archive");
  s.file->Close();
}
} // namespace supera::io::larcv2
