#include "Deposit.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TNamed.h"
#include "TTree.h"

#include <map>

namespace supera::io::larcv2 {
DepositProducts project_deposits(const DepositCollection &input) {
  if (!input.grid || input.kind == DepositKind::ReadoutIDE || input.particle_count < 0 ||
      input.particle_count == std::numeric_limits<RecordIndex>::max())
    throw std::invalid_argument(
        "Energy voxel output requires a physical-energy grid and inventory");
  const auto &g = *input.grid;
  std::uint64_t size = 1;
  for (std::size_t axis = 0; axis < 3; ++axis) {
    if (!std::isfinite(g.minimum_cm[axis]) || !std::isfinite(g.maximum_cm[axis]) ||
        g.maximum_cm[axis] <= g.minimum_cm[axis] ||
        !std::isfinite(g.maximum_cm[axis] - g.minimum_cm[axis]) || g.counts[axis] == 0 ||
        g.counts[axis] > ((std::uint64_t{1} << 53) - 1) ||
        g.counts[axis] > static_cast<std::uint64_t>(std::numeric_limits<RecordIndex>::max()) / size)
      throw std::invalid_argument("Invalid deposit grid");
    size *= g.counts[axis];
  }
  larcv::Voxel3DMeta meta;
  meta.set(g.minimum_cm[0], g.minimum_cm[1], g.minimum_cm[2], g.maximum_cm[0], g.maximum_cm[1],
           g.maximum_cm[2], g.counts[0], g.counts[1], g.counts[2]);
  DepositProducts result;
  result.clusters.resize(input.particle_count + 1);
  result.clusters.meta(meta);
  result.tensor.meta(meta);
  std::map<RecordIndex, double> total;
  std::optional<std::pair<RecordIndex, RecordIndex>> previous;
  for (const auto &voxel : input.voxels) {
    const auto key = std::make_pair(voxel.particle, voxel.voxel);
    if (previous && key <= *previous)
      throw std::invalid_argument("Deposit voxels must have sorted, unique particle/voxel keys");
    previous = key;
    require_index(voxel.voxel, size);
    const auto particle = voxel.particle == missing_index ? input.particle_count : voxel.particle;
    require_index(particle, static_cast<std::uint64_t>(input.particle_count) + 1);
    if (!std::isfinite(voxel.energy_mev) || voxel.energy_mev < 0. ||
        voxel.energy_mev > std::numeric_limits<float>::max())
      throw std::invalid_argument("Unrepresentable deposit energy");
    result.clusters.writeable_voxel_set(particle).emplace(voxel.voxel, float(voxel.energy_mev),
                                                          false);
    total[voxel.voxel] += voxel.energy_mev;
  }
  for (const auto &[id, energy] : total) {
    if (!std::isfinite(energy) || energy > std::numeric_limits<float>::max())
      throw std::overflow_error("Deposit voxel energy exceeds float range");
    result.tensor.emplace(id, float(energy), false);
  }
  return result;
}

class DepositArchive::Impl {
public:
  std::unique_ptr<TFile> file;
  TTree *rows{}, *events{};
  unsigned long long run{}, subrun{}, event{}, count{};
  std::string collection, product, transport;
  int kind{}, attribution{}, raw_track{}, raw_original{}, track{}, original{};
  bool has_raw_original{}, has_track{}, has_original{}, reassigned{};
  long long source_index{}, channel_index{}, ide_index{}, channel{}, tdc{};
  long long original_particle{}, assigned_particle{}, attributed_particle{}, primary{};
  double x{}, y{}, z{}, time_ns{}, energy_mev{}, electrons{};
  bool has_segment{}, has_pdg{};
  int pdg{};
  double start_x{}, start_y{}, start_z{}, end_x{}, end_y{}, end_z{};
  double start_ns{}, end_ns{}, photons{}, yield_ratio{};
  std::vector<unsigned int> channels;
  std::vector<unsigned long long> counts;
  std::vector<double> energies;
  unsigned long long invalid_energy{}, invalid_position{}, outside_grid{};
  double outside_energy{}, invalid_position_energy{};
};

DepositArchive::DepositArchive(const std::string &filename) : impl_(std::make_unique<Impl>()) {
  TDirectory::TContext context;
  auto &s = *impl_;
  s.file.reset(TFile::Open(filename.c_str(), "RECREATE"));
  if (!s.file || s.file->IsZombie())
    throw std::runtime_error("Cannot create deposit archive");
  s.rows = new TTree("deposits", "Native contributions and attribution, version 1");
  s.events = new TTree("events", "Complete event/collection inventory and accounting");
#define FIELD(name) s.rows->Branch(#name, &s.name)
  FIELD(run);
  FIELD(subrun);
  FIELD(event);
  FIELD(collection);
  FIELD(product);
  FIELD(transport);
  FIELD(kind);
  FIELD(source_index);
  FIELD(channel_index);
  FIELD(ide_index);
  FIELD(channel);
  FIELD(tdc);
  FIELD(raw_track);
  FIELD(has_raw_original);
  FIELD(raw_original);
  FIELD(has_track);
  FIELD(track);
  FIELD(has_original);
  FIELD(original);
  FIELD(reassigned);
  FIELD(x);
  FIELD(y);
  FIELD(z);
  FIELD(time_ns);
  FIELD(energy_mev);
  FIELD(electrons);
  FIELD(has_segment);
  FIELD(has_pdg);
  FIELD(pdg);
  FIELD(start_x);
  FIELD(start_y);
  FIELD(start_z);
  FIELD(end_x);
  FIELD(end_y);
  FIELD(end_z);
  FIELD(start_ns);
  FIELD(end_ns);
  FIELD(photons);
  FIELD(yield_ratio);

  FIELD(original_particle);
  FIELD(assigned_particle);
  FIELD(attributed_particle);
  FIELD(primary);
  FIELD(attribution);
#undef FIELD
#define FIELD(name) s.events->Branch(#name, &s.name)
  FIELD(run);
  FIELD(subrun);
  FIELD(event);
  FIELD(collection);
  FIELD(product);
  FIELD(transport);
  FIELD(kind);
  FIELD(count);
  FIELD(counts);
  FIELD(channels);
  FIELD(energies);
  FIELD(invalid_energy);
  FIELD(invalid_position);
  FIELD(outside_grid);
  FIELD(outside_energy);
  FIELD(invalid_position_energy);
#undef FIELD
  TNamed("schema", "deposit-v1; cm,ns,MeV,electrons; TDC is unconverted; "
                   "ReadoutIDE energies are channel sums, not plane-independent physical energy")
      .Write();
}

DepositArchive::~DepositArchive() = default;

void DepositArchive::write(const EventRecord &event) {
  auto &s = *impl_;
  if (!s.file->IsOpen())
    throw std::logic_error("Cannot write closed deposit archive");
  s.run = event.id.run;
  s.subrun = event.id.subrun;
  s.event = event.id.event;
  for (const auto &collection : event.deposits) {
    s.collection = collection.name;
    s.product = collection.product;
    s.transport = collection.transport;
    s.kind = static_cast<int>(collection.kind);
    s.count = collection.deposits.size();
    s.channels.assign(collection.channel_inventory.begin(), collection.channel_inventory.end());
    const auto &report = collection.report;
    s.counts.assign(report.counts.begin(), report.counts.end());
    s.energies.assign(report.energy_mev.begin(), report.energy_mev.end());
    s.invalid_energy = report.invalid_energy;
    s.invalid_position = report.invalid_position;
    s.outside_grid = report.outside_grid;
    s.outside_energy = report.outside_energy_mev;
    s.invalid_position_energy = report.invalid_position_energy_mev;
    if (s.events->Fill() < 0)
      throw std::runtime_error("Cannot write deposit accounting");
    for (const auto &d : collection.deposits) {
      s.source_index = d.source_index;
      s.channel_index = d.channel_index;
      s.ide_index = d.ide_index;
      s.channel = d.channel ? static_cast<long long>(*d.channel) : -1;
      s.tdc = d.tdc ? static_cast<long long>(*d.tdc) : -1;
      s.raw_track = d.raw_track_id;
      s.has_raw_original = d.raw_original_id.has_value();
      s.raw_original = d.raw_original_id.value_or(0);
      s.has_track = d.track_id.has_value();
      s.track = d.track_id.value_or(0);
      s.has_original = d.original_id.has_value();
      s.original = d.original_id.value_or(0);
      s.reassigned = d.reassigned;
      s.x = d.position_cm[0];
      s.y = d.position_cm[1];
      s.z = d.position_cm[2];
      s.time_ns = d.time_ns;
      s.energy_mev = d.energy_mev;
      s.electrons = d.electrons;
      s.has_segment = d.segment.has_value();
      const auto segment = d.segment.value_or(DepositSegment{});
      s.has_pdg = segment.pdg.has_value();
      s.pdg = segment.pdg.value_or(0);
      s.start_x = segment.start_cm[0];
      s.start_y = segment.start_cm[1];
      s.start_z = segment.start_cm[2];
      s.end_x = segment.end_cm[0];
      s.end_y = segment.end_cm[1];
      s.end_z = segment.end_cm[2];
      s.start_ns = segment.start_ns;
      s.end_ns = segment.end_ns;
      s.photons = segment.photons;
      s.yield_ratio = segment.scintillation_yield_ratio;
      s.original_particle = d.original_particle;
      s.assigned_particle = d.assigned_particle;
      s.attributed_particle = d.attributed_particle;
      s.primary = d.primary_ancestor;
      s.attribution = static_cast<int>(d.attribution);
      if (s.rows->Fill() < 0)
        throw std::runtime_error("Cannot write deposit contribution");
    }
  }
}

void DepositArchive::close() {
  auto &s = *impl_;
  if (!s.file->IsOpen())
    return;
  s.file->cd();
  if (s.file->Write() <= 0)
    throw std::runtime_error("Cannot finalize deposit archive");
  s.file->Close();
}
} // namespace supera::io::larcv2
