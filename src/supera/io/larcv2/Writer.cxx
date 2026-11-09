#include "supera/io/larcv2/Writer.h"
#include "AuxDet.h"
#include "Deposit.h"
#include "Generator.h"
#include "SpacePoint.h"
#include "Transport.h"
#include "larcv/core/DataFormat/EventCRTHit.h"
#include "larcv/core/DataFormat/EventNeutrino.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

#include "larcv/core/DataFormat/EventFlash.h"
#include "larcv/core/DataFormat/IOManager.h"

#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>

namespace supera::io::larcv2 {
class Writer::Impl {
public:
  larcv::IOManager manager{larcv::IOManager::kWRITE};
  std::vector<std::string> collections;
  std::vector<std::string> auxiliary_collections;
  std::vector<std::string> generator_collections;
  std::vector<SpacePointOutput> spacepoint_collections;
  std::vector<std::string> transport_collections;
  std::unique_ptr<TransportArchive> archive;
  std::vector<DepositOutput> deposit_collections;
  std::unique_ptr<DepositArchive> deposit_archive;
  bool open{false};
};

Writer::Writer(const std::string &filename, const std::vector<std::string> &collections,
               const std::vector<std::string> &auxiliary_collections,
               const std::vector<std::string> &generator_collections,
               const std::vector<SpacePointOutput> &spacepoint_collections,
               const std::vector<std::string> &transport_collections,
               const std::vector<DepositOutput> &deposit_collections)
    : impl_(std::make_unique<Impl>()) {
  // Validate the complete output inventory before opening a file.
  if (filename.empty() || (collections.empty() && auxiliary_collections.empty() &&
                           generator_collections.empty() && spacepoint_collections.empty() &&
                           transport_collections.empty() && deposit_collections.empty()))
    throw std::invalid_argument("Output filename and collections must not be empty");

  std::set<std::string> names;
  for (const auto &name : collections) {
    if (name.empty() || !names.insert(name).second)
      throw std::invalid_argument("Output collection names must be nonempty and unique");
  }

  names.clear();
  for (const auto &name : auxiliary_collections) {
    if (name.empty() || !names.insert(name).second)
      throw std::invalid_argument("Auxiliary output names must be nonempty and unique");
  }

  names.clear();
  for (const auto &name : generator_collections) {
    if (name.empty() || !names.insert(name).second)
      throw std::invalid_argument("Generator output names must be nonempty and unique");
  }

  for (const auto &name : transport_collections) {
    if (name.empty() || !names.insert(name).second)
      throw std::invalid_argument(
          "Transport output names must be unique, including generator outputs");
  }
  impl_->transport_collections = transport_collections;

  names.clear();
  for (const auto &configuration : spacepoint_collections)
    for (const auto &label : spacepoint_labels(configuration))
      if (!names.insert(label).second)
        throw std::invalid_argument("Duplicate sparse output label");
  impl_->spacepoint_collections = spacepoint_collections;

  std::set<std::string> deposit_names;
  for (const auto &definition : deposit_collections) {
    if (definition.name.empty() || !deposit_names.insert(definition.name).second)
      throw std::invalid_argument("Deposit output names must be nonempty and unique");
    if (definition.voxelize && !names.insert(definition.name).second)
      throw std::invalid_argument("Deposit tensor conflicts with another output label");
  }
  impl_->deposit_collections = deposit_collections;

  // Open the sink only after configuration validation succeeds.
  impl_->auxiliary_collections = auxiliary_collections;
  impl_->collections = collections;
  impl_->generator_collections = generator_collections;
  impl_->manager.set_out_file(filename);
  if (!impl_->manager.initialize())
    throw std::runtime_error("Cannot initialize LArCV2 output");
  impl_->open = true;

  if (!transport_collections.empty())
    impl_->archive = std::make_unique<TransportArchive>(filename + ".transport.root");
  for (const auto &name : transport_collections)
    impl_->manager.get_data<larcv::EventParticle>(name);

  if (!deposit_collections.empty())
    impl_->deposit_archive = std::make_unique<DepositArchive>(filename + ".deposits.root");
  for (const auto &definition : deposit_collections)
    if (definition.voxelize)
      impl_->manager.get_data<larcv::EventClusterVoxel3D>(definition.name);

  // Register every product, including collections that may stay empty.
  for (const auto &label : names)
    impl_->manager.get_data<larcv::EventSparseTensor3D>(label);
  for (const auto &name : collections)
    impl_->manager.get_data<larcv::EventFlash>(name);
  for (const auto &name : auxiliary_collections)
    impl_->manager.get_data<larcv::EventCRTHit>(name);
  for (const auto &name : generator_collections) {
    impl_->manager.get_data<larcv::EventParticle>(name);
    impl_->manager.get_data<larcv::EventNeutrino>(name);
  }
}

Writer::~Writer() {
  // Explicit close() is the checked path; cleanup during unwinding is best effort.
  try {
    close();
  } catch (...) {
  }
}

void Writer::write(const EventRecord &event) {
  if (!impl_->open)
    throw std::logic_error("Cannot write after close");
  if (event.optical.size() != impl_->collections.size())
    throw std::invalid_argument("Event optical collection inventory differs from configuration");

  // Validate and translate the entire event before modifying the output manager.
  std::vector<std::vector<larcv::Flash>> converted;
  converted.reserve(event.optical.size());
  for (std::size_t i = 0; i < event.optical.size(); ++i) {
    const auto &collection = event.optical[i];
    if (collection.name != impl_->collections[i])
      throw std::invalid_argument("Event optical collection order/name differs from configuration");

    auto &flashes = converted.emplace_back();
    flashes.reserve(collection.flashes.size());

    for (const auto &flash : collection.flashes) {
      if (flash.id < 0 || static_cast<std::uint64_t>(flash.id) >= larcv::kINVALID_INSTANCEID)
        throw std::invalid_argument("Flash ID cannot be represented in LArCV2");
      // Validate canonical values before any output conversion or event mutation.
      for (const auto value : {flash.time_us, flash.time_width_us, flash.absolute_time,
                               flash.fast_to_total, flash.x_cm, flash.x_width_cm, flash.y_cm,
                               flash.y_width_cm, flash.z_cm, flash.z_width_cm})
        require_measurement(value);
      for (const auto *values : {&flash.pe, &flash.wire_centers, &flash.wire_widths})
        for (const auto value : *values)
          require_measurement(value);

      flashes.emplace_back(flash.time_us, flash.time_width_us, flash.absolute_time, flash.frame,
                           flash.pe, flash.in_beam_frame, flash.on_beam_time, flash.fast_to_total,
                           std::isnan(flash.x_cm) ? std::numeric_limits<double>::max() : flash.x_cm,
                           std::isnan(flash.x_width_cm) ? std::numeric_limits<double>::max()
                                                        : flash.x_width_cm,
                           flash.y_cm, flash.y_width_cm, flash.z_cm, flash.z_width_cm,
                           flash.wire_centers, flash.wire_widths, flash.id);
    }
  }

  // Prepare auxiliary products without modifying the output manager.
  if (event.auxiliary.size() != impl_->auxiliary_collections.size())
    throw std::invalid_argument("Event auxiliary inventory differs from configuration");

  std::vector<std::vector<larcv::CRTHit>> auxiliary;
  for (std::size_t i = 0; i < event.auxiliary.size(); ++i) {
    if (event.auxiliary[i].name != impl_->auxiliary_collections[i])
      throw std::invalid_argument("Event auxiliary collection name/order mismatch");
    auxiliary.push_back(project_auxiliary(event.auxiliary[i]));
  }

  // Project selected generator records before committing any product.
  if (event.generator.size() != impl_->generator_collections.size())
    throw std::invalid_argument("Event generator inventory differs from configuration");

  std::vector<LegacyGeneratorProducts> generator;
  for (std::size_t i = 0; i < event.generator.size(); ++i) {
    if (event.generator[i].name != impl_->generator_collections[i])
      throw std::invalid_argument("Event generator collection name/order mismatch");
    generator.push_back(project_legacy_generator(event.generator[i]));
  }

  if (event.spacepoints.size() != impl_->spacepoint_collections.size())
    throw std::invalid_argument("Space-point output inventory mismatch");
  std::vector<std::vector<larcv::SparseTensor3D>> spacepoints;
  for (std::size_t i = 0; i < event.spacepoints.size(); ++i)
    spacepoints.push_back(
        project_spacepoints(event.spacepoints[i], impl_->spacepoint_collections[i]));

  if (event.transport.size() != impl_->transport_collections.size())
    throw std::invalid_argument("Transport output inventory mismatch");
  std::vector<std::vector<larcv::Particle>> transport;
  for (std::size_t i = 0; i < event.transport.size(); ++i) {
    if (event.transport[i].name != impl_->transport_collections[i])
      throw std::invalid_argument("Transport output name/order mismatch");
    transport.push_back(project_transport(event.transport[i]));
  }

  if (event.deposits.size() != impl_->deposit_collections.size())
    throw std::invalid_argument("Deposit output inventory mismatch");
  std::vector<std::optional<DepositProducts>> deposits;
  for (std::size_t i = 0; i < event.deposits.size(); ++i) {
    const auto &definition = impl_->deposit_collections[i];
    if (event.deposits[i].name != definition.name)
      throw std::invalid_argument("Deposit output name/order mismatch");
    if (definition.voxelize)
      deposits.emplace_back(project_deposits(event.deposits[i]));
    else
      deposits.emplace_back();
  }

  // LArCV2 stores event IDs as size_t. Do not silently narrow on other platforms.
  for (const auto value : {event.id.run, event.id.subrun, event.id.event})
    if (value >= std::numeric_limits<std::size_t>::max())
      throw std::invalid_argument("Event ID cannot be represented in LArCV2");

  // Commit the fully validated event, replacing even empty product contents.
  impl_->manager.set_id(event.id.run, event.id.subrun, event.id.event);
  for (std::size_t i = 0; i < converted.size(); ++i)
    impl_->manager.get_data<larcv::EventFlash>(impl_->collections[i])
        .emplace(std::move(converted[i]));
  for (std::size_t i = 0; i < auxiliary.size(); ++i)
    impl_->manager.get_data<larcv::EventCRTHit>(impl_->auxiliary_collections[i])
        .emplace(std::move(auxiliary[i]));

  for (std::size_t i = 0; i < generator.size(); ++i) {
    impl_->manager.get_data<larcv::EventParticle>(impl_->generator_collections[i])
        .emplace(std::move(generator[i].particles));
    impl_->manager.get_data<larcv::EventNeutrino>(impl_->generator_collections[i])
        .emplace(std::move(generator[i].neutrinos));
  }

  for (std::size_t i = 0; i < spacepoints.size(); ++i) {
    const auto labels = spacepoint_labels(impl_->spacepoint_collections[i]);
    for (std::size_t j = 0; j < labels.size(); ++j) {
      auto &tensor = spacepoints[i][j];
      const auto meta = tensor.meta();
      impl_->manager.get_data<larcv::EventSparseTensor3D>(labels[j]).emplace(std::move(tensor),
                                                                             meta);
    }
  }

  for (std::size_t i = 0; i < transport.size(); ++i)
    impl_->manager.get_data<larcv::EventParticle>(impl_->transport_collections[i])
        .emplace(std::move(transport[i]));
  if (impl_->archive)
    impl_->archive->write(event);

  for (std::size_t i = 0; i < deposits.size(); ++i) {
    if (!deposits[i])
      continue;
    auto &products = *deposits[i];
    const auto &name = impl_->deposit_collections[i].name;
    const auto meta = products.clusters.meta();
    impl_->manager.get_data<larcv::EventClusterVoxel3D>(name).emplace(std::move(products.clusters),
                                                                      meta);
    impl_->manager.get_data<larcv::EventSparseTensor3D>(name).emplace(std::move(products.tensor),
                                                                      meta);
  }
  if (impl_->deposit_archive)
    impl_->deposit_archive->write(event);

  if (!impl_->manager.save_entry())
    throw std::runtime_error("Cannot save LArCV2 event");
}

void Writer::close() {
  if (impl_->open) {
    impl_->open = false;
    impl_->manager.finalize();
    if (impl_->deposit_archive)
      impl_->deposit_archive->close();
    if (impl_->archive)
      impl_->archive->close();
  }
}
} // namespace supera::io::larcv2
