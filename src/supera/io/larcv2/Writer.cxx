#include "supera/io/larcv2/Writer.h"
#include "AuxDet.h"
#include "Generator.h"
#include "larcv/core/DataFormat/EventCRTHit.h"
#include "larcv/core/DataFormat/EventNeutrino.h"
#include "larcv/core/DataFormat/EventParticle.h"

#include "larcv/core/DataFormat/EventFlash.h"
#include "larcv/core/DataFormat/IOManager.h"

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
  bool open{false};
};

Writer::Writer(const std::string &filename, const std::vector<std::string> &collections,
               const std::vector<std::string> &auxiliary_collections,
               const std::vector<std::string> &generator_collections)
    : impl_(std::make_unique<Impl>()) {
  // Validate the complete output inventory before opening a file.
  if (filename.empty() ||
      (collections.empty() && auxiliary_collections.empty() && generator_collections.empty()))
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

  // Open the sink only after configuration validation succeeds.
  impl_->auxiliary_collections = auxiliary_collections;
  impl_->collections = collections;
  impl_->generator_collections = generator_collections;
  impl_->manager.set_out_file(filename);
  if (!impl_->manager.initialize())
    throw std::runtime_error("Cannot initialize LArCV2 output");
  impl_->open = true;

  // Register every product, including collections that may stay empty.
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
      if (flash.id >= larcv::kINVALID_INSTANCEID)
        throw std::invalid_argument("Flash ID cannot be represented in LArCV2");
      flashes.emplace_back(flash.time_us, flash.time_width_us, flash.absolute_time, flash.frame,
                           flash.pe, flash.in_beam_frame, flash.on_beam_time, flash.fast_to_total,
                           flash.x_cm, flash.x_width_cm, flash.y_cm, flash.y_width_cm, flash.z_cm,
                           flash.z_width_cm, flash.wire_centers, flash.wire_widths, flash.id);
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

  if (!impl_->manager.save_entry())
    throw std::runtime_error("Cannot save LArCV2 event");
}

void Writer::close() {
  if (impl_->open) {
    impl_->open = false;
    impl_->manager.finalize();
  }
}
} // namespace supera::io::larcv2
