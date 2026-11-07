#include "supera/io/larcv2/Writer.h"

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
  bool open{false};
};

Writer::Writer(const std::string &filename, const std::vector<std::string> &collections)
    : impl_(std::make_unique<Impl>()) {
  if (filename.empty() || collections.empty())
    throw std::invalid_argument("Output filename and collections must not be empty");
  std::set<std::string> names;
  for (const auto &name : collections) {
    if (name.empty() || !names.insert(name).second)
      throw std::invalid_argument("Output collection names must be nonempty and unique");
  }
  impl_->collections = collections;
  impl_->manager.set_out_file(filename);
  if (!impl_->manager.initialize())
    throw std::runtime_error("Cannot initialize LArCV2 output");
  impl_->open = true;
  for (const auto &name : collections)
    impl_->manager.get_data<larcv::EventFlash>(name);
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
  // LArCV2 stores event IDs as size_t. Do not silently narrow on other platforms.
  for (const auto value : {event.id.run, event.id.subrun, event.id.event})
    if (value >= std::numeric_limits<std::size_t>::max())
      throw std::invalid_argument("Event ID cannot be represented in LArCV2");
  impl_->manager.set_id(event.id.run, event.id.subrun, event.id.event);
  for (std::size_t i = 0; i < converted.size(); ++i)
    impl_->manager.get_data<larcv::EventFlash>(impl_->collections[i])
        .emplace(std::move(converted[i]));
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
