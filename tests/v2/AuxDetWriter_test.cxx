#include "TFile.h"
#include "TTree.h"
#include "larcv/core/DataFormat/EventCRTHit.h"
#include "supera/io/larcv2/Writer.h"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool value) {
  if (!value)
    throw std::runtime_error("Auxiliary writer contract failed");
}

template <class F> void rejects(F &&action) {
  try {
    action();
  } catch (const std::invalid_argument &) {
    return;
  }
  throw std::runtime_error("Accepted an unsupported auxiliary projection");
}

supera::AuxDetSpacePoint point() {
  using namespace supera;
  AuxDetSpacePoint result;
  result.id = 7;
  result.position = AuxDetPosition{{-10., 20., 30.}, std::array<double, 3>{1., 2., 3.}};
  result.signal = AuxDetSignal{42.5, SignalUnit::Photoelectron};
  result.plane = -1;
  result.region = "panel";
  result.readout_ids = {255, 0, 255};
  result.readouts = {
      {0, {}},
      {255, {{-1, {-0., SignalUnit::Photoelectron}}, {31, {2.5, SignalUnit::Photoelectron}}}}};
  result.times = {
      {"ts0_s",
       (std::uint64_t{1} << 53) + 1,
       TimeUnit::Second,
       TimeReference::SourceDefined,
       "sbn_crt_source",
       {}},
      {"ts0_s_corr", 1.25, TimeUnit::Second, TimeReference::SourceDefined, "sbn_crt_source", {}},
      {"ts0_ns", -2.5, TimeUnit::Nanosecond, TimeReference::SourceDefined, "sbn_crt_source", {}},
      {"ts0_ns_corr",
       3.75,
       TimeUnit::Nanosecond,
       TimeReference::SourceDefined,
       "sbn_crt_source",
       {}},
      {"ts1_ns", -4.5, TimeUnit::Nanosecond, TimeReference::SourceDefined, "sbn_crt_source", {}}};
  return result;
}
} // namespace

int main() try {
  using namespace supera;
  using supera::io::larcv2::Writer;
  const std::string filename = "auxiliary-writer-contract.root";
  rejects([&] { Writer writer(filename, {}, {""}); });
  rejects([&] { Writer writer(filename, {}, {"aux", "aux"}); });

  // Arrange a supported record with a non-cosmic subsystem role.
  EventRecord event;
  event.id = {1, 2, 3};
  event.auxiliary.push_back({"aux", "test", {}, {point()}});
  event.auxiliary[0].conventions.source_schema = AuxDetSourceSchema::SbnCrtHit;
  // The projection is a measurement contract, not a cosmic-role restriction.
  event.auxiliary[0].conventions.roles = {AuxDetRole::BeamTiming};

  // Rejected projections must not leave partial events in the output.
  {
    Writer writer(filename, {}, {"aux"});
    auto reject_change = [&](auto change) {
      auto invalid = event;
      change(invalid.auxiliary.front().points.front());
      rejects([&] { writer.write(invalid); });
    };

    // Reject unsupported units and numeric narrowing.
    reject_change([](auto &p) { p.signal->unit = SignalUnit::ADC; });
    reject_change([](auto &p) { p.position->cm[0] = 0.1; });
    reject_change([](auto &p) { p.position->cm[0] = std::numeric_limits<double>::max(); });

    // Missing estimates and unsupported completeness cannot be fabricated or dropped.

    reject_change([](auto &p) { p.complete = false; });

    // Hardware mappings must fit the legacy identity and signal contract.
    reject_change([](auto &p) { p.readout_ids.push_back(256); });
    reject_change([](auto &p) { p.readouts.push_back(p.readouts.front()); });
    reject_change([](auto &p) {
      p.readouts.back().channels[0].channel = std::numeric_limits<std::int64_t>::max();
    });
    reject_change([](auto &p) { p.readouts.back().channels[0].signal.unit = SignalUnit::MeV; });
    reject_change([](auto &p) { p.id = supera::missing_index; });

    // Clock names, types, units, and references form one explicit contract.
    reject_change([](auto &p) { p.times.pop_back(); });
    reject_change([](auto &p) { p.times[1] = p.times[0]; });
    reject_change([](auto &p) { p.times[0].uncertainty = 1.; });
    reject_change([](auto &p) { p.times[0].value = 1.; });
    reject_change([](auto &p) { p.times[0].unit = TimeUnit::Tick; });
    reject_change([](auto &p) { p.times[0].clock = "other"; });
    reject_change([](auto &p) { p.times[0].reference = TimeReference::Trigger; });

    // Validate collection inventory independently of individual point contents.
    auto invalid = event;
    invalid.auxiliary.clear();
    rejects([&] { writer.write(invalid); });

    invalid = event;
    invalid.auxiliary[0].name = "other";
    rejects([&] { writer.write(invalid); });

    invalid = event;
    invalid.auxiliary[0].conventions.source_schema = AuxDetSourceSchema::Generic;
    rejects([&] { writer.write(invalid); });

    // Valid events and canonical missing estimates must both be representable.
    writer.write(event);

    event.id.event = 4;
    event.auxiliary.front().points.clear();
    writer.write(event);
    event.id.event = 5;
    auto missing = point();
    missing.position.reset();
    missing.signal.reset();
    missing.plane.reset();
    missing.times[0].value = std::monostate{};
    missing.times[1].value = missing_value;
    event.auxiliary.front().points = {missing};
    writer.write(event);
    writer.close();
  }

  // Read back through LArCV directly, without using the conversion under test.
  {
    TFile file(filename.c_str(), "READ");
    auto *tree = file.Get<TTree>("crthit_aux_tree");
    require(tree && tree->GetEntries() == 3);
    larcv::EventCRTHit storage;
    auto *record = &storage;
    tree->SetBranchAddress("crthit_aux_branch", &record);
    require(tree->GetEntry(0) > 0 && record->as_vector().size() == 1);
    const auto &hit = record->as_vector().front();
    require(record->run() == 1 && record->subrun() == 2 && record->event() == 3);
    require(hit.id() == 7 && hit.peshit() == 42.5F && hit.plane() == -1 && hit.tagger() == "panel");
    require(hit.ts0_s() == (std::uint64_t{1} << 53) + 1 && hit.ts0_s_corr() == 1.25 &&
            hit.ts0_ns() == -2.5 && hit.ts0_ns_corr() == 3.75 && hit.ts1_ns() == -4.5);
    require(hit.x_pos() == -10.F && hit.y_pos() == 20.F && hit.z_pos() == 30.F &&
            hit.x_err() == 1.F && hit.y_err() == 2.F && hit.z_err() == 3.F);
    require(hit.feb_id() == std::vector<std::uint8_t>{255, 0, 255});
    const auto &signals = hit.pesmap();
    require(signals.size() == 2 && signals.at(0).empty() && signals.at(255).size() == 2 &&
            signals.at(255)[0].first == -1 && std::signbit(signals.at(255)[0].second) &&
            signals.at(255)[1] == std::make_pair(31, 2.5F));

    // The second event must clear the previous point.
    require(tree->GetEntry(1) > 0 && record->event() == 4 && record->as_vector().empty());
    require(tree->GetEntry(2) > 0 && record->as_vector().size() == 1);
    const auto &missing = record->as_vector().front();
    require(missing.x_pos() == std::numeric_limits<float>::lowest() &&
            missing.y_pos() == std::numeric_limits<float>::lowest() &&
            missing.x_err() == std::numeric_limits<float>::lowest() &&
            missing.peshit() == std::numeric_limits<float>::lowest() && missing.plane() == -1 &&
            missing.ts0_s() == std::numeric_limits<std::uint64_t>::max() &&
            missing.ts0_s_corr() == std::numeric_limits<double>::lowest());
    tree->ResetBranchAddresses();
  }

  std::filesystem::remove(filename);
  return 0;

} catch (const std::exception &error) {

  std::cerr << error.what() << '\n';
  return 1;
}
