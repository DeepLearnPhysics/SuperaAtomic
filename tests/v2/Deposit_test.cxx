#include "supera/process/AttributeDeposits.h"
#include "supera/process/ResolveTransport.h"

#include <iostream>
#include <stdexcept>

namespace {
void check(bool value) {
  if (!value)
    throw std::runtime_error("Deposit attribution contract failed");
}
} // namespace

int main() try {
  using namespace supera;
  TransportCollection transport;
  transport.name = "truth";
  TransportParticle primary;
  primary.track_id = 1;
  primary.mother_id = 0;
  primary.process = "primary";
  TransportParticle child;
  child.track_id = 2;
  child.mother_id = 1;
  transport.particles = {primary, child};
  transport.suppression[1] = {99};
  resolve_transport(transport);
  DepositCollection input;
  input.transport = "truth";
  input.grid = Grid3D{{0., 0., 0.}, {2., 2., 2.}, {2, 2, 2}};
  Deposit exact;
  exact.track_id = 1;
  exact.original_id = 2;
  exact.energy_mev = 3.;
  exact.position_cm = {0.5, 0.5, 0.5};
  auto reassigned = exact;
  reassigned.original_id.reset();
  reassigned.reassigned = true;
  auto ancestor = exact;
  ancestor.track_id = 99;
  ancestor.original_id = 99;
  auto missing = exact;
  missing.track_id.reset();
  missing.original_id.reset();
  input.deposits = {exact, reassigned, ancestor, missing};
  attribute_deposits(input, transport);
  check(input.deposits[0].attribution == DepositAttribution::Exact);
  check(input.deposits[0].original_particle == 1 && input.deposits[0].assigned_particle == 0);
  check(input.deposits[1].attribution == DepositAttribution::Reassigned);
  check(input.deposits[2].attribution == DepositAttribution::Ancestor);
  check(input.deposits[2].original_particle == missing_index);
  check(input.deposits[3].attribution == DepositAttribution::Unresolved);
  check(input.voxels.size() == 3 && input.report.energy_mev[0] == 3.);
  double energy = 0.;
  for (const auto &v : input.voxels)
    energy += v.energy_mev;
  check(energy == 12.);

  // Readout can supply explicit reassignment evidence absent from the suppression map.
  DepositCollection readout;
  readout.transport = "truth";
  readout.kind = DepositKind::ReadoutIDE;
  auto signal = exact;
  signal.original_id = 100;
  signal.track_id = 2;
  signal.reassigned = true;
  readout.deposits = {signal, signal};
  const auto evidence = reassignment_evidence({readout}, "truth");
  check(evidence.at(100) == std::set<std::int32_t>{2});
  check(reassignment_evidence({readout}, "other").empty());
  auto full_steps = readout;
  full_steps.kind = DepositKind::EnergyStep;
  check(reassignment_evidence({full_steps}, "truth") == evidence);
  full_steps.kind = DepositKind::EnergyPoint;
  check(reassignment_evidence({full_steps}, "truth").empty());
  auto suppressed = exact;
  suppressed.original_id = 100;
  suppressed.track_id = 100;
  input.deposits = {suppressed};
  attribute_deposits(input, transport, evidence);
  check(input.deposits[0].attribution == DepositAttribution::EvidenceAncestor);
  check(input.deposits[0].original_particle == missing_index &&
        input.deposits[0].attributed_particle == 1);
  check(input.deposits[0].primary_ancestor == 0);
  // Multiple supported ancestors in one family fall back to the common primary.
  attribute_deposits(input, transport, {{100, {1, 2}}});
  check(input.deposits[0].attribution == DepositAttribution::EvidenceAncestor);
  check(input.deposits[0].attributed_particle == 0);
  input.deposits = {exact, reassigned, ancestor, missing};

  // Conflicting original and upstream assigned identities must not be silently reconciled.
  primary.track_id = 3;
  transport.particles.push_back(primary);
  resolve_transport(transport);
  input.deposits = {suppressed};
  attribute_deposits(input, transport, {{100, {2, 3}}});
  check(input.deposits[0].attribution == DepositAttribution::Conflict);
  check(input.deposits[0].attributed_particle == missing_index);
  input.deposits = {exact};
  input.deposits[0].track_id = 3;
  attribute_deposits(input, transport);
  check(input.deposits[0].attribution == DepositAttribution::Conflict);
  check(input.deposits[0].attributed_particle == missing_index);

  // Invalid measurements and geometric exclusions remain explicit accounting categories.
  input.deposits = {exact, exact, exact, exact};
  input.deposits[0].energy_mev = NAN;
  input.deposits[1].position_cm[0] = NAN;
  input.deposits[2].position_cm[0] = 3.;
  input.deposits[3].energy_mev = -1.;
  attribute_deposits(input, transport);
  check(input.report.invalid_energy == 2 && input.report.invalid_position == 1);
  check(input.report.outside_grid == 1 && input.report.outside_energy_mev == 3.);
  check(input.voxels.empty());
  input.kind = DepositKind::ReadoutIDE;
  bool rejected = false;
  try {
    attribute_deposits(input, transport);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  check(rejected);
  input.grid.reset();
  attribute_deposits(input, transport);
  input.deposits.clear();
  attribute_deposits(input, transport);
  check(input.report.counts[0] == 0 && input.voxels.empty());
  input.transport = "wrong";
  rejected = false;
  try {
    attribute_deposits(input, transport);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  check(rejected);

} catch (const std::exception &error) {

  std::cerr << error.what() << '\n';
  return 1;
}
