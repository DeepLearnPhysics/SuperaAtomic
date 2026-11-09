#ifndef SUPERA_PROCESS_ATTRIBUTEDEPOSITS_H
#define SUPERA_PROCESS_ATTRIBUTEDEPOSITS_H

#include "supera/data/Deposit.h"
#include "supera/data/Transport.h"

namespace supera {
using ReassignmentEvidence = std::map<std::int32_t, std::set<std::int32_t>>;

/// Gather explicit negative-ID step/readout reassignments, scoped to one transport namespace.
ReassignmentEvidence reassignment_evidence(const std::vector<DepositCollection> &sources,
                                           const std::string &transport);

/// Attribute without dropping deposits. A recovered ancestor is never called an exact creator.
/// Optional rasterization is only valid for physical energy records, not plane-repeated IDEs.
void attribute_deposits(DepositCollection &deposits, const TransportCollection &transport,
                        const ReassignmentEvidence &readout = {});
} // namespace supera
#endif
