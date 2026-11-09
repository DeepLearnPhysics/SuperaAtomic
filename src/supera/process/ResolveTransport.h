#ifndef SUPERA_PROCESS_RESOLVETRANSPORT_H
#define SUPERA_PROCESS_RESOLVETRANSPORT_H

#include "supera/data/Transport.h"

namespace supera {

/// Resolve source relationships without removing records or inventing missing parents.
/// Duplicate/absent particle identities are rejected; ambiguous links remain explicit.
/// A source primary requires both the primary process and a configured root mother ID.
void resolve_transport(TransportCollection &collection);

} // namespace supera
#endif
