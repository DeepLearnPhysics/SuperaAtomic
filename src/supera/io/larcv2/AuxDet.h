#ifndef SUPERA_IO_LARCV2_AUXDET_H
#define SUPERA_IO_LARCV2_AUXDET_H
// Private output implementation: never installed as part of the shared model.
#include "larcv/core/DataFormat/CRTHit.h"
#include "supera/data/AuxDet.h"

namespace supera::io::larcv2 {
std::vector<larcv::CRTHit> project_auxiliary(const AuxDetSpacePointCollection &collection);
}
#endif
