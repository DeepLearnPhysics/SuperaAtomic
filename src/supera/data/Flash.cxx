#ifndef __SUPERA_FLASH_CXX__
#define __SUPERA_FLASH_CXX__

#include "Flash.h"
#include "supera/base/meatloaf.h"
#include <numeric> // std::accumulate()

namespace supera {

  std::string Flash::dump() const
  {
      std::stringstream ss;
      std::stringstream buf;
      ss  << "      \033[95m" << "Flash id = (" << id << ")\033[00m " << std::endl;
      buf << "      ";

      ss << buf.str() << "time = (" << time << ")" << std::endl
         << buf.str() << "tpc ID = (" << tpc ")" << std::endl;

      return ss.str();
    }

    double Flash::TotalPE() const {
        return std::accumulate(_PEPerOpDet.begin(), _PEPerOpDet.end(), 0.0);
    }


}
#endif
