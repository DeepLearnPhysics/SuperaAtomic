/**
 * \file TriggerMeta3D.h
 *
 * \ingroup base
 *
 * \brief Storage for metadata about trigger(s) that should be associated to image(s)
 *
 * @author J. Wolcott <jwolcott@fnal.gov>
 */

/** \addtogroup base
    @{*/
#ifndef SUPERA_TRIGGERMETA_H
#define SUPERA_TRIGGERMETA_H

#include <cstddef>

#include "supera/base/SuperaType.h"

namespace supera
{

  /**
     \class TriggerMeta
     @brief Metadata for triggers (id, time, type)
  */
  class TriggerMeta
  {
    public:
      /// Default constructor
      TriggerMeta(InstanceID_t id=kINVALID_INSTANCEID, unsigned long long time_s=kINVALID_ULONGLONG, unsigned int time_ns=kINVALID_UINT, int type=kINVALID_INT)
      : _id(id), _time_s(time_s), _time_ns(time_ns), _type(type)
      {}

      void clear()   { _id = kINVALID_INSTANCEID; _time_s = kINVALID_ULONGLONG; _time_ns = kINVALID_UINT; _type = kINVALID_INT; }

      // accessors
      InstanceID_t       id()      const { return _id; };
      unsigned long long time_s()  const { return _time_s; }
      unsigned int       time_ns() const { return _time_ns; }
      int                type()    const { return _type; }

      // setters
      void id(InstanceID_t id)               { _id = id; };
      void time_s(unsigned long long time_s) { _time_s = time_s; }
      void time_ns(unsigned int time_ns)     { _time_ns = time_ns; }
      void type(unsigned int typ)            { _type = typ; }

    private:
      InstanceID_t       _id;       ///< Trigger ID, if relevant
      unsigned long long _time_s;   ///< Integer seconds part of UNIX time this trigger was initiated
      unsigned int       _time_ns;  ///< Fractional part of UNIX time this trigger was initiated, measured in ns
      int                _type;     ///< DAQ-specific trigger type
  };

} // supera

#endif //SUPERA_TRIGGERMETA_H
/** @} */ // end of doxygen group

