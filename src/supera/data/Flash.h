/**
 * \file Flash.h
 *
 * \ingroup base
 *
 * \brief Class for light information, based on larcv::Flash
 *
 * @author Sindhu
 */

/** \addtogroup base

    @{*/
#ifndef __SUPERA_FLASH_H__
#define __SUPERA_FLASH_H__

#include "supera/base/SuperaType.h"
#include <vector>

namespace supera {
    /**
       \class Flash
       \brief Optical Flash (like recob::OpFlash)
    */
    class Flash {

    public:

        /// Default constructor
        Flash(double time=0, double timeWidth=0, double absTime=0, unsigned int frame=0,
            std::vector<double> PEPerOpDet=std::vector<double>(0),
            bool inBeamFrame=0, int onBeamTime=0, double fastToTotal=1,
            double xCenter=0, double xWidth=0,
            double yCenter=0, double yWidth=0,
            double zCenter=0, double zWidth=0,
            std::vector<double> wireCenters = std::vector<double>(0),
            std::vector<double> wireWidths = std::vector<double>(0),
            InstanceID_t id = kINVALID_INDEX,
            unsigned int tpc = 111);

        /// Default destructor
        virtual ~Flash() {}

        double TotalPE() const;

    public:
        // See also https://internal.dunescience.org/doxygen/OpFlash_8h_source.html
        // Mirroring it to facilitate Larsoft recob::OpFlash -> larcv -> CAF translation

        InstanceID_t   _id; ///< "ID" of this flash in FlashSet collection
        double _time; ///< Time on @ref DetectorClocksHardwareTrigger "trigger time scale" [us]
        unsigned int _tpc; ///< ID of the hit TPC
        double _timeWidth; ///< Width of the flash in time [us]
        double _absTime; ///< Time by PMT readout clock
        unsigned int _frame; ///< Frame number
        std::vector<double> _PEPerOpDet; ///< Number of PE on each PMT
        std::vector<double> _wireCenters; ///< Geometric center in each view
        std::vector<double> _wireWidths; ///< Geometric width in each view
        double _xCenter;
        double _xWidth;
        double _yCenter;
        double _yWidth;
        double _zCenter;
        double _zWidth;
        double _fastToTotal; ///< Fast to total light ratio
        bool _inBeamFrame; ///< Is this in the beam frame?
        int _onBeamTime; ///< Is this in time with beam?
    };

}

#endif
/** @} */ // end of doxygen group
