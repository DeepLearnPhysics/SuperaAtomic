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
            std::vector<double> PEPerOpDet=std::vector<double>(0), int tpc=-1,
            bool inBeamFrame=0, int onBeamTime=0, double fastToTotal=1,
            double xCenter=0, double xWidth=0,
            double yCenter=0, double yWidth=0,
            double zCenter=0, double zWidth=0,
            std::vector<double> wireCenters = std::vector<double>(0),
            std::vector<double> wireWidths = std::vector<double>(0),
            InstanceID_t id = kINVALID_INDEX);

        /// Default destructor
        virtual ~Flash() {}

        double TotalPE() const;
        
    public:
        
        /// Getters
        InstanceID_t id; ///< "ID" of this flash in FlashSet collection
        double time; ///< Time on @ref DetectorClocksHardwareTrigger "trigger time scale" [us]
        int tpc;
        double timeWidth; ///< Width of the flash in time [us]
        double absTime; ///< Time by PMT readout clock
        unsigned int frame;  ///< Frame number
        const std::vector<double>& PEPerOpDet; ///< Number of PE on each PMT
        const std::vector<double>& wireCenters; ///< Geometric center in each view
        const std::vector<double>& wireWidths; ///< Geometric width in each view
        double xCenter;
        double xWidth;
        double yCenter;
        double yWidth;
        double zCenter;
        double zWidth;
        double fastToTotal; ///< Fast to total light ratio
        bool inBeamFrame; ///< Is this in the beam frame?
        int onBeamTime; ///< Is this in time with beam?

    };


#endif
/** @} */ // end of doxygen group
