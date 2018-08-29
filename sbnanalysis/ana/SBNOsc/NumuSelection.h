#ifndef __sbnanalysis_ana_SBNOsc_NumuSelection__
#define __sbnanalysis_ana_SBNOsc_NumuSelection__

/**
 * \file NumuSelection.h
 *
 * SBN nue selection.
 *
 * Author: 
 */

#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"

class TH2D;

namespace ana {
namespace SBNOsc {

/**
 * \class NumuSelection
 * \brief Muon neutrino event selection
 */

class NumuSelection : public core::SelectionBase {
    public:
        /** Constructor. */
        NumuSelection();

        /**
         * Initialization.
         *
         * \param config A configuration, as a JSON object
         */
        void Initialize(Json::Value* config=NULL);

        /** Finalize and write objects to the output file. */
        void Finalize();

        /**
         * Process one event.
         *
         * \param ev A single event, as a gallery::Event
         * \param Reconstructed interactions
         * \return True to keep event
         */
        bool ProcessEvent(const gallery::Event& ev, std::vector<Event::Interaction>& reco);

    protected:
        unsigned fEventCounter; //!< Count processed events
        unsigned fNuAll;        //!< Count all events
        unsigned fNuCount;      //!< Count selected events
        unsigned fNuinFid;      //!< Count nus that reacted in the fiducial volume

        /** Configuration parameters */
        art::InputTag fTruthTag;    //!< art tag for MCTruth information
        art::InputTag fPartTag;     // for use with particle (non-nu) information
        art::InputTag fEventWgtTag;
        art::InputTag fFluxWgtTag;
    
        /** My parameters */
        std::string fDet;
        double ml, mn, mp, Eb;
    
};

} // namespace SBNOsc
} // namespace ana

#endif // __sbnanalysis_ana_SBNOsc_NumuSelection__
