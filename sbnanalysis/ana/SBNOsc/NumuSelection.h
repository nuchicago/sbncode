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

// take the geobox stuff from uboonecode
#include "uboone/LLBasicTool/GeoAlgo/GeoAABox.h"

class TH2D;

namespace ana {
  namespace SBNOsc {

/**
 * \class NumuSelection
 * \brief Electron neutrino event selection
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
  /** Configuration parameters */
  struct Config {
    art::InputTag truthTag; //!< art tag for MCTruth information
    bool doFVCut; //!< Whether to apply fiducial volume cut
    std::vector<geoalgo::AABox> aaBoxes; //!< List of FV containers -- set by "fiducial_volumes"
  };

  /** Returns whether to apply FV cut on neutrino */
  bool passFV(double x, double y, double z);
  /** Run Selection on a neutrino */
  bool Select(const simb::MCNeutrino& nu);

  unsigned fEventCounter;  //!< Count processed events
  unsigned fNuCount;  //!< Count selected events

  Config _config; //!< The config
};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NumuSelection__

