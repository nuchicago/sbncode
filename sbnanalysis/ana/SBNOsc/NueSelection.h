#ifndef __sbnanalysis_ana_SBNOsc_NueSelection__
#define __sbnanalysis_ana_SBNOsc_NueSelection__

/**
 * \file NueSelection.h
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
class TH1D;

namespace ana {
  namespace SBNOsc {

/**
 * \class NueSelection
 * \brief Electron neutrino event selection
 */
class NueSelection : public core::SelectionBase {
public:
  /** Constructor. */
  NueSelection();

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
  unsigned fEventCounter;  //!< Count processed events
  unsigned fNuCount;  //!< Count selected events



  TH1D* fDiffLength;
  TH1D* fTrackLength;
  TH1D* fGenNueHist;
  TH1D* fGenNueFidVolHist;
  TH1D* fSelectedNuHist;
  TH1D* fShowerEnergy;
  TH1D* fEnergeticShowerHist;
  TH1D* fVisibleVertexNuEHist;
  TH1D* fCGSelectionHist;
  TH1D* fRecoSelectionHist;
  TH1D* fShowerCutSelectionHist;
  TH1D* fSelectedTrueNue;

  //shower true type
  TH1D* fMuShowerSelectedNu;
  TH1D* fEShowerSelectedNu;
  TH1D* fGammaShowerSelectedNu;
  TH1D* fOtherShowerSelectedNu;
  TH1D* fShowerdEdx;


  /** Configuration parameters */
  art::InputTag fTruthTag;  //!< art tag for MCTruth information
  art::InputTag fTrackTag; //! <art tag for MCTrack information
  art::InputTag fShowerTag; //! <art tag for MCShower information
  double fEnergyThreshold; //configurable parameter
};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NueSelection__
