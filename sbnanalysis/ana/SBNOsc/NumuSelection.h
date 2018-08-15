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

#include "TH1D.h"

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
    geoalgo::AABox active_volume; //!< Active volume
    double vertexDistanceCut; //!< Value of max distance [cm] between truth and reconstructed vertex. Will not apply cut if value is negative.
    bool verbose; //!< Whether to print out info associated w/ selection.
    double minLength; //!< Minimum length [cm] of contained leptons. Will not apply cut if value is negative.
  };

  /** Histograms made for output */
  struct RootOut {
    TH1D *h_numu_ccqe; //!< histogram w/ CCQE energy veriable
    TH1D *h_numu_trueE; //!< histogram w/ truth energy variable
    TH2D *h_numu_Vxy; //!< 2D x-y vertex histogram
    TH2D *h_numu_Vxz; //!< 2D x-z vertex histogram
    TH2D *h_numu_Vyz; //!< 2D y-z vertex histogram
  };

  /** Returns whether to apply FV cut on neutrino */
  bool passFV(double x, double y, double z);
  /** Applies reco-truth vertex matching cut */
  bool passRecoVertex(double truth_v[3], double reco_v[3]);
  /** Applies truth length cut */
  bool passMinLength(double length, bool stop_in_tpc);
  /** Run Selection on a neutrino */
  std::vector<bool> Select(const gallery::Event& ev, const simb::MCTruth& mctruth, unsigned truth_ind);

  unsigned fEventCounter;  //!< Count processed events
  unsigned fNuCount;  //!< Count selected events

  static const unsigned nCuts = 4; //!< number of cuts
  /** Names of cuts */
  static const std::vector<std::string> cutNames() {
    return {"CC", "FV", "min_L", "reco_V"};
  }

  Config _config; //!< The config
  RootOut _root_out[nCuts];
};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NumuSelection__

