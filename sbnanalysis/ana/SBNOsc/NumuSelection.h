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
#include "TDatabasePDG.h"

#include "nusimdata/SimulationBase/MCTruth.h"

// take the geobox stuff from uboonecode
#include "uboone/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoAlgo.h"

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

  /** Additional information used by the selection per neutrino interaction */
  struct NuMuInteraction {
    bool l_is_contained; //!< whether the lepton track is totally contained in the fiducial volume
    double l_contained_length; //!< the length of the lepton track contained in the fiducial volume
    double l_length; //!< total length of lepton track
  };

protected:
  /** Configuration parameters */
  struct Config {
    art::InputTag truthTag; //!< art tag for MCTruth
    art::InputTag mctTag; //!< art tag for MCTrack
    art::InputTag mcsTag; //!< art tag for MCShower
    art::InputTag mcpTag; //!< art tag for MCParticle
    bool doFVCut; //!< Whether to apply fiducial volume cut
    std::vector<geoalgo::AABox> fiducial_volumes; //!< List of FV containers -- set by "fiducial_volumes"
    geoalgo::AABox active_volume; //!< Active volume
    double vertexDistanceCut; //!< Value of max distance [cm] between truth and reconstructed vertex. Will not apply cut if value is negative.
    bool verbose; //!< Whether to print out info associated w/ selection.
    double minLengthContainedLepton; //!< Minimum length [cm] of contained leptons. Will not apply cut if value is negative.
    double minLengthExitingLepton; //!< Minimum length [cm] of exiting leptons.  Will not apply cut if value is negative.
    double containedLeptonESmear; //!< % to smear truth energy of contained CC leptons.
    double exitingLeptonESmear; //!< % to smear truth energy of exiting CC leptons.
    double hadronESmear; //!< % to smear truth energy of hadrons.
  };

  /** Histograms made for output */
  struct RootHistos {
    TH1D *h_numu_ccqe; //!< histogram w/ CCQE energy veriable
    TH1D *h_numu_trueE; //!< histogram w/ truth energy variable
    TH1D *h_numu_visibleE; //!< histogram w/ visible energy variable (total muon momentum + kinetic hadron energy)
    TH1D *h_numu_true_v_visibleE; //!< histogram w/ difference of visible and truth energy
    TH1D *h_numu_l_is_contained; //!< histogram w/ whether associated lepton is contained in FV 
    TH1D *h_numu_contained_L; //!< histogram w/ FV contained length of lepton in CC event
    TH1D *h_numu_l_length; //!< histogram w/ total length of associated lepton
    TH2D *h_numu_Vxy; //!< 2D x-y vertex histogram
    TH2D *h_numu_Vxz; //!< 2D x-z vertex histogram
    TH2D *h_numu_Vyz; //!< 2D y-z vertex histogram
  };

  /* Applies FV cut 
  * \param v The neutrino interaction vertex
  * \returns Whether to apply FV cut on neutrino 
  *
  * */
  bool passFV(const TVector3 &v) { return containedInFV(v); }
  /** Applies reco-truth vertex matching cut 
 * \param truth_v Truth vertex vector
 * \param reco_v Reconstructed vertex vector
 * \returns Whether to apply reco-truth vertex matching cut: true == passed cut
 * */
  bool passRecoVertex(const TVector3 &truth_v, const TVector3 &reco_v);
  /** Applies truth length cut 
 *  \param length Distance travelled by lepton
 *  \param stop_in_tpc Whether the lepton stopped in the TPC volume
 *  \returns Whether to apply length cut: true == passed cut
 * */
  bool passMinLength(double length, bool stop_in_tpc);
  /** Run Selection on a neutrino 
 *  \param ev The gallery event.
 *  \param mctruth The mctruth object associated with the currently considered interaction.
 *  \param truth_ind The index into the vector of MCTruth objects in the gallery event from which the "mctruth" object is
 *  \param intInfo The interaction info object associated with this mctruth object.
 *
 *  \returns A list containing whether the interaction passed each cut (ordered in the same way as the "cutNames()")
 * */
  std::vector<bool> Select(const gallery::Event& ev, const simb::MCTruth& mctruth, unsigned truth_ind, const NumuSelection::NuMuInteraction &intInfo);
  /** Get associated interaction information from monte carlo 
 * \param ev The gallery event.
 * \param mctruth The mctruth object associated with the currently considered interaction.
 *
 * \return NuMuInteraction object containing information from the mctruth object
 * */
  NuMuInteraction interactionInfo(const gallery::Event& ev, const simb::MCTruth &mctruth);
  /** Helper function -- whether point is contained in fiducial volume list 
 * \param v The point vector.
 *
 * \returns Whether the point is contained in the configured list of fiducial volumes.
 * */
  bool containedInFV(const TVector3 &v);

  unsigned _event_counter;  //!< Count processed events
  unsigned _nu_count;  //!< Count selected events

  static const unsigned nCuts = 4; //!< number of cuts
  /** Names of cuts 
 * \returns List of names of cuts (for histogram names)
 * */
  static const std::vector<std::string> cutNames() {
    return {"CC", "FV", "min_L", "reco_V"};
  }

  Config _config; //!< The config

  std::vector<NuMuInteraction> *_interactionInfo; //!< Branch holder

  RootHistos _root_histos[nCuts]; //!< Histos (one group per cut)
};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NumuSelection__

