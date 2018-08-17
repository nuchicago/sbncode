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
#include <random>

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
    double visible_energy; //!< sum of kinetic energies of particles produced directly in interaction
    double smeared_visible_energy; //!< visible energy with component particle energies smeared
    double smeared_eccqe; //!< CCQE energy w/ lepton energy smeared
  };

protected:
  /** Configuration parameters */
  struct Config {
    art::InputTag truthTag; //!< art tag for MCTruth information
    art::InputTag mctTag;
    art::InputTag mcsTag;
    art::InputTag mcpTag;
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

  /** cpp distributions for smearing energies */
  struct Smearing {
    std::mt19937 _gen;
    std::normal_distribution<double> contained_lepton_E;
    std::normal_distribution<double> exiting_lepton_E;
    std::normal_distribution<double> hadron_E;

    Smearing(double contained_l_smear, double exiting_l_smear, double hadron_smear):
      _gen( time(0) ),
      contained_lepton_E(0., contained_l_smear),
      exiting_lepton_E(0., exiting_l_smear),
      hadron_E(0., hadron_smear) {}

    double Smear(double energy, int pdg, bool is_contained=false) {
      if (pdg == 13 /* muon */) {
        if (is_contained) return energy * contained_lepton_E(_gen);
        else return energy * exiting_lepton_E(_gen);
      }
      else {
        return hadron_E(_gen);
      }
    }
  };


  /** Returns whether to apply FV cut on neutrino */
  bool passFV(const TVector3 &v) { return containedInFV(v); }
  /** Applies reco-truth vertex matching cut */
  bool passRecoVertex(const TVector3 &truth_v, const TVector3 &reco_v);
  /** Applies truth length cut */
  bool passMinLength(double length, bool stop_in_tpc);
  /** Run Selection on a neutrino */
  std::vector<bool> Select(const gallery::Event& ev, const simb::MCTruth& mctruth, unsigned truth_ind, const NumuSelection::NuMuInteraction &intInfo);
  /** Get associated interaction information from monte carlo */
  NuMuInteraction interactionInfo(const gallery::Event& ev, const simb::MCTruth &mctruth);
  /** Helper function -- whether point is contained in fiducial volume list */
  bool containedInFV(const TVector3 &v);

  unsigned fEventCounter;  //!< Count processed events
  unsigned fNuCount;  //!< Count selected events

  static const unsigned nCuts = 4; //!< number of cuts
  /** Names of cuts */
  static const std::vector<std::string> cutNames() {
    return {"CC", "FV", "min_L", "reco_V"};
  }

  Config _config; //!< The config

  std::vector<NuMuInteraction> *_interactionInfo; //!< Branch holder

  RootHistos _root_histos[nCuts]; //!< Histos (one group per cut)

  Smearing *_smear; //!< Instance for smearing FSP energies
};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NumuSelection__

