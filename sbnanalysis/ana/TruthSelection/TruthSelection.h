#ifndef __sbnanalysis_ana_TruthSelection_TruthSelection__
#define __sbnanalysis_ana_TruthSelection_TruthSelection__

/**
 * \file TruthSelection.h
 *
 * A truth-based selection based on the SBN proposal, arxiv:1503.01520,
 * c. page I-20.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>
 */

#include <string>
#include <vector>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"

class TH2D;

namespace ana {
  namespace TruthSelection {

/**
 * \class TruthSelection
 * \brief Truth-based selection per the SBN proposal
 */
class TruthSelection : public core::SelectionBase {
public:
  /** Constructor. */
  TruthSelection();

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
   * \return True to keep event
   */
  bool ProcessEvent(gallery::Event& ev);

protected:
  unsigned fEventCounter;  //!< Count processed events
  unsigned fSelectedCounter;  //!< Count selected events

  /// Configuration parameters
  art::InputTag fTruthTag;  //!< art tag for MCTruth information
  art::InputTag fTrackTag;  //!< art tag for MCTrack information
  art::InputTag fShowerTag;  //!< art tag for MCShower information
  std::string fSelectionType;  //!< Selection type, from configuration parameter

  /// Custom data branches
  double fWeight;  //!< Efficiency weight
  double fRecoEnergy;  //!< Reconstructed neutrino energy
  int fRecoPDG;  //!< Selection PDG
};

  }  // namespace TruthSelection
}  // namespace ana

#endif  // __sbnanalysis_ana_TruthSelection_TruthSelection__

