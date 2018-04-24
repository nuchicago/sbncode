#include <iostream>
#include <vector>
#include <TH2D.h>
#include <json/json.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "TruthSelection.h"
#include "Selections.h"

namespace ana {
  namespace TruthSelection {

TruthSelection::TruthSelection()
    : SelectionBase(), fEventCounter(0), fSelectedCounter(0) {}


void TruthSelection::Initialize(Json::Value* config) {
  // Load configuration parameters
  fSelectionType = "ccnue_true";
  fTruthTag = { "generator" };
  fTrackTag = { "mcreco" };
  fShowerTag = { "mcreco" };

  if (config) {
    fSelectionType = (*config)["TruthSelection"].get("SelectionType", fSelectionType).asString();
    fTruthTag = { (*config)["TruthSelection"].get("MCTrithTag", "generator").asString() };
    fTrackTag = { (*config)["TruthSelection"].get("MCTrackTag", "mcreco").asString() };
    fShowerTag = { (*config)["TruthSelection"].get("MCShowerTag", "mcreco").asString() };
  }

  // Add custom branches
  AddBranch("reco_wgh", &fWeight);
  AddBranch("reco_e", &fRecoEnergy);
  AddBranch("reco_pdg", &fRecoPDG);
}


void TruthSelection::Finalize() {}


bool TruthSelection::ProcessEvent(gallery::Event& ev) {
  if (fEventCounter % 10 == 0) {
    std::cout << "TruthSelection: Processing event " << fEventCounter << " "
              << "(" << fSelectedCounter << " neutrinos selected)"
              << std::endl;
  }
  fEventCounter++;

  // Grab a data product from the event
  auto const& mctruths = *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);
  auto const& mctracks = *ev.getValidHandle<std::vector<sim::MCTrack> >(fTrackTag);
  auto const& mcshowers = *ev.getValidHandle<std::vector<sim::MCShower> >(fShowerTag);

  // Apply selection using tracks and showers
  bool pass = false;
  double fWeight = 1.0;

  if (fSelectionType == "ccnue_true") {
    pass = selections::CCNueTrue(mctruths, mctracks, mcshowers, fRecoEnergy, fWeight);
    fRecoPDG = 11;
  }
  else if (fSelectionType == "ccnum_true") {
    pass = selections::CCNumuTrue(mctruths, mctracks, mcshowers, fRecoEnergy, fWeight);
    fRecoPDG = 13;
  }
  else if (fSelectionType == "ccnue") {
    pass = selections::CCNue(mctruths, mctracks, mcshowers, fRecoEnergy, fWeight);
    fRecoPDG = 11;
  }
  else if (fSelectionType == "ccnumu") {
    pass = selections::CCNumu(mctruths, mctracks, mcshowers, fRecoEnergy, fWeight);
    fRecoPDG = 13;
  }
  else if (fSelectionType == "1e1p") {
    pass = selections::True1l1p0pi0(mctruths, mctracks, mcshowers, 11, fRecoEnergy, fWeight);
    fRecoPDG = 11;
  }
  else if (fSelectionType == "1m1p") {
    pass = selections::True1l1p0pi0(mctruths, mctracks, mcshowers, 13, fRecoEnergy, fWeight);
    fRecoPDG = 13;
  }
  else if (fSelectionType == "ccpi0") {
    pass = selections::CCPi0(mctruths, mctracks, mcshowers, fRecoEnergy, fWeight);
    fRecoPDG = 111;
  }
  else {
    std::cerr << "TruthSelection: Unknown selection type \""
              << fSelectionType << "\"" << std::endl;
    assert(false);
  }

  if (pass) {
    fSelectedCounter++;
    return true;
  }

  return false; 
}

  }  // namespace TruthSelection
}  // namespace ana


// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::TruthSelection::TruthSelection)

