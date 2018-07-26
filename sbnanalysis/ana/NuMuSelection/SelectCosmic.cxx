#include <iostream>
#include <vector>
#include <TH2D.h>
#include <json/json.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "SelectCosmic.h"

namespace ana {
  namespace NuMuSelection {

SelectCosmic::SelectCosmic() : SelectionBase() {}


void SelectCosmic::Initialize(Json::Value* config) {
}


void SelectCosmic::Finalize() {
}


bool SelectCosmic::ProcessEvent(gallery::Event& ev) {
  // Get MCTruth information
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);  
  std::cout << "MCTruth size: " << mctruths.size() << std::endl;
  for (auto const& truth: mctruths) {
    std::cout << truth << std::endl;
  }
  // get Pandora PFParticles
  return true;
}

  }  // namespace NuMuSelection
}  // namespace ana


// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::NuMuSelection::SelectCosmic)

