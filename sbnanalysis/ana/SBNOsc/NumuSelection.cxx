#include <iostream>
#include <vector>
#include <TH2D.h>
#include <json/json.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "core/Event.hh"
#include "NumuSelection.h"
#include "Utilities.h"


namespace ana {
  namespace SBNOsc {

NumuSelection::NumuSelection() : SelectionBase(), fEventCounter(0), fNuCount(0) {}


void NumuSelection::Initialize(Json::Value* config) {
  // Load configuration parameters
  _config.truthTag = { "generator" };

  if (config) {
    _config.truthTag = { (*config)["SBNOsc"].get("MCTruthTag", "generator").asString() };

    // allow multiple fiducial volumes (accomodate for uboone data channels and icarus 2 TPC's)
    auto FVs = (*config)["NumuSelection"]["fiducial_volumes"];
    for (auto FV: FVs) {
      _config.aaBoxes.emplace_back(FV["xmin"].asDouble(), FV["ymin"].asDouble(), FV["zmin"].asDouble(), FV["xmax"].asDouble(), FV["ymax"].asDouble(), FV["zmax"].asDouble());
    }
    _config.doFVCut = (*config)["NumuSelection"].get("doFVcut", true).asBool();
  }

  hello();
}


void NumuSelection::Finalize() {}


bool NumuSelection::ProcessEvent(const gallery::Event& ev, std::vector<Event::Interaction>& reco) {
  if (fEventCounter % 10 == 0) {
    std::cout << "NumuSelection: Processing event " << fEventCounter << " "
              << "(" << fNuCount << " neutrinos selected)"
              << std::endl;
  }
  fEventCounter++;

  // Grab a data product from the event
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth> >(_config.truthTag);

  // Iterate through the neutrinos
  for (size_t i=0; i<mctruths.size(); i++) {
    Event::Interaction interaction;
    auto const& mctruth = mctruths.at(i);
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();

    if (Select(nu)) {
      Event::Interaction interaction = TruthReco(mctruth);
      
      reco.push_back(interaction);

    }
  }

  bool selected = !reco.empty();

  if (selected) {
    fNuCount++;
  }

  return selected;
}

bool NumuSelection::Select(const simb::MCNeutrino& nu) {
  // select true CC events
  bool pass_true_CC = (nu.CCNC() == simb::kCC && nu.Mode() == 0 && nu.Nu().PdgCode() == 14);

  // pass fiducial volume cut
  bool pass_FV = passFV(nu.Nu().Vx(), nu.Nu().Vy(), nu.Nu().Vz());

  return pass_true_CC && pass_FV;
}

bool NumuSelection::passFV(double x, double y, double z) {
  if (!_config.doFVCut) return true;
  for (auto const& aaBox: _config.aaBoxes) {
    geoalgo::Point_t p(x, y, z); 
    if (aaBox.Contain(p)) return true;
  }
  return false;

}

  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NumuSelection)

