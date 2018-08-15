#include <iostream>
#include <vector>
#include <TH2D.h>
#include <json/json.h>
#include <cmath>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include "gallery/Event.h"
#include <TLorentzVector.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <random>
#include "core/Event.hh"
#include "NueSelection.h"
#include "Utilities.h"

namespace ana {
  namespace SBNOsc {

NueSelection::NueSelection() : SelectionBase(), fEventCounter(0), fNuCount(0),fMatchedNuCount(0){}


void NueSelection::Initialize(Json::Value* config) {

  fDiffLength = new TH1D ("diff_length","",60,0,6);

  // Load configuration parameters
  fTruthTag = { "generator" };
  fTrackTag = { "mcreco" };
  fShowerTag = { "mcreco" };

  if (config) {
    fTruthTag = { (*config)["SBNOsc"].get("MCTruthTag", "generator").asString() };
    fTrackTag = { (*config)["SBNOsc"].get("MCTrackTag", "mcreco").asString() };
    fShowerTag = { (*config)["SBNOsc"].get("MCShowerTag","mcreco").asString() };
  }

  AddBranch("nucount",&fNuCount);
  AddBranch("matched_nucount",&fMatchedNuCount);

  hello();
}


void NueSelection::Finalize() {}
  fOutputFile->cd();

bool NueSelection::ProcessEvent(const gallery::Event& ev, std::vector<Event::Interaction>& reco) {
  if (fEventCounter % 10 == 0) {
    std::cout << "NueSelection: Processing event " << fEventCounter << " "
              << "(" << fNuCount << " neutrinos selected)"
              << std::endl;
  }
  fEventCounter++;

  // Grab data products from the event
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);
  auto const& mctracks = \
    *ev.getValidHandle<std::vector<sim::MCTrack> >(fTrackTag);
  auto const& mcshowers = \
    *ev.getValidHandle<std::vector<sim::MCShower> >(fShowerTag);

  // Iterate through the showers
  for (size_t k=0;k<mcshowers.size();k++) {
    auto const& mcshower = mcshowers.at(k);
    auto shower_energy = mcshower.DetProfile().E();
    auto shower_pos = mcshower.DetProfile().Position();
    if (shower_energy >= 0.2) {
      std::vector<const simb::MCNeutrino&> nearby_nus;
      for (size_t i=0; i<mctruths.size();i++) {
        auto const& mctruth = mctruths.at(i);
        const simb::MCNeutrino& nu = mctruth.GetNeutrino();
        const simb::MCParticle& nu_particle = nu.Nu();
        auto nu_pos = nu_particle.Position();
        auto distance = (nu_pos.Vect() - shower_pos.Vect()).Mag();
        fDiffLength->Fill(distance);
        if (distance <= 5.) nearby_nus.push_back(nu);
      }
    }
  }

  // Iterate through the neutrinos/MCTruth
  for (size_t i=0; i<mctruths.size(); i++) {
    Event::Interaction interaction;
    auto const& mctruth = mctruths.at(i);
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();

    if (nu.CCNC() == simb::kCC && nu.Mode() == 0 && nu.Nu().PdgCode() == 12) {
      Event::Interaction interaction = TruthReco(mctruth);
      reco.push_back(interaction);
    }
  }

  bool selected = !reco.empty(); // true if reco info is not empty

  if (selected) {
    fNuCount++;
  }

  return selected;
}

  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NueSelection)
