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
//#include "core/Event.hh"
#include "NueSelection.h"
#include "Utilities.h"

namespace ana {
  namespace SBNOsc {

NueSelection::NueSelection() : SelectionBase(), fEventCounter(0), fNuCount(0){}


void NueSelection::Initialize(Json::Value* config) {

  fDiffLength = new TH1D ("diff_length","",200,0,200);

  fTrackLength = new TH1D ("diff_length","",200,0,200);
  fGenNueHist = new TH1D ("generated_nue_hist","",60,0,6);
  fGenNueFidVolHist = new TH1D ("generated_nue_in_fiducial_volume","",60,0,6);
  fSelectedNuHist = new TH1D ("selected_nu_hist","",60,0,6);

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


  hello();
}


void NueSelection::Finalize() {
  fOutputFile->cd();
  fDiffLength->Write();
  fTrackLength->Write();
  fGenNueHist->Write();
  fGenNueFidVolHist->Write();
  fSelectedNuHist->Write();
}



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

  std::vector<bool> matchedness;
  // Iterate through the neutrinos
  for (size_t i=0;i<mctruths.size();i++) {
    auto const& mctruth = mctruths.at(i);
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();
    auto nu_E = nu.Nu().E();
    if (nu.Nu().PdgCode() == 12) fGenNueHist->Fill(nu_E);
    auto vx = nu.Nu().Vx();
    auto vy = nu.Nu().Vy();
    auto vz = nu.Nu().Vz();
    if ((nu.Nu().PdgCode() ==12)&&(((-174.15 < vx && vx < -27.65) || (27.65 < vx && vx < 174.15)) && (-175 < vy && vy < 175) && (25 < vz && vz < 475))) fGenNueFidVolHist->Fill(nu_E);

    auto nu_pos = nu.Nu().Position();
    int matched_shower_count = 0;
    for (size_t j=0;j<mcshowers.size();j++) {
      auto const& mcshower = mcshowers.at(j);
      auto shower_pos = mcshower.DetProfile().Position();
      double distance = (nu_pos.Vect()-shower_pos.Vect()).Mag();
      fDiffLength->Fill(distance);
      if (distance <= 5.) {
        matched_shower_count++;
      }
    }
    if (matched_shower_count>0) matchedness.push_back(true);
    else matchedness.push_back(false);
  }

  // Iterate through the tracks, compute the active length
  std::vector<bool> NotMuTrackness;
  int MuTrackCount =0;
  for (size_t i=0;i<mctracks.size();i++) {
    auto const& mctrack = mctracks.at(i);
    if (IsUsableTrack(mctrack)) {
      double this_active_length;
      this_active_length = GetActiveLength(mctrack);
      fTrackLength->Fill(this_active_length);
      if (this_active_length>=100.) {
        MuTrackCount++;
      }
    }
  }
  if (MuTrackCount>0) {
    return false; //questionable
  }


  // Iterate through the neutrinos/MCTruth
  for (size_t i=0; i<mctruths.size(); i++) {
    Event::Interaction interaction;
    auto const& mctruth = mctruths.at(i);
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();
    auto nu_E = nu.Nu().E();
    auto vx = nu.Nu().Vx();
    auto vy = nu.Nu().Vy();
    auto vz = nu.Nu().Vz();
    bool IsFid = (((-174.15 < vx && vx < -27.65) || (27.65 < vx && vx < 174.15)) && (-175 < vy && vy < 175) && (25 < vz && vz < 475));
    if (matchedness[i]&&IsFid) {
      fSelectedNuHist->Fill(nu_E);
      Event::Interaction interaction = TruthReco(mctruth);
      reco.push_back(interaction);
    }
    /*
    if (nu.CCNC() == simb::kCC && nu.Mode() == 0 && nu.Nu().PdgCode() == 12) {
      Event::Interaction interaction = TruthReco(mctruth);
      reco.push_back(interaction);
    }
    */
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
