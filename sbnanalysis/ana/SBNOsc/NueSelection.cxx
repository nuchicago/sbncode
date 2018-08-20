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

  fTrackLength = new TH1D ("track_length","",200,0,200);
  fGenNueHist = new TH1D ("generated_nue_hist","",60,0,6);
  fGenHist = new TH1D ("generated_particles","",60,0,6);
  fGenNumuHist = new TH1D("generated_numu","",60,0,6);
  fGenBarNueHist = new TH1D ("gen_bar_nue","",60,0,6);
  fGenOtherHist = new TH1D ("gen_other","",60,0,6);
  fGenNueFidVolHist = new TH1D ("generated_nue_in_fiducial_volume","",60,0,6);
  fSelectedNuHist = new TH1D ("selected_nu_hist","",60,0,6);
  fNodEdxNuHist = new TH1D ("no_dEdx","",60,0,6);

  fShowerEnergy = new TH1D ("shower_energy","",1000,0,1000);
  fEnergeticShowerHist = new TH1D("energetic_shower_energy","",1000,0,100);
  fVisibleVertexNuEHist = new TH1D("visible_vertex_nu_energy","",60,0,6);
  fCGSelectionHist = new TH1D("final_selected_nu","",60,0,6);
  fRecoSelectionHist = new TH1D("final_selected_nu_w_reco_efficiency","",60,0,6);
  fShowerCutSelectionHist = new TH1D("shower_cut_fid_only_nu_energy","",60,0,6);
  fSelectedTrueNue = new TH1D ("true_nue_selected","",60,0,6);

  //shower true type
  fMuShowerSelectedNu = new TH1D("mu_shower","",60,0,6);
  fEShowerSelectedNu = new TH1D ("e_shower","",60,0,6);
  fGammaShowerSelectedNu = new TH1D ("gamma_shower","",60,0,6);
  fOtherShowerSelectedNu = new TH1D("other_shower","",60,0,6);

  //shower dE/dx
  fShowerdEdx = new TH1D("shower_dEdx","",60,0,6);
  fEShowerdEdx = new TH1D("electron_shower_dEdx","",60,0,6);
  fGammaShowerdEdx = new TH1D("gamma_shower_dEdx","",60,0,6);
  fMuShowerdEdx = new TH1D("mu_shower_dEdx","",60,0,6);
  fPositronShowerdEdx = new TH1D("positron_shower_dEdx","",60,0,6);
  fOtherShowerdEdx = new TH1D("other_shower_dEdx","",60,0,6);




  // Load configuration parameters
  fEnergyThreshold =0.;
  fNuCount=0;
  fTruthTag = { "generator" };
  fTrackTag = { "mcreco" };
  fShowerTag = { "mcreco" };

  if (config) {
    fEnergyThreshold = (*config)["SBNOsc"].get("energy_threshold",123.0).asDouble();
    //fEnergyThreshold = { (*config)["SBNOsc"].get("energy_threshold",12.34).asDouble() };
    fTruthTag = { (*config)["SBNOsc"].get("MCTruthTag", "generator").asString() };
    fTrackTag = { (*config)["SBNOsc"].get("MCTrackTag", "mcreco").asString() };
    fShowerTag = { (*config)["SBNOsc"].get("MCShowerTag","mcreco").asString() };
  }
  AddBranch("energy_threshold",&fEnergyThreshold);
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
  fShowerEnergy->Write();
  fEnergeticShowerHist->Write();
  fVisibleVertexNuEHist->Write();
  fCGSelectionHist->Write();
  fRecoSelectionHist->Write();
  fShowerCutSelectionHist->Write();
  fSelectedTrueNue->Write();
  fNodEdxNuHist->Write();

  fMuShowerSelectedNu->Write();
  fEShowerSelectedNu->Write();
  fGammaShowerSelectedNu->Write();
  fOtherShowerSelectedNu->Write();

  fShowerdEdx->Write();
  fEShowerdEdx->Write();
  fGammaShowerdEdx->Write();
  fOtherShowerdEdx->Write();
  fMuShowerdEdx->Write();
  fPositronShowerdEdx->Write();

  fGenHist->Write();
  fGenBarNueHist->Write();
  fGenNumuHist->Write();
  fGenOtherHist->Write();
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

 //Conversion gap cut: mark visibility of vertices
 std::vector<bool> IsVisibleVertex;
  for (size_t i=0;i<mctruths.size();i++) {
    auto const& mctruth = mctruths.at(i);
    auto const& nu = mctruth.GetNeutrino();
    // condition 1: nu vertex within active volume
    double vx = nu.Nu().Vx();
    double vy = nu.Nu().Vy();
    double vz = nu.Nu().Vz();
    bool InActiveVol = (((-199.15 < vx && vx < -2.65) || (2.65 < vx && vx < 199.15)) && (-200 < vy && vy < 200) && (0 < vz && vz < 500));

    // condition 2: charged hadron activity > 50 MeV

    auto nu_pos = nu.Nu().Position();
    double ChargedHadronEnergy=0.;
    for (size_t j=0; j<mctracks.size();j++) {
      auto const& mctrack = mctracks.at(j);
      // 2a. loop through all tracks, and find ones that are within 5 cm to the vertex
      auto track_start_pos = mctrack.Start().Position();
      double dis = (nu_pos.Vect() - track_start_pos.Vect()).Mag(); //compute the distance between track starting position and vertex position
      bool AssnVertex = (dis<=5.); //vertex association marker

      //2b. is proton track
      bool IsChargedHadron = (mctrack.PdgCode() ==2212); //truth-level proton pdg code

      double track_energy = mctrack.Start().E();

      if (AssnVertex&&IsChargedHadron) ChargedHadronEnergy+=track_energy;
    }
    bool Visible = (ChargedHadronEnergy>=0.05);
    IsVisibleVertex.push_back(Visible);
    if (Visible) fVisibleVertexNuEHist->Fill(nu.Nu().E());
  }
  assert (IsVisibleVertex.size() == mctruths.size());

  //Conversion gap cut: look at only visible vertices and apply the cut
  std::vector<bool> PassConversionGap; //bool vector storing passing info
  for (size_t i=0;i<mctruths.size();i++) {
    auto const& mctruth = mctruths.at(i);
    auto const& nu = mctruth.GetNeutrino();
    auto nu_pos = nu.Nu().Position();
    if (IsVisibleVertex[i]) { //if visible, check number of assn showers
      int AssnShowerCount=0;
      for (size_t j=0;j<mcshowers.size();j++) {
        auto const& mcshower = mcshowers.at(j);
        auto shower_pos = mcshower.DetProfile().Position();
        double dis = (nu_pos.Vect() - shower_pos.Vect()).Mag();
        if (dis<=3) AssnShowerCount++;
      }
      PassConversionGap.push_back(AssnShowerCount>0); //pass if the number of assn showers is nonzero
    }
    else PassConversionGap.push_back(true); //if not visible, automatically pass
  }
  assert (PassConversionGap.size()==mctruths.size()); //sanity check

  // shower energy cut
  std::vector<int> EnergeticShowersIndices;
  // mark dEdx
  std::vector<bool> HasGooddEdx;
  for (size_t i=0;i<mcshowers.size();i++) {
    auto const& mcshower = mcshowers.at(i);
    double shower_E = mcshower.DetProfile().E();
    //fill in the dEdx hist
    double shower_dEdx = mcshower.dEdx();
    bool dEdxQuality = (shower_dEdx<=1.5);
    HasGooddEdx.push_back(dEdxQuality);
    int showerPDG = mcshower.PdgCode();
    fShowerdEdx->Fill(shower_dEdx);
    if (showerPDG == 11) fEShowerdEdx->Fill(shower_dEdx);
    if (showerPDG == 22) fGammaShowerdEdx->Fill(shower_dEdx);
    if (showerPDG == 13) fMuShowerdEdx->Fill(shower_dEdx);
    if (showerPDG == -11) fPositronShowerdEdx->Fill(shower_dEdx);
    if ((showerPDG!=11)&&(showerPDG!=22)&&(showerPDG!=-11)&&(showerPDG!=13)) fOtherShowerdEdx->Fill(shower_dEdx);
    fShowerEnergy->Fill(shower_E);
    if (shower_E >= fEnergyThreshold) {
      EnergeticShowersIndices.push_back(i); //have yet to implement the configurable energy threshold parameter
    // if (Shower_E > fEnergyThreshold) EnergeticShowers.push_back(mcshower);
      fEnergeticShowerHist->Fill(shower_E);
      }
  }
  assert (HasGooddEdx.size()==mcshowers.size());



  //matching
  std::vector<bool> matchedness;

  std::vector<int> ShowerPDG; //shower true pdg code

  // Iterate through the neutrinos
  std::vector<bool> AssnShowerGooddEdx;
  for (size_t i=0;i<mctruths.size();i++) {
    auto const& mctruth = mctruths.at(i);
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();
    auto nu_E = nu.Nu().E();
    fGenHist->Fill(nu_E);
    if (nu.Nu().PdgCode() == 12) fGenNueHist->Fill(nu_E);
    if (nu.Nu().PdgCode() == 14) fGenNumuHist->Fill(nu_E);
    if (nu.Nu().PdgCode() == -12) fGenBarNueHist->Fill(nu_E);
    if ((nu.Nu().PdgCode() != 12)&&(nu.Nu().PdgCode() != 14)&&(nu.Nu().PdgCode() != -12)) fGenOtherHist->Fill(nu_E);
    auto vx = nu.Nu().Vx();
    auto vy = nu.Nu().Vy();
    auto vz = nu.Nu().Vz();
    if ((nu.Nu().PdgCode() ==12)&&(((-174.15 < vx && vx < -27.65) || (27.65 < vx && vx < 174.15)) && (-175 < vy && vy < 175) && (25 < vz && vz < 475))) fGenNueFidVolHist->Fill(nu_E);
    auto nu_pos = nu.Nu().Position();
    int matched_shower_count = 0;
    std::vector<int> assn_showers_pdg;
    std::vector<bool> assn_showers_quality; //dEdx
    // loop through only energetic showers
    for (auto j : EnergeticShowersIndices) {
      auto const& shower = mcshowers.at(j);
      auto shower_pos = shower.DetProfile().Position();
      double distance = (nu_pos.Vect()-shower_pos.Vect()).Mag();
      fDiffLength->Fill(distance);
      if (distance <= 5.) {
        matched_shower_count++;
        assn_showers_pdg.push_back(shower.PdgCode());
        assn_showers_quality.push_back(HasGooddEdx[j]);
      }
    }
    if (!assn_showers_pdg.empty()) ShowerPDG.push_back(assn_showers_pdg[0]);
    else ShowerPDG.push_back(0);
    if (!assn_showers_quality.empty()) AssnShowerGooddEdx.push_back(assn_showers_quality[0]);
    else AssnShowerGooddEdx.push_back(false);
    if (matched_shower_count>0) matchedness.push_back(true);
    else matchedness.push_back(false);
  }
  assert(ShowerPDG.size()==matchedness.size());
  assert(AssnShowerGooddEdx.size()==matchedness.size());

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
    int MuTrackCount=0;
    auto nu_pos = nu.Nu().Position();
    for (size_t j=0;j<mctracks.size();j++) {
      auto const& mctrack = mctracks.at(j);
      auto track_start_pos = mctrack.Start().Position();
      auto nu_track_dis = (track_start_pos.Vect() - nu_pos.Vect()).Mag();
      //do total length cut first
      auto total_length = (mctrack.Start().Position().Vect() - mctrack.End().Position().Vect()).Mag();
      if ((total_length >= 100.) && (nu_track_dis<=5.)) MuTrackCount++;
    }
    bool NotMuTrackness = (MuTrackCount==0);

    if (matchedness[i]&&IsFid) fShowerCutSelectionHist->Fill(nu_E);
    if (matchedness[i]&&IsFid&&NotMuTrackness) fSelectedNuHist->Fill(nu_E);
    if (matchedness[i]&&IsFid&&NotMuTrackness&&PassConversionGap[i]) fNodEdxNuHist->Fill(nu_E);
    if (matchedness[i]&&IsFid&&NotMuTrackness&&PassConversionGap[i]&&AssnShowerGooddEdx[i]) {
      fCGSelectionHist->Fill(nu_E);
      if (nu.Nu().PdgCode() ==12) fSelectedTrueNue->Fill(nu_E);

      //fill in shower pdg histograms
      if (ShowerPDG[i]==13) fMuShowerSelectedNu->Fill(nu_E);
      if (ShowerPDG[i]==11) fEShowerSelectedNu->Fill(nu_E);
      if (ShowerPDG[i]==22) fGammaShowerSelectedNu->Fill(nu_E);
      if ((ShowerPDG[i]!=13)&&(ShowerPDG[i]!=11)&&(ShowerPDG[i]!=22))
        fOtherShowerSelectedNu->Fill(nu_E);

      double lower_bound=0.;
      double upper_bound=1.;
      std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
      std::random_device r;
      std::default_random_engine e1(r());
      double rnd_double = unif(e1);
      if (rnd_double <= 0.8) {
        fRecoSelectionHist->Fill(nu_E);
        Event::Interaction interaction = TruthReco(mctruth);
        reco.push_back(interaction);
      }
      else continue;
    }

    /*
    if (nu.CCNC() == simb::kCC && nu.Mode() == 0 && nu.Nu().PdgCode() == 12) {
      Event::Interaction interaction = TruthReco(mctruth);
      reco.push_back(interaction);
    }
    */
  }
    /*
    if (nu.CCNC() == simb::kCC && nu.Mode() == 0 && nu.Nu().PdgCode() == 12) {
      Event::Interaction interaction = TruthReco(mctruth);
      reco.push_back(interaction);
    }
    */


  bool selected = !reco.empty(); // true if reco info is not empty

  if (selected) {
    fNuCount++;
  }

  return selected;
}

  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NueSelection)
