#include <TH1D.h>

#include "../sbnanalysis/core/Event.hh"
#include "../DataTypes/RecobInteraction.h"
#include "MakeHistos.h"
#include "../NuMuUtil.hh"

#include "TFile.h"
#include "TTree.h"

NuMuSelection::MakeHistos::MakeHistos(const char *input_fname, const char *output_fname): 
  fin(new TFile(input_fname)),
  fout(new TFile(output_fname, "RECREATE")),
  ev(0),
  reco_interactions(0)   
{
  // load libraries for branch classes
  //gROOT->ProcessLine(".L $SBN_LIB_DIR/libsbnanalysis_Event.so"); 
  //gROOT->ProcessLine(".L $SBN_LIB_DIR/libNuMuSelection_classes.so"); 
  
  // hook up branches
  TTree *event_tree = (TTree*)fin->Get("sbnana");
  // hook up TTree w/ branch class containers
  //event_tree->SetBranchAddress("reco_interactions", &reco_interactions);
  event_tree->SetBranchAddress("events", &ev);
  
  fout->cd();
  // setup histos 
  h_numu_ccqe = new TH1D("numu_ccqe", "numu_ccqe", 50, 0, 100);
  h_numu_trueE = new TH1D("numu_trueE", "numu_trueE", 100, 0 , 10);
}


void NuMuSelection::MakeHistos::Run() {
  // iterate over events
  long nEntries = event_tree->GetEntries();
  for (int i = 0; i < nEntries; i++) {
    event_tree->GetEntry(i);
    ProcessEvent(i);
  }
  // end
  FinishAnalysis();
}

void NuMuSelection::MakeHistos::ProcessEvent(int entry) {
  for (auto const &interaction: ev->interactions) {
    h_numu_ccqe->Fill(NuMuUtil::ECCQE(interaction.lepton.momentum, interaction.lepton.energy)); 
    h_numu_trueE->Fill(interaction.neutrino.energy);
  }
}

void NuMuSelection::MakeHistos::FinishAnalysis() {
  fout->Write();
  fout->Close();
}

