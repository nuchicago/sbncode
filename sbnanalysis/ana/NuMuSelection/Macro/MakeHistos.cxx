#include "../sbnanalysis/core/Event.hh"
#include "../RecobInteraction.h"
#include "../NuMuUtil.hh"

// class for holding information
class Histos {
public:
  explicit Histos(const char *output_fname):
    fout(new TFile(output_fname, "RECREATE")),
    ev(0),
    reco_interactions(0)   
  {
    fout->cd();
    h_numu_ccqe = new TH1D("numu_ccqe", "numu_ccqe", 50, 0, 100);
    h_numu_trueE = new TH1D("numu_trueE", "numu_trueE", 100, 0 , 10);
  }

  TFile *fout;
  TH1D *h_numu_ccqe;
  TH1D *h_numu_trueE;
  // class containers
  Event *ev;
  std::vector<RecobInteraction> *reco_interactions;
};

// function per event
void ProcessEvent(int entry, Histos& histos);
// function at beginning of analysis
Histos BeginAnalysis(const char *output_fname);
// function at end of analysis
void FinishAnalysis(Histos& histos);

void MakeHistos(const char *input_fname="output.root", const char *output_fname="histos.root") {
  // load libraries for branch classes
  gROOT->ProcessLine(".L $SBN_LIB_DIR/libsbnanalysis_Event.so"); 
  gROOT->ProcessLine(".L $SBN_LIB_DIR/libNuMuSelection_classes.so"); 

  // open file and load TTree
  TFile signal_file(input_fname);
  TTree *event_tree = (TTree*)signal_file.Get("sbnana");

  // start
  auto histos = BeginAnalysis(output_fname);

  // hook up TTree w/ branch class containers
  //event_tree->SetBranchAddress("reco_interactions", &histos.reco_interactions);
  event_tree->SetBranchAddress("events", &histos.ev);

  // iterate over events
  long nEntries = event_tree->GetEntries();
  for (int i = 0; i < nEntries; i++) {
    event_tree->GetEntry(i);
    ProcessEvent(i, histos);
  }

  // end
  FinishAnalysis(histos);
}

void ProcessEvent(int entry, Histos& histos) {
  for (auto const &interaction: histos.ev->interactions) {
    histos.h_numu_ccqe->Fill(NuMuUtil::ECCQE(interaction.lepton.momentum, interaction.lepton.energy)); 
    histos.h_numu_trueE->Fill(interaction.neutrino.energy);
  }
}

Histos BeginAnalysis(const char *output_fname) {
  return Histos(output_fname);
}

void FinishAnalysis(Histos& histos) {
  histos.fout->Write();
  histos.fout->Close();
}
