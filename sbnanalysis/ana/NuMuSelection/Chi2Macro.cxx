#include "sbnanalysis/core/Event.hh"
#include "RecobInteraction.h"

// function per event
void ProcessEvent();
// function at beginning of analysis
void BeginAnalysis();
// function at end of analysis
void FinishAnalysis();

void Chi2Macro(const char *file_name="output.root") {
  // load libraries for branch classes
  gROOT->ProcessLine(".L $SBN_LIB_DIR/libsbnanalysis_Event.so"); 
  gROOT->ProcessLine(".L $SBN_LIB_DIR/libNuMuSelection_classes.so"); 

  // define class containers
  Event *ev = 0;
  std::vector<RecobInteraction> *reco_interactions = 0;

  // open file and load TTree
  TFile signal_file(file_name);
  TTree *event_tree = (TTree*)signal_file.Get("sbnana");

  // hook up TTree w/ branch class containers
  event_tree->SetBranchAddress("reco_interactions", &reco_interactions);
  event_tree->SetBranchAddress("events", &ev);

  // start
  BeginAnalysis();

  // iterate over events
  long nEntries = event_tree->GetEntries();
  for (int i = 0; i < nEntries; i++) {
    event_tree->GetEntry(i);
    ProcessEvent();
  }

  // end
  FinishAnalysis();
}

void ProcessEvent() {
  std::cout << "Event!\n";
}

void BeginAnalysis() {
  std::cout << "Begin!\n";
}

void EndAnalysis() {
  std::cout << "End!\n";
}

