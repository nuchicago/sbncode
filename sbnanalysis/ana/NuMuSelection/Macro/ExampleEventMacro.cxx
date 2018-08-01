#include "sbnanalysis/core/Event.hh"

// function per event
void ProcessEvent(int event_no);
// function at beginning of analysis
void BeginAnalysis();
// function at end of analysis
void FinishAnalysis();

void ExampleEventMacro(const char *file_name="output.root") {
  // load libraries for branch classes
  gROOT->ProcessLine(".L $SBN_LIB_DIR/libsbnanalysis_Event.so"); 

  // define class containers
  Event *ev = 0;

  // open file and load TTree
  TFile signal_file(file_name);
  TTree *event_tree = (TTree*)signal_file.Get("sbnana");

  // hook up TTree w/ branch class containers
  event_tree->SetBranchAddress("events", &ev);

  // start
  BeginAnalysis();

  // iterate over events
  long nEntries = event_tree->GetEntries();
  for (int i = 0; i < nEntries; i++) {
    event_tree->GetEntry(i);
    ProcessEvent(i);
  }

  // end
  FinishAnalysis();
}

void ProcessEvent(int event_no) {
  std::cout << "Event " << event_no << "!\n";
}

void BeginAnalysis() {
  std::cout << "Begin!\n";
}

void FinishAnalysis() {
  std::cout << "End!\n";
}

