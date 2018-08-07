#ifndef __NuMuSelection_MakeHistos_hh
#define __NuMuSelection_MakeHistos_hh
#include "../sbnanalysis/core/Event.hh"
#include "../DataTypes/RecobInteraction.h"

#include "TFile.h"
#include "TTree.h"

class TH1D;

namespace NuMuSelection {
class MakeHistos {
public:
  explicit MakeHistos(const char *input_fname="signal.root", const char *output_fname="histos.root");
  void Run();

private:
  // function per event
  void ProcessEvent(int entry);
  // function at end of analysis
  void FinishAnalysis();

  // input file
  TFile *fin;
  // input tree
  TTree *event_tree;
  // output stuff
  TFile *fout;
  TH1D *h_numu_ccqe;
  TH1D *h_numu_trueE;
  // class containers
  Event *ev;
  std::vector<RecobInteraction> *reco_interactions;
};
}
#endif
