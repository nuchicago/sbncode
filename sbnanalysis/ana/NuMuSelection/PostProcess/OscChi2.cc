#include <TH1D.h>

#include "../sbnanalysis/core/Event.hh"
#include "../DataTypes/RecobInteraction.h"
#include "OscChi2.h"
#include "../NuMuUtil.hh"

#include "TFile.h"
#include "TTree.h"

NuMuSelection::OscChi2::OscChi2(double osc_dm, double osc_angle, double osc_distance, const char *input_fname, const char *output_fname):
  signal_file(new TFile(input_fname)),
  signal_histo((TH1D*)signal_file->Get("numu_trueE")),
  output_file(new TFile(output_fname, "RECREATE")),
  _config({osc_dm, osc_angle, osc_distance})
{
  output_file->cd();
  // make your own histogram
  osc_histo = new TH1D(*signal_histo);
}

void NuMuSelection::OscChi2::Run() {
  int n_bins = osc_histo->GetNbinsX();
  for (int i = 0; i < n_bins; i++) {
    double energy = osc_histo->GetBinCenter(i);
    // scale by oscillation probability
    osc_histo->Scale(NuMuUtil::OscillationWeight(energy, _config.osc_dm, _config.osc_angle, _config.osc_distance));
  }
  // print the chi2
  double chi2 = NuMuUtil::Chi2Stat(signal_histo, osc_histo);
  std::cout << "CHI2: " << chi2 << std::endl;

  // write the histo and close the file
  output_file->Write();
  output_file->Close();
}

