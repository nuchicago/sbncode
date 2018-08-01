#include "sbnanalysis/core/Event.hh"
#include "../RecobInteraction.h"
#include "../NuMuUtil.hh"

void OscChi2(double osc_dm, double osc_angle, double osc_distance, const char *input_fname="histos.root", const char *output_fname="scaled_histos.root") {
  // open file and load TTree
  TFile signal_file(input_fname);
  TH1D *signal_histo = (TH1D*)signal_file.Get("numu_trueE");

  TFile output_file(output_fname, "RECREATE");
  output_file.cd();
  // make your own histogram
  TH1D osc_histo(*signal_histo);
  int n_bins = osc_histo.GetNbinsX();
  for (int i = 0; i < n_bins; i++) {
    double energy = osc_histo.GetBinCenter(i);
    // scale by oscillation probability
    osc_histo.Scale(NuMuUtil::OscillationWeight(energy, osc_dm, osc_angle, osc_distance));
  }
  // print the chi2
  double chi2 = NuMuUtil::Chi2Stat(signal_histo, &osc_histo);
  std::cout << "CHI2: " << chi2 << std::endl;

  // write the histo and close the file
  output_file.Write();
  output_file.Close();
}

