#ifndef __NuMuSelection_OscChi2_hh
#define __NuMuSelection_OscChi2_hh

class TTree;
class TFile;
class TH1D;

namespace NuMuSelection {
class OscChi2 {
  OscChi2(double osc_dm, double osc_angle, double osc_distance, const char *input_fname="histos.root", const char *output_fname="scaled_histos.root");
  void Run();

  struct Config {
    double osc_dm;
    double osc_angle;
    double osc_distance;
  };

private:
  // input
  TFile *signal_file;
  TH1D *signal_histo;
  // output
  TFile *output_file;
  TH1D *osc_histo;
  Config _config;

};
}

#endif
