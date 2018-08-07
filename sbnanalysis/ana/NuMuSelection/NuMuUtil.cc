#include <TH1D.h>
#include <TVector3.h>

#include "NuMuUtil.hh"

namespace NuMuUtil {
double ECCQE(const TVector3& inp_v, double inp_energy, double energy_distortion, double angle_distortion) {
  // Based on D. Kaleko, LowEnergyExcess LArLite module ECCQECalculator
  double M_n = 0.939565; // GeV/c^2
  double M_p = 0.938272; // GeV/c^2
  double M_e = 0.000511; // GeV/c^2
  double bindingE = 0.0300; // GeV
  
  double mp2 = M_p * M_p;
  double me2 = M_e * M_e;
  double mnb = M_n - bindingE;

  // mess around with lorentz vector
  TVector3 v(inp_v);
  v.SetTheta( v.Theta() + angle_distortion);
  double l_energy = inp_energy + energy_distortion;
  double l_mom = sqrt(l_energy*l_energy - me2);
  double l_theta = \
    acos(v.Pz() / sqrt(v.Px()*v.Px() + v.Py()*v.Py() + v.Pz()*v.Pz()));
  double enu_top = mp2 - mnb*mnb - me2 + 2.0 * mnb * l_energy;
  double enu_bot = 2.0 * (mnb - l_energy + l_mom * cos(l_theta));
  return enu_top / enu_bot;
}

double OscillationWeight(double numu_energy, double osc_dm, double osc_angle, double osc_dist) {
  // if any parameters weren't configured, just return 1
  if (osc_dm < 0 || osc_dist < 0 || osc_angle < 0) return 1;  

  double overlap = sin(2*osc_angle);
  double energy_factor = sin(1.27 * osc_dm * osc_dist / numu_energy);

  return 1 - overlap * overlap * energy_factor * energy_factor;
}

double Chi2Stat(TH1D *baseline_histo, TH1D *osc_histo) {
  double chi2 = 0;
  int n_bins = osc_histo->GetNbinsX();
  for (int i = 0; i < n_bins; i++) {
    double baseline = baseline_histo->GetBinContent(i);
    double osc = osc_histo->GetBinContent(i);
    if (baseline < 1e-4) continue;
     
    chi2 += (baseline - osc)*(baseline - osc) / baseline;   
  }
  return chi2;

}

} // namespace NuMuUtil
