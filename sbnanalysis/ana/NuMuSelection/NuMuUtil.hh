#ifndef __NuMuUTIL_hh_
#define __NuMuUTIL_hh_

class TH1D;
class TVector3;

namespace NuMuUtil {

double ECCQE(const TVector3& inp_v, double inp_energy, double energy_distortion=0., double angle_distortion=0.);

// return probability of muon neutrino being around at a given energy
// 
// if given a sample of non-oscilalted numu's, we can add in the osccilation
// by adding this in as an event weight
double OscillationWeight(double numu_energy, double osc_dm, double osc_angle, double osc_dist);

// calculates chi2 given a oscillated/non-oscillated histograms
// uses only statistical errors
double Chi2Stat(TH1D *baseline_histo, TH1D *osc_histo);
} // namespace NuMuUtil
#endif 
