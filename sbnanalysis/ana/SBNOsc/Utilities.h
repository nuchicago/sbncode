#ifndef __sbnanalysis_ana_SBNOsc_Utilities__
#define __sbnanalysis_ana_SBNOsc_Utilities__

/**
 * \file Utilities.h
 *
 * Common utilties
 *
 * This is some auxiliary code that is not a selection, but does a piece
 * of the analysis. We can define any number of other functions, classes,
 * etc. which we use in the selection.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>
 */

#include "nusimdata/SimulationBase/MCTruth.h"
#include "core/Event.hh"

namespace ana {
  namespace SBNOsc {

/** A function that says hello. */
void hello();


/** Extract truth information to approximate reconstruction. */
Event::Interaction TruthReco(const simb::MCTruth& mctruth);

/** Calculate CCQE energy from associated lepton information (and optional distortion). Energy in GeV. 
 *
 * \param l_momentum Lepton momentum (in any units -- used only to get angle info)
 * \param l_energy Lepton energy in GeV
 * \param energy_distortion Optional energy distortion in GeV
 * \param angle_distortion Optiona langle distortion 
 *
 * */
double ECCQE(const TVector3& l_momentum, double l_energy, double energy_distortion=0., double angle_distortion=0.);

/** Get oscillation probability of muon neutrino in a 3+1 model. I.e. probability that the numu will stay a numu.
 *
 * \param numu_energy Energy of incident muon neutrino in GeV
 * \param numu_dist Distance travelled by muon neutrino in km
 * \param osc_dm2 dm^2 of sterile netrino in eV^2
 * \param osc_angle Sterile neutrino mixing angle
 *
 * */
double NuMuOscillation(double numu_energy, double numu_dis, double osc_dm2, double osc_angle);


  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_Utilities__

