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

#include "TDatabasePDG.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"

#include "core/Event.hh"

#include "uboone/LLBasicTool/GeoAlgo/GeoAABox.h"

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

/** Get mass from PDGID of particle.
 *
 * \param pdg The Particle Data Group ID of the particle (as returned by i.e. an MCTruth object)
 *
 * */
double PDGMass(int pdg);

/** Returns whether track/shower object is from the neutrino vertex
 *
 * \param mc MCTruth corresponding to neutrino interaction
 * \param show The object to be matched
 * \param distance between shower start and interaction vertex
 * */
bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCShower& show,
                           float distance=5.0);
/** Returns whether track/shower object is from the neutrino vertex
 *
 * \param mc MCTruth corresponding to neutrino interaction
 * \param track The object to be matched
 * \param distance between track start and interaction vertex
 * */
bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track,
                            float distance=5.0);

/* Finds length of line segment contained inside AABox
 *
 * \param v0 the first point of the line segment
 * \param v1 the second point of the line segment
 * \param boxes a list of fiducial volumes instantiated as AABoxes
 * 
 * */
double containedLength(const TVector3 &v0, const TVector3 &v1, const std::vector<geoalgo::AABox> &boxes);


  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_Utilities__

