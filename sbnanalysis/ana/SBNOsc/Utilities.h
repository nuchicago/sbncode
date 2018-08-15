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

/** Calculate CCQE energy from associated lepton information (and optional distortion). Energy in GeV. */
double ECCQE(const TVector3& l_momentum, double l_energy, double energy_distortion=0., double angle_distortion=0.);


  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_Utilities__

