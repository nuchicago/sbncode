#ifndef __UTILITY_INTERACTION_HH
#define __UTILITY_INTERACTION_HH

#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "gallery/Event.h"
#include "canvas/Utilities/InputTag.h"

#include "core/Event.hh"

namespace util {
/* Get the "visible" energy from a neutrino interaction. Is equal to sum of non-neutral hadronic kinetic energies and lepton total energies. 
 * Does not account for containment or tracking.
 *
 * \param ev The gallery event.
 * \param mctruth The MCTruth object corresponding to the interaction.
 * \param mctTag Art tag to get a handle to MCTruth objects.
 * \param mcsTag Art tag to get a handle to MCShower objects.
 * \returns Visble energy in MeV.
 */
double visibleEnergy(const gallery::Event &ev, const simb::MCTruth &mctruth, const art::InputTag &mctTag={ "mcreco" }, const art::InputTag &mcsTag={ "mcreco" });

/** Get mass from PDGID of particle in MeV/c^2.
 *
 * \param pdg The Particle Data Group ID of the particle (as returned by i.e. an MCTruth object)
 *
 * \returns Mass of particle in MeV/c^2
 *
 * */
double PDGMass(int pdg);

/** Get charge from PDGID of particle in |e|/3.
 * \param pdg The Particle Data Group ID of the particle (as returned by i.e. an MCTruth object)
 *
 * \returns Charge of particle in |e|/3.
 */
double PDGCharge(int pdg);

/** Returns whether track/shower object is from the neutrino vertex
 *
 * \param mc MCTruth corresponding to neutrino interaction
 * \param show The object to be matched
 * \param distance between shower start and interaction vertex
 * \returns Whether track/shower object is from neutrino vertex
 * */
bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCShower& show,
                           float distance=5.0);
/** Returns whether track/shower object is from the neutrino vertex
 *
 * \param mc MCTruth corresponding to neutrino interaction
 * \param track The object to be matched
 * \param distance between track start and interaction vertex
 * \returns Whether track/shower object is from neutrino vertex
 * */
bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track,
                            float distance=5.0);

}  // namespace util
#endif
