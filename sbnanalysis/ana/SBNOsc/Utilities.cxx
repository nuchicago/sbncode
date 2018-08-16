#include <iostream>
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include "core/Event.hh"

namespace ana {
  namespace SBNOsc {

void hello() {
  std::cout << "Hello SBNOsc!" << std::endl;
}


Event::Interaction TruthReco(const simb::MCTruth& mctruth) {
  Event::Interaction interaction;

  // Neutrino
  const simb::MCNeutrino& nu = mctruth.GetNeutrino();
  interaction.neutrino.energy = nu.Nu().EndMomentum().Energy();

  // Primary lepton
  const simb::MCParticle& lepton = nu.Lepton();
  interaction.lepton.pdg = lepton.PdgCode();
  interaction.lepton.energy = lepton.Momentum(0).Energy();
  interaction.lepton.momentum = lepton.Momentum(0).Vect();

  // Hadronic system
  for (int iparticle=0; iparticle<interaction.finalstate.size(); iparticle++) {
    Event::FinalStateParticle fsp;
    const simb::MCParticle& particle = mctruth.GetParticle(iparticle);

    if (particle.Process() != "primary") {
      continue;
    }

    fsp.pdg = particle.PdgCode();
    fsp.energy = particle.Momentum(0).Energy();
    fsp.momentum = particle.Momentum(0).Vect();

    interaction.finalstate.push_back(fsp);
  }

  return interaction;
}

bool IsUsableTrack(const sim::MCTrack& mctrack) {
  int GoodStepCount =0;
  for (size_t i=0;i<mctrack.size();i++) {
    auto const& mcstep = mctrack.at(i);
    auto stepX = mcstep.X();
    auto stepY = mcstep.Y();
    auto stepZ = mcstep.Z();
    if(((-199.15 < stepX && stepX < -2.65) || (2.65 < stepX && stepX < 199.15)) && (-200 < stepY && stepY < 200) && (0 < stepZ && stepZ < 500)) GoodStepCount++;
  }
  if (GoodStepCount>1) return true;
}

double GetActiveLength(const sim::MCTrack& mctrack) {
  const sim::MCStep entering_step_obj;
  const sim::MCStep exiting_step_obj;
  const sim::MCStep& entering_step = entering_step_obj;
  const sim::MCStep& exiting_step = exiting_step_obj;
  for (size_t i=0;i<mctrack.size();i++) {
    auto const& mcstep = mctrack.at(i);
    auto stepX = mcstep.X();
    auto stepY = mcstep.Y();
    auto stepZ = mcstep.Z();
    if(((-199.15 < stepX && stepX < -2.65) || (2.65 < stepX && stepX < 199.15)) && (-200 < stepY && stepY < 200) && (0 < stepZ && stepZ < 500)) {
      entering step = mcstep;
      break;
    }
    else continue;
  }

  for (size_t i=0;i<mctrack.size();i++){
    auto const& reverse_mcstep = mctrack.at(mctrack.size()-1-i);
    auto reverse_mcstepX = reverse_mcstep.X();
    auto reverse_mcstepY = reverse_mcstep.Y();
    auto reverse_mcstepZ = reverse_mcstep.Z();
    bool IsOut = !(((-199.15 < reverse_mcstepX && reverse_mcstepX < -2.65) || (2.65 < reverse_mcstepX && reverse_mcstepX < 199.15)) && (-200 < reverse_mcstepY && reverse_mcstepY < 200) && (0 < reverse_mcstepZ && reverse_mcstepZ < 500));
    auto const& reverse_mcstep_former = mctrack.at(mctrack.size()-2-i);
    auto reverse_mcstepX_former = reverse_mcstep_former.X();
    auto reverse_mcstepY_former = reverse_mcstep_former.Y();
    auto reverse_mcstepZ_former = reverse_mcstep_former.Z();
    bool FormerOut = !(((-199.15 < reverse_mcstepX_former && reverse_mcstepX_former < -2.65) || (2.65 < reverse_mcstepX_former && reverse_mcstepX_former < 199.15)) && (-200 < reverse_mcstepY_former && reverse_mcstepY_former < 200) && (0 < reverse_mcstepZ_former && reverse_mcstepZ_former < 500));
    if (IsOut && (!FormerOut)) {
      exiting_step = reverse_mcstep_former;
      break;
    }
    if ((!IsOut)&&(!FormerOut)) {
      auto const& exiting_step = reverse_mcstep;
      break;
    }
    else continue;
  }
  double active_length = (exiting_step.Position().Vect()-entering_step.Position().Vect()).Mag();
  return active_length;
  }


  }  // namespace SBNOsc
}  // namespace ana
