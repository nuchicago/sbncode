#include <iostream>

#include "TDatabasePDG.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "core/Event.hh"

#include "Utilities.h"

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

double ECCQE(const TVector3& l_momentum, double l_energy, double energy_distortion, double angle_distortion) {
  // Based on D. Kaleko, LowEnergyExcess LArLite module ECCQECalculator
  double M_n = 0.939565; // GeV/c^2
  double M_p = 0.938272; // GeV/c^2
  double M_e = 0.000511; // GeV/c^2
  double bindingE = 0.0300; // GeV
  
  double mp2 = M_p * M_p;
  double me2 = M_e * M_e;
  double mnb = M_n - bindingE;

  // mess around with lorentz vector
  TVector3 v(l_momentum);
  v.SetTheta( v.Theta() + angle_distortion);
  l_energy = l_energy + energy_distortion;
  double l_mom = sqrt(l_energy*l_energy - me2);
  double l_theta = \
    acos(v.Pz() / sqrt(v.Px()*v.Px() + v.Py()*v.Py() + v.Pz()*v.Pz()));
  double enu_top = mp2 - mnb*mnb - me2 + 2.0 * mnb * l_energy;
  double enu_bot = 2.0 * (mnb - l_energy + l_mom * cos(l_theta));
  return enu_top / enu_bot;
}

double NuMuOscillation(double numu_energy, double numu_dist, double osc_dm2, double osc_angle) {
  double overlap = sin(2*osc_angle);
  double energy_factor = sin(1.27 * osc_dm2 * numu_dist / numu_energy);
  return 1 - overlap * overlap * energy_factor * energy_factor;

}

double PDGMass(int pdg) {
  TDatabasePDG *PDGTable = new TDatabasePDG;
  // regular particle
  if (pdg < 1000000000) {
    TParticlePDG* ple = PDGTable->GetParticle(pdg);
    return ple->Mass() * 1000.0;
  }
  // ion
  else {
    int p = (pdg % 10000000) / 10000;
    int n = (pdg % 10000) / 10 - p;
    return (PDGTable->GetParticle(2212)->Mass() * p +
            PDGTable->GetParticle(2112)->Mass() * n) * 1000.0;
  }
}

bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCShower& show, float distance)  {
  TLorentzVector nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0);
  TLorentzVector showStart = show.Start().Position();
  return (showStart - nuVtx).Mag() < distance;
}

bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track, float distance) {
  TLorentzVector nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0);
  TLorentzVector trkStart = track.Start().Position();
  return (trkStart - nuVtx).Mag() < distance;
}

  }  // namespace SBNOsc
}  // namespace ana

