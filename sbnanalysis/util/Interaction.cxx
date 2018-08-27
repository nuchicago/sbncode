#include "TDatabasePDG.h"

#include "Interaction.hh"

namespace util {

double visibleEnergy(const gallery::Event &ev, const simb::MCTruth &mctruth, const art::InputTag &mctTag, const art::InputTag &mcsTag) {
  // pull up mctracks and showers
  auto const& mctrack_list = \
    *ev.getValidHandle<std::vector<sim::MCTrack> >(mctTag);
  auto const& mcshower_list = \
    *ev.getValidHandle<std::vector<sim::MCShower> >(mcsTag);

  double visible_E = 0;

  // total up visible energy from tracks...
  for (auto const &mct: mctrack_list) {
    // ignore neutral particles
    if (isFromNuVertex(mctruth, mct) && abs(PDGCharge(mct.PdgCode())) > 1e-4) {
      double mass = (mct.PdgCode() == 13 || mct.PdgCode() == 11) ? 0:PDGMass(mct.PdgCode());
      double this_visible_energy = mct.Start().E() - mass;
      visible_E += this_visible_energy;
    }
  }

  // ...and showers
  for (auto const &mcs: mcshower_list) {
    if (isFromNuVertex(mctruth, mcs) && abs(PDGCharge(mcs.PdgCode())) > 1e-4) {
      double mass = PDGMass(mcs.PdgCode());
      double this_visible_energy = mcs.Start().E() - mass;
      visible_E += this_visible_energy;
    }
  }

  return visible_E;
}

// define global static const PDGTable to be used by helper functions
static const TDatabasePDG *PDGTable(new TDatabasePDG);

double PDGMass(int pdg) {
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

double PDGCharge(int pdg) {
  // regular particle
  if (pdg < 1000000000) {
    TParticlePDG* ple = PDGTable->GetParticle(pdg);
    return ple->Charge();
  }
  // ion
  else {
    int p = (pdg % 10000000) / 10000;
    return p * 3;
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


}
