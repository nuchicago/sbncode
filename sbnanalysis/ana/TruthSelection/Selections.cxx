#include <cassert>
#include <iostream>
#include <numeric>
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

namespace ana {
  namespace TruthSelection {
    namespace selections {

bool CCNueTrue(const std::vector<simb::MCTruth>& mctruths,
               const std::vector<sim::MCTrack>& mctracks,
               const std::vector<sim::MCShower>& mcshowers,
               double& ereco, double& weight) {
  weight = 0.8;

  // Pass event if 1 or more true nue CC interactions
  bool ccnue = false;
  for (auto const& t : mctruths) {
    if (std::abs(t.GetNeutrino().Nu().PdgCode()) == 12 &&
        t.GetNeutrino().CCNC() == simb::kCC) {
      ccnue = true;
      break;
    }
  }

  if (!ccnue) {
    return false;
  }

  // Require one electron shower over threshold
  for (auto const& s : mcshowers) {
    if (s.Process() == "primary" &&
        std::abs(s.PdgCode()) == 11 &&
        s.Start().E() - 511e-3 > 200) {
      ereco = s.Start().E();
      return true;
    }
  }

  return false;
}


bool CCNumuTrue(const std::vector<simb::MCTruth>& mctruths,
                const std::vector<sim::MCTrack>& mctracks,
                const std::vector<sim::MCShower>& mcshowers,
                double& ereco, double& weight) {
  weight = 1.0;
  // Pass event if 1 or more true numu CC interactions
  bool ccnum = false;
  for (auto const& t : mctruths) {
    if (std::abs(t.GetNeutrino().Nu().PdgCode()) == 14 &&
        t.GetNeutrino().CCNC() == simb::kCC) {
      ccnum = true;
      break;
    }
  }

  if (!ccnum) {
    return false;
  }

  // Require one muon track over threshold
  for (auto const& t : mctracks) {
    if (t.Process() == "primary" &&
        std::abs(t.PdgCode()) == 13 &&
        t.Start().E() - 105.6 > 60) {
      ereco = t.Start().E();
      return true;
    }
  }

  return false;
}


bool True1l1p0pi0(const std::vector<simb::MCTruth>& mctruths,
                  const std::vector<sim::MCTrack>& mctracks,
                  const std::vector<sim::MCShower>& mcshowers,
                  int leptonPDG,
                  double& ereco, double& weight) {
  weight = 1.0;

  unsigned np = 0;  // Protons
  unsigned nl = 0;  // Leptons
  unsigned nt = 0;  // Other tracks
  unsigned ns = 0;  // Other showers

  for (auto const& t : mctracks) {
    if (t.Process() != "primary" || t.Origin() != simb::kBeamNeutrino) {
      continue;
    }

    if (t.PdgCode() == leptonPDG) {
      if (t.Start().E() - 105.6 > 60) {
        ereco = t.Start().E();
        nl++;
      }
    }
    else if (abs(t.PdgCode()) == 2212) {
      if (t.Start().E() - 938.3 > 30) {
        np++;
      }
    }
    else if (abs(t.PdgCode()) == 211) {
      if (t.Start().E() - 139.5 > 35) { 
        nt++;
      }
    }
    else if (t.PdgCode() != 2112 || t.PdgCode() < 100000) {
      nt++;  // Some other track with unknown PDG, cut it
    }
  }

  for (auto const& t : mcshowers) {
    if (t.Process() != "primary" || t.Origin() != simb::kBeamNeutrino) {
      continue;
    }

    if (t.PdgCode() == leptonPDG) {
      if (t.Start().E() - 0.511 > 30) {
        ereco = t.Start().E();
        nl++;
      }
    }
    else if ((t.PdgCode() == 11  && t.Start().E() - 0.511 > 30) ||
             (t.PdgCode() == 22  && t.Start().E() > 30) ||
             (t.PdgCode() == 111)) {
      ns++;
    }
    else {
      ns++;  // Some other shower with unknown PDG, cut it
    }
  }

  return (np == 1 && nl == 1 && nt == 0 && ns == 0);
}


bool CCNue(const std::vector<simb::MCTruth>& mctruths,
           const std::vector<sim::MCTrack>& mctracks,
           const std::vector<sim::MCShower>& mcshowers,
           double& ereco, double& weight) {
  
  assert(false);
  for (auto const& track : mctracks) {
    // No tracks over 1 m
    double d = (track.End().Position().Vect() -
                track.Start().Position().Vect()).Mag();
    if (d > 100) {
      return false;
    }
  }

  return true;
}


bool CCNumu(const std::vector<simb::MCTruth>& mctruths,
            const std::vector<sim::MCTrack>& mctracks,
            const std::vector<sim::MCShower>& mcshowers,
            double& ereco, double& weight) {
  assert(false);
}


bool CCPi0(const std::vector<simb::MCTruth>& mctruths,
           const std::vector<sim::MCTrack>& mctracks,
           const std::vector<sim::MCShower>& mcshowers,
           double& ereco, double& weight) {
  bool cc = false;
  for (auto const& t : mctruths) {
    if (t.GetNeutrino().CCNC() == simb::kCC) {
      cc = true;
      break;
    }
  }

  if (!cc) {
    return false;
  }

  for (auto const& t : mctracks) {
    if (t.Process() != "primary" || t.Origin() != simb::kBeamNeutrino) {
      continue;
    }

    if (t.PdgCode() == 111) {
      ereco = t.Start().E();
      return true;
    }
  }

  return false;
}

    }  // namespace selections
  }  // namespace TruthSelection
}  // namespace ana



