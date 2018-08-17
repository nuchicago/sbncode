#include <iostream>
#include <vector>
#include <TH2D.h>
#include <json/json.h>
#include <cmath>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "gallery/Event.h"
#include <TLorentzVector.h>
#include "ExampleSelection.h"
#include "ExampleTools.h"
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <random>
#include "lardataobj/Simulation/SimPhotons.h"

void ActExp() {
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
}
