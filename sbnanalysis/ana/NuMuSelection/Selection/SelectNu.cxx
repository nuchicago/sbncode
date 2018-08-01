#include <iostream>
#include <vector>
#include <TH2D.h>
#include <json/json.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "SelectNu.h"

namespace ana {
  namespace NuMuSelection {

SelectNu::SelectNu() : SelectionBase() {}


void SelectNu::Initialize(Json::Value* config) {
  // tree filled per truth neutrino
  _neutrinoTree = new TTree("neutrino", "Truth Neutrino Tree");
  // distance from MCNeutrino to closest reconstructed neutrino
  _neutrinoTree->Branch("min_reconstructed_R", &_min_reconstructed_R);

  // add in custom TBranches
  // distance associated w/ event tree
  _min_reconstructed_Rs = new std::vector<double>();
  AddBranch("min_reconstructed_R", &_min_reconstructed_Rs);

  // setup config
  _config.dist_cut = config->get("distance_cut", -1.).asDouble();
}


void SelectNu::Finalize() {
  // write the TTree
  fOutputFile->cd();
  _neutrinoTree->Write("neutrino", TObject::kOverwrite);
}


bool SelectNu::ProcessEvent(gallery::Event& ev) {
  std::cout << "\nNew Event\n";

  // clear branch containers
  _min_reconstructed_Rs->clear();

  // Get MCTruth information
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);  
  // get Pandora PFParticles
  auto const pfp_handle = ev.getValidHandle<std::vector<recob::PFParticle>>("pandoraNu");

  // get associations to tracks and vertcies
  art::FindMany<recob::Track> ftrk(pfp_handle, ev, "pandoraNu");
  art::FindMany<recob::Vertex> fvtx(pfp_handle, ev, "pandoraNu");
  
  // loop over truth and find associated PFParticle
  for (auto const& truth: mctruths) {
    // closest distance
    _min_reconstructed_R = -1;

    std::cout << "MCTruth" << std::endl;
    auto const& neutrino = truth.GetNeutrino().Nu();
    std::cout << "Pos: " << neutrino.EndX() << " " << neutrino.EndY() << " " << neutrino.EndZ() << std::endl;

    int i = 0;
    for (auto const& pfp: *pfp_handle) {
      // Get the primary PFParticle (the neutrino)
      if (pfp.IsPrimary() && pfp.PdgCode() == 14) {
        std::cout << "PFParticle!" << std::endl;
        std::cout << pfp;

	// get vertex and track
	std::vector<const recob::Track*> tracks = ftrk.at(i);
	std::vector<const recob::Vertex*> vertices = fvtx.at(i);

        // Print positions
        for (unsigned j = 0; j < vertices.size(); j++) {
          double xyz[3];
          vertices[j]->XYZ(xyz);
          std::cout << "Pos: " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " " << std::endl;

          TVector3 pfp_pos(xyz);
          double R = (neutrino.EndPosition().Vect() - pfp_pos).Mag();
          if (_min_reconstructed_R < 0 || R < _min_reconstructed_R) _min_reconstructed_R = R; 
        } 
      }
      i ++;
    }
    std::cout << "Closest reconstructed: " << _min_reconstructed_R << std::endl;

    // fill per-neutrino containers
    _min_reconstructed_Rs->push_back(_min_reconstructed_R);

    // fill the per-neutrino tree
    _neutrinoTree->Fill();
  } 
   
  // how to deal with multiple neutrino events?
  // for now, only accept event if all neutrinos were reconstructed
  for (double r: *_min_reconstructed_Rs) {
    if (_config.dist_cut > 0 && r > _config.dist_cut) return false;
  }
  return true;
}


  }  // namespace NuMuSelection
}  // namespace ana


// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::NuMuSelection::SelectNu)

