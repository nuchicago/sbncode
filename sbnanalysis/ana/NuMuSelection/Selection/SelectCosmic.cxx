#include <iostream>
#include <vector>
#include <TH2D.h>
#include <json/json.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/MCSFitResult.h"

#include "SelectCosmic.h"
#include "../DataTypes/RecobInteraction.h"


namespace ana {
  namespace NuMuSelection {

SelectCosmic::SelectCosmic() : SelectionBase() {}


void SelectCosmic::Initialize(Json::Value* config) {
  _reco_interactions = new std::vector<RecobInteraction>();
  AddBranch("reco_interactions", _reco_interactions);
}


void SelectCosmic::Finalize() {
}


bool SelectCosmic::ProcessEvent(gallery::Event& ev) {
  std::cout << "\nNEW EVENT!\n";

  // clear out branch containers
  _reco_interactions->clear();

  // Get MCTruth information
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);  
  // get Pandora PFParticles
  auto const pfp_handle = ev.getValidHandle<std::vector<recob::PFParticle>>("pandoraNu");
  // get MCS reco info
  auto const mcs_handle = ev.getValidHandle<std::vector<recob::MCSFitResult>>("pandoraNuMCSMu"); 

  // get associations to tracks and vertcies
  art::FindManyP<recob::Track> ftrk(pfp_handle, ev, "pandoraNu");
  art::FindManyP<recob::Vertex> fvtx(pfp_handle, ev, "pandoraNu");

  // loop over PFP particles and get neutrino candidates
  int i = 0;
  for (auto const& pfp: *pfp_handle) {
    // interaction object
    RecobInteraction interaction;

    // Get the primary PFParticle (the neutrino)
    if (pfp.IsPrimary() && pfp.PdgCode() == 14) {
      std::cout << "PFParticle!" << std::endl;
      std::cout << pfp;
    
      // get vertex and track and momentum
      std::vector<art::Ptr<recob::Track>> tracks = ftrk.at(i);
      std::vector<art::Ptr<recob::Vertex>> vertices = fvtx.at(i);
    
      // Print positions
      for (unsigned j = 0; j < vertices.size(); j++) {
        double xyz[3];
        vertices[j]->XYZ(xyz);
        std::cout << "Pos: " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " " << std::endl;
        //TVector3 pfp_pos(xyz);
      } 

      // get daughter particles
      std::vector<size_t> const& daughter_particles = pfp.Daughters();
      std::vector<double> d_mcs_momenta;
      std::vector<double> d_truth_momenta;
      for (size_t daughter_ind: daughter_particles) {
        auto const& daughter = pfp_handle->at(daughter_ind);
        std::cout << "Daughter: " << daughter;
        std::cout << "At: " << daughter_ind << std::endl;

        // get tracks and vertices of daughter
        std::vector<art::Ptr<recob::Track>> d_tracks = ftrk.at(daughter_ind);
        std::vector<art::Ptr<recob::Vertex>> d_vertices = fvtx.at(daughter_ind);

        // if tracks, get MCS momenta
        if (d_tracks.size() != 0) {
	  // get the mcs from the daughter track
	  auto const& d_mcs_momentum = mcs_handle->at(d_tracks[0].key());
	  
	  std::cout << "Forward P: " << d_mcs_momentum.fwdMomentum() << std::endl;
	  std::cout << "Backward P: " << d_mcs_momentum.bwdMomentum() << std::endl;
	  // and push it into the vector   
	  d_mcs_momenta.push_back(d_mcs_momentum.fwdMomentum());
        }
      }

      // use largest momentum of muon candidates?
      interaction.l_mcs_momentum = *std::max_element(d_mcs_momenta.begin(), d_mcs_momenta.end());
      // all muons
      interaction.l_pandora_pdgid = 13; 

      _reco_interactions->push_back(interaction);     

    }
    i ++;
  }
  return true;
}

  }  // namespace NuMuSelection
}  // namespace ana


// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::NuMuSelection::SelectCosmic)

