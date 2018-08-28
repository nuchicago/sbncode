#include <iostream>
#include <cmath> 
#include <vector>

#include <json/json.h>

#include <TH2D.h>
#include <TH1D.h>

#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "core/Event.hh"
#include "NumuSelection.h"
#include "Utilities.h"
#include "util/Interaction.hh"

using namespace util;

namespace ana {
  namespace SBNOsc {

NumuSelection::NumuSelection() : 
  SelectionBase(), 
  _event_counter(0), 
  _nu_count(0), 
  _interactionInfo(new std::vector<NuMuInteraction>) {}


void NumuSelection::Initialize(Json::Value* config) {
  if (config) {
    // setup active volume bounding box
    {
      double xmin, xmax, ymin, ymax, zmin, zmax;
      xmin = xmax = ymin = ymax = zmin = zmax = 0;
      if ((*config)["NumuSelection"].isMember("active_volume")) {
        xmin = (*config)["NumuSelection"]["active_volume"]["xmin"].asDouble();
        xmax = (*config)["NumuSelection"]["active_volume"]["xmax"].asDouble();
        ymin = (*config)["NumuSelection"]["active_volume"]["ymin"].asDouble();
        ymax = (*config)["NumuSelection"]["active_volume"]["ymax"].asDouble();
        zmin = (*config)["NumuSelection"]["active_volume"]["zmin"].asDouble();
        zmax = (*config)["NumuSelection"]["active_volume"]["zmax"].asDouble();
      }
      else if ((*config)["NumuSelection"].isMember("fiducial_volumes")) {
        xmin = std::min_element(_config.fiducial_volumes.begin(), _config.fiducial_volumes.end(), 
          [](const auto& lhs, const auto &rhs) { return lhs.Min()[0] < rhs.Min()[0];})->Min()[0];
        xmax = std::max_element(_config.fiducial_volumes.begin(), _config.fiducial_volumes.end(), 
          [](const auto& lhs, const auto &rhs) { return lhs.Max()[0] < rhs.Max()[0];})->Max()[0];
        ymin = std::min_element(_config.fiducial_volumes.begin(), _config.fiducial_volumes.end(), 
          [](const auto& lhs, const auto &rhs) { return lhs.Min()[1] < rhs.Min()[1];})->Min()[1];
        ymax = std::max_element(_config.fiducial_volumes.begin(), _config.fiducial_volumes.end(), 
          [](const auto& lhs, const auto &rhs) { return lhs.Max()[1] < rhs.Max()[1];})->Max()[1];
        zmin = std::min_element(_config.fiducial_volumes.begin(), _config.fiducial_volumes.end(), 
          [](const auto& lhs, const auto &rhs) { return lhs.Min()[2] < rhs.Min()[2];})->Min()[2];
        zmax = std::max_element(_config.fiducial_volumes.begin(), _config.fiducial_volumes.end(), 
          [](const auto& lhs, const auto &rhs) { return lhs.Max()[2] < rhs.Max()[2];})->Max()[2];
      }
      _config.active_volume = geoalgo::AABox(xmin, ymin, zmin, xmax, ymax, zmax);
    }

    // allow multiple fiducial volumes (accomodate for uboone data channels and icarus 2 TPC's)
    auto FVs = (*config)["NumuSelection"]["fiducial_volumes"];
    for (auto FV: FVs) {
      _config.fiducial_volumes.emplace_back(FV["xmin"].asDouble(), FV["ymin"].asDouble(), FV["zmin"].asDouble(), FV["xmax"].asDouble(), FV["ymax"].asDouble(), FV["zmax"].asDouble());
    }
    _config.doFVCut = (*config)["NumuSelection"].get("doFVcut", true).asBool();
    _config.vertexDistanceCut = (*config)["NumuSelection"].get("vertexDistance", -1).asDouble();
    _config.minLengthContainedLepton = (*config)["NumuSelection"].get("minLengthContainedLepton", -1).asDouble();
    _config.minLengthExitingLepton = (*config)["NumuSelection"].get("minLengthExitingLepton", -1).asDouble();
    _config.verbose = (*config)["NumuSelection"].get("verbose", false).asBool();
  }

  // Setup histo's for root output
  fOutputFile->cd();
  auto cut_names = cutNames();
  for (unsigned i = 0; i < nCuts; i++) {
    _root_histos[i].h_numu_ccqe = new TH1D(("numu_ccqe_" + cut_names[i]).c_str(), "numu_ccqe", 100, 0, 10);
    _root_histos[i].h_numu_trueE = new TH1D(("numu_trueE_" + cut_names[i]).c_str(), "numu_trueE", 100, 0 , 10);
    _root_histos[i].h_numu_visibleE = new TH1D(("numu_visibleE_" + cut_names[i]).c_str(), "numu_visibleE", 100, 0, 10);
    _root_histos[i].h_numu_true_v_visibleE = new TH1D(("numu_true_v_visibleE_" + cut_names[i]).c_str(), "numu_true_v_visibleE", 100, -10, 10);
    _root_histos[i].h_numu_l_length = new TH1D(("numu_l_length_" + cut_names[i]).c_str(), "numu_l_length", 101, -10, 1000);
    _root_histos[i].h_numu_contained_L = new TH1D(("numu_contained_L_" + cut_names[i]).c_str(), "numu_contained_L", 101, -10 , 1000);
    _root_histos[i].h_numu_l_is_contained = new TH1D(("l_is_contained_" + cut_names[i]).c_str(), "l_is_contained", 3, -1.5, 1.5); 
    _root_histos[i].h_numu_Vxy = new TH2D(("numu_Vxy_" + cut_names[i]).c_str(), "numu_Vxy", 
      20, _config.active_volume.Min()[0], _config.active_volume.Max()[0], 
      20, _config.active_volume.Min()[1], _config.active_volume.Max()[1]);
    _root_histos[i].h_numu_Vxz = new TH2D(("numu_Vxz_" + cut_names[i]).c_str(), "numu_Vxz", 
      20, _config.active_volume.Min()[0], _config.active_volume.Max()[0], 
      20, _config.active_volume.Min()[2], _config.active_volume.Min()[2]); 
    _root_histos[i].h_numu_Vyz = new TH2D(("numu_Vyz_" + cut_names[i]).c_str(), "numu_Vyz", 
      20, _config.active_volume.Min()[1], _config.active_volume.Max()[1], 
      20, _config.active_volume.Min()[2], _config.active_volume.Min()[2]);
  }

  // add branches
  fTree->Branch("numu_interaction", &_interactionInfo);

  hello();
}


void NumuSelection::Finalize() {
  // write out histos
  fOutputFile->cd();
  for (unsigned i = 0; i < nCuts; i++) {
    _root_histos[i].h_numu_ccqe->Write();
    _root_histos[i].h_numu_trueE->Write();
    _root_histos[i].h_numu_visibleE->Write();
    _root_histos[i].h_numu_true_v_visibleE->Write();
    _root_histos[i].h_numu_l_length->Write();
    _root_histos[i].h_numu_contained_L->Write();
    _root_histos[i].h_numu_l_is_contained->Write();
    _root_histos[i].h_numu_Vxy->Write();
    _root_histos[i].h_numu_Vxz->Write();
    _root_histos[i].h_numu_Vyz->Write();
  }
}


bool NumuSelection::ProcessEvent(const gallery::Event& ev, std::vector<Event::Interaction>& reco) {
  if (_event_counter % 10 == 0) {
    std::cout << "NumuSelection: Processing event " << _event_counter << " "
              << "(" << _nu_count << " neutrinos selected)"
              << std::endl;
  }

  // clean up containers
  _event_counter++;
  _interactionInfo->clear();

  // Get truth
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);

  // Iterate through the neutrinos
  bool selected = false;
  for (size_t i=0; i<mctruths.size(); i++) {
    auto const& mctruth = mctruths.at(i);
    // get the neutrino
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();

    // build the interaction
    Event::Interaction interaction = TruthReco(ev, mctruth);

    // Get selection-specific info
    //
    // Start with the interaction stuff
    NuMuInteraction intInfo = interactionInfo(ev, mctruth);

    // run selection
    std::array<bool, NumuSelection::nCuts> selection = Select(ev, mctruth, i, intInfo);

    // pass iff pass each cut
    bool pass_selection = std::find(selection.begin(), selection.end(), false) == selection.end();
    if (pass_selection) { 
      // store interaction in reco tree
      reco.push_back(interaction);
      // store local info
      _interactionInfo->push_back(intInfo);

      _nu_count++;
      selected = true;
    }
    
    // fill histos
    for (size_t select_i=0; select_i < selection.size(); select_i++) {
      if (selection[select_i]) {
        _root_histos[select_i].h_numu_trueE->Fill(interaction.neutrino.energy);
        _root_histos[select_i].h_numu_ccqe->Fill(ECCQE(interaction.lepton.momentum, interaction.lepton.energy));
        _root_histos[select_i].h_numu_visibleE->Fill(interaction.neutrino.visible_energy);
        _root_histos[select_i].h_numu_true_v_visibleE->Fill(interaction.neutrino.visible_energy - interaction.neutrino.energy);
        _root_histos[select_i].h_numu_l_length->Fill(intInfo.t_length);
        _root_histos[select_i].h_numu_contained_L->Fill(intInfo.t_contained_length);
        _root_histos[select_i].h_numu_l_is_contained->Fill(intInfo.t_is_contained);
        _root_histos[select_i].h_numu_Vxy->Fill(nu.Nu().Vx(), nu.Nu().Vy());
        _root_histos[select_i].h_numu_Vxz->Fill(nu.Nu().Vx(), nu.Nu().Vz());
        _root_histos[select_i].h_numu_Vyz->Fill(nu.Nu().Vy(), nu.Nu().Vz());
      }
    }
  }
  return selected;
}

NumuSelection::NuMuInteraction NumuSelection::trackInfo(const sim::MCTrack &track) {
  double contained_length = 0;
  double length = 0;
  // If the interaction is outside the active volume, then g4 won't generate positions for the track.
  // So size == 0 => outside FV
  //
  // If size != 0, then we have to check fiducial volume
  bool contained_in_FV = track.size() > 0;
  int pdgid = track.PdgCode(); 

  // Get the length and determine if any point leaves the fiducial volume
  TLorentzVector pos = track.Start().Position();
  for (int i = 1; i < track.size(); i++) {
    // update if track is contained
    if (contained_in_FV) contained_in_FV = containedInFV(pos.Vect());
      
    // update length
    contained_length += containedLength(track[i].Position().Vect(), pos.Vect(), _config.fiducial_volumes);
    length += (track[i].Position().Vect() - pos.Vect()).Mag();
      
    pos = track[i].Position();
  }
  
  return NumuSelection::NuMuInteraction({contained_in_FV, contained_length, length, pdgid});
}

NumuSelection::NuMuInteraction NumuSelection::interactionInfo(const gallery::Event &ev, const simb::MCTruth &mctruth) {
  // get handle to tracks and showers
  auto const& mctrack_list = \
    *ev.getValidHandle<std::vector<sim::MCTrack> >(fMCTrackTag);
  auto const& mcshower_list = \
    *ev.getValidHandle<std::vector<sim::MCShower> >(fMCShowerTag);

  // print out track/shower info
  if (_config.verbose) {
    std::cout << "\n\nINTERACTION:\n";
    std::cout << "MODE: " << mctruth.GetNeutrino().Mode() << std::endl;
    std::cout << "CC: " << mctruth.GetNeutrino().CCNC() << std::endl;
    std::cout << "\nTRACKS:\n";
    for (auto const &mct: mctrack_list) {
      std::cout << "TRACK:\n";
      std::cout << "PDG: " << mct.PdgCode() << std::endl;
      std::cout << "Vertex: " << isFromNuVertex(mctruth, mct) << std::endl;
      std::cout << "Energy: " << mct.Start().E() << std::endl;
      std::cout << "Kinetic: " << mct.Start().E() - PDGMass(mct.PdgCode()) << std::endl;
      std::cout << "Length: " << (mct.Start().Position().Vect() - mct.End().Position().Vect()).Mag() << std::endl;
      std::cout << "Process: " << mct.Process() << std::endl;
    }
    std::cout << "\nSHOWERS:\n";
    for (auto const &mcs: mcshower_list) { 
      std::cout << "SHOWER:\n";
      std::cout << "PDG: " << mcs.PdgCode() << std::endl;
      std::cout << "Vertex: " << isFromNuVertex(mctruth, mcs) << std::endl;
      std::cout << "Energy: " << mcs.Start().E() << std::endl;
      std::cout << "Kinetic: " << mcs.Start().E() - PDGMass(mcs.PdgCode()) << std::endl;
      std::cout << "Length: " << (mcs.Start().Position().Vect() - mcs.End().Position().Vect()).Mag() << std::endl;
      std::cout << "Process: " << mcs.Process() << std::endl;
    }
  }

  // and particles
  auto const& mcparticle_list = \
    *ev.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);

  // Get the length and determine if any point leaves the fiducial volume
  //
  // get lepton track
  int track_ind = -1;
  for (int i = 0; i < mctrack_list.size(); i++) {
    if (isFromNuVertex(mctruth, mctrack_list[i]) && mctrack_list[i].PdgCode() == 13 && mctrack_list[i].Process() == "primary") {
      track_ind = i;
      break;
    }
  } 
  // if there's no lepton, look for a pi+ that can "fake" a muon 
  // if there's multiple, get the longest one
  if (track_ind == -1) {
    double track_contained_length = -1;
    for (int i = 0; i < mctrack_list.size(); i++) {
      if (isFromNuVertex(mctruth, mctrack_list[i]) && mctrack_list[i].PdgCode() == 211 && mctrack_list[i].Process() == "primary") {
        double this_contained_length = trackInfo(mctrack_list[i]).t_contained_length; 
        if (track_contained_length < 0 || this_contained_length > track_contained_length) {
          track_ind = i;
          track_contained_length = this_contained_length;
        }
      }
    }
  }

  return (track_ind != -1) ? trackInfo(mctrack_list[track_ind]) : NumuSelection::NuMuInteraction({false, -1, -1, -1});
}

std::array<bool, NumuSelection::nCuts> NumuSelection::Select(const gallery::Event& ev, const simb::MCTruth& mctruth, 
      unsigned truth_ind, const NumuSelection::NuMuInteraction &intInfo) {
  // get the neutrino
  const simb::MCNeutrino& nu = mctruth.GetNeutrino();

  // has valid track
  bool pass_valid_track = intInfo.t_pdgid != -1;

  // pass fiducial volume cut
  bool pass_FV = passFV(nu.Nu().Position().Vect());

  // min length cut
  bool pass_min_length = passMinLength(intInfo.t_contained_length, intInfo.t_is_contained);

  // pass vertex reconstruction cut
  bool pass_reco_vertex = true;
  if (_config.vertexDistanceCut > 0) {
    double truth_v[3];
    truth_v[0] = nu.Nu().Vx();
    truth_v[1] = nu.Nu().Vy();
    truth_v[2] = nu.Nu().Vz();
  
    // get the reco vertex information
    auto const pfp_handle = ev.getValidHandle<std::vector<recob::PFParticle>>("pandoraNu");
    art::FindMany<recob::Vertex> fvtx(pfp_handle, ev, "pandoraNu");
    // index of truth info is same as index of vertex info
    std::vector<const recob::Vertex*> vertices = fvtx.at(truth_ind);
    // neutrino will only have one associated vetex
    auto vertex = vertices[0]->position();
    TVector3 reco_v(vertex.X(), vertex.Y(), vertex.Z());
    
    pass_reco_vertex = passRecoVertex(nu.Nu().Position().Vect(), reco_v);
  }

  // print selection information
  if (_config.verbose) {
    std::cout << "NEW EVENT" << std::endl;
    std::cout << "CCNC: " << nu.CCNC() << " MODE: " << nu.Mode() << " PDG: " << nu.Nu().PdgCode() << std::endl;
    std::cout << "Track PDG: " << intInfo.t_pdgid <<std::endl;
    std::cout << "pass Valid Track: " << pass_valid_track << std::endl;
    std::cout << "Pos: " << nu.Nu().Vx() << " " << nu.Nu().Vy() << " " << nu.Nu().Vz() << std::endl;
    std::cout << "pass FV: " << pass_FV << std::endl;
    std::cout << "pass Reco: " << pass_reco_vertex << std::endl;
    std::cout << "Length: " << intInfo.t_contained_length << " Contained: " << intInfo.t_is_contained << std::endl;
    std::cout << "pass Length: " << pass_min_length << std::endl;
  }

  // retrun list of cuts
  return {pass_valid_track, pass_FV, pass_valid_track && pass_FV && pass_min_length, pass_valid_track && pass_FV && pass_reco_vertex};
}

bool NumuSelection::containedInFV(const TVector3 &v) {
  if (!_config.doFVCut) return true;

  geoalgo::Point_t p(v); 
  for (auto const& FV: _config.fiducial_volumes) {
    if (FV.Contain(p)) return true;
  }
  return false;
}

bool NumuSelection::passRecoVertex(const TVector3 &truth_v, const TVector3 &reco_v) {
  if (_config.vertexDistanceCut < 0) return true;

  return (truth_v - reco_v).Mag() < _config.vertexDistanceCut;
}

bool NumuSelection::passMinLength(double length, bool stop_in_tpc) {
  if (!stop_in_tpc)
    return _config.minLengthExitingLepton < 0 || length > _config.minLengthExitingLepton;
  else
    return _config.minLengthContainedLepton < 0 || length > _config.minLengthContainedLepton;
}

  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NumuSelection)

