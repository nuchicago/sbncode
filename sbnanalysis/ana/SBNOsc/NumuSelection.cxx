#include <iostream>
#include <math.h> 
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

namespace ana {
  namespace SBNOsc {

NumuSelection::NumuSelection() : 
  SelectionBase(), 
  fEventCounter(0), 
  fNuCount(0), 
  _interactionInfo(new std::vector<NuMuInteraction>) {}


void NumuSelection::Initialize(Json::Value* config) {
  // Load configuration parameters
  _config.truthTag = { "generator" };
  _config.mctTag = { "mcreco" };
  _config.mcsTag = { "mcreco" };
  _config.mcpTag = { "largeant" };

  if (config) {
    _config.truthTag = { (*config)["SBNOsc"].get("MCTruthTag", "generator").asString() };
    _config.mctTag = { (*config)["SBNOsc"].get("MCTrackTag", "mcreco").asString() };
    _config.mcsTag = { (*config)["SBNOsc"].get("MCShowerTag", "mcreco").asString() };
    _config.mcpTag = { (*config)["SBNOsc"].get("MCParticleTag", "largeant").asString() };

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
    _config.minLength = (*config)["NumuSelection"].get("minLength", -1).asDouble();
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
    _root_histos[i].h_numu_Vxy->Write();
    _root_histos[i].h_numu_Vxz->Write();
    _root_histos[i].h_numu_Vyz->Write();
  }
}


bool NumuSelection::ProcessEvent(const gallery::Event& ev, std::vector<Event::Interaction>& reco) {
  if (fEventCounter % 10 == 0) {
    std::cout << "NumuSelection: Processing event " << fEventCounter << " "
              << "(" << fNuCount << " neutrinos selected)"
              << std::endl;
  }

  // clean up containers
  fEventCounter++;
  _interactionInfo->clear();

  // Get truth
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth> >(_config.truthTag);

  // Iterate through the neutrinos
  bool selected = false;
  for (size_t i=0; i<mctruths.size(); i++) {
    auto const& mctruth = mctruths.at(i);
    // get the neutrino
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();

    // build the interaction
    Event::Interaction interaction = TruthReco(mctruth);

    // Get selection-specific info
    //
    // Start with the interaction stuff
    auto intInfo = interactionInfo(ev, mctruth);

    // run selection
    std::vector<bool> selection = Select(ev, mctruth, i, intInfo);

    // pass iff pass each cut
    bool pass_selection = std::find(selection.begin(), selection.end(), false) == selection.end();
    if (pass_selection) { 
      // store interaction in reco tree
      reco.push_back(interaction);
      // store local info
      _interactionInfo->push_back(intInfo);

      fNuCount++;
      selected = true;
    }
    
    // fill histos
    unsigned select_i = 0; 
    for (bool pass: selection) {
      if (pass) {
        _root_histos[select_i].h_numu_trueE->Fill(interaction.neutrino.energy);
        _root_histos[select_i].h_numu_ccqe->Fill(ECCQE(interaction.lepton.momentum, interaction.lepton.energy));
        _root_histos[i].h_numu_visibleE->Fill(intInfo.visible_energy);
        _root_histos[i].h_numu_true_v_visibleE->Fill(intInfo.visible_energy - interaction.neutrino.energy);
        _root_histos[select_i].h_numu_Vxy->Fill(nu.Nu().Vx(), nu.Nu().Vy());
        _root_histos[select_i].h_numu_Vxz->Fill(nu.Nu().Vx(), nu.Nu().Vz());
        _root_histos[select_i].h_numu_Vyz->Fill(nu.Nu().Vy(), nu.Nu().Vz());
      }
      select_i ++;
    }
  }
  return selected;
}

NumuSelection::NuMuInteraction NumuSelection::interactionInfo(const gallery::Event &ev, const simb::MCTruth &mctruth) {
  // get handle to tracks and showers
  auto const& mctrack_list = \
    *ev.getValidHandle<std::vector<sim::MCTrack> >(_config.mctTag);
  auto const& mcshower_list = \
    *ev.getValidHandle<std::vector<sim::MCShower> >(_config.mcsTag);

  // and particles
  auto const& mcparticle_list = \
    *ev.getValidHandle<std::vector<simb::MCParticle>>(_config.mcpTag);

  // get visible energy
  // "Reconstruct" visible energy from tracks + showers
  double visible_E = 0.;
  // first the tracks
  for (auto const &mct: mctrack_list) {
    if (isFromNuVertex(mctruth, mct)) {
      double mass = PDGMass(mct.PdgCode());
      visible_E += mct.Start().E() - mass;
    }
  }
  // now the showers
  for (auto const &mcs: mcshower_list) {
    if (isFromNuVertex(mctruth, mcs)) {
      double mass = PDGMass(mcs.PdgCode());
      visible_E += mcs.Start().E() - mass;
    }
  }

  // get lepton track
  int lepton_ind = -1;
  for (int i = 0; i < mctrack_list.size(); i++) {
    if (isFromNuVertex(mctruth, mctrack_list[i]) && mctrack_list[i].PdgCode() == 13 && mctrack_list[i].size() > 0) {
      lepton_ind = i;
      break;
    }
  } 

  // return failure if there was no lepton
  if (lepton_ind == -1) {
    return { false, -1, visible_E };
  }

  auto const& lepton_track = mctrack_list.at(lepton_ind);

  // Get the length and determine if any point leaves the fiducial volume
  bool contained_in_FV = true;
  double l_contained_length = 0;
  TLorentzVector pos = lepton_track.Start().Position();
  for (int i = 1; i < lepton_track.size(); i++) {
    // update if track is contained
    if (contained_in_FV) contained_in_FV = containedInFV(pos.Vect());

    // update length
    l_contained_length += containedLength(lepton_track[i].Position().Vect(), pos.Vect(), _config.fiducial_volumes);

    pos = lepton_track[i].Position();
  }

  // determine if stopped in TPC
  // sim::MCStep end = lepton_track.End();
  // geoalgo::Point_t l_end_pos(end.Position().Vect());
  // bool stop_in_tpc = _config.active_volume.Contain(l_end_pos);

  return {contained_in_FV, l_contained_length, visible_E}; 
}

std::vector<bool> NumuSelection::Select(const gallery::Event& ev, const simb::MCTruth& mctruth, unsigned truth_ind, const NumuSelection::NuMuInteraction &intInfo) {
  // get the neutrino
  const simb::MCNeutrino& nu = mctruth.GetNeutrino();

  // select true CC events
  bool pass_true_CC = (nu.CCNC() == simb::kCC && (nu.Mode() == 0 || nu.Mode() == 10) && nu.Nu().PdgCode() == 14);

  // pass fiducial volume cut
  bool pass_FV = passFV(nu.Nu().Position().Vect());

  // min length cut
  bool pass_min_length = passMinLength(intInfo.l_contained_length, intInfo.l_is_contained);

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
    double reco_v[3];
    vertices[0]->XYZ(reco_v);
    
    pass_reco_vertex = passRecoVertex(truth_v, reco_v);
  }

  // print selection information
  if (_config.verbose) {
    std::cout << "NEW EVENT" << std::endl;
    std::cout << "CCNC: " << nu.CCNC() << " MODE: " << nu.Mode() << " PDG: " << nu.Nu().PdgCode() << std::endl;
    std::cout << "pass CC: " << pass_true_CC << std::endl;
    std::cout << "Pos: " << nu.Nu().Vx() << " " << nu.Nu().Vy() << " " << nu.Nu().Vz() << std::endl;
    std::cout << "pass FV: " << pass_FV << std::endl;
    std::cout << "pass Reco: " << pass_reco_vertex << std::endl;
  }

  // retrun list of cuts
  return {pass_true_CC, pass_true_CC && pass_FV, pass_true_CC && pass_FV && pass_min_length, pass_true_CC && pass_FV && pass_reco_vertex};
}

bool NumuSelection::containedInFV(const TVector3 &v) {
  if (!_config.doFVCut) return true;

  geoalgo::Point_t p(v); 
  for (auto const& FV: _config.fiducial_volumes) {
    if (FV.Contain(p)) return true;
  }
  return false;
}

bool NumuSelection::passRecoVertex(double truth_v[3], double reco_v[3]) {
  if (_config.vertexDistanceCut < 0) return true;

  double R2 =  
    (truth_v[0] - reco_v[0]) * (truth_v[0] - reco_v[0]) +
    (truth_v[1] - reco_v[1]) * (truth_v[1] - reco_v[1]) +
    (truth_v[2] - reco_v[2]) * (truth_v[2] - reco_v[2]);
 
  double R = sqrt(R2); 
  return R < _config.vertexDistanceCut;
}

bool NumuSelection::passMinLength(double length, bool stop_in_tpc) {
  if (!stop_in_tpc) return true;
  if (_config.minLength < 0) return true;
  return length > _config.minLength;
}

  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NumuSelection)

