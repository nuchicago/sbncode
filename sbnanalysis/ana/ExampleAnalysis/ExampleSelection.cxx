#include <iostream>
#include <vector>
#include <string>
#include <TMath.h>
#include <TH2D.h>
#include <THStack.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <json/json.h>
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "ExampleSelection.h"
#include "ExampleTools.h"

namespace ana {
namespace ExampleAnalysis {

/* Function to check if vector is within 10cm of SBND active volume bounadries */
bool in_SBND(std::vector <double> coords) {

    return (coords[0] > -174.15)*(coords[0] < -27.65)*(coords[1] > -175)*(coords[1] < 175)*(coords[2] > 30)*(coords[2] < 450) + (coords[0] > 27.65)*(coords[0] < 174.15)*(coords[1] > -175)*(coords[1] < 175)*(coords[2] > 30)*(coords[2] < 450);

}

bool in_MicroBooNE(std::vector <double> coords) {
    
    return (coords[0] > 23.45)*(coords[0] < 229.8)*(coords[1] > -90.53)*(coords[1] < 92.47)*(coords[2] > 30.1)*(coords[2] < 986.9);
    
}

bool in_ICARUS(std::vector <double> coords) {
    
    return (coords[0] > -339.49)*(coords[0] < -241.29)*(coords[1] > -158.41)*(coords[1] < 118.41)*(coords[2] > -884.950652)*(coords[2] < 854.950652) + (coords[0] > -191.14)*(coords[0] < -92.94)*(coords[1] > -163.41)*(coords[1] < 133.41)*(coords[2] > -899.950652)*(coords[2] < 869.950652) + (coords[0] > -42.94)*(coords[0] < -191.14)*(coords[1] > -158.41)*(coords[1] < 118.41)*(coords[2] > -884.950652)*(coords[2] < 854.950652) + (coords[0] > 241.29)*(coords[0] < 339.49)*(coords[1] > -163.41)*(coords[1] < 133.41)*(coords[2] > -899.950652)*(coords[2] < 869.950652);
    
}

bool in_detector(std::vector <double> coords, std::string det) {
    
    if (det == (std::string)"SBND") {
        return in_SBND(coords);
    } else if (det == (std::string)"MicroBooNE") {
        return in_MicroBooNE(coords);
    } else if (det == (std::string)"ICARUS") {
        return in_ICARUS(coords);
    } else {
        std::cout << "Not a valid detector!" << std::endl;
        return 0;
    }
    
}

/* Function that gets length travelled inside the fiducial volume (assuming the particle's trajectory starts there) */
std::vector <double> get_lens(auto mcpart, std::string det) {
    
    /* Loop over positions until one is outside of the detector. Add up the distance travelled as we go. 
       Check if the particle leaves detector. If so then retun in second position of return vector. Else, 
       return in first position. */
    
    int retind = 0;                     // Index of return vector to return length in
    double len = 0;
    std::vector <double> oldcoords;
    
    for (int i = 0; i < mcpart->NumberTrajectoryPoints(); i++) {
        
        TLorentzVector event = mcpart->Position(i);
        std::vector <double> newcoords = {event.X(), event.Y(), event.Z()};
        if (i == 0) { oldcoords = newcoords; }
        
        if (!in_detector(newcoords, det)) {
            retind = 1;
            break;
        }
        
        len += TMath::Sqrt(TMath::Power(newcoords[0] - oldcoords[0], 2) +
                           TMath::Power(newcoords[1] - oldcoords[1], 2) +
                           TMath::Power(newcoords[2] - oldcoords[2], 2));
        
        oldcoords = newcoords;
        
    }
    
    std::vector <double> ret_vec(2);
    ret_vec[retind] = len;
    
    return ret_vec;
    
}

ExampleSelection::ExampleSelection() : SelectionBase(), /*fNuCount(0),*/ fEventCounter(0) {}

void ExampleSelection::Initialize(Json::Value* config) {
    
    /* Initialise variables */
    
    det = "ICARUS";
    
    n_unis = 100;
    Double_t nuE_bins[] = { 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5, 3 }; // { 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.25, 1.5, 2, 2.5, 3 };
    Int_t num_nuE_bins = sizeof(nuE_bins)/sizeof(Double_t) - 1; // 19;
    //min_nuE = 0.2; max_nuE = 3;
    //nuE_interval = 0.1;
    //num_nuE_bins = (int)((max_nuE-min_nuE)/nuE_interval);
    
    nucount = 0; nuindet = 0;
    mu_passes = 0; pi_passes = 0; double_pass = 0;
    
    mn = 0.93956541; mp = 0.93827208; ml = 0.10565837; Eb = 0.030;
    
    unis1 = {};
    for (int u = 0; u < n_unis+1; u++) { unis1.push_back(" (" + std::to_string(u) + ")"); }
    unis1[0] = ""; unis1[5] = " ";
    
    unis2 = {};
    for (int u = 0; u < n_unis+1; u++) { unis2.push_back("Alt. Uni. " + std::to_string(u)); }
    unis2[0] = "Base Uni."; unis2[5] = "Alt. Uni. 5";
    
    /* Initialise Histograms */
    
    // Pre-cut µ and π track lengths
    len_stack = new THStack("len_stack", (det + (std::string)": All Pions and Muons (Pre-Cut); Track Length (cm); Counts").c_str());
    len_hists.push_back(new TH1D("CC -> #mu^{#pm}  ", "", 50, 0, 150));
    len_hists[0]->SetFillColor(38); len_hists[0]->SetLineColor(38);
    len_hists.push_back(new TH1D("NC -> #pi^{#pm}   ", "", 50, 0, 150));
    len_hists[1]->SetFillColor(30); len_hists[1]->SetLineColor(30);
    
    // Neutrino energies. Vector of vector of histograms; two for each universe
    for (int u = 0; u < n_unis+1; u++) {
        
        nuE_hists.push_back({
            new TH1D(((std::string)"CC -> #mu^{#pm}" + unis1[u]).c_str(), "", num_nuE_bins, nuE_bins), 
            new TH1D(((std::string)"NC -> #pi^{#pm}" + unis1[u]).c_str(), "", num_nuE_bins, nuE_bins)
        });
        
        nuE_hists[u][0]->SetFillColor(38); nuE_hists[u][0]->SetLineColor(38);
        nuE_hists[u][1]->SetFillColor(30); nuE_hists[u][1]->SetLineColor(30);
        
    }
    
    /* Load configuration parameters */
    fTruthTag = { "generator" };
    fPartTag = { "largeant" };
    fEventWgtTag = { "genieeventweight" };
    fFluxWgtTag = { "fluxeventweight" };

    if (config) {
        fTruthTag = { (*config)["ExampleAnalysis"].get("MCTruthTag", "generator").asString() };
        fPartTag = { (*config)["ExampleAnalysis"].get("MCParticleTag", "largeant").asString() };
        fEventWgtTag = { (*config)["ExampleAnalysis"].get("MCEventWeightTag", "genieeventweight").asString() };
        fFluxWgtTag = { (*config)["ExampleAnalysis"].get("MCEventWeightTag", "fluxeventweight").asString() };
    }
    
    /* Add custom branches */
    // AddBranch("nucount", &fNuCount);
    // AddBranch("myvar", &fMyVar);
    
}

bool ExampleSelection::ProcessEvent(gallery::Event& ev) {

    if (fEventCounter%50 == 0) {
        std::cout << "Event " << fEventCounter << std::endl;
    }
    fEventCounter++;

    /* Grab a data product from the event */
    auto const& mctruths = *ev.getValidHandle<std::vector<simb::MCTruth>>(fTruthTag);
    auto const& mcevweights = *ev.getValidHandle<std::vector<evwgh::MCEventWeight>>(fEventWgtTag);
    auto const& mcfluxweights = *ev.getValidHandle<std::vector<evwgh::MCEventWeight>>(fFluxWgtTag);

    auto const& mctruth_handle = ev.getValidHandle<std::vector<simb::MCTruth>>(fTruthTag);
    auto const& mcparticle_handle = ev.getValidHandle<std::vector<simb::MCParticle>>(fPartTag);
    
    /* Associations b/w MCTruths (neutrinos) and MCParticles */
    const art::FindManyP <simb::MCParticle, sim::GeneratedParticleInfo> find_many_mcparticle(mctruth_handle, ev, fPartTag);

    /* Fill in the custom branches */
    //fNuCount = mctruths.size();         // Number of neutrinos in this event
    
    /* Iterate through the neutrinos */
    for (size_t i=0; i<mctruths.size(); i++) {
        
        nucount++; nuindet++;
        
        /* Definitions */
        auto const& mctruth = mctruths.at(i);
        auto const& mcevweight = mcevweights.at(i);
        auto const& mcfluxweight = mcfluxweights.at(i);
        
        /* Skip if vertex is outside fiducial volume */
        std::vector <double> vertex = {mctruth.GetNeutrino().Nu().Vx(), mctruth.GetNeutrino().Nu().Vy(), mctruth.GetNeutrino().Nu().Vz()};
        if (!in_detector(vertex, det)) { nuindet--; continue; }
        
        /* Loop over associated MCParticles to see if this neutrino passes the cuts */
        std::vector <art::Ptr <simb::MCParticle> > const& mcparticle = find_many_mcparticle.at(i);
        
        std::vector <int> pass_cut(2);
        double nuE;
        
        for (auto const& mcpart : mcparticle) {
            
            // Look only at muons and pions
            if ((mcpart->PdgCode()*mcpart->PdgCode() == 13*13 && nu.CCNC() == simb::kCC) || 
                (mcpart->PdgCode()*mcpart->PdgCode() == 211*211 && nu.CCNC() == simb::kNC)) {
                
                // Cut off those that start outside fiducial volume
                std::vector <double> startcoords = {mcpart->Vx(), mcpart->Vy(), mcpart->Vz()};
                if (!in_detector(startcoords, det)) { continue; }
                
                // Cut off those whose total track length is less than 50cm
                if (mcpart->Trajectory().TotalLength() < 50) { continue; }
                
                // Muon (0) or pion (1)?
                int part_ind = 0;
                if (mcpart->PdgCode()*mcpart->PdgCode() == 211*211) { part_ind = 1; }
                
                // Store lens
                std::vector <double> lens_forcuts = get_lens(mcpart, det);
                len_hists[part_ind]->Fill(lens_forcuts[0] + lens_forcuts[1]);
                
                // Pass cuts?
                if (lens_forcuts[0] > 50 || lens_forcuts[1] > 100) {
                    
                    pass_cut[part_ind] += 1;
                    
                    // Get neutrino energy
                    double px = mcpart->Px(), py = mcpart->Py(), pz = mcpart->Pz(), E = mcpart->E(),
                        p = TMath::Sqrt(px*px + py*py + pz*pz), pxy = TMath::Sqrt(px*px + py*py),
                        theta = TMath::ATan(pxy/pz);
                    nuE = 0.5 * (mp*mp - (mn - Eb)*(mn - Eb) - ml*ml + 2*(mn - Eb)*E) / ((mn - Eb) - E + p*TMath::Cos(theta));
                    
                    // nuE = mctruth.GetNeutrino().Nu().E();
                    
                }
                
            }
                
        } // End loop over MCParticles
        
        // pass_cut = {1, 0}; nuE = mctruth.GetNeutrino().Nu().E();
        
        /* If passed cut, get (weighted) count(s) for this event */
        if (pass_cut[0] + pass_cut[1] == 1) {
            
            // Which particle passed the cut?
            int pi_pass = 0;
            if (pass_cut[0] == 0) { pi_pass = 1; }
            mu_passes += (1-pi_pass); pi_passes += pi_pass;
            
            // Base counts (neutrino energies)
            nuE_hists[0][pi_pass]->Fill(nuE);
            
            // Weights for the diff 'universes'. Can increase n_unis w/out problem till n_unis = 1000 (size of longest weight param vector)
            double weight;
            for (int u = 0; u < n_unis; u++) {
                
                weight = 1;
                int wind;
                
                /*
                // Loop over event weight parameters
                for (auto it : mcevweight.fWeight) {
                    if (it.first == "genie_IntraNukeNel_Genie" || it.first == "genie_IntraNukePIinel_Genie") { continue; }
                    wind = u;
                    while (wind > it.second.size()-1) { wind -= it.second.size(); }
                    weight *= it.second.at(wind);
                }
                */
                
                // Loop over flux weight parameters
                for (auto it : mcfluxweight.fWeight) {
                    wind = u;
                    while (wind > it.second.size()-1) { wind -= it.second.size(); }
                    weight *= it.second.at(wind);
                }
                
                nuE_hists[u+1][pi_pass]->Fill(nuE, weight);
                
            }
                        
        } else if (pass_cut[0] + pass_cut[1] > 1) {
            
            std::cout << "Multiple pass!" << std::endl;
            double_pass++;
        
        }
        
    } // End iteration over neutrinos

    return true;

}

void ExampleSelection::Finalize() {

    /* Set up output of files (I think...) */
    fOutputFile->cd();
    /* to in fact output something: thing_to_output->Write(); */
    
    std::cout << std::endl << "Out of " << nucount << " neutrinos, we had " << nuindet << " react inside the fiducial volume of the detector." << std::endl << 
    "   Of these, " << mu_passes + pi_passes << " passed the cut." << std::endl;
    
    std::cout << std::endl << "Got " << mu_passes << " muons and " << pi_passes << " pions passing the cut." << std::endl;
    std::cout << "And got " << double_pass << " events with more than one particle passing the cut." << std::endl;
    
    std::cout << "This is equivalent to a proportion " << (mu_passes + pi_passes)/(double)nucount << " of total neutrinos." << std::endl << 
                 "                  and a proportion " << (mu_passes + pi_passes)/(double)nuindet << " of those that reacted in the fiducial volume." << std::endl << std::endl;
    
    /* Pre-Cut Track Lengths Histogram */
    len_stack->Add(len_hists[0]);
    len_stack->Add(len_hists[1]);
    len_stack->Write();
    
    /* Fill out stacked histograms and output them */
    std::vector <double> nuE_bins = { 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5, 3 } /* { 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.25, 1.5, 2, 2.5, 3 } */, nuE_binwidths;
    int num_nuE_bins = nuE_bins.size()-1;
    for (int i = 0; i < num_nuE_bins; i++) {
        nuE_binwidths.push_back(nuE_bins[i+1] - nuE_bins[i]);
    }
    
    for (int u = 0; u < n_unis+1; u++) {
        
        count_stack = new THStack(((std::string)"count_stack_" + std::to_string(u)).c_str(),
                                  (det + (std::string)": " + unis2[u] + (std::string)" Neutrino Counts; Energy (GeV); Counts (GeV^{-1})").c_str());
        
        // Rescale bins and add to stack (pi first)
        for (int p = 1; p > -1; p--) {
            for (int i = 0; i < num_nuE_bins; i++) {
                nuE_hists[u][p]->SetBinContent(i+1, (nuE_hists[u][p]->GetBinContent(i+1)/nuE_binwidths[i]));
            }
            count_stack->Add(nuE_hists[u][p]);
        }
        
        count_stack->Write();
        
    }
    
}

}// namespace ExampleAnalysis
}// namespace ana

/* This line must be included for all selections! */
DECLARE_SBN_PROCESSOR(ana::ExampleAnalysis::ExampleSelection)
