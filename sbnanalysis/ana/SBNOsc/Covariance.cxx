#include <string>
#include <vector>
#include <TChain.h>
#include <TTree.h>
#include "Covariance.h"

// My includes
//#include <vector>
//#include <string>
#include <map>
#include <cmath>
#include <iostream>
#include <ctime>
#include <cassert>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
// #include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <core/Event.hh>

namespace ana {
namespace SBNOsc {

EventSample::EventSample(std::vector<std::string> filenames, float scaleFactor) : fScaleFactor(scaleFactor) {

    TChain* t = new TChain("sbnana");

    for (size_t i=0; i<filenames.size(); i++) {
        t->Add(filenames[i].c_str());
    }

    tree = dynamic_cast<TTree*>(t);

};

// Function that gets scale factor for different universes
std::vector <double> get_uni_weights(std::map <std::string, std::vector <double> > weights, int n_unis) {
    
    // Tentative format: universe u scale factor is the product of the u-th entries on each vector 
    // inside the map. For vectors with less than u entries, use the (u - vec_size)-th entry
    
    std::vector <double> uweights;
    
    for (int u = 0; u < n_unis; u++) {
                
        double weight = 1;
        int wind;
        
        for (auto it : weights) {
            wind = u;
            while (wind > it.second.size()-1) { wind -= it.second.size(); }
            weight *= it.second.at(wind);
        }
        
        uweights.push_back(weight);
        
    }
    
    return uweights;
    
    // Probem: this method may impose a limit on number of possible universes that is less than the 
    // absolute limit of (if p is number of parameters (weights.size()) and s the number of weights total
    // in the map (sum weights[i].size() for all i)): pC1 + pC2 + ... pCs
    
    // Solution: use random number generators (with set seed to make it reproducible on an 
    // event-by-event basis).
    
}

// Function that gets statistical errors
std::vector <double> stat_err(TH1D *counts) {
    
    std::vector <double> errors;
    for (int b = 0; b < counts->GetNbinsX(); b++) {
        errors.push_back(TMath::Sqrt(counts->GetBinContent(b+1)));
    }
    
    return errors;
    
}

// Function that gets systematic errors
std::vector <double> syst_err(TH1D *counts, TH2D *cov) {
    
    std::vector <double> errors;
    for (int b = 0; b < counts->GetNbinsX(); b++) {
        errors.push_back(TMath::Sqrt(cov->GetBinContent(b+1, b+1)));
    }
    
    return errors;
}

Covariance::Covariance(std::vector<EventSample> samples) {
    
    //// Get counts on each (base and alternative) universe
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Params (may in future come from config file)
    int n_alt_unis = 100;
    
    Double_t nuE_bins[] = { 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5, 3 }; 
    // Double_t nuE_bins[] = { 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.25, 1.5, 2, 2.5, 3 };
    Int_t num_nuE_bins = sizeof(nuE_bins)/sizeof(Double_t) - 1;
    
    std::vector <float> scale_tgt = {6.6e20, 13.2e21, 6.6e20};
    
    // Large (meaningless x-axis) histograms for cov
    std::vector <TH1D*> count_hists = {new TH1D("base", "Base Uni. Counts; Bin; Counts", 
                                                num_nuE_bins*samples.size(), 0, num_nuE_bins*samples.size())};
    for (int u = 0; u < n_alt_unis; u++) {
        std::string name = "alt" + std::to_string(u+1),
                title = "Alt. Uni. " + std::to_string(u+1) + " Counts; Bin; Counts";
            count_hists.push_back(new TH1D(name.c_str(), title.c_str(), 
                                           num_nuE_bins*samples.size(), 0, num_nuE_bins*samples.size()));
    }
    
    // Get counts
    int done_base_nuEs = 0;
    for (int s = 0; s < samples.size(); s++) {
        
        EventSample sample = samples[s];
        
        // Hists to store counts
        std::string basename = "Base" + sample.sDet + sample.sDesc,
            basetitle = sample.sDet + " " + sample.sDesc + " Counts; Energy (GeV); Counts";
        std::vector <TH1D*> temp_count_hists = {new TH1D("tempbase", " Counts; Energy (GeV); Counts", num_nuE_bins, nuE_bins)};
        for (int u = 0; u < n_alt_unis; u++) {
            std::string name = "tempalt" + std::to_string(u+1),
                title = "Alt. Uni. " + std::to_string(u+1) + " Counts; Energy (GeV); Counts";
            temp_count_hists.push_back(new TH1D(name.c_str(), title.c_str(), num_nuE_bins, nuE_bins));
        }
        
        // Tree stuff
        Event *event = new Event;
        sample.tree->SetBranchAddress("events", &event);
        
        // Loop over neutrinos (events)
        std::cout << std::endl << "For " << sample.sDet << ", sample.tree->GetEntries() = " << sample.tree->GetEntries() << std::endl;
        int nucount = 0;
        for (int e = 0; e < sample.tree->GetEntries(); e++) {
            
            sample.tree->GetEntry(e);
            if (event->truth.size() == 0) { std::cout << "EMPTY truth IN EVENT" << e << "!!!" << std::endl; }
            
            for (int t = 0; t < event->truth.size(); t++) {
                
                nucount++;
                
                // Add (reconstructed) energy to base universe histogram
                double nuE = event->reco[t].neutrino.energy;
                temp_count_hists[0]->Fill(nuE);
                
                // Get weights for each universe
                std::vector <double> uweights = get_uni_weights(event->truth[t].weights, n_alt_unis);
                for (int u = 0; u < uweights.size(); u++) {
                    temp_count_hists[u+1]->Fill(nuE*uweights[u]);
                }
                
            }
        }
        std::cout << std::endl << "For sample: " << sample.sDet << ", " << sample.sDesc << ", there were " << nucount << " neutrinos." << std::endl;
        
        // Rescale each bin's counts to /GeV
        for (int u = 0; u < temp_count_hists.size(); u++) {
            for (int b = 0; b < temp_count_hists[u]->GetNbinsX(); b++) {
                double bincontent = temp_count_hists[u]->GetBinContent(b+1),
                         binwidth = temp_count_hists[u]->GetBinWidth(b+1);
                temp_count_hists[u]->SetBinContent(b+1, bincontent/binwidth);
            }
        }
        
        // Put into (vector of) big histogram(s)
        int offset = num_nuE_bins*(3 + s - done_base_nuEs);
        if (sample.sDesc == "Neutrino") {
            
            if (sample.sDet == "SBND") { offset = 0;
            } else if (sample.sDet == "MicroBooNE") { offset = 2*num_nuE_bins;
            } else if (sample.sDet == "ICARUS") { offset = 3*num_nuE_bins; }
            
            done_base_nuEs++;
            
        }
        
        for (int h = 0; h < temp_count_hists.size(); h++) {
            for (int bin = 0; bin < num_nuE_bins; bin++) {
                count_hists[h]->SetBinContent(offset+bin+1, temp_count_hists[h]->GetBinContent(bin+1));
            }
        }
        
        // Extension: for when theres more types of samples than just 'nu'/'Neutrino', we should have each 'type'
        // go together and in SBND-MicroBooNE-ICARUS order.
        
    }
    
    
    //// Get covariances, fractional covariances and correlation coefficients
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Getting covs..." << std::endl;
    
    // Covariance and fractional covariance
    int nbins = count_hists[0]->GetNbinsX();
    TH2D *cov = new TH2D("cov", "Covariance Matrix", nbins, 0, nbins, nbins, 0, nbins),
         *fcov = new TH2D("fcov", "Fractional Covariance Matrix", nbins, 0, nbins, nbins, 0, nbins);
    
    for (int i = 0; i < cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov->GetNbinsY(); j++) {
            
            double covij = 0;
            for (int u = 0; u < n_alt_unis; u++) {
                covij += (count_hists[0]->GetBinContent(i+1) - count_hists[u]->GetBinContent(i+1)) * 
                         (count_hists[0]->GetBinContent(j+1) - count_hists[u]->GetBinContent(j+1));
            }
            covij /= n_alt_unis;
            cov->SetBinContent(i+1, j+1, covij);
            
            double fcovij = covij / (count_hists[0]->GetBinContent(i+1) * count_hists[0]->GetBinContent(j+1));
            fcov->SetBinContent(i+1, j+1, fcovij);
            
        }
    }
    
    covmat = cov; fcovmat = fcov;
    
    // Pearson Correlation Coefficients
    TH2D *corr = new TH2D("corr", "Correlation Matrix", nbins, 0, nbins, nbins, 0, nbins);
    for (int i = 0; i < cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov->GetNbinsY(); j++) {
            
            double corrij = cov->GetBinContent(i+1, j+1) / TMath::Sqrt(cov->GetBinContent(i+1, i+1) * cov->GetBinContent(j+1, j+1));
            corr->SetBinContent(i+1, j+1, corrij);
            
        }
    }
    
    corrmat = corr;
    
    std::cout << std::endl << "Done with that." << std::endl;
    
}

/*Covariance::Save(std::string directory) {
    
    cov->Write("/sbnd/data/users/gavarela/selection/new/cov_output/root/cov_output_cov.root");
    fcov->Write("/sbnd/data/users/gavarela/selection/new/cov_output/root/cov_output_fcov.root");
    corr->Write("/sbnd/data/users/gavarela/selection/new/cov_output/root/cov_output_corr.root");
    
    TCanvas *canvas = new TCanvas();
    cov->Draw("colz"); canvas->Write((directory + (std::string)"/cov_plot.png").c_str());
    fcov->Draw("colz"); canvas->Write((directory + (std::string)"/fcov_plot.png").c_str());
    corr->Draw("colz"); canvas->Write((directory + (std::string)"/corr_plot.png").c_str());
    
}*/

}   // namespace SBNOsc
}   // namespace ana