#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
#include <ctime>
#include <cassert>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include "core/Event.hh"

/*

DESCRIPTION:

Point of this is just to get us a covariance matrix and plot it. No inverting it or anything.
So just need to:
- Get data and put it all into one single histogram
- Get cov mat, frac cov mat, corr mat; plot (with labels somehow... new input into EventSample? som string with description of sample?)
- Draw event dist with (alpha != 0) pink boxes that show total, systematic and statistical uncertainties

Inputs we'll need from config file:
- Binning (standard is proposal, non-standard can be given)
- Additional cuts on other vars (how to format these?)
- What ‘samples’ to use (why not just input only the ones we want to use rather than make this a config input?)
- What to plot (thought it was set...)
- POT to normalise each det to (would have to then have some sort of indication of what det each sample corresponds         to and this would have to be in EventSample)
- Which weights to use (standard is all, also make two options – "flux" and "event" – that, predictably, only 
    use flux and event weights, respectively)
- (Tentative:) How many universes to create?
    
Extra inputs we'll need in EventSample:
- some sort of (std::string) label – very concise, I'll put it into the plot too
- which det it came from (std::string) from "sbnd", "uboone" or "icarus"
    (to know what to normalise it to and for label in plot)

*/


// Simple version of class we'll get, just so I can get started here...
class EventSample {
    
    public:
        
    EventSample(TTree* _tree, float ScaleFactor, std::string Det, std::string Desc) : tree(_tree), fScaleFactor(ScaleFactor), sDet(Det), sDesc(Desc) {
        
        // Detector
        if (Det == "sbnd") {
            sDet = "SBND";
        } else if (Det == "uboone") {
            sDet = "MicroBooNE";
        } else if (Det == "icarus") {
            sDet = "ICARUS";
        } else {
            std::cout << std::endl << "ERROR: " << Det << " not valid detector." << std::endl << std::endl;
            assert(false);
        }
            
        // Description
        if (Desc == "nu") {
            sDesc = "Neutrino";
        } else {
            std::cout << std::endl << "ERROR: Sample type" << Desc << " not supported." << std::endl << std::endl;
            assert(false);
        }
        
    }
    
    TTree* tree;            // Should contain events (counts) and weights we'll analyse
    float fScaleFactor;     // Some factor to scale this by, eg. exposure
    std::string sDet;       // What detector it comes from
    std::string sDesc;      // (Very concise) Description of sample
    
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


void testCov() {
    
    clock_t start = clock(); // n of ticks since start of program
    gROOT->ProcessLine(".L $SBN_LIB_DIR/libsbnanalysis_Event.so");
    
    
    //// Prepare (imitation) input
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // For now I'll assume all samples are neutrino counts and that fScaleFactor is the total scale factor (i.e.
    // contains what to normalise TO, not just what to normalise FROM). This way I can just use the current format
    // of EventSample. Later I'll add more features.
    
    std::vector <std::string> DETLIST = {"SBND", "MicroBooNE", "ICARUS"},
        detlist = {"sbnd", "uboone", "icarus"},
        desclist = {"nu", "nu", "nu"};
    std::vector <float> scalelist = {3.0958e18, 8.87435e19, 6.59165e18};
    std::vector <EventSample> samples;
    
    for (int d = 0; d < detlist.size(); d++) {
        
        std::string det = detlist[d], desc = desclist[d];
        float scale = scalelist[d];
        
        TFile f(((std::string)"/sbnd/data/users/gavarela/selection/new/output_" + det + (std::string)".root").c_str());
        
        samples.push_back(EventSample((TTree*)f.Get("sbnana"), scale, det, desc));
        
    }
    
    // Now we can start: ...
    
    
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
        for (int e = 0; e > sample.tree->GetEntries(); e++) {
            for (int t = 0; t > event->truth.size(); t++) {
                
                tree->GetEntry(e);
                
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
        if (sample.desc == "Neutrino") {
            
            if (sample.det == "SBND") { offset = 0;
            } else if (sample.det == "MicroBooNE") { offset = 2*num_nuE_bins;
            } else if (sample.det == "ICARUS") { offset = 3*num_nuE_bins; }
            
            done_base_nuEs++;
            
        }
        
        for (int h = 0; h < temp_count_hists.size(); h++) {
            for (int bin = 0; bin < num_nuE_bins; bin++) {
                count_hists[h]->SetBinContent(offset+b+1, temp_count_hists[h]->GetBinContent(b+1));
            }
        }
        
        // Extension: for when theres more types of samples than just 'nu'/'Neutrino', we should have each 'type'
        // go together and in SBND-MicroBooNE-ICARUS order.
        
    }
    
    
    //// Get covariances, fractional covariances and correlation coefficients
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Covariance and fractional covariance
    nbins = count_hists[0]->GetNbinsX();
    TH2D *cov = new TH2D("cov", "Covariance Matrix", nbins, 0, nbins, nbins, 0, nbins),
         *fcov = new TH2D("fcov", "Fractional Covariance Matrix", nbins, 0, nbins, nbins, 0, nbins);
    
    for (int i = 0; i < cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov->GetNbinsY(); j++) {
            
            double covij = 0
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
    
    cov->Write("/sbnd/data/users/gavarela/selection/new/cov_output/root/cov_output_cov.root");
    fcov->Write("/sbnd/data/users/gavarela/selection/new/cov_output/root/cov_output_fcov.root");
    
    TCanvas *canvas = new TCanvas();
    cov->Draw("colz"); canvas->Write("/sbnd/data/users/gavarela/selection/new/cov_output/png/cov_plot.png");
    fcov->Draw("colz"); canvas->Write("/sbnd/data/users/gavarela/selection/new/cov_output/png/fcov_plot.png");
    
    
    // Pearson Correlation Coefficients
    TH2D *corr = new TH2D("corr", "Correlation Matrix", nbins, 0, nbins, nbins, 0, nbins);
    for (int i = 0; i < cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov->GetNbinsY(); j++) {
            
            double corrij = cov->GetBinContent(i+1, j+1) / TMath::Sqrt(cov->GetBinContent(i+1, i+1) * cov->GetBinContent(j+1, j+1));
            corr->SetBinContent(i+1, j+1, corrij);
            
        }
    }
    
    corr->Write("/sbnd/data/users/gavarela/selection/new/cov_output/root/cov_output_corr.root");
    corr->Draw("colz"); canvas->Write("/sbnd/data/users/gavarela/selection/new/cov_output/png/corr_plot.png");
    
    
    //// Pink box plot
    //// ~~~~~~~~~~~~~
    
    /*
    // Get errors
    std::vector <double> stat_errors = stat_err(count_hists[0], cov),
                         syst_errors = syst_err(count_hists[0], cov),
                          tot_errors = stat_errors;
    for (int i = 0; i < stat_errors.size(); i++) {
        tot_errors[i] += syst_errors[i];
    }
    std::vector <std::vector <double> > errors = {stat_errors, syst_errors, tot_errors};
    
    // Make all hists (in vectors)
    std::vector <std::vector <TH1D*> > pink_box_hists;
    std::vector <std::string> errs = {"stat", "syst", "tot"}, ERRS = {"Statistical", "Systematic", "Total"},
        dets = {"sbnd", "uboone", "icarus"}, DETS = {"SBND", "MicroBooNE", "ICARUS"};
    
    for (int d = 0; d < dets.size(); d++) {
        
        pink_box_hists.push_back({});
        for (int e = 0; e < errs.size(); e++) {
            
            std::string name = ERRS[e] + " errors" + " "*(d==1) + "  "*(d==2);
            TH1D *temphist = new TH1D(name, "", )
            
        }
        
    }
    */  
    
    
    
    
    // What's left?
    // 1. DONE: Calculate covariances, frac cov and corr (in TH2D*s)
    // 2. Plot of pink boxes (TCanvas of 3 TH1Ds)
    // 3. Save those in a tree in a .root file
                          
    
    
    
    
    
}

int main(int argc, char* argv[]) {
    
    testCov();
    
    return 0;
    
}