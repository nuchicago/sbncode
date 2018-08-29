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
#include <TCanvas.h>
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
    
// My own constructor that deals with the extra public members I added.
EventSample::EventSample(TFile* _file, TTree* _tree, float ScaleFactor, std::string Det, std::string Desc) : file(_file), tree(_tree), fScaleFactor(ScaleFactor), fDet(Det), fDesc(Desc) {

    // Detector
    std::map <std::string, std::string> det_to_DET = {{"sbnd", "SBND"}, {"uboone", "MicroBooNE"}, {"icarus", "ICARUS"}};

    if (det_to_DET.find(Det) != det_to_DET.end()) {
        fDet = det_to_DET[Det];
    } else {
        std::cout << std::endl << "ERROR: " << Det << " not valid detector." << std::endl << std::endl;
        assert(false);
    }

    // Description (can add more)
    std::map <std::string, std::string> desc_to_DESC = {{"numu", "#nu_{#mu}"}, {"nue", "#nu_{e}"}, {"pi", "#pi"}};

    if (desc_to_DESC.find(Desc) != desc_to_DESC.end()) {
        fDesc = desc_to_DESC[Desc];
    } else {
        std::cout << std::endl << "ERROR: Sample type " << Desc << " not supported." << std::endl << std::endl;
        assert(false);
    }

};

// Function that gets scale factors (weights) for different universes
std::vector <double> get_uni_weights(std::map <std::string, std::vector <double> > weights, int n_unis) {
    
    // Tentative format: universe u scale factor is the product of the u-th entries on each vector 
    // inside the map. For vectors with less than u entries, use the (u - vec_size)-th entry
    
    std::vector <double> uweights;
    
    for (int u = 0; u < n_unis; u++) {
                
        double weight = 1;
        int wind;
        
        for (auto it : weights) {
            wind = u;
            while (wind >= it.second.size()) { wind -= it.second.size(); }
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
    
// Function that gets plotting order (SBND-ICARUS-MicroBooNE and, within each det, nue-numu-rest)
std::vector <std::string> get_plot_order(std::vector<EventSample> samples) {
    
    // Create map of samples: detector -> description
    std::map <std::string, std::vector <std::string> > sample_descs;
    std::vector <std::string> all_descs, all_dets;
    for (int s = 0; s < samples.size(); s++) {
        
        // If this is a new sample detector, make new key and add detector 
        if (sample_descs.find(samples[s].fDet) == sample_descs.end()) {
        
            sample_descs.insert({samples[s].fDet, {samples[s].fDesc}});
        
        // Else, add detector to existing key
        } else {
            
            sample_descs[samples[s].fDet].push_back(samples[s].fDesc);
        
        }
        
        // If new sample description, add to vector
        int newdesc = 1;
        for (int d = 0; d < all_descs.size(); d++) {
            if (all_descs[d] == samples[s].fDesc) { newdesc = 0; }
        }
        if (newdesc == 1) { all_descs.push_back(samples[s].fDesc); }
        
        // If new sample detector, add to vector
        int newdet = 1;
        for (int d = 0; d < all_dets.size(); d++) {
            if (all_dets[d] == samples[s].fDet) { newdet = 0; }
        }
        if (newdet == 1) { all_dets.push_back(samples[s].fDet); }
        
    }
    
    // Get order descriptions will be plotted in (nue, numu, rest)
    std::vector <std::string> desc_order;
    if (std::find(all_descs.begin(), all_descs.end(), "#nu_{e}") != all_descs.end()) {
        desc_order.push_back("#nu_{e}");
    }
    if (std::find(all_descs.begin(), all_descs.end(), "#nu_{#mu}") != all_descs.end()) {
        desc_order.push_back("#nu_{#mu}");
    }
    for (int d = 0; d < all_descs.size(); d++) {
        if (all_descs[d] != "#nu_{e}" && all_descs[d] != "#nu_{#mu}") {
            desc_order.push_back(all_descs[d]);
        }
    }
    
    // Get order detectors will be plotted in
    std::vector <std::string> det_order;
    if (std::find(all_dets.begin(), all_dets.end(), "SBND") != all_dets.end()) {
        det_order.push_back("SBND");
    }
    if (std::find(all_dets.begin(), all_dets.end(), "ICARUS") != all_dets.end()) {
        det_order.push_back("ICARUS");
    }
    if (std::find(all_dets.begin(), all_dets.end(), "MicroBooNE") != all_dets.end()) {
        det_order.push_back("MicroBooNE");
    }
    
    // Get actual order all samples will be plotted in 
    // (this method allows for uneven samples (eg. SBND-nue and ICARUS-numu and no others))
    std::vector <std::string> plot_order;
    for (std::string det : all_dets) {
        for (std::string desc : all_descs) {
            
            if (std::find(sample_descs[det].begin(), sample_descs[det].end(), desc) != sample_descs[det].end()) {
                plot_order.push_back(det+"_"+desc);
            }
            
        }
    }
    
    return plot_order;
    
}

Covariance::Covariance(std::vector<EventSample> samples) {
    
    //// Get counts on each (base and alternative) universe
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Params (may in future come from config file)
    int n_alt_unis = 100;
    
    Double_t nue_bins[] = { 0.2, 0.35, 0.5, 0.65, 0.8, 0.95, 1.1, 1.3, 1.5, 1.75, 2, 3 };
    Double_t numu_bins[] = { 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5, 3 }; 
        // Proposal uses below but I don't have enough events for it to look good:
        // { 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.25, 1.5, 2, 2.5, 3 };
    
    std::map <std::string, float> scale_tgt = {{"SBND", 6.6e20}, {"MicroBooNE", 13.2e21}, {"ICARUS", 6.6e20}};
    
    // Some stuff related to binning and plotting
    Int_t num_nue_bins = sizeof(nue_bins)/sizeof(Double_t) - 1,
          num_numu_bins = sizeof(numu_bins)/sizeof(Double_t) - 1;
        
        // Get plotting order of samples
    std::vector <std::string> plot_order = get_plot_order(samples);
    
        // Number of bins needed in big histogram (for covariance)
    int num_bins = 0;
    std::vector <int> offset;
    for (int o = 0; o < plot_order.size(); o++) {
        
        if (plot_order[o].find("#nu_{e}") != std::string::npos) {
            offset.push_back(num_bins); num_bins += num_nue_bins;
        } else if (plot_order[o].find("#nu_{#mu}") != std::string::npos) {
            offset.push_back(num_bins); num_bins += num_numu_bins;
        } else {
            offset.push_back(num_bins); num_bins += num_numu_bins;
        }
        
        // Note: for now am assuming that nue and numu get plotted as above and all the rest get plotted as numu does
        //       if this changes, need to change the else right above here and the first else in sample loop below.
        
    }
    offset.push_back(num_bins); // for consistency of the vector[ind] call at the end, setting labels to covs
    
    // Large (meaningless x-axis) histograms for cov
    std::vector <TH1D*> count_hists = {new TH1D("base", "Base Uni. Counts; Bin; Counts", num_bins, 0, num_bins)};
    
    for (int u = 0; u < n_alt_unis; u++) {
        
        std::string name = "alt" + std::to_string(u+1), 
                    title = "Alt. Uni. " + std::to_string(u+1) + " Counts; Bin; Counts";
        
        count_hists.push_back(new TH1D(name.c_str(), title.c_str(), num_bins, 0, num_bins));
    
    }
    
    // Get counts
    std::cout << std::endl << "Getting counts for each sample..." << std::endl;
    
    TCanvas *base_canvas = new TCanvas("base_canvas", "Base Universe", 950, 345), 
            *alt5_canvas = new TCanvas("alt5_canvas", "Fifth Alternative Universe", 950, 345);
    base_canvas->Divide(3, 1); alt5_canvas->Divide(3, 1);
    
    std::vector <double> energies;
    
    for (int o = 0; o < plot_order.size(); o++) {
        
        // Get relevant sample
        EventSample sample = samples[o];
        
        // Initialise temp hists to store counts
        Double_t *bins; Int_t nbins;
        if (sample.fDesc == "#nu_{e}") {
            bins = nue_bins; nbins = num_nue_bins;
        } else if (sample.fDesc == "#nu_{#mu}") {
            bins = numu_bins; nbins = num_numu_bins;
        } else {
            bins = numu_bins; nbins = num_numu_bins;
        }
        
        std::string title = sample.fDet+"; Energy (GeV); Counts";
        std::vector <TH1D*> temp_count_hists = {new TH1D("tempbase", title.c_str(), nbins, bins)};
        
        for (int u = 0; u < n_alt_unis; u++) {
            
            std::string name = "tempalt" + std::to_string(u+1);
            temp_count_hists.push_back(new TH1D(name.c_str(), title.c_str(), nbins, bins));
            
        }
        
        // Tree stuff
        Event *event = new Event;
        sample.tree->SetBranchAddress("events", &event);
        
        // Loop over neutrinos (events)
        int nucount = 0;
        for (int e = 0; e < sample.tree->GetEntries(); e++) {
            
            sample.tree->GetEntry(e);
            
            for (int t = 0; t < event->truth.size(); t++) {
                
                nucount++;
                
                // Add (reconstructed) energy to base universe histogram
                double nuE = event->reco[t].neutrino.energy;
                temp_count_hists[0]->Fill(nuE, 1 + 1*(sample.fDet=="ICARUS"));
                
                // Get weights for each universe
                std::vector <double> uweights = get_uni_weights(event->truth[t].weights, n_alt_unis);
                for (int u = 0; u < uweights.size(); u++) {
                    temp_count_hists[u+1]->Fill(nuE, (1 + 1*(sample.fDet=="ICARUS"))*uweights[u]);
                }
                
            }
        }
        std::cout << std::endl << "For sample: " << sample.fDet << ", " << sample.fDesc << ", there were " << nucount << " neutrinos." << std::endl;
        
        // Rescale each bin's counts to /GeV and to desired POT
        for (int u = 0; u < temp_count_hists.size(); u++) {
            
            for (int b = 0; b < temp_count_hists[u]->GetNbinsX(); b++) {
                
                double binwidth = temp_count_hists[u]->GetBinWidth(b+1), 
                     bincontent = temp_count_hists[u]->GetBinContent(b+1);
                temp_count_hists[u]->SetBinContent(b+1, bincontent / binwidth * scale_tgt[sample.fDet] / sample.fScaleFactor);
                
            }
            
        }
        
        // Add base and alt5 hists to a canvas with three pads
        if (sample.fDesc == "#nu_{#mu}") {
            if (sample.fDet == "SBND") { 
                TH1D *SBND_0 = temp_count_hists[0], *SBND_5 = temp_count_hists[5];
                base_canvas->cd(1); SBND_0->SetStats(kFALSE); SBND_0->Draw("hist");
                alt5_canvas->cd(1); SBND_5->SetStats(kFALSE); SBND_5->Draw("hist");
            } else if (sample.fDet == "MicroBooNE") { 
                TH1D *MicroBooNE_0 = temp_count_hists[0], *MicroBooNE_5 = temp_count_hists[5];
                base_canvas->cd(2); MicroBooNE_0->SetStats(kFALSE); MicroBooNE_0->Draw("hist");
                alt5_canvas->cd(2); MicroBooNE_5->SetStats(kFALSE); MicroBooNE_5->Draw("hist");
            } if (sample.fDet == "ICARUS") { 
                TH1D *ICARUS_0 = temp_count_hists[0], *ICARUS_5 = temp_count_hists[5];
                base_canvas->cd(3); ICARUS_0->SetStats(kFALSE); ICARUS_0->Draw("hist");
                alt5_canvas->cd(3); ICARUS_5->SetStats(kFALSE); ICARUS_5->Draw("hist");
            }
        }
        
        // Pass onto the big histograms and get energies
        for (int h = 0; h < temp_count_hists.size(); h++) {
            
            for (int bin = 0; bin < nbins; bin++) {
                
                if (h == 0) { energies.push_back(temp_count_hists[h]->GetBinCenter(bin+1)); }
                count_hists[h]->SetBinContent(1+offset[o]+bin, temp_count_hists[h]->GetBinContent(bin+1));
                
            }
            
            count_hists[h]->GetXaxis()->SetBinLabel(offset[o]+nbins/2, (sample.fDet+" "+sample.fDesc).c_str());
            count_hists[h]->GetXaxis()->LabelsOption("h");
            // temp_count_hists[h]->Delete();
            
        }
        
    }
    
    std::cout << std::endl << "Created energies in Covariance. Length " << energies.size() << ". Here it is: ";
    for (int i = 0; i < energies.size(); i++) {
        std::cout << energies[i] << ", ";
    }
    std::endl;
    
    TCanvas *tempcanvas = new TCanvas();
    count_hists[0]->Draw(); count_hists[0]->SetStats(kFALSE);
    tempcanvas->SaveAs("/sbnd/data/users/gavarela/selection/new/cov_output/base.pdf");
    count_hists[5]->Draw(); count_hists[5]->SetStats(kFALSE);
    tempcanvas->SaveAs("/sbnd/data/users/gavarela/selection/new/cov_output/alt5.pdf");
    
    base_canvas->SaveAs("/sbnd/data/users/gavarela/selection/new/cov_output/nicebase.pdf");
    alt5_canvas->SaveAs("/sbnd/data/users/gavarela/selection/new/cov_output/nicealt5.pdf");
    
    
    //// Get covariances, fractional covariances and correlation coefficients
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Getting covs..." << std::endl;
    
    // Covariance and fractional covariance
    TH2D *cov = new TH2D("cov", "Covariance Matrix", num_bins, 0, num_bins, num_bins, 0, num_bins),
         *fcov = new TH2D("fcov", "Fractional Covariance Matrix", num_bins, 0, num_bins, num_bins, 0, num_bins);
    
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
    
    // Pearson Correlation Coefficients
    TH2D *corr = new TH2D("corr", "Correlation Matrix", num_bins, 0, num_bins, num_bins, 0, num_bins);
    for (int i = 0; i < cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov->GetNbinsY(); j++) {
            
            double corrij = cov->GetBinContent(i+1, j+1) / TMath::Sqrt(cov->GetBinContent(i+1, i+1) * cov->GetBinContent(j+1, j+1));
            corr->SetBinContent(i+1, j+1, corrij);
            
        }
    }
    
    // Add bin labels
    for (int o = 0; o < plot_order.size(); o++) {
        
        // Get label and position
        std::string label = plot_order[o].replace(plot_order[o].find("_"), 1, " ");
        int pos = 1 + offset[o] + (offset[o+1] - offset[o])/2;
        
        // Set label
        cov->GetXaxis()->SetBinLabel(pos, label.c_str());
        cov->GetYaxis()->SetBinLabel(pos, label.c_str());
        
        fcov->GetXaxis()->SetBinLabel(pos, label.c_str());
        fcov->GetYaxis()->SetBinLabel(pos, label.c_str());
        
        corr->GetXaxis()->SetBinLabel(pos, label.c_str());
        corr->GetYaxis()->SetBinLabel(pos, label.c_str());
        
    }
    
    cov->GetXaxis()->LabelsOption("h"); cov->GetYaxis()->LabelsOption("v");
    fcov->GetXaxis()->LabelsOption("h"); fcov->GetYaxis()->LabelsOption("v");
    corr->GetXaxis()->LabelsOption("h"); corr->GetYaxis()->LabelsOption("v");
    
    std::cout << std::endl << "Done with that." << std::endl;
    
    
    //// Output relevant objects
    //// ~~~~~~~~~~~~~~~~~~~~~~~
    
    CV_counts = count_hists[0];
    energies = energies;
    
    sample_order = plot_order;
    sample_bins = offset;
    
    covmat = cov;
    fcovmat = fcov;
    corrmat = corr;
    
    
    
}

}   // namespace SBNOsc
}   // namespace ana