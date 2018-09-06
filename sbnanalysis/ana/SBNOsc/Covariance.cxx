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
EventSample::EventSample(TFile* _file, float ScaleFactor, std::string Det, std::string Desc) : file(_file), fScaleFactor(ScaleFactor), fDet(Det), fDesc(Desc) {

    // Tree
    tree = (TTree*) file->Get("sbnana");
    
    /*
    // Detector
    std::map <std::string, std::string> det_to_DET = {{"sbnd", "SBND"}, {"uboone", "MicroBooNE"}, {"icarus", "ICARUS"}};

    if (det_to_DET.find(Det) != det_to_DET.end()) {
        fDet = det_to_DET[Det];
    } else {
        std::cout << std::endl << "ERROR: " << Det << " not valid detector." << std::endl << std::endl;
        assert(false);
    }
    */

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
        if (std::find(all_descs.begin(), all_descs.end(), samples[s].fDesc) == all_descs.end()) {
            all_descs.push_back(samples[s].fDesc);
        }
        
        // If new sample detector, add to vector
        if (std::find(all_dets.begin(), all_dets.end(), samples[s].fDet) == all_dets.end()) {
            all_dets.push_back(samples[s].fDet);
        }
        
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

// Order dets in sample SBND-MicroBooNE-ICARUS
std::vector <std::string> get_dets_inorder(std::vector<EventSample> samples) {
    
    std::vector <std::string> dets, dets_inorder;
    
    for (EventSample sample : samples) {
        if (std::find(dets.begin(), dets.end(), sample.fDet) == dets.end()) {
            dets.push_back(sample.fDet);
        }
    }
    
    if (std::find(dets.begin(), dets.end(), "SBND") != dets.end()) {
        dets_inorder.push_back("SBND");
    }
    if (std::find(dets.begin(), dets.end(), "MicroBooNE") != dets.end()) {
        dets_inorder.push_back("MicroBooNE");
    }
    if (std::find(dets.begin(), dets.end(), "ICARUS") != dets.end()) {
        dets_inorder.push_back("ICARUS");
    }
    
    return dets_inorder;
    
}

Covariance::Covariance(std::vector<EventSample> samples, char *configFileName) {
    
    //// Getting parameters from config file
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // What's in our samples?
    std::vector <std::string> descs;
    for (EventSample sample : samples) {
        if (std::find(descs.begin(), descs.end(), sample.fDesc) == descs.end()) {
            descs.push_back(sample.fDesc);
        }
    }
    std::vector <std::string> dets = get_dets_inorder(samples);
    
    // Get configuration file parameters
    Json::Value* config = core::LoadConfig(configFileName);
    
    if (config != NULL) {
        
        // Weight and universe stuff
        fWeightKey = (*config)["Covariance"].get("WeightKey", "").asString();
        if (fWeightKey == "GetWeights") {
            fNumAltUnis = (*config)["Covariance"].get("NumAltUnis", -7).asInt();
        } else {
            Event *tempev = new Event;
            samples[0].tree->SetBranchAddress("events", &tempev);
            samples[0].tree->GetEntry(0);
            fNumAltUnis = tempev->truth[0].weights[fWeightKey].size();
                // Relies on all 'weights' branches being the same for all samples – should be true
        }
        
        // Type of energy
        fEnergyType = (*config)["Covariance"].get("EnergyType", "").asString();
        
        // Histogram bins
        std::vector <double> default_bins = { 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5, 3 };
        for (std::string desc : descs) {
            
            fBins.insert({desc, {}});
            
            if ((*config)["Covariance"]["BinDefs"].isMember(desc)) {
                for (auto binlim : (*config)["Covariance"]["BinDefs"][desc]) {
                    fBins[desc].push_back(binlim.asDouble());
                }
            } else {
                fBins[desc] = default_bins;
            }
            
        }
        
        // Exposure normalisation
        for (std::string det : dets) {
            fScaleTargets.insert({det, (*config)["Covariance"]["ScaleTargets"].get(det, -1).asFloat()});
        }
        
        
        // more documentation is at: https://open-source-parsers.github.io/jsoncpp-docs/doxygen/index.html#_example
    
    }
    
    std::cout << std::endl << "Got all parameters from the config file." << std::endl;
    
    
    //// Get counts on each (base and alternative) universe
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Some stuff related to binning and plotting
        // If nue appearance, set binning for all numus to same as nues
    int nue_appearance = (fBins.find("#nu_{e}") != fBins.end());
    if (nue_appearance == 1) {
        fBins["#nu_{#mu}"] = fBins["#nu_{e}"];
    }
    
        // Get plotting order of samples
    std::vector <std::string> plot_order = get_plot_order(samples);
    
        // Number of bins needed in big histogram (for covariance)
        // And bin 'boundaries' between each separate sample
    int num_bins = 0;
    std::vector <int> offset;
    for (int o = 0; o < plot_order.size(); o++) {
        
        offset.push_back(num_bins);
        
        std::string binkey = plot_order[o].substr(plot_order[o].find("_")+1, plot_order[o].length());
        num_bins += fBins[binkey].size() - 1;
        
    }
    offset.push_back(num_bins); // for later when setting labels to cov plots
    
    // Large (meaningless x-axis) histograms for cov
    std::vector <TH1D*> count_hists = {new TH1D("base", "Base Uni. Counts; Bin; Counts", num_bins, 0, num_bins)};
    // count_hists[0]->GetXaxis()->LabelsOption("h");
    
    for (int u = 0; u < fNumAltUnis; u++) {
        
        std::string name = "alt" + std::to_string(u+1), 
                    title = "Alt. Uni. " + std::to_string(u+1) + " Counts; Bin; Counts";
        
        count_hists.push_back(new TH1D(name.c_str(), title.c_str(), num_bins, 0, num_bins));
        
        // count_hists[u]->GetXaxis()->LabelsOption("h");
    
    }
    
    // Canvases for nice histograms
    TCanvas *nue_canvas = new TCanvas("nue_canvas", "#nu_{e} Distribution", 950/3*dets.size(), 345),
           *numu_canvas = new TCanvas("numu_canvas", "#nu_{#mu} Distribution", 950/3*dets.size(), 345);
    numu_canvas->Divide(dets.size(), 1); nue_canvas->Divide(dets.size(), 1);
    
    // Get counts
    std::cout << std::endl << "Getting counts for each sample..." << std::endl;
    for (int o = 0; o < plot_order.size(); o++) {
        
        // Get relevant sample
        EventSample sample = samples[o];
        
        // Initialise temp hists to store counts
        std::string title = sample.fDet+"; Reconstructed Energy (GeV); Counts";
        std::vector <TH1D*> temp_count_hists = {new TH1D("tempbase", title.c_str(), fBins[sample.fDesc].size() - 1, &fBins[sample.fDesc][0])};
        
        for (int u = 0; u < fNumAltUnis; u++) {
            
            std::string name = sample.fDet + "tempalt" + std::to_string(u+1);
            temp_count_hists.push_back(new TH1D(name.c_str(), title.c_str(), fBins[sample.fDesc].size() - 1, &fBins[sample.fDesc][0]));
            
        }
        
        // Tree stuff
        Event *event = new Event;
        sample.tree->SetBranchAddress("events", &event);
        
        // Loop over neutrinos (events)
        int nucount = 0;
        for (int e = 0; e < sample.tree->GetEntries(); e++) {
            
            sample.tree->GetEntry(e);
            
            for (int n = 0; n < event->reco.size(); n++) {
                
                nucount++;
                
                // Add energy to base universe histogram
                double nuE;
                if (fEnergyType == "CCQE") {
                    nuE = event->reco[n].truth.neutrino.energy;
                } else if (fEnergyType == "True") {
                    nuE = event->reco[n].truth.neutrino.energy;
                }
                std::cout << "Neutrino " << nucount << " in " << sample.fDet << " has " << fEnergyType << " energy " << nuE << std::endl;
                temp_count_hists[0]->Fill(nuE);
                
                // Get weights for each alternative universe and fill
                std::vector <double> uweights;
                unsigned truth_ind = event->reco[n].truth_index;
                if (fWeightKey == "GetWeights") {
                    uweights = get_uni_weights(event->truth[truth_ind].weights, fNumAltUnis);
                } else {
                    uweights = event->truth[truth_ind].weights[fWeightKey];
                }
                
                for (int u = 0; u < uweights.size(); u++) {
                    temp_count_hists[u+1]->Fill(nuE, uweights[u]);
                }
                
            }
        }
        
        std::cout << std::endl << "For sample: " << sample.fDet << ", " << sample.fDesc << ", there were " << nucount << " neutrinos." << std::endl;
        
        // Rescale to desired POT
        for (int u = 0; u < temp_count_hists.size(); u++) {
            for (int b = 0; b < temp_count_hists[u]->GetNbinsX(); b++) {
                
                double bincontent = temp_count_hists[u]->GetBinContent(b+1);
                temp_count_hists[u]->SetBinContent(b+1, bincontent * fScaleTargets[sample.fDet] 
                                                            / sample.fScaleFactor);
                
            }
        }
        
        // Pass onto the big histograms and get energies
        for (int h = 0; h < temp_count_hists.size(); h++) {
            
            for (int bin = 0; bin < temp_count_hists[h]->GetNbinsX(); bin++) {
                
                if (h == 0) energies.push_back(temp_count_hists[h]->GetBinCenter(bin+1));
                
                count_hists[h]->SetBinContent(1+offset[o]+bin, temp_count_hists[h]->GetBinContent(1+bin));
                
            }
            
            std::string label = sample.fDet + " " + sample.fDesc;
            count_hists[h]->GetXaxis()->SetBinLabel((offset[o]+offset[o+1])/2, label.c_str());
            
        }
        
        // Add numu and nue hists to canvases
        std::vector <std::string> dets_inorder = get_dets_inorder(samples);
        for (int d = 0; d < dets_inorder.size(); d++) {
            if (sample.fDet == dets_inorder[d] && (sample.fDesc == "#nu_{#mu}" || sample.fDesc == "#nu_{e}")) {
                
                for (int b = 0; b < temp_count_hists[0]->GetNbinsX(); b++) {
                    double binwidth = temp_count_hists[0]->GetBinWidth(b+1), 
                         bincontent = temp_count_hists[0]->GetBinContent(b+1);
                    temp_count_hists[0]->SetBinContent(b+1, bincontent/binwidth);
                }
                
                numu_canvas->cd(d+1);
                if (sample.fDesc == "#nu_{e}") nue_canvas->cd(d+1);
            
            }
        }
        
        temp_count_hists[0]->SetStats(kFALSE);
        temp_count_hists[0]->Draw("hist");
        
    }
    
    TCanvas *tempcanvas = new TCanvas();
    count_hists[0]->Draw(); count_hists[0]->SetStats(kFALSE);
    tempcanvas->SaveAs("/sbnd/data/users/gavarela/selection/new/cov_output/base.pdf");
    count_hists[5]->Draw(); count_hists[5]->SetStats(kFALSE);
    tempcanvas->SaveAs("/sbnd/data/users/gavarela/selection/new/cov_output/alt5.pdf");
    
    numu_canvas->SaveAs("/sbnd/data/users/gavarela/selection/new/cov_output/numu_cts.pdf");
    if (nue_appearance == 1) {
        nue_canvas->SaveAs("/sbnd/data/users/gavarela/selection/new/cov_output/nue_cts.pdf");
    }
    
    
    //// Get covariances, fractional covariances and correlation coefficients
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Getting covs..." << std::endl;
    
    // Bins for cov, fcov and corr
    std::vector <double> covbins = {};
    for (int o = 0; o < plot_order.size(); o++) {
        std::string desc = plot_order[o].substr(plot_order[o].find("_")+1, plot_order[o].length());
        for (int i = 0; i < fBins[desc].size(); i++) {
            if (o == 0 && i == 0) {
                covbins.push_back(fBins[desc][i]);
            } else if (i > 0) {
                covbins.push_back(covbins[covbins.size()-1] + fBins[desc][i] - fBins[desc][i-1]);
            }
        }
    }
    covbins.pop_back();
    
    if (num_bins == covbins.size() - 1) {
        std::cout << std::endl << "NOT SAME SIZE COVBINS!!!" << std::endl; 
    } else {
        std::cout << std::endl << "covbins same size :)" << std::endl;
    }
    
    int do_varied_bins = 0;
    
    // Covariance and fractional covariance
    TH2D *cov, *fcov;
    if (do_varied_bins == 0) {
        cov = new TH2D("cov", "Covariance Matrix", num_bins, 0, num_bins, num_bins, 0, num_bins);
        fcov = new TH2D("fcov", "Fractional Covariance Matrix", num_bins, 0, num_bins, num_bins, 0, num_bins);
    } else {
        cov = new TH2D("cov", "Covariance Matrix", num_bins, &covbins[0], num_bins, &covbins[0]);
        fcov = new TH2D("fcov", "Fractional Covariance Matrix", num_bins, &covbins[0], num_bins, &covbins[0]);
    }
    
    for (int i = 0; i < cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov->GetNbinsY(); j++) {
            
            double covij = 0;
            for (int u = 0; u < fNumAltUnis; u++) {
                covij += (count_hists[0]->GetBinContent(i+1) - count_hists[u+1]->GetBinContent(i+1)) * 
                         (count_hists[0]->GetBinContent(j+1) - count_hists[u+1]->GetBinContent(j+1));
            }
            covij /= fNumAltUnis;
            cov->SetBinContent(i+1, j+1, covij);
            
            double fcovij = covij / (count_hists[0]->GetBinContent(i+1) * count_hists[0]->GetBinContent(j+1));
            fcov->SetBinContent(i+1, j+1, fcovij);
            
        }
    }
    
    // Pearson Correlation Coefficients
    TH2D *corr;
    if (do_varied_bins == 0) {
        corr = new TH2D("corr", "Correlation Matrix", num_bins, 0, num_bins, num_bins, 0, num_bins);
    } else {
        corr = new TH2D("corr", "Correlation Matrix", num_bins, &covbins[0], num_bins, &covbins[0]);
    }
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
    //energies was already done directly;
    
    sample_order = plot_order;
    sample_bins = offset;
    
    covmat = cov;
    fcovmat = fcov;
    corrmat = corr;
    
    
    
}

}   // namespace SBNOsc
}   // namespace ana
