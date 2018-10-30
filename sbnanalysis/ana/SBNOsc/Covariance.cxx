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
#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THStack.h>
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
EventSample::EventSample(TFile* _file, float ScaleFactor, std::string Det, std::string Desc, std::vector <double> bins, int scale_sample) : file(_file), fScaleFactor(ScaleFactor), fDet(Det), fDesc(Desc), fBins(bins), fScaleSample(scale_sample) {

    // Tree
    tree = (TTree*) file->Get("sbnana");

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
    
    // Get configuration file parameters
    Json::Value* config = core::LoadConfig(configFileName);
    
    if (config != NULL) {
        
        // Output directory
        fOutputDirectory = (*config).get("OutputDirectory", "./").asString();
        
        // Weight and universe stuff
        fWeightKey = (*config)["Covariance"].get("WeightKey", "").asString();
        if (fWeightKey == "GetWeights" || fWeightKey == "Flux" || fWeightKey == "Cross-Section") {
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
        
        // Further selection and rejection 'efficiencies'
        fSelectionEfficiency = (*config)["Covariance"].get("SelectionEfficiency", -1e99).asDouble();
        fRejectionEfficiency = (*config)["Covariance"].get("RejectionEfficiency", -1e99).asDouble();
        
        // Exposure normalisation
        for (auto sample : samples) {
            if (fScaleTargets.find(sample.fDet) == fScaleTargets.end()) {
                fScaleTargets.insert({sample.fDet, (*config)["Covariance"]["ScaleTargets"].get(sample.fDet, -1).asFloat()});
            }
        }
        
        // more documentation is at: https://open-source-parsers.github.io/jsoncpp-docs/doxygen/index.html#_example
    
    }
    
    std::cout << "Got all parameters from the configuration file." << std::endl;
    
    
    //// Get counts on each (base and alternative) universe
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Some stuff related to binning and plotting
    
        // Number of bins needed in big histogram (for covariance)
        // And bin 'boundaries' between each separate sample
    int num_bins = 0;
    std::vector <int> sample_bins;
    for (auto sample : samples) {
        
        sample_bins.push_back(num_bins);
        num_bins += sample.fBins.size() - 1;
        
    }
    sample_bins.push_back(num_bins);
    
    // Large (meaningless x-axis) histograms for cov
    std::vector <TH1D*> count_hists = {new TH1D("base", "Base Uni. Counts; Bin; Counts", num_bins, 0, num_bins)};
    
    for (int u = 0; u < fNumAltUnis; u++) {
        
        std::string name = "alt" + std::to_string(u+1), 
                    title = "Alt. Uni. " + std::to_string(u+1) + " Counts; Bin; Counts";
        
        count_hists.push_back(new TH1D(name.c_str(), title.c_str(), num_bins, 0, num_bins));
    
    }
    
    // Vectors to hold nice histograms
    std::vector <TH1D*> numu_cts(fScaleTargets.size(), new TH1D()), numu_bkg(fScaleTargets.size(), new TH1D()),
                        nue_cts(fScaleTargets.size(), new TH1D()), nue_bkg(fScaleTargets.size(), new TH1D());
    
    // Get counts
    std::cout << std::endl << "Getting counts for each sample..." << std::endl;
    
    int o = 0;
    for (auto sample : samples) {
        
        // What sample did we get?
        std::cout << "Doing " << sample.fDesc << " sample in " << sample.fDet << std::endl;
        
        // Initialise temp hists to store counts
            // Base
        std::string title = sample.fDet+"; Reconstructed Energy (GeV); Counts";
        std::vector <TH1D*> temp_count_hists = {new TH1D((sample.fDet+"tempbase").c_str(), title.c_str(), sample.fBins.size() - 1, &sample.fBins[0])};
        
            // Alt unis
        for (int u = 0; u < fNumAltUnis; u++) {
            
            std::string name = sample.fDet + "tempalt" + std::to_string(u+1);
            temp_count_hists.push_back(new TH1D(name.c_str(), title.c_str(), sample.fBins.size() - 1, &sample.fBins[0]));
            
        }
        
        // Loop over neutrinos (events in tree)
        Event *event = new Event;
        sample.tree->SetBranchAddress("events", &event);
        
        int nucount = 0;
        for (int e = 0; e < sample.tree->GetEntries(); e++) {
            
            sample.tree->GetEntry(e);
            
            for (int n = 0; n < event->reco.size(); n++) {
                
                nucount++;
                
                unsigned truth_ind = event->reco[n].truth_index;
                
                // Get energy
                double nuE, true_nuE = event->reco[n].truth.neutrino.energy;
                if (fEnergyType == "CCQE") {
                    nuE = event->truth[truth_ind].neutrino.eccqe;
                } else if (fEnergyType == "True") {
                    nuE = true_nuE;
                } else if (fEnergyType == "Reco") {
                    nuE = event->reco[n].reco_energy;
                }
                
                // Apply selection (or rejection) efficiencies
                int isCC = event->truth[truth_ind].neutrino.iscc;
                double wgt = isCC*(fSelectionEfficiency) + (1-isCC)*(1 - fRejectionEfficiency);
                
                // Add to base count histogram
                temp_count_hists[0]->Fill(nuE, wgt);
                
                // Get weights for each alternative universe
                std::vector <double> uweights;
                if (fWeightKey == "GetWeights") {
                    uweights = get_uni_weights(event->truth[truth_ind].weights, fNumAltUnis);
                } else if (fWeightKey == "Flux") {
                    std::map <std::string, std::vector<double> > tempweights;
                    for (auto it : event->truth[truth_ind].weights) {
                        if (it.first.find("genie") > it.first.size()) tempweights.insert(it);
                    }
                    uweights = get_uni_weights(tempweights, fNumAltUnis);
                } else if (fWeightKey == "Cross-Section") {
                    std::map <std::string, std::vector<double> > tempweights;
                    for (auto it : event->truth[truth_ind].weights) {
                        if (it.first.find("genie") <= it.first.size()) tempweights.insert(it);
                    }
                    uweights = get_uni_weights(tempweights, fNumAltUnis);
                } else {
                    uweights = event->truth[truth_ind].weights[fWeightKey];
                }
                
                // Fill alternative universe histograms
                for (int u = 0; u < uweights.size(); u++) {
                    temp_count_hists[u+1]->Fill(nuE, wgt*uweights[u]);
                }
                
            }
        }
        
        // Rescale to desired POT
        for (int h = 0; h < temp_count_hists.size(); h++) {
            temp_count_hists[h]->Scale(fScaleTargets[sample.fDet] / sample.fScaleFactor);
        }
        
        // Pass onto the big histograms
            // base and alt uni counts
        for (int h = 0; h < temp_count_hists.size(); h++) {
            
            for (int bin = 0; bin < temp_count_hists[h]->GetNbinsX(); bin++) {
                count_hists[h]->SetBinContent(1+sample_bins[o]+bin, temp_count_hists[h]->GetBinContent(1+bin));
            }
            
        }
        
        /*
        // Add numu and nue hists to vectors
        if (sample.fDesc == "#nu_{#mu}" ) {
            
            numu_bkg[detind] = temp_bkg_counts;
            numu_bkg[detind]->SetName((sample.fDet+"_numu_bkg").c_str());
            
            numu_cts[detind] = (TH1D*) temp_nu_counts->ProjectionY();
            numu_cts[detind]->SetName((sample.fDet+"_numu_cts").c_str());
            
        } else if (sample.fDesc == "#nu_{e}") {
            
            nue_bkg[detind] = temp_bkg_counts;
            nue_bkg[detind]->SetName((sample.fDet+"_nue_bkg").c_str());
            
            nue_cts[detind] = (TH1D*) temp_nu_counts->ProjectionY();
            nue_cts[detind]->SetName((sample.fDet+"_nue_cts").c_str());
            
        }
        */
        
        o++;
        
    }
    
    
    //// Get covariances, fractional covariances and correlation coefficients
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Getting covs..." << std::endl;
    
    std::string pre_title = "";
    if (fWeightKey != "GetWeights") pre_title = fWeightKey + " ";
    
    // Covariance and fractional covariance
    cov = new TH2D("cov", (pre_title+"Covariance Matrix").c_str(), num_bins, 0, num_bins, num_bins, 0, num_bins);
    fcov = new TH2D("fcov", ("Fractional "+pre_title+"Covariance Matrix").c_str(), num_bins, 0, num_bins, num_bins, 0, num_bins);
    
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
    corr = new TH2D("corr", (pre_title+"Correlation Matrix").c_str(), num_bins, 0, num_bins, num_bins, 0, num_bins);
    
    for (int i = 0; i < cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov->GetNbinsY(); j++) {
            
            double corrij = cov->GetBinContent(i+1, j+1) / TMath::Sqrt(cov->GetBinContent(i+1, i+1) * cov->GetBinContent(j+1, j+1));
            corr->SetBinContent(i+1, j+1, corrij);
            
        }
    }
    
    // Add bin labels
    std::vector <TH2D*> hists = {cov, fcov, corr};
    for (int o = 0; o < samples.size(); o++) {

        // Get label and position
        std::string label = samples[o].fDet + " " + samples[o].fDesc;
        int bin = (sample_bins[o] + sample_bins[o+1])/2;

        // Set label
        for (TH2D* hist : hists) {
            hist->GetXaxis()->SetBinLabel(1+bin, label.c_str());
            hist->GetYaxis()->SetBinLabel(1+bin, label.c_str());
        }

    }

    for (TH2D* hist : hists) { hist->GetXaxis()->LabelsOption("h"); hist->GetYaxis()->LabelsOption("v"); }
    
    std::cout << std::endl << "  Got covs." << std::endl;
    
    
    //// Output relevant objects
    //// ~~~~~~~~~~~~~~~~~~~~~~~
    
    /*
    numu_counts = numu_cts;
    numu_bkgs = numu_bkg;
    if (nue_appearance) {
        nue_counts = nue_cts;
        nue_bkgs = nue_bkg;
    }
    */
    
    
}

}   // namespace SBNOsc
}   // namespace ana
