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
EventSample::EventSample(TFile* _file, float ScaleFactor, std::string Det, std::string Desc, std::vector <double> bins, int scale_sample, std::string nutype) : file(_file), fScaleFactor(ScaleFactor), fDet(Det), fDesc(Desc), fBins(bins), fScaleSample(scale_sample), fNuType(nutype) {

    // Tree
    tree = (TTree*) file->Get("sbnana");

};   
    
// Gets scale factors (weights) for different universes
GetUniWeights(std::map <std::string, std::vector <double> > weights, int n_unis) {
    
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

Covariance::Covariance(std::vector<EventSample> samples, char *configFileName) {
    
    //// Gets parameters from config file
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Get configuration file parameters
    Json::Value* config = core::LoadConfig(configFileName);
    
    if (config != NULL) {
        
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
        
        // Use only signal events in cov?
        fSignalOnly = (*config)["Covariance"].get("SignalOnly", 7).asBool();
        
        // more documentation is at: https://open-source-parsers.github.io/jsoncpp-docs/doxygen/index.html#_example
    
    }
    
    ev_samples = samples;
    
    std::cout << "Got all parameters from the configuration file." << std::endl;
    
}

ScanEvents() {
    
    //// Geta counts for each (base and alternative) universe
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Bin 'boundaries'
    num_bins = 0;
    for (auto sample : samples) {
        sample_bins.push_back(num_bins);
        num_bins += sample.fBins.size() - 1;   
    }
    sample_bins.push_back(num_bins);
    
    // Vector to hold counts in all 'universes'
    for (int u = 0; u < fNumAltUnis; u++) nu_counts.push_back({});
    
    // Get counts
    std::cout << std::endl << "Getting counts for each sample..." << std::endl;
    
    int o = 0;
    for (auto sample : ev_samples) {
        
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
        
            // Background counts
        TH1D *temp_bkg_counts = new TH1D((sample.fDet+"tempbkg").c_str(), title.c_str(), sample.fBins.size()-1, &sample.fBins[0]);
        
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
                
                // Add to base and bkg count histogram
                temp_count_hists[0]->Fill(nuE, wgt);
                if ((sample.fNuType == "numu" || sample.fNuType == "nue") && !isCC) {
                    temp_bkg_counts->Fill(nuE, wgt);
                }
                
                if (fSignalOnly && 
                    ((sample.fNuType == "numu" || sample.fNuType == "nue") || isCC)) {
                    
                    // Get weights for each alternative universe
                    std::vector <double> uweights;
                    if (fWeightKey == "GetWeights") {
                        uweights = GetUniWeights(event->truth[truth_ind].weights, fNumAltUnis);
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
        }
        
        // Add to vector and scale to desired POT
        for (int h = 0; h < temp_count_hists.size(); h++) {
            for (int bin = 0; bin < temp_count_hists[h]->GetNbinsX(); bin++) {

                double scaled_counts = temp_count_hists[h]->GetBinContent(1+bin) * fScaleTargets[sample.fDet] / sample.fScaleFactor;
                nu_counts[h].push_back(scaled_counts);
                
                if (h == 0) {
                    double scaled_bkg_counts = temp_bkg_hist->GetBinContent(1+bin) * fScaleTargets[sample.fDet] / sample.fScaleFactor;
                    bkg_counts.push_back(scaled_bkg_counts);
                }
            
            }
        }
        
        o++;
        
    }
    
}

GetCovs() {
    
    //// Get covariances, fractional covariances and correlation coefficients
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Getting covs..." << std::endl;
    
    // Plot data
    std::vector <double> cov_counts = nu_counts[0];
    if (fSignalOnly) {
        for (int i = 0; i < cov_counts.size(); i++) cov_counts[i] -= bkg_counts[i];
    }
    
    // Plot name
    std::string pre_title = "";
    if (fWeightKey != "GetWeights") pre_title = fWeightKey + " ";
    
    // Covariance and fractional covariance
    cov = new TH2D("cov", (pre_title+"Covariance Matrix").c_str(), num_bins, 0, num_bins, num_bins, 0, num_bins);
    fcov = new TH2D("fcov", ("Fractional "+pre_title+"Covariance Matrix").c_str(), num_bins, 0, num_bins, num_bins, 0, num_bins);
    
    for (int i = 0; i < cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov->GetNbinsY(); j++) {
            
            double covij = 0;
            for (int u = 0; u < fNumAltUnis; u++) {
                covij += (cov_counts[i] - nu_counts[u+1][i]) * (cov_counts[j] - nu_counts[u+1][j]);
            }
            covij /= fNumAltUnis;
            cov->SetBinContent(i+1, j+1, covij);
            
            double fcovij = covij / (cov_counts[i] * cov_counts[j]);
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
    for (int o = 0; o < ev_samples.size(); o++) {

        // Get label and position
        std::string label = ev_samples[o].fDet + " " + ev_samples[o].fDesc;
        int bin = (sample_bins[o] + sample_bins[o+1])/2;

        // Set label
        for (TH2D* hist : hists) {
            hist->GetXaxis()->SetBinLabel(1+bin, label.c_str());
            hist->GetYaxis()->SetBinLabel(1+bin, label.c_str());
        }

    }

    for (TH2D* hist : hists) { hist->GetXaxis()->LabelsOption("h"); hist->GetYaxis()->LabelsOption("v"); }
    
    std::cout << std::endl << "  Got covs." << std::endl;
    
}
    
GetCounts() {
    
    // Vectors to hold nice histograms
    std::vector <TH1D*> numu_cts(fScaleTargets.size(), new TH1D()), numu_bkg(fScaleTargets.size(), new TH1D()),
                        nue_cts(fScaleTargets.size(), new TH1D()), nue_bkg(fScaleTargets.size(), new TH1D());
    
    // Loop
    int numu_ind = 0, nue_ind = 0, sample_ind = 0;
    for (auto sample: ev_samples) {
        
        if (sample.fNuType == "numu") {
            
            numu_bkg[numu_ind] = new TH1D((sample.fDet+"_numu_bkg").c_str(), "",
                                          sample.fBins.size() - 1, &sample.fBins[0]);
            numu_cts[numu_ind] = new TH1D((sample.fDet+"_numu_cts").c_str(), "",
                                          sample.fBins.size() - 1, &sample.fBins[0]);
            
            for (int o = sample_bins[sample_ind]; o < sample_bins[sample_ind+1]; o++) {
                
                numu_bkg[numu_ind]->SetBinContent(bkg_counts[o]);
                numu_cts[numu_ind]->SetBinContent(nu_counts[o] - bkg_counts[o]);
                
            }
            
            numu_ind++;
            
        } else if (sample.fNuType == "nue") {
            
            nue_bkg[nue_ind] = new TH1D((sample.fDet+"_nue_bkg").c_str(), "",
                                        sample.fBins.size() - 1, &sample.fBins[0]);
            nue_cts[nue_ind] = new TH1D((sample.fDet+"_nue_cts").c_str(), "",
                                        sample.fBins.size() - 1, &sample.fBins[0]);
            
            for (int o = sample_bins[sample_ind]; o < sample_bins[sample_ind+1]; o++) {
                
                nue_bkg[nue_ind]->SetBinContent(bkg_counts[o]);
                nue_cts[nue_ind]->SetBinContent(nu_counts[o] - bkg_counts[o]);
                
            }
            
            nue_ind++;
            
        }
        
        sample_ind++;
        
    }
    
    // Output relevant objects
    numu_counts = numu_cts; numu_bkgs = numu_bkg;
    nue_counts = nue_cts; nue_bkgs = nue_bkg;
    
}

Write(std::string directory) {
    
    // Write Covariances
    TFile* covfile = TFile::Open((directory + "cov.root").c_str(), "recreate");
    assert(covfile && covfile->IsOpen());
    
    cov.cov->Write();
    cov.fcov->Write();
    cov.corr->Write();
    
    // Write Counts
    if (numu_counts.size() + nue_counts.size() > 0) {
        
        TFile* countfile = TFile::Open((directory + "counts.root").c_str(), "recreate");
        assert(countfile && countfile->IsOpen());

        std::vector <std::vector <TH1D*> > hist_vecs = {cov.numu_counts, cov.numu_bkgs};
        if (cov.nue_counts.size() > 0) { hist_vecs.push_back(cov.nue_counts); hist_vecs.push_back(cov.nue_bkgs); }

        for (std::vector <TH1D*> hist_vec : hist_vecs) {
            for (TH1D* hist : hist_vec) hist->Write();
        }
        
    }
    
}

}   // namespace SBNOsc
}   // namespace ana
