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
EventSample::EventSample(TFile* _file, float ScaleFactor, std::string Det, std::string Desc) : file(_file), fScaleFactor(ScaleFactor), fDet(Det), fDesc(Desc) {

    // Tree
    tree = (TTree*) file->Get("sbnana");

};

// Function that gets scale factors (weights) for different universes
std::vector <double> get_uni_weights(std::map <std::string, std::vector <double> > weights, std::vector<std::string> keys, int n_unis) {
    
    // Tentative format: universe u scale factor is the product of the u-th entries on each vector 
    // inside the map. For vectors with less than u entries, use the (u - vec_size)-th entry
    
    std::vector <double> uweights;
    
    for (int u = 0; u < n_unis; u++) {
                
        double weight = 1;
        int wind;
        
        for (auto const &key: keys) {
            std::vector<double>& this_weights = weights[key];
            wind = u % this_weights.size();
            weight *= this_weights.at(wind);
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
        
        // Output directory
        fOutputDirectory = (*config).get("OutputDirectory", "./").asString();
        
        // Weight and universe stuff
        // WeightKey can either be the name of a single weight, or a name signifying a list of weights,
        // or a list of weights to be oscillated

        Json::Value configWeightKey = (*config)["Covariance"].get("WeightKey", ""); 
        if (configWeightKey.type() == Json::arrayValue) {
          for (auto const& keyName: configWeightKey) {
            fWeightKeys.push_back(keyName.asString());
          }
        }
        else {
          std::string WeightKey = configWeightKey.asString();

          // either way, we'll need an event to get some info about the weights
          // Relies on all 'weights' branches being the same for all samples – should be true
          Event *tempev = new Event;
          samples[0].tree->SetBranchAddress("events", &tempev);
          samples[0].tree->GetEntry(0);

          // get the number of universes
          if (WeightKey == "GetWeights" || WeightKey == "Flux" || WeightKey == "Cross-Section") {
              fNumAltUnis = (*config)["Covariance"].get("NumAltUnis", -7).asInt();
          } else {
              fNumAltUnis = tempev->truth[0].weights[WeightKey].size();
          }
          // get the list of weights
          if (WeightKey == "GetWeights") {
              for (auto it : tempev->truth[0].weights) {
                  fWeightKeys.push_back(it.first);
              }
          } else if (WeightKey == "Flux") {
              for (auto it : tempev->truth[0].weights) {
                  if (it.first.find("genie") > it.first.size()) fWeightKeys.push_back(it.first);
              }
          } else if (WeightKey == "Cross-Section") {
              for (auto it : tempev->truth[0].weights) {
                  if (it.first.find("genie") <= it.first.size()) fWeightKeys.push_back(it.first);
              }
          } else {
              fWeightKeys.push_back(WeightKey);
          }
        }
        
        // Type of energy
        fEnergyType = (*config)["Covariance"].get("EnergyType", "").asString();
        
        // Further selection and rejection 'efficiencies'
        fSelectionEfficiency = (*config)["Covariance"].get("SelectionEfficiency", -1e99).asDouble();
        fRejectionEfficiency = (*config)["Covariance"].get("RejectionEfficiency", -1e99).asDouble();


        // what to include in covariance calculation
        fIncludeBackground = (*config)["Covariance"].get("IncludeBackground", true).asBool();
        
        // True energy binning
        fTrueELims = {};
        for (auto binlim : (*config)["Covariance"]["TrueELims"]) {
            fTrueELims.push_back(binlim.asDouble());
        }
        fNumTrueEBins = (*config)["Covariance"].get("NumTrueEBins", -1).asInt();
        
        // Detector dimensions and distances
        fNumDistBinsPerMeter = (*config)["Covariance"].get("NumDistanceBinsPerMeter", -1).asInt();
        fDetDists = {}; fDetDims = {};
        for (auto det : (*config)["Covariance"]["DetectorDimensions"]) {
            
            std::string detname = det["Detector"].asString();
            
            std::vector <std::vector <double > > xyz_lims = {{}, {}, {}};
            for (auto lim : det["X"]) xyz_lims[0].push_back(lim.asDouble());
            for (auto lim : det["Y"]) xyz_lims[1].push_back(lim.asDouble());
            for (auto lim : det["Z"]) xyz_lims[2].push_back(lim.asDouble());
            
            fDetDims.insert({detname, xyz_lims});

            float distance = det["Distance"].asFloat();
            fDetDists.insert({detname, distance});
            
        }
        
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
    
    std::cout << "Got all parameters from the configuration file." << std::endl;
    
    
    //// Get counts on each (base and alternative) universe
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Some stuff related to binning and plotting
        // If nue appearance, set binning for all numus to same as nues
    int nue_appearance = (fBins.find("#nu_{e}") != fBins.end());
    if (nue_appearance == 1) {
        fBins["#nu_{#mu}"] = fBins["#nu_{e}"];
    }
    
        // Get plotting order of samples
    sample_order = get_plot_order(samples);
    
        // Number of bins needed in big histogram (for covariance)
        // And bin 'boundaries' between each separate sample
    int num_bins = 0;
    sample_bins;
    for (int o = 0; o < sample_order.size(); o++) {
        
        sample_bins.push_back(num_bins);
        
        std::string binkey = sample_order[o].substr(sample_order[o].find("_")+1, sample_order[o].length());
        num_bins += fBins[binkey].size() - 1;
        
    }
    sample_bins.push_back(num_bins);
    
        // Number of distance bins needed in 3D histogram of true CC interactions
    std::map <std::string, std::vector <double> > dist_bin_lims;
    std::map <std::string, int> dist_bin_nums;
    for (auto it : fDetDims) {
        
        double min_dist, max_dist; 
        min_dist = fDetDists[it.first];
       
std::cout << "Detector " << it.first << " gave min_dist" << min_dist << std::endl;
 
        float xlen = (it.second[0][1] - it.second[0][0])/100000 /* cm -> km */,
              ylen = (it.second[1][1] - it.second[1][0])/100000 /* cm -> km */,
              zlen = (it.second[2][1] - it.second[2][0])/100000 /* cm -> km */;
        max_dist = TMath::Sqrt( xlen*xlen + ylen*ylen + (min_dist+zlen)*(min_dist+zlen) );
        
std::cout << "    and max_dist " << max_dist << std::endl;
std::cout << "    and number of bins - (max_dist - min_dist)*(fNumDistBinsPerMeter*1000 /* 1/m -> 1/km */) - : " << (max_dist - min_dist)/(fNumDistBinsPerMeter*1000 /* 1/m -> 1/km */) << std::endl;

std::cout << fNumDistBinsPerMeter << " is num bins per meter" << std::endl;

        dist_bin_lims.insert({it.first, {min_dist, max_dist}});
        dist_bin_nums.insert({it.first, (max_dist - min_dist)*(fNumDistBinsPerMeter*1000 /* 1/m -> 1/km */)});
    }
    
std::cout << std::endl << "Dist bin Lims:";
std::cout << std::endl << "uBooNE: ";
for (auto x : dist_bin_lims["MicroBooNE"]) std::cout << x << "  ";
std::cout << std::endl << "SBND: ";
for (auto x : dist_bin_lims["SBND"]) std::cout << x << "  ";
std::cout << std::endl << "ICARUS: ";
for (auto x : dist_bin_lims["ICARUS"]) std::cout << x << "  ";

    std::cout << std::endl << "Doing distance bins:" << std::endl;
    
    int num_dist_bins = 0;
    dist_bins = {}; sample_dist_bins = {0};
    for (std::string sample : sample_order) {
        
        std::cout << "   Sample = " << sample << " with ";
        
        for (auto it : dist_bin_nums) {
            if (sample.find(it.first) != std::string::npos) {
                
                num_dist_bins += it.second;
                sample_dist_bins.push_back(sample_dist_bins[sample_dist_bins.size()-1] + it.second);
                
                std::cout << it.second << " bins:" << std::endl;
                
                for (int i = 0; i < it.second + 1; i++) {
                    dist_bins.push_back(dist_bin_lims[it.first][0] + i*(dist_bin_lims[it.first][1] - dist_bin_lims[it.first][0])/(it.second));
                    
                    std::cout << dist_bins[dist_bins.size()-1] << ", ";
                }
                
                std::cout << std::endl;
                
                break;
            }
        }
        
    }
    
    std::cout << "And sample_dist_bins came out as ";
    for (double bin : sample_dist_bins) {
        std::cout << bin << ", ";
    }
    std::cout << std::endl;
    
std::cout << std::endl << "With the bins themselves:" << std::endl;
for (auto bin : dist_bins) std::cout << bin << ", ";
std::cout << std::endl << std::endl;


    // Large (meaningless x-axis) histograms for cov
    std::vector <TH1D*> count_hists = {new TH1D("base", "Base Uni. Counts; Bin; Counts", num_bins, 0, num_bins)};
    
    for (int u = 0; u < fNumAltUnis; u++) {
        
        std::string name = "alt" + std::to_string(u+1), 
                    title = "Alt. Uni. " + std::to_string(u+1) + " Counts; Bin; Counts";
        
        count_hists.push_back(new TH1D(name.c_str(), title.c_str(), num_bins, 0, num_bins));
    
    }
    
    // Large (meaningless x-axis) histograms (one 3d) for oscillation calculations later on
    bkg_counts = new TH1D("Background", "Background Counts; Reconstructed Energy Bin; Counts", num_bins, 0, num_bins);
    
    nu_counts = new TH3D("Neutrinos", "Neutrino Counts; True Energy Bin; Reconstructed Energy Bin; Distance Bin", fNumTrueEBins*sample_order.size(), 0, fNumTrueEBins*sample_order.size(), num_bins, 0, num_bins, num_dist_bins, 0, num_dist_bins);
    
    // Vectors to hold nice histograms
    std::vector <TH1D*> numu_cts(dets.size(), new TH1D()), numu_bkg(dets.size(), new TH1D()),
                        nue_cts(dets.size(), new TH1D()), nue_bkg(dets.size(), new TH1D());
    
    // Get counts
    std::cout << std::endl << "Getting counts for each sample..." << std::endl;
    for (int o = 0; o < sample_order.size(); o++) {
        
        // Get relevant sample
        EventSample sample = samples[0];
        int s = 0;
        while (sample.fDet+"_"+sample.fDesc != sample_order[o]) { s++; sample = samples[s]; }
        
        // Initialise temp hists to store counts
            // Base
        std::string title = sample.fDet+"; Reconstructed Energy (GeV); Counts";
        std::vector <TH1D*> temp_count_hists = {new TH1D((sample.fDet+"tempbase").c_str(), title.c_str(), fBins[sample.fDesc].size() - 1, &fBins[sample.fDesc][0])};
        
            // Alt unis
        for (int u = 0; u < fNumAltUnis; u++) {
            
            std::string name = sample.fDet + "tempalt" + std::to_string(u+1);
            temp_count_hists.push_back(new TH1D(name.c_str(), title.c_str(), fBins[sample.fDesc].size() - 1, &fBins[sample.fDesc][0]));
            
        }
        
            // Bkg
        TH1D *temp_bkg_counts = new TH1D((sample.fDet+"tempbkg").c_str(), "", fBins[sample.fDesc].size() - 1, &fBins[sample.fDesc][0]);
        
            // Neutrinos
        std::cout << std::endl << "For sample " << sample.fDet << " " << sample.fDesc << " there are:" << std::endl;
        
        std::cout << "  " << fNumTrueEBins << " true E bins: ";
        std::vector <double> temp_trueE_bins = {};
        for (int i = 0; i < fNumTrueEBins + 1; i++) {
            temp_trueE_bins.push_back(fTrueELims[0] + i*(fTrueELims[1]-fTrueELims[0])/(fNumTrueEBins));
            std::cout << fTrueELims[0] + i*(fTrueELims[1]-fTrueELims[0])/(fNumTrueEBins) << ", ";
        }
        std::cout << std::endl;
        
        std::cout << "  " << dist_bin_nums[sample.fDet] << " distance bins: ";
        std::vector <double> temp_dist_bins = {};
        for (int i = 0; i < dist_bin_nums[sample.fDet]+1; i++) {
            temp_dist_bins.push_back(dist_bin_lims[sample.fDet][0] + i*(dist_bin_lims[sample.fDet][1]-dist_bin_lims[sample.fDet][0])/(dist_bin_nums[sample.fDet]));
            std::cout << dist_bin_lims[sample.fDet][0] + i*(dist_bin_lims[sample.fDet][1]-dist_bin_lims[sample.fDet][0])/(dist_bin_nums[sample.fDet]) << ", ";
        }
        std::cout << std::endl;
        
        std::cout << std::endl << "Will now build the TH3D with axes:";
        std::cout << std::endl << "X (size = " << fNumTrueEBins << "): ";
        for (auto a : temp_trueE_bins) std::cout << a << ", ";
        std::cout << std::endl << "Y (size = " << fBins[sample.fDesc].size() - 1 << "): ";
        for (auto a : fBins[sample.fDesc]) std::cout << a << ", ";
        std::cout << std::endl << "Z (size = " << dist_bin_nums[sample.fDet] << "): ";
        for (auto a : temp_dist_bins) std::cout << a << ", ";
         
        TH3D *temp_nu_counts = new TH3D((sample.fDet+"tempnu").c_str(), "", fNumTrueEBins, &temp_trueE_bins[0], fBins[sample.fDesc].size() - 1, &fBins[sample.fDesc][0], dist_bin_nums[sample.fDet], &temp_dist_bins[0]);
        
        std::cout << "Defined the TH3D in the line above." << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl;

std::cout << "Will loop over nus now. The distances are (in m, for select events)";

        // Loop over neutrinos (events in tree)
        Event *event = new Event;
        sample.tree->SetBranchAddress("events", &event);
        
        // Number of neutrinos in sample -- NOT necessarily number of neutrinos included in covariance matrix
        int nucount = 0;
        for (int e = 0; e < sample.tree->GetEntries(); e++) {
            
            sample.tree->GetEntry(e);
            
            for (int n = 0; n < event->reco.size(); n++) {
                nucount++;

                unsigned truth_ind = event->reco[n].truth_index;
                int isCC = event->truth[truth_ind].neutrino.iscc;
 
                
                // Get energy
                double nuE, true_nuE = event->reco[n].truth.neutrino.energy;
                if (fEnergyType == "CCQE") {
                    nuE = event->truth[truth_ind].neutrino.eccqe;
                } else if (fEnergyType == "True") {
                    nuE = true_nuE;
                } else if (fEnergyType == "Reco") {
                    nuE = event->reco[n].reco_energy;
                }
                if (nuE < fBins[sample.fDesc][0] || nuE > fBins[sample.fDesc][fBins[sample.fDesc].size()-1]) {
                    continue;
                } else if (true_nuE < fTrueELims[0] || true_nuE > fTrueELims[1]) {
                    std::cout << std::endl << "NUE IN RANGE, TRUE E NOT!!!   nuE = " << nuE << " and true_nuE = " << true_nuE << std::endl;
                    continue;
                }

                // Apply selection (or rejection) efficiencies
                double wgt = isCC*(fSelectionEfficiency) + (1-isCC)*(1 - fRejectionEfficiency);
                
                // if not including background sample in covariance, skip this one if bkg
                if (!fIncludeBackground && !isCC) {
                    // Add to base count histogram
                    temp_count_hists[0]->Fill(nuE, wgt);
                
                    // Get weights for each alternative universe
                    std::vector <double> uweights = uweights = get_uni_weights(event->truth[truth_ind].weights, fWeightKeys, fNumAltUnis);
                
                    // Fill alternative universe histograms
                    for (int u = 0; u < uweights.size(); u++) {
                        temp_count_hists[u+1]->Fill(nuE, wgt*uweights[u]);
                    }
                }
                // always fill the CV counts histo for output
                CV_counts->Fill(nuE, wgt);
                
                // Get distance travelled (assuming nu started at (x, y, z) = (0, 0, min_det_zdim - det_dist))
                double dx = (event->truth[truth_ind].neutrino.position.X() - (fDetDims[sample.fDet][0][1] + fDetDims[sample.fDet][0][0])/2) / 100000 /* cm -> km */,
                       dy = (event->truth[truth_ind].neutrino.position.Y() - (fDetDims[sample.fDet][1][1] + fDetDims[sample.fDet][1][0])/2) / 100000 /* cm -> km */,
                       dz = (event->truth[truth_ind].neutrino.position.Z() - fDetDims[sample.fDet][2][0]) / 100000 /* cm -> km */;
                double dist = TMath::Sqrt( dx*dx + dy*dy + (fDetDists[sample.fDet] + dz)*(fDetDists[sample.fDet] + dz) );

if (e%200 == 0) {
std::cout << std::endl << "For event " << e << std::endl
<< "X = " << event->truth[truth_ind].neutrino.position.X() << " and dx = " << dx << std::endl
<< "Y = " << event->truth[truth_ind].neutrino.position.Y() << " and dy = " << dy << std::endl
<< "Z = " << event->truth[truth_ind].neutrino.position.Z() << " and dz = " << dz << std::endl
<< "dist = " << dist * 1000;
}

if ((dist < temp_dist_bins[0]) | (dist > temp_dist_bins[temp_dist_bins.size()-1])) {
std::cout << std::endl << std::endl << " DIST OUTSIDE OF BOUNDS!!! " << dist << " !in (" << temp_dist_bins[0] << ", " << temp_dist_bins[temp_dist_bins.size()-1] << ")" << std::endl;
} 

                // Fill chisq histograms
                bool isnu = (sample.fDesc.find("#nu") != std::string::npos);
                if (isnu && isCC) {
                    temp_nu_counts->Fill(true_nuE, nuE, dist, wgt);
                } else {
                    temp_bkg_counts->Fill(nuE, wgt);
                }
                
            }
        }
        
std::cout << std::endl << std::endl << std::endl << std::endl;

        std::cout << "Sample " << sample.fDet << " - " << sample.fDesc << " had " << nucount << " neutrinos." << std::endl;
        
        // Rescale to desired POT
        for (int h = 0; h < temp_count_hists.size(); h++) {
            temp_count_hists[h]->Scale(fScaleTargets[sample.fDet] / sample.fScaleFactor);
        }
        temp_bkg_counts->Scale(fScaleTargets[sample.fDet] / sample.fScaleFactor);
        temp_nu_counts->Scale(fScaleTargets[sample.fDet] / sample.fScaleFactor);
        
        
        // Pass onto the big histograms and get energies
            // base and alt uni counts
        for (int h = 0; h < temp_count_hists.size(); h++) {
            
            for (int bin = 0; bin < temp_count_hists[h]->GetNbinsX(); bin++) {
                
                if (h == 0) energies.push_back(temp_count_hists[h]->GetBinCenter(bin+1));
                
                count_hists[h]->SetBinContent(1+sample_bins[o]+bin, temp_count_hists[h]->GetBinContent(1+bin));
                
            }
        }
        
        for (int rb = 0; rb < temp_nu_counts->GetNbinsY(); rb++) {
            
            // bkg_counts
            bkg_counts->SetBinContent(1+sample_bins[o]+rb, temp_bkg_counts->GetBinContent(1+rb));
            
            // nu_counts
            for (int tb = 0; tb < temp_nu_counts->GetNbinsX(); tb++) {
                for (int db = 0; db < temp_nu_counts->GetNbinsZ(); db++) {
                    
                    nu_counts->SetBinContent(1 + o*fNumTrueEBins + tb, 1 + sample_bins[o] + rb, 
                                             1 + sample_dist_bins[o] + db, 
                                             temp_nu_counts->GetBinContent(1+tb, 1+rb, 1+db));
                    
                }
            }
            
        }
        
        
        // Add numu and nue hists to vectors
        int detind = std::find(dets.begin(), dets.end(), sample.fDet) - dets.begin();
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
        
        
    }
    
    
    //// Get covariances, fractional covariances and correlation coefficients
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Getting covs..." << std::endl;
    
    // Covariance and fractional covariance
    cov = new TH2D("cov", "Covariance Matrix", num_bins, 0, num_bins, num_bins, 0, num_bins);
    fcov = new TH2D("fcov", "Fractional Covariance Matrix", num_bins, 0, num_bins, num_bins, 0, num_bins);
    
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
    corr = new TH2D("corr", "Correlation Matrix", num_bins, 0, num_bins, num_bins, 0, num_bins);
    
    for (int i = 0; i < cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov->GetNbinsY(); j++) {
            
            double corrij = cov->GetBinContent(i+1, j+1) / TMath::Sqrt(cov->GetBinContent(i+1, i+1) * cov->GetBinContent(j+1, j+1));
            corr->SetBinContent(i+1, j+1, corrij);
            
        }
    }
    
    // Add bin labels
    std::vector <TH2D*> hists = {cov, fcov, corr};
    for (int o = 0; o < sample_order.size(); o++) {

        // Get label and position
        std::string label = sample_order[o].replace(sample_order[o].find("_"), 1, " ");
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
    trueEs = {};
    double trueE_binwidth = (fTrueELims[1] - fTrueELims[0])/fNumTrueEBins;
    for (int o = 0; o < sample_order.size(); o++) {
        for (int i = 0; i < fNumTrueEBins; i++) {
            trueEs.push_back(fTrueELims[0] + (i+0.5)*trueE_binwidth);
        }
    }
    
    numu_counts = numu_cts;
    numu_bkgs = numu_bkg;
    if (nue_appearance) {
        nue_counts = nue_cts;
        nue_bkgs = nue_bkg;
    }
    
    
}

}   // namespace SBNOsc
}   // namespace ana
