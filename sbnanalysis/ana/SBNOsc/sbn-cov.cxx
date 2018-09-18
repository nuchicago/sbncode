/**
 * Generate covariance matrices.
 *
 * Document...
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "Covariance.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>

int main(int argc, char* argv[]) {
    
    std::cout << std::endl << "Hello!" << std::endl << std::endl;
    
    //// Build sample
    //// ~~~~~~~~~~~~
    
    char *configFileName = argv[1];
    
    Json::Value* config = core::LoadConfig(configFileName);
    assert(config);
    
    std::vector <ana::SBNOsc::EventSample> samples;
    
    for (auto sample : (*config)["EventSamples"]) {
        
        TFile *file = new TFile((sample["path"].asString()).c_str());
        float scalefactor = sample["scalefactor"].asFloat();
        std::string det = sample["det"].asString(),
                        desc = sample["desc"].asString();
        
        samples.push_back(ana::SBNOsc::EventSample(file, scalefactor, det, desc));
        
    }
    
    assert(!samples.empty());
    
    
    //// Get covariances and write outputs to ROOT file(s)
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ana::SBNOsc::Covariance cov(samples, configFileName);
    
    std::string directory = (*config).get("OutputDirectory", "./").asString();

    TFile* newfile = TFile::Open((directory + "cov_output.root").c_str(), "recreate");
    assert(newfile && newfile->IsOpen());
    
    cov.cov->Write();
    cov.fcov->Write();
    cov.corr->Write();
    
    // Write counts to file
    TFile* countfile = TFile::Open((directory + "counts.root").c_str(), "recreate");
    assert(countfile && countfile->IsOpen());
    
    std::vector <std::vector <TH1D*> > hist_vecs = {cov.numu_counts, cov.numu_bkgs};
    if (cov.nue_counts.size() > 0) { hist_vecs.push_back(cov.nue_counts); hist_vecs.push_back(cov.nue_bkgs); }
    
    for (std::vector <TH1D*> hist_vec : hist_vecs) {
        for (TH1D* hist : hist_vec) hist->Write();
    }
    
    
    
    return 0;
    
}
