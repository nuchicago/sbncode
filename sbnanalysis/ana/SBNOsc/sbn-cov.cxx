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

int main(int argc, char* argv[]) {
    
    std::cout << std::endl << "Hello!" << std::endl << std::endl;
    
    //// Build sample (temporary)
    //// ~~~~~~~~~~~~~~~~~~~~~~~~
    
    char *configFileName = argv[1];
    
    Json::Value* config = core::LoadConfig(configFileName);
    assert(config != NULL);
    
    std::vector <ana::SBNOsc::EventSample> samples;
    
    for (auto sample : (*config)["EventSamples"]) {
        
        TFile *file = new TFile((sample["path"].asString()).c_str());
        float scalefactor = sample["scalefactor"].asFloat();
        std::string det = sample["det"].asString(),
                        desc = sample["desc"].asString();
        
        samples.push_back(ana::SBNOsc::EventSample(filename, scalefactor, det, desc));
        
    }
    
    assert(!samples.empty());
    
    
    //// Get covariances and write outputs to ROOT file(s)
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ana::SBNOsc::Covariance cov(samples, configFileName);
    
    std::string directory = "/sbnd/data/users/gavarela/selection/new/cov_output/";
    TFile* newfile = TFile::Open((directory + "cov_output.root").c_str(), "recreate");
    assert(newfile && newfile->IsOpen());
    
    cov.covmat->Write();
    cov.fcovmat->Write();
    cov.corrmat->Write();
    
    TCanvas *canvas = new TCanvas();
    cov.covmat->Draw("colz"); cov.covmat->SetStats(kFALSE); canvas->SaveAs((directory + "cov_plot.pdf").c_str());
    cov.fcovmat->Draw("colz"); cov.fcovmat->SetStats(kFALSE); canvas->SaveAs((directory + "fcov_plot.pdf").c_str());
    cov.corrmat->Draw("colz"); cov.corrmat->SetStats(kFALSE); canvas->SaveAs((directory + "corr_plot.pdf").c_str());
    
    
    
    return 0;
    
}