/**
 * Run sensitivity calculation.
 *
 * Document...
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "Covariance.h"
#include "Chi2Sensitivity.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>

int main(int argc, char* argv[]) {

    std::cout << std::endl << "Hello!" << std::endl << std::endl;
    
    //// Build sample
    //// ~~~~~~~~~~~~
    
    char *configFileName = argv[1];
    
    Json::Value* config = core::LoadConfig(configFileName);
    assert(config);
    
    // Output directory
    std::string fOutputDirectory = (*config).get("OutputDirectory", "./").asString();
    
    // Samples
    std::vector <ana::SBNOsc::EventSample> samples;
    
    for (auto sample : (*config)["EventSamples"]) {
        
        TFile *file = new TFile((sample["path"].asString()).c_str());
        float scalefactor = sample["scalefactor"].asFloat();
        std::string det = sample["det"].asString(),
                        desc = sample["desc"].asString();
        std::vector <double> bins = {};
        
        for (auto binlim : sample["binlims"]) {
            bins.push_back(binlim.asDouble());
        }
        int scale_sample = sample.get("scalesample", 0).asInt();
        std::string nutype = sample.get("nutype", "").asString();
        
        samples.push_back(ana::SBNOsc::EventSample(file, scalefactor, det, desc, bins, scale_sample, nutype));
        
    }
    
    assert(!samples.empty());
    
    
    //// Get covariances
    //// ~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Starting covariance procedure..." << std::endl << std::endl;
    
    ana::SBNOsc::Covariance cov(samples, configFileName);
    
    cov.ScanEvents();
    cov.GetCovs();
    cov.GetCounts();
    cov.Write(fOutputDirectory);
    
    
    //// Get sensitivity contours
    //// ~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Starting chisq/sensitivity procedure..." << std::endl << std::endl;
    
    ana::SBNOsc::Chi2Sensitivity chi2(samples, cov, configFileName);
    
    chi2.ScanEvents();
    chi2.GetChi2();
    chi2.GetContours();
    chi2.Write();
    
    
    //// Save plots
    //// ~~~~~~~~~~
    
    int SavePDFs = (*config).get("SavePDFs", 0).asInt();
    if (SavePDFs == 1) /* call python function from here... */;
    
    
    return 0;

}
