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
    
    std::vector <ana::SBNOsc::EventSample> samples;
    
    for (auto sample : (*config)["EventSamples"]) {
        
        TFile *file = new TFile((sample["path"].asString()).c_str());
        float scalefactor = sample["scalefactor"].asFloat();
        std::string det = sample["det"].asString(),
                        desc = sample["desc"].asString();
        
        samples.push_back(ana::SBNOsc::EventSample(file, scalefactor, det, desc));
        
    }
    
    assert(!samples.empty());
    
    
    //// Get covariances
    //// ~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Starting covariance procedure..." << std::endl << std::endl;
    
    ana::SBNOsc::Covariance cov(samples, configFileName);
    
    // Write matrices to file
    std::string directory = (*config).get("OutputDirectory", "./").asString();
    
    TFile* covfile = TFile::Open((directory + "cov.root").c_str(), "recreate");
    assert(covfile && covfile->IsOpen());
    
    cov.cov->Write();
    cov.fcov->Write();
    cov.corr->Write();
    
    /*
    // Write counts to file
    TFile* countfile = TFile::Open((directory + "counts.root").c_str(), "recreate");
    assert(countfile && countfile->IsOpen());
    
    std::vector <std::vector <TH1D*> > forloop = {numu_counts, numu_bkgs};
    if (cov.nue_counts.size() > 0) { forloop.push_back(nue_counts); forloop.pushback(nue_bkgs); }
    
    for (int v = 0; v < forloop.size(); v++) {
        
        std::string toappend = "_numu";
        if (v > 1) toappend = "_nue";
        hist->SetName(((std::string)hist->GetName + "_numu").c_str());
        hist->Write();
        
    }
    cov.numu_counts->Write(); cov.numu_bkgs->Write();
    if (cov.nue_counts.size() > 0) { cov.nue_counts->Write(); cov.nue_bkgs->Write(); }
    */
    
    //// Get sensitivity contours
    //// ~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Starting chisq/sensitivity procedure..." << std::endl << std::endl;
    
    ana::SBNOsc::Chi2Sensitivity chi2(cov, configFileName);
    
    // Wite to file
    TFile* chi2file = TFile::Open((directory + "chi2.root").c_str(), "recreate");
    assert(chi2file && chi2file->IsOpen());
    
    chi2.contour_90pct->SetName("90pct"); chi2.contour_90pct->Write();
    chi2.contour_3sigma->SetName("3sigma"); chi2.contour_3sigma->Write();
    chi2.contour_5sigma->SetName("5sigma"); chi2.contour_5sigma->Write();
    
    chi2.chisqplot->SetName("chisq"); chi2.chisqplot->Write();
    
    
    //// Save plots
    int SavePDFs = (*config).get("SavePDFs", 0).asInt();
    if (SavePDFs == 1) /* call python function from here... */;
    
    
    return 0;

}
