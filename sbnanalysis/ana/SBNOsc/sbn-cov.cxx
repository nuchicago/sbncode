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

    /* Andy's sample-building code...
    std::vector<ana::SBNOsc::EventSample> samples;
    for (int i=2; i<argc; i++) {
    // Build sample list
    }
    */
    
    std::cout << std::endl << "Hello!" << std::endl << std::endl;
    
    // My sample-building code
    std::vector <std::string> DETLIST = {"SBND", "MicroBooNE", "ICARUS"},
        detlist = {"sbnd", "uboone", "icarus"},
        desclist = {"nu", "nu", "nu"};
    std::vector <float> scalelist = {3.0958e18, 8.87435e19, 6.59165e18};
    std::vector <ana::SBNOsc::EventSample> samples;
    
    for (int d = 0; d < detlist.size(); d++) {
        
        std::string det = detlist[d], DET = DETLIST[d], desc = desclist[d];
        float scale = scalelist[d];
        
        TFile f(((std::string)"/sbnd/data/users/gavarela/selection/new/output_" + DET + (std::string)".root").c_str());
        
        samples.push_back(ana::SBNOsc::EventSample((TTree*)f.Get("sbnana"), scale, det, desc));
        
    }
    
    assert(!samples.empty());
    
    std::cout << std::endl << "Doing cov now:" << std::endl << std::endl;
    
    ana::SBNOsc::Covariance cov(samples);
    
    std::cout << std::endl << "Did cov." << std::endl << std::endl;
    
    // Write matrix out to a ROOT file...
    
    std::string directory = "/sbnd/data/users/gavarela/selection/new/cov_output";
    
    cov.covmat->Write((directory + (std::string)"/cov_output_cov.root").c_str());
    cov.fcovmat->Write((directory + (std::string)"/cov_output_fcov.root").c_str());
    cov.corrmat->Write((directory + (std::string)"/cov_output_corr.root").c_str());
    
    TCanvas *canvas = new TCanvas();
    cov.covmat->Draw("colz"); canvas->Write((directory + (std::string)"/cov_plot.png").c_str());
    cov.fcovmat->Draw("colz"); canvas->Write((directory + (std::string)"/fcov_plot.png").c_str());
    cov.corrmat->Draw("colz"); canvas->Write((directory + (std::string)"/corr_plot.png").c_str());
    
    return 0;
    
}