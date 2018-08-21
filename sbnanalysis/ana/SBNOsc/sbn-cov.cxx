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
    
    std::vector <TFile*> tfiles;
    for (int d = 0; d < detlist.size(); d++) {
        
        std::string det = detlist[d], DET = DETLIST[d], desc = desclist[d];
        float scale = scalelist[d];
        
        tfiles.push_back(new TFile(((std::string)"/sbnd/data/users/gavarela/selection/new/output_" + DET + (std::string)".root").c_str()));
        
        samples.push_back(ana::SBNOsc::EventSample(tfiles[d], (TTree*)tfiles[d]->Get("sbnana"), scale, det, desc));
        
    }
    
    std::cout << std::endl << "Done. samples.size() = " << samples.size() << " with samples: " << std::endl;
    for (int s = 0; s < samples.size(); s++) {
        
//        Event *event = new Event;
//        samples[s].tree->SetBranchAddress("events", &event);
//        for (int e = 0; e < tree->GetEntries())
        
        std::cout << std::endl << "(" << samples[s].sDet << " " << samples[s].sDesc << ")" << std::endl;
    }
    
    assert(!samples.empty());
    
    std::cout << std::endl << "Doing cov now:" << std::endl << std::endl;
    
    ana::SBNOsc::Covariance cov(samples);
    
    std::cout << std::endl << "Did cov." << std::endl << std::endl;
    
    // Write matrix out to a ROOT file...
    
    std::string directory = "/sbnd/data/users/gavarela/selection/new/cov_output/";
    TFile* newfile = TFile::Open((directory + "cov_output.root").c_str(), "recreate");
    assert(newfile && newfile->IsOpen());
    
    cov.covmat->Write();
    cov.fcovmat->Write();
    cov.corrmat->Write();
    
    TCanvas *canvas = new TCanvas();
    cov.covmat->Draw("colz"); canvas->SaveAs((directory + "cov_plot.pdf").c_str());
    cov.fcovmat->Draw("colz"); canvas->SaveAs((directory + "fcov_plot.pdf").c_str());
    cov.corrmat->Draw("colz"); canvas->SaveAs((directory + "corr_plot.pdf").c_str());
    
    return 0;
    
}