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

#include <TChain.h>
#include <TFileCollection.h>
#include <TTreePlayer.h>

double getPOTs(std::string filelist) {
    
    /* Stolen from Swapnil */
    
    TChain chain("SubRuns");
    TFileCollection f("dum");
    f.AddFromFile(filelist.c_str());
    chain.AddFileInfoList(f.GetList());
    ((TTreePlayer*)(chain.GetPlayer()))->SetScanRedirect(true);
    ((TTreePlayer*)(chain.GetPlayer()))->SetScanFileName("total_pot.list");
    chain.Scan("sumdata::POTSummary_generator__GenieGen.obj.totpot");
    
    
    std::ifstream input("total_pot.list");
    double total = 0.0;
    
    std::string line;
    while(std::getline(input,line)) {
        
        if (line.empty()) continue;
        std::istringstream is(line);
        
        std::string dim1, dim2, dim3;
        int row;
        double pot;
        if (is >> dim1 >> row >> dim2 >> pot >> dim3) { total += pot; }
        
    }
    
    //remove((std::string)"total_pot.list");
    
    //std::cout << "Total Exposure in ICARUS = " << total << std::endl;
    return total;
    
}

int main(int argc, char* argv[]) {

    std::cout << std::endl << "Hello!" << std::endl << std::endl;
    
    //// Build sample (temporary)
    //// ~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::vector <std::string> DETLIST = {"SBND", "MicroBooNE", "ICARUS"},
                              detlist = {"sbnd", "uboone", "icarus"},
                             desclist = {"numu", "numu", "numu"};
    
    //std::string prelist = "/sbnd/data/users/gavarela/selection/", postlist = "/spatel_output/spatel.list";
    //std::vector <float> scalelist = {getPOTs(prelist+"sbnd"+postlist), getPOTs(prelist+"uboone"+postlist), getPOTs(prelist+"icarus"+postlist)};
    
    std::vector <float> scalelist = {3.0958e18, 8.87435e19, 6.59165e18};
    std::vector <ana::SBNOsc::EventSample> samples;
    
    std::vector <TFile*> tfiles;
    for (int d = 0; d < detlist.size(); d++) {
        
        std::string det = detlist[d], DET = DETLIST[d], desc = desclist[d];
        float scale = scalelist[d];
        
        tfiles.push_back(new TFile(((std::string)"/sbnd/data/users/gavarela/selection/new/output_" + DET + (std::string)"_new.root").c_str()));
        
        samples.push_back(ana::SBNOsc::EventSample(tfiles[d], (TTree*)tfiles[d]->Get("sbnana"), scale, det, desc));
        
    }
    
    assert(!samples.empty());
    
    
    //// Get covariances and contours
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ana::SBNOsc::Covariance cov(samples);
    
    ana::SBNOsc::Chi2Sensitivity chi2(cov.covmat, cov.CV_counts, cov.sample_order, cov.sample_bins);

    // Write out contours to a file
    
    

    return 0;

}