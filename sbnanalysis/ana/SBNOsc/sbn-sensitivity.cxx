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
#include <fstream>

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>

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
    //std::cout << std::endl << "SBND's POT was " << scalelist[0] << ", MicroBooNE's was " << scalelist[1] << " and ICARUS' was " << scalelist[2] << std::endl;
    
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
    
    
    //// Get covariances
    //// ~~~~~~~~~~~~~~~
    
    ana::SBNOsc::Covariance cov(samples);
    
    // Write to file
    std::string directory = "/sbnd/data/users/gavarela/selection/new/cov_output/";
    
    TFile* covfile = TFile::Open((directory + "cov.root").c_str(), "recreate");
    assert(covfile && covfile->IsOpen());
    
    cov.covmat->Write();
    cov.fcovmat->Write();
    cov.corrmat->Write();
    
    // Save plots
    TCanvas *canvas = new TCanvas();
    cov.covmat->Draw("colz"); cov.covmat->SetStats(kFALSE); canvas->SaveAs((directory + "cov_plot.pdf").c_str());
    cov.fcovmat->Draw("colz"); cov.fcovmat->SetStats(kFALSE); canvas->SaveAs((directory + "fcov_plot.pdf").c_str());
    cov.corrmat->Draw("colz"); cov.corrmat->SetStats(kFALSE); canvas->SaveAs((directory + "corr_plot.pdf").c_str());
    
    
    //// Get sensitivity contours
    //// ~~~~~~~~~~~~~~~~~~~~~~~~
    
    ana::SBNOsc::Chi2Sensitivity chi2(cov);
    
    // Wite to file
    TFile* chi2file = TFile::Open((directory + "chi2.root").c_str(), "recreate");
    assert(chi2file && chi2file->IsOpen());
    
    chi2.contour_90pct->Write();
    chi2.contour_3sigma->Write();
    chi2.contour_5sigma->Write();
    
    // Save plot
    TCanvas *contour_canvas = new TCanvas();
    
    std::vector <int> colours = {30, 38, 46};
    std::vector <TGraph*> contour_graphs = {chi2.contour_90pct, chi2.contour_3sigma, chi2.contour_5sigma};
    
    for (int g = 0; g < contour_graphs.size(); g++) {
        
        contour_graphs[g]->SetMarkerStyle(20);
        contour_graphs[g]->SetMarkerSize(0.1);
        contour_graphs[g]->SetMarkerColor(colours[g]);
        contour_graphs[g]->SetLineColor(colours[g]);
    
    }
    
    TGraph *range = new TGraph();
    range->SetPoint(0, 0.001, 0.01);
    range->SetPoint(1, 1, 100);
    range->SetMarkerColor(0);
    
    TGraph *gr_bestfit = new TGraph();
    gr_bestfit->SetPoint(0, 0.062, 1.7);
    gr_bestfit->SetMarkerStyle(29);
    gr_bestfit->SetMarkerSize(1.6);
    gr_bestfit->SetMarkerColor(40);
    
    contour_graphs[0]->SetTitle("SBN Sensitivity; sin^{2}(2#theta); #Delta m^{2} (eV^{2})");
    
    TLegend *legend = new TLegend();
    legend->AddEntry(contour_graphs[0], "90% CL", "l");
    legend->AddEntry(contour_graphs[1], "3#sigma CL", "l");
    legend->AddEntry(contour_graphs[2], "5#sigma CL", "l");
    legend->AddEntry(gr_bestfit, "Best Fit Point", "p");
    
    contour_canvas->SetLogy();
    contour_canvas->SetLogx();
    
    contour_graphs[0]->Draw("AP");
    contour_graphs[0]->GetXaxis()->SetRangeUser(0.001, 1);
    contour_graphs[0]->GetYaxis()->SetRangeUser(0.01, 100);
    
    contour_graphs[1]->Draw("P same");
    contour_graphs[2]->Draw("P same");
    legend->Draw();
    range->Draw("P same");
    gr_bestfit->Draw("P same");
    
    contour_canvas->SaveAs("PNGs/Sensitivity.pdf");
    

    return 0;

}