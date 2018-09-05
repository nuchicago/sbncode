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
        
        std::cout << desc << " " << det << std::endl;
        
        samples.push_back(ana::SBNOsc::EventSample(file, scalefactor, det, desc));
        
    }
    
    assert(!samples.empty());
    
    
    //// Get covariances
    //// ~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Starting covariance procedure..." << std::endl << std::endl;
    
    ana::SBNOsc::Covariance cov(samples, configFileName);
    
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
    
    std::cout << std::endl << "Starting chisq/sensitivity procedure..." << std::endl << std::endl;
    
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
        contour_graphs[g]->SetMarkerSize(0.25);
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
    
    range->SetTitle("SBN Sensitivity (Reconstructed Energy); sin^{2}(2#theta); #Delta m^{2} (eV^{2})");
    
    TLegend *legend = new TLegend();
    legend->AddEntry(contour_graphs[0], "90% CL", "l");
    legend->AddEntry(contour_graphs[1], "3#sigma CL", "l");
    legend->AddEntry(contour_graphs[2], "5#sigma CL", "l");
    legend->AddEntry(gr_bestfit, "Best Fit Point", "p");
    
    contour_canvas->SetLogy();
    contour_canvas->SetLogx();
    
    range->Draw("AP");
    range->GetXaxis()->SetRangeUser(0.001, 1);
    range->GetYaxis()->SetRangeUser(0.01, 100);
    
    contour_graphs[0]->Draw("P same");
    contour_graphs[1]->Draw("P same");
    contour_graphs[2]->Draw("P same");
    
    legend->Draw();
    gr_bestfit->Draw("P same");
    
    contour_canvas->SaveAs((directory + "Sensitivity.pdf").c_str());
    

    return 0;

}