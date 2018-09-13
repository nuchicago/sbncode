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
    
    // Write to file
    std::string directory = (*config).get("OutputDirectory", "./").asString();
    
    TFile* covfile = TFile::Open((directory + "cov.root").c_str(), "recreate");
    assert(covfile && covfile->IsOpen());
    
    cov.cov->Write();
    cov.fcov->Write();
    cov.corr->Write();
    
    // Save plots
    int savePDFs = (*config).get("SavePDFs", 0).asInt();
    if (savePDFs == 1) {
        
        TCanvas *canvas = new TCanvas();
        gStyle->SetPadLeftMargin(0.15); gStyle->SetPadRightMargin(0.15);
        
        cov.cov->GetYaxis()->SetLabelSize(0.07); cov.cov->GetXaxis()->SetLabelSize(0.07);
        cov.cov->Draw("colz"); cov.cov->SetStats(kFALSE); canvas->SaveAs((directory + "cov_plot.pdf").c_str());
        
        cov.fcov->GetYaxis()->SetLabelSize(0.07); cov.fcov->GetXaxis()->SetLabelSize(0.07);
        cov.fcov->Draw("colz"); cov.fcov->SetStats(kFALSE); canvas->SaveAs((directory + "fcov_plot.pdf").c_str());
        
        cov.corr->GetYaxis()->SetLabelSize(0.07); cov.corr->GetXaxis()->SetLabelSize(0.07);
        cov.corr->Draw("colz"); cov.corr->SetStats(kFALSE); canvas->SaveAs((directory + "corr_plot.pdf").c_str());
        
    }
    
    
    //// Get sensitivity contours
    //// ~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Starting chisq/sensitivity procedure..." << std::endl << std::endl;
    
    ana::SBNOsc::Chi2Sensitivity chi2(cov, configFileName);
    
    // Wite to file
    TFile* chi2file = TFile::Open((directory + "chi2.root").c_str(), "recreate");
    assert(chi2file && chi2file->IsOpen());
    
    chi2.contour_90pct->Write();
    chi2.contour_3sigma->Write();
    chi2.contour_5sigma->Write();
    
    // Plot chi squareds
    TCanvas *chisqcanvas = new TCanvas();
    
    chi2.chisqplot->SetTitle("#chi^{2}; log_{10}(sin^{2}(2#theta)); log_{10}(#Delta m^{2}); #chi^{2}");
    gStyle->SetPalette(1);
    chi2.chisqplot->Draw("surf1");
    if (savePDFs == 1) chisqcanvas->SaveAs((directory + "chisq.pdf").c_str());
    
    // Plot contours
    TCanvas *contour_canvas = new TCanvas("cont_canvas", "", 1020, 990);
    
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
    
    range->SetTitle("SBN Sensitivity; sin^{2}(2#theta); #Delta m^{2} (eV^{2})");
    
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
    
    if (savePDFs == 1) contour_canvas->SaveAs((directory + "Sensitivity.pdf").c_str());
    
    // Save as .root
    TFile* contourfile = TFile::Open((directory + "contours.root").c_str(), "recreate");
    assert(contourfile && contourfile->IsOpen());
    
    TGraph *contour_90pct = contour_graphs[0], 
          *contour_3sigma = contour_graphs[1],
          *contour_5sigma = contour_graphs[2];
    contour_90pct->Write();
    contour_3sigma->Write();
    contour_5sigma->Write();
    
    chi2.chisqplot->Write();
    
    
    
    return 0;

}
