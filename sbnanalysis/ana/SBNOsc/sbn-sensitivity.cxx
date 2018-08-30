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

#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <core/Event.hh>
#include <TFile.h>
#include <TGraph2D.h>
#include <TMatrixD.h>
#include <TDecompLU.h>
#include <TMatrixDSym.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>

int main(int argc, char* argv[]) {

    std::cout << std::endl << "Hello!" << std::endl << std::endl;

    /* From scratch */
    
    std::string dir = "/sbnd/data/users/gavarela/selection/new/";

    /* Get counts */

    int n_unis = 100;
    Double_t bins[] = { 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5, 3 };
    Int_t nbins = sizeof(bins)/sizeof(Double_t) - 1;

    std::vector <std::string> dets = {"SBND", "MicroBooNE", "ICARUS"};
    std::vector <double> scalefactor = {6.6e20/3.0958e18, 1.32e21/8.87435e19, 6.6e20/6.59165e18};

    std::vector <TH1D*> hists = {new TH1D("base", "Base; GeV;", 3*nbins, 0, 3*nbins)}, basehists;
    for (int i = 0; i < n_unis; i++) {
        hists.push_back(new TH1D(("alt" + std::to_string(i+1)).c_str(), "Alt; GeV;", 3*nbins, 0, 3*nbins));
    }

    std::vector <double> energies;

    for (int d = 0; d < dets.size(); d++) {

        std::string det = dets[d];

        TFile f((dir + "output_" + det + "_new.root").c_str());
        TTree *tree = (TTree*) f.Get("sbnana");

        Event *ev = new Event;
        tree->SetBranchAddress("events", &ev);

        std::vector <TH1D*> temphists = {new TH1D((det+"tempbase").c_str(), (det + "; GeV;").c_str(), nbins, bins)};
        for (int i = 0; i < n_unis; i++) {
            temphists.push_back(new TH1D((det+"tempalt"+std::to_string(i+1)).c_str(), "Alt; GeV;", nbins, bins));
        }

        for (int e = 0; e < tree->GetEntries(); e++) {

            tree->GetEntry(e);
            if (ev->reco.size() != ev->truth.size()) { continue; }

            for (int n = 0; n < ev->reco.size(); n++) {

                temphists[0]->Fill(ev->reco[n].neutrino.energy);
                
                std::cout << "Det " << det << " and neutrino " << e << std::endl;
                
                for (int u = 0; u < n_unis; u++) {

                    double weight = 1;
                    int wind;
                    for (auto it : ev->truth[n].weights) {
                        wind = u;
                        while (wind >= it.second.size()) { wind -= it.second.size(); }
                        weight *= it.second.at(wind);
                    }
                    
                    temphists[u+1]->Fill(ev->reco[n].neutrino.energy, weight);

                }

            }
            
        }

        for (int h = 0; h < temphists.size(); h++) {
            for (int b = 0; b < temphists[h]->GetNbinsX(); b++) {

                int offset = 0 + nbins*(det == "MicroBooNE") + 2*nbins*(det == "ICARUS");

                hists[h]->SetBinContent(offset+b+1, temphists[h]->GetBinContent(b+1) * scalefactor[d] / temphists[h]->GetBinWidth(b+1));
            
                if (h == 0) {

                    temphists[h]->SetBinContent(b+1, temphists[h]->GetBinContent(b+1) * scalefactor[d] / temphists[h]->GetBinWidth(b+1));

                    energies.push_back(temphists[h]->GetBinCenter(b+1));

                }
                
            }
        
            basehists.push_back(temphists[0]);
        
        }

    }

    /* Plot base */
    /*
    std::cout << "Size of basehists " << basehists.size() << std::endl;
    
    TCanvas *basec = new TCanvas();
    basec->Divide(3, 1);
    for (int h = 0; h < basehists.size(); h++) {
        
        std::cout << "Doing basehist " << h << std::endl;
        
        basec->cd(h+1);
        //basehists[h]->Draw("hist");

    }

    //basec->SaveAs((dir+"test/basecounts.png").c_str());
    */
    
    TCanvas *can = new TCanvas();
    hists[0]->Draw("hist");
    can->SaveAs((dir+"test/basecounts.png").c_str());
    
    /* Get cov, fcov and corr */

    TH2D *cov = new TH2D("cov", "Covariance Matrix", nbins, 0, nbins, nbins, 0, nbins),
         *fcov = new TH2D("fcov", "Fractional Covariance Matrix", nbins, 0, nbins, nbins, 0, nbins);

    for (int i = 0; i < cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov->GetNbinsY(); j++) {

            double covij = 0;
            for (int u = 0; u < n_unis; u++) {
                covij += (hists[0]->GetBinContent(i+1) - hists[u]->GetBinContent(i+1)) * 
                         (hists[0]->GetBinContent(j+1) - hists[u]->GetBinContent(j+1));
            }
            covij /= n_unis;
            cov->SetBinContent(i+1, j+1, covij);

            double fcovij = covij / (hists[0]->GetBinContent(i+1) * hists[0]->GetBinContent(j+1));
            fcov->SetBinContent(i+1, j+1, fcovij);

        }
    }

    TH2D *corr = new TH2D("corr", "Correlation Matrix", nbins, 0, nbins, nbins, 0, nbins);
    for (int i = 0; i < cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov->GetNbinsY(); j++) {

            double corrij = cov->GetBinContent(i+1, j+1) / TMath::Sqrt(cov->GetBinContent(i+1, i+1) * cov->GetBinContent(j+1, j+1));
            corr->SetBinContent(i+1, j+1, corrij);

        }
    }


    /* Plot */

    TCanvas *c = new TCanvas();
    cov->Draw("colz"); c->SaveAs((dir + "test/cov.png").c_str());
    fcov->Draw("colz"); c->SaveAs((dir+"test/fcov.png").c_str());
    corr->Draw("colz"); c->SaveAs((dir+"test/corr.png").c_str());


    /* Invert Error */

    TMatrixDSym E_mat(cov->GetNbinsX());

    for (int i = 0; i < cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov->GetNbinsY(); j++) {

            E_mat[i][j] = cov->GetBinContent(i+1, j+1);
            if (i == j) { E_mat[i][i] += hists[0]->GetBinContent(i+1); }

        }
    }

    TMatrixD E_inv = E_mat.Invert();

    /* Get chisq */

    TF1 numu_to_numu("numu_to_numu", "1 - [0] * (sin(1.27 * [1] * x))^2", 0, 25);

    std::vector <double> distance, detdist = {0.1, 0.47, 0.6};
    for (int d = 0; d < detdist.size(); d++){
        for (int i = 0; i < nbins; i++){
            distance.push_back(detdist[d]);
        }
    }

    int np = 500;
    std::vector <double> dm2(np), sin2theta(np);
    for (int i = 0; i < np; i++) {
        dm2[i] = TMath::Power(10, -2.0 + i*4.0/(np-1));
        sin2theta[i] = TMath::Power(10, -3.0 + i*3.0/(np-1));
    }

    clock_t startchi = clock();
    std::cout << std::endl << "Calculating chi squareds..." << std::endl;

    double minchisq = 1e99;
    std::vector <double> npzeros(np, 0);
    std::vector <std::vector <double> > chisq(np, npzeros);
    for (int i = 0; i < np; i++){
        for (int j = 0; j < np; j++) {

            numu_to_numu.SetParameters(sin2theta[i], dm2[j]);

            for (int k = 0; k < hists[0]->GetNbinsX(); k++) {
                for (int l = 0; l < hists[0]->GetNbinsX(); l++) {

                    if (E_inv[k][l] != 0) {

                        chisq[i][j] += (hists[0]->GetBinContent(k+1) * (1 - numu_to_numu(distance[k]/energies[k])));

                        chisq[i][j] *= E_inv[k][l];

                        chisq[i][j] *= (hists[0]->GetBinContent(l+1) * (1 - numu_to_numu(distance[l]/energies[l])));

                    }
                }

            }

            if (chisq[i][j] < minchisq) {
                minchisq = chisq[i][j];
            }

        }
    }

    clock_t endchi = clock();
    clock_t tickschi = endchi - startchi;                    // in n of ticks
    double timechi = tickschi / (double) CLOCKS_PER_SEC;     // make into secs

    std::cout << "   Done in " << timechi << "s. " << std::endl;

    std::vector <std::vector <double> > chisq_diffs(np, npzeros);
    for (int i = 0; i < chisq.size(); i++) {
        for (int j = 0; j < chisq[0].size(); j++) {
            chisq_diffs[i][j] = chisq[i][j] - minchisq;
        }
    }


    /* Plot */

    TCanvas *chisqcanvas = new TCanvas();

    TGraph2D *logchisqplot = new TGraph2D();
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < np; j++) {
            logchisqplot->SetPoint(i*np + j, TMath::Log10(sin2theta[i]), TMath::Log10(dm2[j]), chisq[i][j]);
        }
    }

    logchisqplot->SetTitle("#chi^{2}; log_{10}(sin^{2}(2#theta)); log_{10}(#Delta m^{2}); #chi^{2}");
    logchisqplot->Draw("surf1");
    chisqcanvas->SaveAs((dir+"test/chisq.png").c_str());


    /* Get contour */
    




    return 0;

}