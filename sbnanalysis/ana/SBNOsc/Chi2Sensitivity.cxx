#include "Covariance.h"
#include "Chi2Sensitivity.h"

#include <TMatrixD.h>
#include <TDecompLU.h>
#include <TMatrixDSym.h>

namespace ana {
  namespace SBNOsc {

Chi2Sensitivity::Chi2Sensitivity(std::vector<EventSample> samples) {
    
    Covariance cov(samples);
    
    // ...
    
}
      
// Personal preference (more explicit about what is used)...
Chi2Sensitivity::Chi2Sensitivity(TH2D* cov, TH1D *counts, std::vector <std::string> sample_order, std::vector <int> sample_bins) {
    
    /*
    Inputs: - cov   , the covariance matrix as a TH2D*
            - counts, a TH1D* containing the CV counts
    
    */
    
    //// Invert Error Matrix
    //// ~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Inverting full error matrix, E_{ij}..." << std::endl;
    
    // Create error (statistical and systematic) matrix
    TMatrixDSym E_mat(cov.GetNbinsX());
    
    for (int i = 0; i < cov.GetNbinsX(); i++) {
        for (int j = 0; j < cov.GetNbinsY(); j++) {
            
            E_mat[i][j] = cov.GetBinContent(i+1, j+1);
            if (i == j) { E_mat[i][i] += counts.GetBinContent(i+1); }
            
        }
    }
    
    // Invert it
    TMatrixD E_inv = E_mat.Invert();
        // Find a way to stop code if the matrix doesn't invert! 
        // Check error message? See if E_inv exists?
    
    std::cout << "   Inverted." << std::endl;
    
    
    //// Get chi squareds
    //// ~~~~~~~~~~~~~~~~
    
    /*
       Note: If the samples given contain nues, the program will automatically assume that we're looking at nue appearance. Else, numu disappearance.
       In both cases, the numu sample is the oscillated one, though with different interpretations of the angle theta.
    */
    
    // Nue appearance or numu disappearance?
    int nue_appearance = 0;
    for (std::string sample : sample_order) {
        if (sample.find("#nu_{e}") != std::string::npos) {
            nue_appearance = 1;
            break;
        }
    }
    
    // Oscillation function
    TF1 numu_to_numu("numu_to_numu", "1 - [0] * (sin(1.27 * [1] * x))^2", 0, 25);
        // [0] is sin^2(2theta), [1] is delta m^2, x is L/E
    TF1 numu_to_nue("numu_to_nue", "[0] * (sin(1.27 * [1] * x))^2", 0, 25);
        // [0] is sin^2(2theta) {diff theta!}, [1] is delta m^2 {same}, x is L/E
    
    // Distances
    double SBND_dist = 0.1, MicroBooNE_dist = 0.47, ICARUS_dist = 0.6; // all in km
    std::vector <double> distance;
    for (int s = 0; s < sample_order.size(); s++) {
        
        // Set dist for the sample
        double tempdist = 0;
        if (std::find(sample_order[s].begin(), sample_order.end(), "SBND") != std::string::npos) {
            tempdist = SBND_dist;
        } else if (std::find(sample_order[s].begin(), sample_order.end(), "MicroBooNE") != std::string::npos) {
            tempdist = MicroBooNE_dist;
        } else if (std::find(sample_order[s].begin(), sample_order.end(), "SBND") != std::string::npos) {
            tempdist = ICARUS_dist;
        }
        
        // Add to distance vector
        for (int d = sample_bins[s]; d < sample_bins[s+1]; d++) {
            distance.push_back(tempdist);
        }
        
    }
    
    // Should we oscillate this index/bin? 0 = no, 1 = numu, 2 = nue.
    std::vector <int> oscillate(counts.GetNbinsX(), 0);
    for (int i = 0; i < oscillate.size(); i++) {
        if (std::find(sample_order[i].begin(), sample_order[i].end(), "#nu_{#mu}") != std::string::npos) {
            oscillate[i] = 1;
        } else if (std::find(sample_order[i].begin(), sample_order[i].end(), "#nu_{e}") != std::string::npos) {
            oscillate[i] = 2;
        }
    }
    
    // Phase space
    int np = 500;
    std::vector <double> dm2(np), sin2theta(np);
    for (int i = 0; i < np; i++) {
        dm2[i] = TMath::Power(10, -2.0 + i*4.0/(np-1));
        sin2theta[i] = TMath::Power(10, -3.0 + i*3.0/(np-1));
    }
    
    // Loop over phase space calculating Chisq
    clock_t startchi = clock();
    std::cout << std::endl << "Calculating chi squareds..." << std::endl;
    
    double minchisq = 1e99, fardetected_osc;
    std::vector <double> npzeros(np, 0);
    std::vector <std::vector <double> > chisq(np, npzeros);
    for (int i = 0; i < np; i++){
        for (int j = 0; j < np; j++) {
        
            // Set function parameters
            numu_to_numu.SetParameters(sin2theta[i], dm2[j]);
            
            // Find null and oscillation fluxes and detections and calculate chisq
            for (int k = 0; k < counts[0].size(); k++) {
                for (int l = 0; l < counts[0].size(); l++) {
                    
                    if ((E_inv[k][l] != 0) && (oscillate[k] == 1 && oscillate[l] == 1)) {
                        chisq[i][j] += (counts[0][k] * (1 - numu_to_numu(distance[k]/energy[k]))) * E_inv[k][l] * (counts[0][l] * (1 - nutonu(distance[l]/energy[l])));
                    }
                }
                
            }
            
            // Check if min chisq
            if (chisq[i][j] < minchisq) {
                minchisq = chisq[i][j];
            }
        
        }
    }
    
    clock_t endchi = clock();
    clock_t tickschi = endchi - startchi;                    // in n of ticks
    double timechi = tickschi / (double) CLOCKS_PER_SEC;     // make into secs
    
    std::cout << "   Done in " << timechi << "s. " << std::endl;
    
    /*
    // Plot
    TCanvas *chisqcanvas = new TCanvas();
    
    TGraph2D *logchisqplot = new TGraph2D();
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < np; j++) {
            logchisqplot->SetPoint(i*np + j, TMath::Log10(sin2theta[i]), TMath::Log10(dm2[j]), chisq[i][j]);
        }
    }
    
    logchisqplot->SetTitle("#chi^{2}; log_{10}(sin^{2}(2#theta)); log_{10}(#Delta m^{2}); #chi^{2}");
    gStyle->SetPalette(1);
    logchisqplot->Draw("surf1");
    chisqcanvas->SaveAs("PNGs/chisq.png");
    */
    
    // Get differences
    std::vector <std::vector <double> > chisq_diffs(np, npzeros);
    for (int i = 0; i < chisq.size(); i++) {
        for (int j = 0; j < chisq[0].size(); j++) {
            chisq_diffs[i][j] = chisq[i][j] - minchisq;
        }
    }
    
    
    //// Get contours
    //// ~~~~~~~~~~~~
    
    /*
       Note: I'm just copying my own code I'd written to get contours. I've heard that there is some function on TH2Ds that does it automatically. This one is pretty fast though so I'm leaving it in. I could try both options/methods... I'm also kind of proud of this one...
       
    */
    
    
    
    
}

  }  // namespace SBNOsc
}  // namespace ana
