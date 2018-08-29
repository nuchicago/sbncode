#include "Covariance.h"
#include "Chi2Sensitivity.h"

#include <TMatrixD.h>
#include <TDecompLU.h>
#include <TMatrixDSym.h>

#include <TF1.h>
#include <TGraph.h>

namespace ana {
  namespace SBNOsc {

Chi2Sensitivity::Chi2Sensitivity(std::vector<EventSample> samples) {
    
    Covariance cov(samples);
    
    // ...
    
}
      
// Personal preference (more explicit about what is used)...
Chi2Sensitivity::Chi2Sensitivity(Covariance cov) {
    
    /*
        Inputs: 
        From the (Covariance::Covariance) cov, we'll use the following arributes:
            - cov         , the covariance matrix as a TH2D*;
            - CV_counts   , a TH1D* containing the CV counts;
            - energies    , a vector <double> with the energies corresponding 
                            to the bin centers in both hists above;
            - sample_order, the order of samples plotted in the two hists above;
            - sample_bins , the bins in the hists above that separate the samples 
                            described in sample_order;
    */
    
    //// Invert Error Matrix
    //// ~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Inverting full error matrix, E_{ij}..." << std::endl;
    
    // Create error (statistical and systematic) matrix
    TMatrixDSym E_mat(cov.covmat->GetNbinsX());
    
    for (int i = 0; i < cov.covmat->GetNbinsX(); i++) {
        for (int j = 0; j < cov.covmat->GetNbinsY(); j++) {
            
            E_mat[i][j] = cov.covmat->GetBinContent(i+1, j+1);
            if (i == j) { E_mat[i][i] += cov.CV_counts->GetBinContent(i+1); }
            
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
    for (std::string sample : cov.sample_order) {
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
    for (int s = 0; s < cov.sample_order.size(); s++) {
        
        // Set dist for the sample
        double tempdist = 0;
        if (cov.sample_order[s].find("SBND") != std::string::npos) {
            tempdist = SBND_dist;
        } else if (cov.sample_order[s].find("MicroBooNE") != std::string::npos) {
            tempdist = MicroBooNE_dist;
        } else if (cov.sample_order[s].find("ICARUS") != std::string::npos) {
            tempdist = ICARUS_dist;
        }
        
        // Add to distance vector
        for (int d = cov.sample_bins[s]; d < cov.sample_bins[s+1]; d++) {
            distance.push_back(tempdist);
        }
        
    }
    
    std::cout << std::endl << "cov.sample_bins: ";
    for (int i = 0; i < cov.sample_bins.size(); i++) {
        std::cout << cov.sample_bins[i] << ", ";
    }
    std::cout << std::endl;
    
    // Should we oscillate this index/bin? 0 = no, 1 = numu, 2 = nue.
    std::vector <int> oscillate(cov.CV_counts->GetNbinsX(), 0);
    for (int i = 0; i < cov.sample_order.size(); i++) {
        
        if (cov.sample_order[i].find("#nu_{#mu}") != std::string::npos) {
            for (int j = cov.sample_bins[i]; j < cov.sample_bins[i+1]; j++) {
                oscillate[j] = 1;
            }
        } else if (cov.sample_order[i].find("#nu_{e}") != std::string::npos) {
            for (int j = cov.sample_bins[i]; j < cov.sample_bins[i+1]; j++) {
                oscillate[j] = 2;
            }
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
            for (int k = 0; k < cov.CV_counts->GetNbinsX(); k++) {
                for (int l = 0; l < cov.CV_counts->GetNbinsX(); l++) {
                    
                    if ((E_inv[k][l] != 0) && (oscillate[k] == 1 && oscillate[l] == 1)) {
                        
                        chisq[i][j] += (cov.CV_counts->GetBinContent(k+1) * (1 - numu_to_numu(distance[k]/cov.energies[k])));
                        
                        chisq[i][j] *= E_inv[k][l];
                            
                        chisq[i][j] *= (cov.CV_counts->GetBinContent(l+1) * (1 - numu_to_numu(distance[l]/cov.energies[l])));
                        
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
       Note: I'm just copying my own code I'd written to get contours. I've heard that there is some function on TH2Ds that does it automatically. This one is pretty fast though so I'm leaving it in. I could try both options/methods... 
    */
    
    clock_t startcont = clock();
    std::cout << std::endl << "Getting contours." << std::endl;
    
    // Get contours
    double target;
    std::vector <double> target_dchisq = {1.64, 7.75, 23.40}, corners(4, 0), twozeros(2, 0);
    std::vector <std::vector <double> > contour, sin_contour(target_dchisq.size()), dm2_contour(target_dchisq.size()), vecof2vecs(np, twozeros);
    std::vector <std::vector <std::vector <double> > > box_minmax(np, vecof2vecs);
    for  (int k = 0; k < target_dchisq.size(); k++){
        
        target = target_dchisq[k];
        
        // Initialise box_minmax
        for (int i = 0; i < np-1; i++) {
            for (int j = 0; j < np-1; j++) {
                box_minmax[i][j] = {0, 0};
            }
        }
        
        // Fill box_minmax out
        for (int i = 0; i < np-1; i++) {
            for (int j = 0; j < np-1; j++) {
                
                corners = {chisq_diffs[i][j], chisq_diffs[i+1][j], chisq_diffs[i][j+1], chisq_diffs[i+1][j+1]};
                box_minmax[i][j] = {corners[0], corners[0]};
                for (int l = 1; l < 4; l++) {
                    if (corners[l] > box_minmax[i][j][1]) {
                        box_minmax[i][j][1] = corners[l];
                    } else if (corners[l] < box_minmax[i][j][0]) {
                        box_minmax[i][j][0] = corners[l];
                    }
                }
                
            }
        }
        
        // Get the contour
        contour.clear();
        for (int i = 0; i < np-1; i++) {
            for (int j = 0; j < np-1; j++) {
                
                if ((target >= box_minmax[i][j][0]) && (target <= box_minmax[i][j][1])) {
                    contour.push_back({(sin2theta[i] + sin2theta[i+1])/2, (dm2[j] + dm2[j+1])/2});
                }
            
            }
        }
        
        // Save the contour
        for (int j = 0; j < contour.size(); j++) {
            sin_contour[k].push_back(contour[j][0]);
            dm2_contour[k].push_back(contour[j][1]);
        }
        
    }
    
    clock_t endcont = clock();
    clock_t tickscont = endcont - startcont;                    // in n of ticks
    double timecont = tickscont / (double) CLOCKS_PER_SEC;      // make into secs
    
    std::cout << "   Done in " << timecont << "s." << std::endl;
    
    std::cout << "sin contour sizes: " << sin_contour[0].size() << ", " << sin_contour[1].size() << ", " << sin_contour[2].size() << std::endl;
    std::cout << "dm2 contour sizes: " << dm2_contour[0].size() << ", " << dm2_contour[1].size() << ", " << dm2_contour[2].size() << std::endl;
    
    // Create TGraphs
    std::vector <TGraph*> graphs;
    for (int i = 0; i < sin_contour.size(); i++) {
        
        graphs.push_back(new TGraph());
        
        for (int j = 0; j < sin_contour[i].size(); j++) {
            graphs[i]->SetPoint(j, sin_contour[i][j], dm2_contour[i][j]);
        }
        
    }
    
    
    //// Output stuff
    //// ~~~~~~~~~~~~
    
    contour_90pct = graphs[0];
    contour_3sigma = graphs[1];
    contour_5sigma = graphs[2];
    
    
}

  }  // namespace SBNOsc
}  // namespace ana
