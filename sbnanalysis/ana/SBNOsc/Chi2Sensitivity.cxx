#include "Covariance.h"
#include "Chi2Sensitivity.h"

#include <TMatrixD.h>
#include <TDecompLU.h>
#include <TMatrixDSym.h>

#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TGraph2D.h>


namespace ana {
  namespace SBNOsc {

Chi2Sensitivity::Chi2Sensitivity(std::vector<EventSample> samples, char *configFileName) {
    
    Covariance cov(samples, configFileName);
    
    // ...
    
}
      
// Personal preference (more explicit about what is used)...
Chi2Sensitivity::Chi2Sensitivity(Covariance cov, char *configFileName) {
    
    /*
        Inputs: 
        From the (Covariance::Covariance) cov, we'll use the following arributes:
            - cov         , the covariance matrix as a TH2D*;
            - bkg_counts  , a TH1D* containing the CV counts of non-oscillating particles
            - nu_counts   , a TH2D* containing the CV counts of oscillating particles, with
                            true E bins on the x-axis and reconstructed E bins on the y-axis
            - CV_counts   , a TH1D* containing the CV counts;
            - energies    , a vector <double> with the energies corresponding 
                            to the bin centers in the three count hists above;
            - sample_order, a vector <string> with the samples plotted in the count hists above;
            - sample_bins , a vector <int> with the bins in the count hists above that separate the
                            samples described in sample_order.
                            
        From the config file:
            - fNP         , number of points in each dimension of the (dm2, sin2(2theta)) phase space.
    */
    
    //// Get parameters from config file
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Json::Value* config = core::LoadConfig(configFileName);
    
    if (config != NULL) {
        
        // Output directory
        fOutputDirectory = (*config).get("OutputDirectory", "./").asString();
        fSavePDFs = (*config).get("SavePDFs", 0).asInt();
        
        // Size of phase space
        fNP = (*config)["Sensitivity"].get("NP", -1).asInt();
        
        // Sample according to which we scale (for shape-only chi squared)
        fScaleSample = (*config)["Sensitivity"].get("ScaleSample", "").asString();
    
    }
    
    
    //// Invert Error Matrix
    //// ~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Inverting full error matrix, E_{ij}..." << std::endl;
    
    // Create error (statistical and systematic) matrix
    TMatrixDSym E_mat(cov.cov->GetNbinsX());
    
    for (int i = 0; i < cov.cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov.cov->GetNbinsY(); j++) {
            
            E_mat[i][j] = cov.cov->GetBinContent(i+1, j+1);
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
    
    // Num true energy bins in each sample
    int num_trueE_bins = cov.trueEs.size()/cov.sample_order.size();
    
    // Phase space
    std::vector <double> dm2(fNP), sin2theta(fNP);
    for (int i = 0; i < fNP; i++) {
        dm2[i] = TMath::Power(10, -2.0 + i*4.0/(fNP-1));
        sin2theta[i] = TMath::Power(10, -3.0 + i*3.0/(fNP-1));
    }
    
    
    
    //// TO-DO:
    
    // Transfer matrices b/w near and far detectors
    // TMatrixD ... 
    
    //
    //
    //   generate transfer matrices
    //
    //
    
    
    
    
    // Loop over phase space calculating Chisq
    clock_t startchi = clock();
    std::cout << std::endl << "Calculating chi squareds..." << std::endl;
    
    double minchisq = 1e99;
    std::vector <double> npzeros(fNP, 0);
    std::vector <std::vector <double> > chisq(fNP, npzeros);
    for (int i = 0; i < fNP; i++){
        for (int j = 0; j < fNP; j++) {
        
            if (j == 0) {
                std::cout << "\r\r\r\r\r\r\r\r\r\r\r";
                int tempi = (int)( (float)i/fNP*100 );
                while (tempi > 0) { tempi /= 10; std::cout << "\r"; }
                std::cout << "Progress: " << (int)( (float)i/fNP*100 ) << "%";
                if (i == fNP-1) std::cout << "\r\r\r\r\r\r\r\r\r\r\r\r" << "Progress: 100%" << std::endl;
            }
            
            // Set function parameters
            numu_to_numu.SetParameters(sin2theta[i], dm2[j]);
            
            // Create and fill hist to hold oscillated counts
            TH1D *osc_counts = new TH1D("temp", "", cov.bkg_counts->GetNbinsX(), 0, cov.bkg_counts->GetNbinsX());
            
            for (int o = 0; o < cov.sample_order.size(); o++) { // For limits on bin loops inside:
                
                for (int rb = cov.sample_bins[o]; rb < cov.sample_bins[o+1]; rb++) {
                    
                    double dosc_counts_rb = 0;
                    for (int tb = o*num_trueE_bins; tb < (o+1)*num_trueE_bins; tb++) {

                        // Numus
                        if (oscillate[rb] == 1) {
                            dosc_counts_rb += cov.nu_counts->GetBinContent(1+tb, 1+rb) * numu_to_numu(distance[rb]/cov.trueEs[tb]);
                        // Nues
                        } else if (oscillate[rb] == 2) {
                            // For the future...
                        }
                        
                    }
                    
                    osc_counts->SetBinContent(1+rb, cov.bkg_counts->GetBinContent(1+rb) + dosc_counts_rb);
                    
                }
                
            }
            
            // If configured to, scale bins based on distribution in given detector for "shape-only" chi2
            if (fScaleSample.size() != 0) /* checks if configured */ {
                
                int fit_sample_index = std::find(cov.sample_order.begin(), cov.sample_order.end(), fScaleSample) - cov.sample_order.begin(); 
                assert(fit_sample_index != cov.sample_order.size());
                
                // Histogram of scale factors (only part with scaled sample will be used)
                TH1D *scale_hist = new TH1D("temp", "", cov.bkg_counts->GetNbinsX(), 0, cov.bkg_counts->GetNbinsX()); 
                scale_hist->Divide(cov.CV_counts, osc_counts);
                
                // Update each value in osc_counts histogram
                for (int o = 0; o < cov.sample_order.size(); o++) {
                    for (int rb = cov.sample_bins[o]; rb < cov.sample_bins[o+1]; rb++) {
                        
                        // Get corresponding bin in the scaled sample histogram
                        int scale_bin = (rb - cov.sample_bins[o]) + cov.sample_bins[fit_sample_index];
                        
                        // Scale accordingly
                        double this_count = osc_counts->GetBinContent(1+rb);
                        double scaled_this_count = scale_hist->GetBinContent(1+scale_bin) * this_count;
                        osc_counts->SetBinContent(1+rb, scaled_this_count);
                        
                    }
                }
                
                scale_hist->Delete();
                
            }
	    
            // Calculate chisq
            for (int k = 0; k < cov.CV_counts->GetNbinsX(); k++) {
                for (int l = 0; l < cov.CV_counts->GetNbinsX(); l++) {
                    
                    double dchisqij = 0;

                    dchisqij += (cov.CV_counts->GetBinContent(k+1) - osc_counts->GetBinContent(k+1));
                    dchisqij *= E_inv[k][l];
                    dchisqij *= (cov.CV_counts->GetBinContent(l+1) - osc_counts->GetBinContent(l+1));

                    chisq[i][j] += dchisqij;

                }
            }
            
            osc_counts->Delete();
            
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
    std::cout << "   Min chisq = " << minchisq << std::endl;
    
    
    // Plot
    TCanvas *chisqcanvas = new TCanvas();
    
    TGraph2D *logchisqplot = new TGraph2D();
    for (int i = 0; i < fNP; i++) {
        for (int j = 0; j < fNP; j++) {
            logchisqplot->SetPoint(i*fNP + j, TMath::Log10(sin2theta[i]), TMath::Log10(dm2[j]), chisq[i][j]);
        }
    }
    
    logchisqplot->SetTitle("#chi^{2}; log_{10}(sin^{2}(2#theta)); log_{10}(#Delta m^{2}); #chi^{2}");
    //gStyle->SetPalette(1);
    logchisqplot->Draw("surf1");
    if (fSavePDFs == 1) chisqcanvas->SaveAs((fOutputDirectory + "chisq.pdf").c_str());
    
    
    // Get differences
    std::vector <std::vector <double> > chisq_diffs(fNP, npzeros);
    for (int i = 0; i < chisq.size(); i++) {
        for (int j = 0; j < chisq[0].size(); j++) {
            chisq_diffs[i][j] = chisq[i][j] - minchisq;
        }
    }
    
    
    //// Get contours
    //// ~~~~~~~~~~~~
    
    /*
       Note: I'm just copying my own code I'd written to get contours. I've heard that there is some 
       function on TH2Ds that does it automatically but haven't used it. This one is pretty fast so 
       I'm leaving it in, at least for now.
    */
    
    clock_t startcont = clock();
    std::cout << std::endl << "Getting contours..." << std::endl;
    
    // Get contours
    double target;
    std::vector <double> target_dchisq = {1.64, 7.75, 23.40}, corners(4, 0), twozeros(2, 0);
    std::vector <std::vector <double> > contour, sin_contour(target_dchisq.size()), dm2_contour(target_dchisq.size()), vecof2vecs(fNP, twozeros);
    std::vector <std::vector <std::vector <double> > > box_minmax(fNP, vecof2vecs);
    for  (int k = 0; k < target_dchisq.size(); k++){
        
        target = target_dchisq[k];
        
        // Initialise box_minmax
        for (int i = 0; i < fNP-1; i++) {
            for (int j = 0; j < fNP-1; j++) {
                box_minmax[i][j] = {0, 0};
            }
        }
        
        // Fill box_minmax out
        for (int i = 0; i < fNP-1; i++) {
            for (int j = 0; j < fNP-1; j++) {
                
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
        for (int i = 0; i < fNP-1; i++) {
            for (int j = 0; j < fNP-1; j++) {
                
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
