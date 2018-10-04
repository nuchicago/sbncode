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
        From (Covariance::Covariance) cov, we'll use the following arributes:
            - cov         , the covariance matrix as a TH2D*;
            - bkg_counts  , a TH1D* containing the CV counts of non-oscillating particles
            - nu_counts   , a TH3D* containing the CV counts of oscillating particles, by
                            true E (x-axis), reconstructed E (y) and distance from target (z)
            - CV_counts   , a TH1D* containing the CV counts;
            - energies    , a vector <double> with the energies corresponding 
                            to the bin centers in the three count hists above;
            - sample_order, a vector <string> with the samples plotted in the count hists above;
            - sample_bins , a vector <int> with the bins in the count hists above that separate the
                            samples described in sample_order;
            - trueEs      , a vector containing the true energy value for each bin of nu_counts;
            - sample_dist_bins, a vector containing the distance bins separating different samples;
            - dist_bins   , a vector containing the distance value for each bin of nu_counts.
                            
        From the config file:
            - fNumDm2     , number of points in dm2 dimension of (dm2, sin2(2theta)) phase space,
            - fLogDm2Lims , limits of the dm2 dimension of the phase space, in log units,
            - fNumSin     , number of points in sin dimension of phase space,
            - fLogSinLims , limits of the sin dimension of the phase space, in log units.
    */
    
    //// Get parameters from config file
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Json::Value* config = core::LoadConfig(configFileName);
    
    if (config != NULL) {
        
        // Output directory
        fOutputDirectory = (*config).get("OutputDirectory", "./").asString();
        
        // Phase space parameters
        fNumDm2 = (*config)["Sensitivity"].get("NumDm2", -1).asInt();
        fLogDm2Lims = {};
        for (auto binlim : (*config)["Sensitivity"]["LogDm2Lims"]) {
            fLogDm2Lims.push_back(binlim.asDouble());
        }
        fNumSin = (*config)["Sensitivity"].get("NumSin", -1).asInt();
        fLogSinLims = {};
        for (auto binlim : (*config)["Sensitivity"]["LogSinLims"]) {
            fLogSinLims.push_back(binlim.asDouble());
        }
        
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
    /*
    double SBND_dist = 0.1, MicroBooNE_dist = 0.47, ICARUS_dist = 0.6; // all in km <- make configurable?
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
    */
    
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
    std::vector <double> sin2theta, dm2;
    for (int i = 0; i < fNumSin; i++) {
        sin2theta.push_back(TMath::Power(10, fLogSinLims[0] + i*(fLogSinLims[1] - fLogSinLims[0])/(fNumSin-1)));
    }
    for (int j = 0; j < fNumDm2; j++) {
        dm2.push_back(TMath::Power(10, fLogDm2Lims[0] + j*(fLogDm2Lims[1] - fLogDm2Lims[0])/(fNumDm2-1)));
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
    std::vector <double> chisq_builder(fNumDm2, 0);
    std::vector <std::vector <double> > chisq(fNumSin, chisq_builder);
    for (int i = 0; i < fNumSin; i++){
        for (int j = 0; j < fNumDm2; j++) {
            
            // Progress counter
            if (j == 0) {
                std::cout << "\r\r\r\r\r\r\r\r\r\r\r";
                int tempi = (int)( (float)(i-1)/fNumSin*100 );
                while (tempi > 0) { tempi /= 10; std::cout << "\r"; }
                std::cout << "Progress: " << (int)( (float)i/fNumSin*100 ) << "%";
                if (i == fNumSin-1) std::cout << "\r\r\r\r\r\r\r\r\r\r\r\r" << "Progress: 100%" << std::endl;
            }
            
            // Set function parameters
            numu_to_numu.SetParameters(sin2theta[i], dm2[j]);
            
            // Create and fill hist to hold oscillated counts
            TH1D *osc_counts = new TH1D("temposc", "", cov.bkg_counts->GetNbinsX(), 0, cov.bkg_counts->GetNbinsX());
            
            for (int o = 0; o < cov.sample_order.size(); o++) { // For limits on bin loops inside:
                
                for (int rb = cov.sample_bins[o]; rb < cov.sample_bins[o+1]; rb++) {
                    
                    double dosc_counts_rb = 0;
                    for (int tb = o*num_trueE_bins; tb < (o+1)*num_trueE_bins; tb++) {
                        for (int db = cov.sample_dist_bins[o]; db < cov.sample_dist)bins[o+1]; db++) {

                            // Numus
                            if (oscillate[rb] == 1) {
                                dosc_counts_rb += cov.nu_counts->GetBinContent(1+tb, 1+rb, 1+db) * numu_to_numu(cov.dist_bins[db]/cov.trueEs[tb]);
                            // Nues
                            } else if (oscillate[rb] == 2) {
                                // For the future...
                            }
                            
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
                TH1D *scale_hist = new TH1D("tempscale", "", cov.bkg_counts->GetNbinsX(), 0, cov.bkg_counts->GetNbinsX());
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
    
    
    // Gen TGraph2D for output
    chisqplot = new TGraph2D();
    for (int i = 0; i < fNumSin; i++) {
        for (int j = 0; j < fNumDm2; j++) {
            chisqplot->SetPoint(i*fNumSin + j, TMath::Log10(sin2theta[i]), TMath::Log10(dm2[j]), chisq[i][j]);
        }
    }
    
    // Get differences
    std::vector <std::vector <double> > chisq_diffs = chisq;
    for (int i = 0; i < fNumSin; i++) {
        for (int j = 0; j < fNumDm2; j++) {
            chisq_diffs[i][j] -= minchisq;
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
    std::vector <std::vector <double> > contour, sin_contour(target_dchisq.size()), dm2_contour(target_dchisq.size()), vecof2vecs(fNumDm2, twozeros);
    std::vector <std::vector <std::vector <double> > > box_minmax(fNumSin, vecof2vecs);
    for  (int k = 0; k < target_dchisq.size(); k++){
        
        target = target_dchisq[k];
        
        // Initialise box_minmax
        for (int i = 0; i < fNumSin-1; i++) {
            for (int j = 0; j < fNumDm2-1; j++) {
                box_minmax[i][j] = {0, 0};
            }
        }
        
        // Fill box_minmax out
        for (int i = 0; i < fNumSin-1; i++) {
            for (int j = 0; j < fNumDm2-1; j++) {
                
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
        for (int i = 0; i < fNumSin-1; i++) {
            for (int j = 0; j < fNumDm2-1; j++) {
                
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
    
    //std::cout << "sin contour sizes: " << sin_contour[0].size() << ", " << sin_contour[1].size() << ", " << sin_contour[2].size() << std::endl;
    //std::cout << "dm2 contour sizes: " << dm2_contour[0].size() << ", " << dm2_contour[1].size() << ", " << dm2_contour[2].size() << std::endl;
    
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
