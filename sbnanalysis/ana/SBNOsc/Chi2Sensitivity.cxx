#include "Covariance.h"
#include "Chi2Sensitivity.h"

#include <TMatrixD.h>
#include <TDecompLU.h>
#include <TMatrixDSym.h>

#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TGraph2D.h>

#include <TFile.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THStack.h>
#include <TROOT.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCanvas.h>
#include <core/Event.hh>

namespace ana {
  namespace SBNOsc {

// Oscillation function
double numu_to_numu(double x, double sin, double dm2) {
    
    // x is L/E, sin is sin^2(2theta) and dm2 is delta m^2
    return 1 - sin * TMath::Power(TMath::Sin(1.27 * dm2 * x), 2);
    
}

double numu_to_nue(double x, double sin, double dm2) {
    
    // x is L/E, sin is sin^2(2theta) and dm2 is delta m^2
    return sin * TMath::Power(TMath::Sin(1.27 * dm2 * x), 2);
}

Chi2Sensitivity::Chi2Sensitivity(std::vector<EventSample> samples, char *configFileName) {
    
    Covariance cov(samples, configFileName);
    
    // ...
    
}
      
// Personal preference (more explicit about what is used)...
Chi2Sensitivity::Chi2Sensitivity(std::vector<EventSample> samples, Covariance cov, char *configFileName) {
    
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
        
        // Type of energy
        fEnergyType = (*config)["Covariance"].get("EnergyType", "").asString();
        
        // Further selection and rejection 'efficiencies'
        fSelectionEfficiency = (*config)["Covariance"].get("SelectionEfficiency", -1e99).asDouble();
        fRejectionEfficiency = (*config)["Covariance"].get("RejectionEfficiency", -1e99).asDouble();
        
        // True energy binning
        fTrueELims = {};
        for (auto binlim : (*config)["Covariance"]["TrueELims"]) {
            fTrueELims.push_back(binlim.asDouble());
        }
        fNumTrueEBins = (*config)["Covariance"].get("NumTrueEBins", -1).asInt();
        
        // Detector dimensions and distances
        fNumDistBinsPerMeter = (*config)["Covariance"].get("NumDistanceBinsPerMeter", -1).asInt();
        fDetDists = {}; fDetDims = {};
        for (auto det : (*config)["Covariance"]["DetectorDimensions"]) {
            
            std::string detname = det["Detector"].asString();
            
            std::vector <std::vector <double > > xyz_lims = {{}, {}, {}};
            for (auto lim : det["X"]) xyz_lims[0].push_back(lim.asDouble());
            for (auto lim : det["Y"]) xyz_lims[1].push_back(lim.asDouble());
            for (auto lim : det["Z"]) xyz_lims[2].push_back(lim.asDouble());
            
            fDetDims.insert({detname, xyz_lims});

            float distance = det["Distance"].asFloat();
            fDetDists.insert({detname, distance});
            
        }
        
        // Exposure normalisation
        for (auto sample : samples) {
            if (fScaleTargets.find(sample.fDet) == fScaleTargets.end()) {
                fScaleTargets.insert({sample.fDet, (*config)["Covariance"]["ScaleTargets"].get(sample.fDet, -1).asFloat()});
            }
        }
               
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
        
        // Shape-only chisq?
        fShapeOnly = (*config)["Sensitivity"].get("ShapeOnly", 0).asInt();
    
    }
    
    //// Get neutrino counts from samples
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Some stuff related to binning and plotting
    
        // Number of bins needed in big histogram (for covariance)
        // And bin 'boundaries' between each separate sample
    int num_bins = 0;
    std::vector <int> sample_bins;
    for (auto sample : samples) {
        
        sample_bins.push_back(num_bins);
        num_bins += sample.fBins.size() - 1;
        
    }
    sample_bins.push_back(num_bins);
    
        // Number of distance bins needed in 3D histogram of true CC interactions
    std::map <std::string, std::vector <double> > dist_bin_lims;
    std::map <std::string, int> dist_bin_nums;
    for (auto it : fDetDims) {
        
        double min_dist, max_dist; 
        min_dist = fDetDists[it.first];
 
        float xlen = (it.second[0][1] - it.second[0][0])/100000 /* cm -> km */,
              ylen = (it.second[1][1] - it.second[1][0])/100000 /* cm -> km */,
              zlen = (it.second[2][1] - it.second[2][0])/100000 /* cm -> km */;
        max_dist = TMath::Sqrt( xlen*xlen + ylen*ylen + (min_dist+zlen)*(min_dist+zlen) );

        dist_bin_lims.insert({it.first, {min_dist, max_dist}});
        dist_bin_nums.insert({it.first, (max_dist - min_dist)*(fNumDistBinsPerMeter*1000 /* 1/m -> 1/km */)});
    
    }
    
    int num_dist_bins = 0;
    std::vector <double> dist_bins = {}, sample_dist_bins = {0};
    for (auto sample : samples) {
        
        for (auto it : dist_bin_nums) {
            if (sample.fDet == it.first) {
                
                num_dist_bins += it.second;
                sample_dist_bins.push_back(sample_dist_bins[sample_dist_bins.size()-1] + it.second);
                
                for (int i = 0; i < it.second + 1; i++) {
                    dist_bins.push_back(dist_bin_lims[it.first][0] + i*(dist_bin_lims[it.first][1] - dist_bin_lims[it.first][0])/(it.second));
                }
                
                break;
                
            }
        }
        
    }
    
    // Large (meaningless x-axis) histograms for cov
    TH1D *CV_counts = {new TH1D("CV", "Base Uni. Counts; Bin; Counts", num_bins, 0, num_bins)};
    
    // Large (meaningless x-axis) histograms (one 3d) for oscillation calculations later on
    TH1D *bkg_counts = new TH1D("Background", "Background Counts; Reconstructed Energy Bin; Counts", num_bins, 0, num_bins);
    
    TH3D *nu_counts = new TH3D("Neutrinos", "Neutrino Counts; True Energy Bin; Reconstructed Energy Bin; Distance Bin", fNumTrueEBins*samples.size(), 0, fNumTrueEBins*samples.size(), num_bins, 0, num_bins, num_dist_bins, 0, num_dist_bins);
    
    
    // Get counts
    std::cout << std::endl << "Getting counts for each sample..." << std::endl;
    
    std::vector <double> energies = {}, trueEs = {};;
    int o = 0;
    for (auto sample : samples) {
        
        // What sample did we get?
        std::cout << "Doing " << sample.fDesc << " sample in " << sample.fDet << std::endl;
        
        // Initialise temp hists to store counts
            // Base
        std::string title = sample.fDet+"; Reconstructed Energy (GeV); Counts";
        TH1D *temp_CV_counts = {new TH1D((sample.fDet+"tempCV").c_str(), title.c_str(), sample.fBins.size() - 1, &sample.fBins[0])};
        
            // Bkg
        TH1D *temp_bkg_counts = new TH1D((sample.fDet+"tempbkg").c_str(), "", sample.fBins.size() - 1, &sample.fBins[0]);
        
            // Neutrinos
        std::vector <double> temp_trueE_bins = {};
        for (int i = 0; i < fNumTrueEBins + 1; i++) {
            temp_trueE_bins.push_back(fTrueELims[0] + i*(fTrueELims[1]-fTrueELims[0])/(fNumTrueEBins));\
        }
        
        std::vector <double> temp_dist_bins = {};
        for (int i = 0; i < dist_bin_nums[sample.fDet]+1; i++) {
            temp_dist_bins.push_back(dist_bin_lims[sample.fDet][0] + i*(dist_bin_lims[sample.fDet][1]-dist_bin_lims[sample.fDet][0])/(dist_bin_nums[sample.fDet]));
        }
        
        TH3D *temp_nu_counts = new TH3D((sample.fDet+"tempnu").c_str(), "", fNumTrueEBins, &temp_trueE_bins[0], sample.fBins.size() - 1, &sample.fBins[0], dist_bin_nums[sample.fDet], &temp_dist_bins[0]);
        
        // Loop over neutrinos (events in tree)
        Event *event = new Event;
        sample.tree->SetBranchAddress("events", &event);
        
        int nucount = 0;
        for (int e = 0; e < sample.tree->GetEntries(); e++) {
            
            sample.tree->GetEntry(e);
            
            for (int n = 0; n < event->reco.size(); n++) {
                
                nucount++;
                
                unsigned truth_ind = event->reco[n].truth_index;
                
                // Get energy
                double nuE, true_nuE = event->reco[n].truth.neutrino.energy;
                if (fEnergyType == "CCQE") {
                    nuE = event->truth[truth_ind].neutrino.eccqe;
                } else if (fEnergyType == "True") {
                    nuE = true_nuE;
                } else if (fEnergyType == "Reco") {
                    nuE = event->reco[n].reco_energy;
                }
                
                if (nuE < sample.fBins[0] || nuE > sample.fBins[sample.fBins.size()-1]) {
                    continue;
                } else if (true_nuE < fTrueELims[0] || true_nuE > fTrueELims[1]) {
                    std::cout << std::endl << "NUE IN RANGE, TRUE E NOT!!!   nuE = " << nuE << " and true_nuE = " << true_nuE << std::endl;
                    continue;
                }
                
                // Apply selection (or rejection) efficiencies
                int isCC = event->truth[truth_ind].neutrino.iscc;
                double wgt = isCC*(fSelectionEfficiency) + (1-isCC)*(1 - fRejectionEfficiency);
                
                // Add to base count histogram
                temp_CV_counts->Fill(nuE, wgt);
                
                // Get distance travelled (assuming nu started at (x, y, z) = (0, 0, min_det_zdim - det_dist))
                double dx = (event->truth[truth_ind].neutrino.position.X() - (fDetDims[sample.fDet][0][1] + fDetDims[sample.fDet][0][0])/2) / 100000 /* cm -> km */,
                       dy = (event->truth[truth_ind].neutrino.position.Y() - (fDetDims[sample.fDet][1][1] + fDetDims[sample.fDet][1][0])/2) / 100000 /* cm -> km */,
                       dz = (event->truth[truth_ind].neutrino.position.Z() - fDetDims[sample.fDet][2][0]) / 100000 /* cm -> km */;
                double dist = TMath::Sqrt( dx*dx + dy*dy + (fDetDists[sample.fDet] + dz)*(fDetDists[sample.fDet] + dz) );
                
                // Fill bkg and nu histograms
                if (sample.fIsNu && isCC) {
                    temp_nu_counts->Fill(true_nuE, nuE, dist, wgt);
                } else {
                    temp_bkg_counts->Fill(nuE, wgt);
                }
                
            }
        }
        
        // Rescale to desired POT
        temp_CV_counts->Scale(fScaleTargets[sample.fDet] / sample.fScaleFactor);
        temp_bkg_counts->Scale(fScaleTargets[sample.fDet] / sample.fScaleFactor);
        temp_nu_counts->Scale(fScaleTargets[sample.fDet] / sample.fScaleFactor);
        
        
        // Pass onto the big histograms and get energies
        for (int bin = 0; bin < temp_CV_counts->GetNbinsX(); bin++) {
            
            CV_counts->SetBinContent(1+sample_bins[o]+bin, temp_CV_counts->GetBinContent(1+bin));
            energies.push_back(temp_CV_counts->GetBinCenter(bin+1));
            
        }
        
        for (int rb = 0; rb < temp_nu_counts->GetNbinsY(); rb++) {
            
            // bkg_counts
            bkg_counts->SetBinContent(1+sample_bins[o]+rb, temp_bkg_counts->GetBinContent(1+rb));
            
            // nu_counts
            for (int tb = 0; tb < temp_nu_counts->GetNbinsX(); tb++) {
                for (int db = 0; db < temp_nu_counts->GetNbinsZ(); db++) {
                    
                    nu_counts->SetBinContent(1 + o*fNumTrueEBins + tb, 1 + sample_bins[o] + rb, 
                                             1 + sample_dist_bins[o] + db, 
                                             temp_nu_counts->GetBinContent(1+tb, 1+rb, 1+db));
                    
                }
            }
            
        }
        
        // Fill out a vector with the true energies
        double trueE_binwidth = (fTrueELims[1] - fTrueELims[0])/fNumTrueEBins;
        for (int i = 0; i < fNumTrueEBins; i++) {
            trueEs.push_back(fTrueELims[0] + (i+0.5)*trueE_binwidth);
        }
        
        o++;
        
    }
    
    
    //// Invert Error Matrix
    //// ~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Inverting full error matrix, E_{ij}..." << std::endl;
    
    // Create error (statistical and systematic) matrix
    TMatrixDSym E_mat(cov.cov->GetNbinsX());
    
    for (int i = 0; i < cov.cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov.cov->GetNbinsY(); j++) {
            
            E_mat[i][j] = cov.cov->GetBinContent(i+1, j+1);
            if (i == j) { E_mat[i][i] += CV_counts->GetBinContent(i+1); }
            
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
    for (auto sample : samples) {
        if (sample.fDesc == "#nu_{e}") {
            nue_appearance = 1;
            break;
        }
    }
    
    // Should we oscillate this index/bin? 0 = no, 1 = numu, 2 = nue.
    std::vector <int> oscillate(CV_counts->GetNbinsX(), 0);
    for (int i = 0; i < samples.size(); i++) {
        
        if (samples[i].fDesc == "#nu_{#mu}") {
            for (int j = sample_bins[i]; j < sample_bins[i+1]; j++) {
                oscillate[j] = 1;
            }
        } else if (samples[i].fDesc == "#nu_{e}") {
            for (int j = sample_bins[i]; j < sample_bins[i+1]; j++) {
                oscillate[j] = 2;
            }
        }
        
    }
    
    // Num true energy bins in each sample
    int num_trueE_bins = trueEs.size()/samples.size();
    
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
    TH1D *osc_counts = new TH1D("temposc", "", bkg_counts->GetNbinsX(), 0, bkg_counts->GetNbinsX());
    for (int i = 0; i < fNumSin; i++){
        for (int j = 0; j < fNumDm2; j++) {
    
            std::cout << "i = " << i << ", j = " << j << std::endl;
       
            /*
            // Progress counter
            if (j == 0) {
                std::cout << "\r\r\r\r\r\r\r\r\r\r\r";
                int tempi = (int)( (float)(i-1)/fNumSin*100 );
                while (tempi > 0) { tempi /= 10; std::cout << "\r"; }
                std::cout << "Progress: " << (int)( (float)i/fNumSin*100 ) << "%";
                if (i == fNumSin-1) std::cout << "\r\r\r\r\r\r\r\r\r\r\r\r" << "Progress: 100%" << std::endl;
            }
            */
            
            // Fill hist that holds oscillated counts
            for (int o = 0; o < samples.size(); o++) { // For limits on bin loops inside:
                
                for (int rb = sample_bins[o]; rb < sample_bins[o+1]; rb++) {
                    
                    double dosc_counts_rb = 0;
                    for (int tb = o*num_trueE_bins; tb < (o+1)*num_trueE_bins; tb++) {
                        for (int db = sample_dist_bins[o]; db < sample_dist_bins[o+1]; db++) {

                            // Numus
                            if (oscillate[rb] == 1) {
                                dosc_counts_rb += nu_counts->GetBinContent(1+tb, 1+rb, 1+db) * numu_to_numu(dist_bins[db]/trueEs[tb], sin2theta[i], dm2[j]);
                            // Nues
                            } else if (oscillate[rb] == 2) {
                                // For the future...
                            }
                            
                        }
                    }
                    
                    osc_counts->SetBinContent(1+rb, bkg_counts->GetBinContent(1+rb) + dosc_counts_rb);
                    
                }
                
            }
            
            // If configured to, scale bins based on distribution in given detector for "shape-only" chi2
            if (fShapeOnly == 1) /* checks if configured */ {
                
                int fit_sample_index = 0;
                for (auto sample : samples) {
                    if (sample.fScaleSample == 1) break;
                    fit_sample_index++;
                }
                assert(fit_sample_index != samples.size());
                
                // Histogram of scale factors (only part with scaled sample will be used)
                TH1D *scale_hist = new TH1D("tempscale", "", bkg_counts->GetNbinsX(), 0, bkg_counts->GetNbinsX());
                scale_hist->Divide(CV_counts, osc_counts);
                
                // Update each value in osc_counts histogram
                for (int o = 0; o < samples.size(); o++) {
                    for (int rb = sample_bins[o]; rb < sample_bins[o+1]; rb++) {
                        
                        // Get corresponding bin in the scaled sample histogram
                        int scale_bin = (rb - sample_bins[o]) + sample_bins[fit_sample_index];
                        
                        // Scale accordingly
                        double this_count = osc_counts->GetBinContent(1+rb);
                        double scaled_this_count = scale_hist->GetBinContent(1+scale_bin) * this_count;
                        osc_counts->SetBinContent(1+rb, scaled_this_count);
                        
                    }
                }
                
                scale_hist->Delete();
                
            }
	    
            // Calculate chisq
            for (int k = 0; k < CV_counts->GetNbinsX(); k++) {
                for (int l = 0; l < CV_counts->GetNbinsX(); l++) {
                    
                    double dchisqij = 0;

                    dchisqij += (CV_counts->GetBinContent(k+1) - osc_counts->GetBinContent(k+1));
                    dchisqij *= E_inv[k][l];
                    dchisqij *= (CV_counts->GetBinContent(l+1) - osc_counts->GetBinContent(l+1));

                    chisq[i][j] += dchisqij;

                }
            }
            
            //osc_counts->Delete();
            
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
