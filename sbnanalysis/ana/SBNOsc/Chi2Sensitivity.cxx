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

// Oscillate Event counts
std::vector<double> Chi2Sensitivity::EventSample::Oscillate(double sinth, double dm2) const {
    // return histogram binned by reco energy
    std::vector<double> ret(fBkgCounts->GetNbinsX(), 0);
    // iterate over signal bins
    for (unsigned Ebin = 0; Ebin < fSignalCounts->GetNbinsY(); Ebin++) {
        for (unsigned Dbin = 0; Dbin < fSignalCounts->GetNbinsZ(); Dbin++) {
            double energy = (fTrueEBins[Ebin] + fTrueEBins[Ebin+1]) / 2.;
            double distance = (fDistBins[Dbin] + fDistBins[Dbin+1]) / 2.;
            for (unsigned bin = 0; bin < fSignalCounts->GetNbinsX(); bin++) {
                double counts = fSignalCounts->GetBinContent(bin+1, Ebin+1, Dbin+1); 

                if (fOscType == 1) counts *= numu_to_nue(distance/energy, sinth, dm2); 
                else if (fOscType == 2) counts *= numu_to_numu(distance/energy, sinth, dm2);

                ret[bin] += counts;  
            }
        }
    }

    return ret;
}
// Signal Event counts
std::vector<double> Chi2Sensitivity::EventSample::Signal() const {
    // return histogram binned by reco energy
    std::vector<double> ret(fBkgCounts->GetNbinsX(), 0);
    // iterate over signal bins
    for (unsigned Ebin = 0; Ebin < fSignalCounts->GetNbinsY(); Ebin++) {
        for (unsigned Dbin = 0; Dbin < fSignalCounts->GetNbinsZ(); Dbin++) {
            for (unsigned bin = 0; bin < fSignalCounts->GetNbinsX(); bin++) {
                double counts = fSignalCounts->GetBinContent(bin+1, Ebin+1, Dbin+1); 
                ret[bin] += counts;  
            }
        }
    }

    return ret;
}

// return vector of background counts
std::vector<double> Chi2Sensitivity::EventSample::Background() const {
    std::vector<double> ret(fBkgCounts->GetNbinsX(), 0);
    for (unsigned bin = 0; bin < fBkgCounts->GetNbinsX(); bin++) {
        ret[bin] = fBkgCounts->GetBinContent(bin+1);
    }
    return ret;
}

// Initialize
void Chi2Sensitivity::Initialize(Json::Value *config) {
    //// Get parameters from config file
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (config != NULL) {

        // initialize covariance
        fCovariance.Initialize(config);
        
        // Output directory
        fOutputFile = (*config)["Sensitivity"].get("OutputFile", "").asString();
        
        // Get energy from covariance
        fEnergyType = (*config)["Sensitivity"].get("EnergyType", "").asString(); 
        
        // Further selection and rejection 'efficiencies'
        fSelectionEfficiency = (*config)["Sensitivity"].get("SelectionEfficiency", 1.0).asDouble();
        fBackgroundRejection = (*config)["Sensitivity"].get("BackgroundRejection", 0.0).asDouble();
        
        // Event Samples
        for (auto const &json: (*config)["EventSamples"]) {
            fEventSamples.emplace_back(json);
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


        // whether to save stuff
        fSavePDFs = (*config)["Sensitivity"].get("SavePDFs", true).asBool();
        fSaveSignal = (*config)["Sensitivity"].get("SaveSignal", true).asBool();
        fSaveBackground = (*config)["Sensitivity"].get("SaveBackground", true).asBool();
        if ((*config)["Sensitivity"].isMember("SaveOscillations") &&
            (*config)["Sensitivity"]["SaveOscillations"].isArray()) {
            for (auto const &osc_pair: (*config)["Sensitivity"]["SaveOscillations"]) {
              fSaveOscillations.push_back({osc_pair[0].asDouble(), osc_pair[1].asDouble()});
            }
        }
    }
    // start at 0th event sample
    fSampleIndex = 0;
}

Chi2Sensitivity::EventSample::EventSample(const Json::Value &config) {
    // scaling stuff
    fName = config.get("name", "").asString();
    fScaleFactor = config.get("scalefactor", 0.).asDouble(); 

    // setup detector stuff
    fDistance = config.get("Distance", "").asDouble();
    fXlim[0] = config["DetX"][0].asDouble();
    fXlim[1] = config["DetX"][1].asDouble();
    fYlim[0] = config["DetY"][0].asDouble();
    fYlim[1] = config["DetY"][1].asDouble();
    fZlim[0] = config["DetZ"][0].asDouble();
    fZlim[1] = config["DetZ"][1].asDouble();

    // oscillation type
    fOscType = config.get("OscType", 0).asInt();

    // bins of things
    //
    // Reco energy
    Json::Value bins = config.get("binlims", "");
    for (auto const& bin: bins) {
        fBins.push_back(bin.asDouble());
    }

    // True energy
    double minTrueE = config["TrueELims"][0].asDouble();
    double maxTrueE = config["TrueELims"][1].asDouble();
    int numTrueEBinLimits = config.get("numTrueEBins", 1).asInt() + 1;
    double trueE_binwidth = (maxTrueE - minTrueE) / numTrueEBinLimits;
    for (int i = 0; i < numTrueEBinLimits; i ++) {
        fTrueEBins.push_back(minTrueE + i * trueE_binwidth);
    }

    // distance
    //
    // Take distance along z-axis as limits
    int numBinsPerMeter = config.get("NumDistanceBinsPerMeter", 1).asInt(); 
    double numBinsPerCM = numBinsPerMeter * 0.01; // unit conversion
    unsigned n_bins = (unsigned) (fZlim[1] - fZlim[0]) / numBinsPerCM; 
    unsigned n_limits = n_bins + 1;
    double dist_binwidth = (fZlim[1] - fZlim[0]) / n_limits;
    for (unsigned i = 0; i < n_limits; i++) {
        fDistBins.push_back(fDistance + (fZlim[0] + i * dist_binwidth) / 100000. /* cm -> km */);
    }

    // setup histograms
    fBkgCounts = new TH1D((fName + " Background").c_str(), (fName + " Background").c_str(), fBins.size(), &fBins[0]);
    fSignalCounts = new TH3D((fName + " Signal").c_str(), (fName + " Signal").c_str(), 
        fBins.size(), &fBins[0], fTrueEBins.size(), &fTrueEBins[0], fDistBins.size(), &fDistBins[0]);

}

void Chi2Sensitivity::FileCleanup(TTree *eventTree) {
    // cleanup for Covariance
    fCovariance.FileCleanup(eventTree);

    // onto the next sample
    fSampleIndex ++;
}
            
void Chi2Sensitivity::ProcessEvent(const Event *event) {
    // have the covariance process the event
    fCovariance.ProcessEvent(event);

    // Iterate over each interaction in the event
    for (int n = 0; n < event->reco.size(); n++) {
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
    
        // check if energy is within bounds
        if (nuE < fEventSamples[fSampleIndex].fBins[0] || nuE > fEventSamples[fSampleIndex].fBins[fEventSamples[fSampleIndex].fBins.size()-1]) {
            continue;
        } 
        else if (nuE < fEventSamples[fSampleIndex].fTrueEBins[0] || 
            nuE > fEventSamples[fSampleIndex].fTrueEBins[fEventSamples[fSampleIndex].fTrueEBins.size()-1]) {
            std::cout << std::endl << "NUE IN RANGE, TRUE E NOT!!!   nuE = " << nuE << " and true_nuE = " << true_nuE << std::endl;
            continue;
        }
    
        // Apply selection (or rejection) efficiencies
        int isCC = event->truth[truth_ind].neutrino.iscc;
        double wgt = isCC*(fSelectionEfficiency) + (1-isCC)*(1 - fBackgroundRejection);
        // and scale weight
        wgt *= fEventSamples[fSampleIndex].fScaleFactor;
    
        // Get distance travelled along z in km
        double dist = fEventSamples[fSampleIndex].fDistance + 
            (event->truth[truth_ind].neutrino.position.Z() - fEventSamples[fSampleIndex].fZlim[0])/100000. /* cm -> km */;
        
        // Get distance travelled (assuming nu started at (x, y, z) = (0, 0, min_det_zdim - det_dist))
        //double dx = (event->truth[truth_ind].neutrino.position.X() - (fDetDims[sample.fDet][0][1] + fDetDims[sample.fDet][0][0])/2) / 100000 /* cm -> km */,
        //dy = (event->truth[truth_ind].neutrino.position.Y() - (fDetDims[sample.fDet][1][1] + fDetDims[sample.fDet][1][0])/2) / 100000 /* cm -> km */,
        //dz = (event->truth[truth_ind].neutrino.position.Z() - fDetDims[sample.fDet][2][0]) / 100000 /* cm -> km */;
        //double dist = TMath::Sqrt( dx*dx + dy*dy + (fDetDists[sample.fDet] + dz)*(fDetDists[sample.fDet] + dz) );

        // fill in hitograms
        //
        // Signal
        if (isCC) {
            fEventSamples[fSampleIndex].fSignalCounts->Fill(nuE, true_nuE, dist, wgt);
        }
        // background
        else {
            fEventSamples[fSampleIndex].fBkgCounts->Fill(nuE, wgt);
        }
    }
}

void Chi2Sensitivity::GetChi2() {
    
    //// Invert Error Matrix
    //// ~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Inverting full error matrix, E_{ij}..." << std::endl;
    
    // Create error (statistical and systematic) matrix
    //
    // Take covariance matrix from Covariance calculator
    TMatrixDSym E_mat = fCovariance.CovarianceMatrix();
    // Invert it
    TMatrixD E_inv = E_mat.Invert();
        // Find a way to stop code if the matrix doesn't invert! 
        // Check error message? See if E_inv exists?
    
    std::cout << "   Inverted." << std::endl;
    
    
    //// Get chi squareds
    //// ~~~~~~~~~~~~~~~~
    
    // Phase space of oscillation parameters
    for (int i = 0; i < fNumSin; i++) {
        sin2theta.push_back(TMath::Power(10, fLogSinLims[0] + i*(fLogSinLims[1] - fLogSinLims[0])/(fNumSin-1)));
    }
    for (int j = 0; j < fNumDm2; j++) {
        dm2.push_back(TMath::Power(10, fLogDm2Lims[0] + j*(fLogDm2Lims[1] - fLogDm2Lims[0])/(fNumDm2-1)));
    }

    // Loop over phase space calculating Chisq
    clock_t startchi = clock();
    std::cout << std::endl << "Calculating chi squareds..." << std::endl;
    
    double minchisq = 1e99;
    std::vector <double> chisq_dm2(fNumDm2, 0);
    std::vector <std::vector <double> > chisq(fNumSin, chisq_dm2);
    for (int j = 0; j < fNumDm2; j++) {
        std::vector<double> signal; // non-oscillated signal
        std::vector<std::vector<double>> oscillated; // oscilalted signal at sin2==1

        // push all event samples into the vector
        for (auto const &sample: fEventSamples) {
            // oscillate each event sample w/ sin == 1
            oscillated.push_back(sample.Oscillate(1., dm2[j]));

            // flatten the signal
            std::vector<double> this_signal = sample.Signal();
            signal.insert(signal.end(), this_signal.begin(), this_signal.end());
        }

        // oscillated signal for general sin2 (set in loop below)
        std::vector<double> sin2_oscillated(signal.size(), 0);

        for (int i = 0; i < fNumSin; i++) {
            // Scale the oscillated samples based on the sin2 value and flatten them
            unsigned count_ind = 0;
            for (unsigned sample_ind = 0; sample_ind < oscillated.size(); sample_ind++) {
                for (unsigned bin_ind = 0; bin_ind < oscillated[sample_ind].size(); bin_ind++) {
                    // No oscillation
                    if (fEventSamples[sample_ind].fOscType == 0) {
                        sin2_oscillated[count_ind] = signal[count_ind];
                    }
                    // numu -> nue
                    else if (fEventSamples[sample_ind].fOscType == 1) {
                        sin2_oscillated[count_ind] = sin2theta[i] * oscillated[sample_ind][bin_ind];
                    }
                    // numu -> numu
                    else if (fEventSamples[sample_ind].fOscType == 2) {
                        sin2_oscillated[count_ind] = signal[count_ind] - 
                            sin2theta[i] * (signal[count_ind] - oscillated[sample_ind][bin_ind]);
                    }

                    // update count ind
                    count_ind ++;
                } 
            }
            assert(count_ind == signal.size());

            // Calculate chisq for sin2th = 1
            for (int k = 0; k < signal.size(); k++) {
                for (int l = 0; l < signal.size(); l++) {
                    chisq[i][j] += (signal[k] - sin2_oscillated[k]) * (signal[l] - sin2_oscillated[l]) * E_inv[k][l];
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
    
    
    // Gen TGraph2D for output
    chisqplot = new TGraph2D();
    for (int i = 0; i < fNumSin; i++) {
        for (int j = 0; j < fNumDm2; j++) {
            chisqplot->SetPoint(i*fNumSin + j, TMath::Log10(sin2theta[i]), TMath::Log10(dm2[j]), chisq[i][j]);
        }
    }
    
    // Get differences
    chisq_diffs = chisq;
    for (int i = 0; i < fNumSin; i++) {
        for (int j = 0; j < fNumDm2; j++) {
            chisq_diffs[i][j] -= minchisq;
        }
    }
    
}

void Chi2Sensitivity::GetContours() {
    
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

void Chi2Sensitivity::Write() {
    
    // Wite to file
    TFile* chi2file = TFile::Open(fOutputFile.c_str(), "recreate");
    assert(chi2file && chi2file->IsOpen());
    
    if (fSavePDFs) {
	contour_90pct->SetName("90pct"); contour_90pct->Write();
	contour_3sigma->SetName("3sigma"); contour_3sigma->Write();
	contour_5sigma->SetName("5sigma"); contour_5sigma->Write();
	
	chisqplot->SetName("chisq"); chisqplot->Write();
    }

    // save histos
    if (fSaveBackground) {
        for (auto const &sample: fEventSamples) {
            sample.fBkgCounts->Write();
        }
    }

    if (fSaveSignal) {
        for (auto const &sample: fEventSamples) {
            sample.fSignalCounts->Write();
        }
    }

    for (auto const &osc_vals: fSaveOscillations) {
        for (auto const &sample: fEventSamples) {
            std::string name = sample.fName + " sin2th: " + std::to_string(osc_vals[0]) + " dm2: " + std::to_string(osc_vals[1]);
            TH1D *hist = new TH1D(name.c_str(), name.c_str(), sample.fBins.size(), &sample.fBins[0]);
            std::vector<double> oscillated = sample.Oscillate(osc_vals[0], osc_vals[1]);
            for (unsigned i = 0; i < oscillated.size(); i++) {
                hist->SetBinContent(i+1, oscillated[i]);
            }
            hist->Write();
        }
    }

    chi2file->Close();
    
}

}  // namespace SBNOsc
}  // namespace ana

DECLARE_SBN_POSTPROCESSOR(ana::SBNOsc::Chi2Sensitivity);

