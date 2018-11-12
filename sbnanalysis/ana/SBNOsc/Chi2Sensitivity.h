#ifndef __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__
#define __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__

/**
 * \file Chi2Sensitivity.h
 */

#include "Covariance.h"

#include <vector>
#include <string>

#include <TGraph2D.h>
#include <TGraph.h>

namespace ana {
namespace SBNOsc {

class Chi2Sensitivity {
    
    public:
        
        // Functions
        
        Chi2Sensitivity(std::vector <EventSample> samples, char *configFileName);
        Chi2Sensitivity(std::vector <EventSample> samples, Covariance cov, char *configFileName);
        
        void ScanEvents(), GetChi2(), GetContours(), Write(std::string directory);
        
        // Output
        
        TGraph2D *chisqplot;
        TGraph *contour_90pct, *contour_3sigma, *contour_5sigma;
    
    private:
        
        // From config file
        
        std::string fEnergyType;
        
        double fSelectionEfficiency, fRejectionEfficiency;
        
        int fNumDistBinsPerMeter;
        std::map <std::string, float> fDetDists;
        std::map <std::string, std::vector <std::vector <double > > > fDetDims;
        
        std::vector <double> fTrueELims;
        int fNumTrueEBins;
        
        std::map <std::string, float> fScaleTargets;
        
        int fNumDm2, fNumSin;
        std::vector <double> fLogDm2Lims, fLogSinLims;
        
        int fShapeOnly;
        
        std::string fOutputDirectory;
        int fSavePDFs;
        
        // Internal
        
        Covariance covar;
        std::vector <EventSample> ev_samples;
        
        int num_bins;
        std::vector <int> sample_bins;
    
        int num_dist_bins;
        std::vector <double> dist_bins, sample_dist_bins;
    
        std::vector <double> trueEs;
        
        std::vector <std::vector <double> > chisq_diffs;
    
};

}  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__
