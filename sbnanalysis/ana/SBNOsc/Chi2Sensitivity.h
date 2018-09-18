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
        Chi2Sensitivity(std::vector<EventSample> samples, char *configFileName);
        Chi2Sensitivity(Covariance cov, char *configFileName);
        
        TGraph2D *chisqplot;
        TGraph *contour_90pct, *contour_3sigma, *contour_5sigma;
    
    private:
        
        std::string fScaleSample;
        int fNumDm2, fNumSin;
        std::vector <double> fLogDm2Lims, fLogSinLims;
    
        std::string fOutputDirectory;
        int fSavePDFs;
    
};

}  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__
