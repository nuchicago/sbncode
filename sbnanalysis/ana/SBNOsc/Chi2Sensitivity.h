#ifndef __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__
#define __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__

/**
 * \file Chi2Sensitivity.h
 */

#include "Covariance.h"

#include <TGraph.h>

namespace ana {
namespace SBNOsc {

class Chi2Sensitivity {
    
    public:
        Chi2Sensitivity(std::vector<EventSample> samples);
        Chi2Sensitivity(Covariance cov);
        
        TGraph *contour_90pct, *contour_3sigma, *contour_5sigma;
    
};

}  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__