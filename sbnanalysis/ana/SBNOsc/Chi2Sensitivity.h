#ifndef __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__
#define __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__

/**
 * \file Chi2Sensitivity.h
 */

#include "Covariance.h"

namespace ana {
namespace SBNOsc {

class Chi2Sensitivity {
    
    public:
        Chi2Sensitivity(std::vector<EventSample> samples);
        Chi2Sensitivity(TH2D* cov, TH1D *counts, std::vector <std::string> sample_order, std::vector <int> sample_bins);
        
        // TGraph
    
    
};

}  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__