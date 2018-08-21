#ifndef __sbnanalysis_ana_SBNOsc_Covariance__
#define __sbnanalysis_ana_SBNOsc_Covariance__

/**
 * \file Covariance.h
 */

#include <string>
#include <vector>
#include <iostream>
#include <cassert>

#include <TH2D.h>

class TTree;

namespace ana {
namespace SBNOsc {

class EventSample {
    
    public:
        /** Constructors. */
        EventSample(TTree* _tree, float scaleFactor) : tree(_tree), fScaleFactor(scaleFactor) {}
        EventSample(std::vector<std::string> filenames, float fScaleFactor);
        
        // My own... (Added two new public members: string for detector and string for description of sample)
        EventSample(TTree* _tree, float ScaleFactor, std::string Det, std::string Desc) : tree(_tree), fScaleFactor(ScaleFactor), sDet(Det), sDesc(Desc) {

            // Detector
            if (Det == "sbnd") {
                sDet = "SBND";
            } else if (Det == "uboone") {
                sDet = "MicroBooNE";
            } else if (Det == "icarus") {
                sDet = "ICARUS";
            } else {
                std::cout << std::endl << "ERROR: " << Det << " not valid detector." << std::endl << std::endl;
                assert(false);
            }

            // Description
            if (Desc == "nu") {
                sDesc = "Neutrino";
            } else {
                std::cout << std::endl << "ERROR: Sample type" << Desc << " not supported." << std::endl << std::endl;
                assert(false);
            }

        }

        TTree* tree;            //!< Event tree
        float fScaleFactor;     //!< Factor for POT (etc.) scaling
        std::string sDet;       // What detector it comes from
        std::string sDesc;      // (Very concise) Description of sample
    
};


class Covariance {
    
    public:
        Covariance(std::vector<EventSample> samples);
        //SavePNGs(std::string directory);
        
        TH2D *covmat, *fcovmat, *corrmat;   // Covariance, fractional covariance and correlation matrices.
    
};

}   // namespace SBNOsc
}   // namespace ana

#endif// __sbnanalysis_ana_SBNOsc_Covariance__
