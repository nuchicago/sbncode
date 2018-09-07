#ifndef __sbnanalysis_ana_SBNOsc_Covariance__
#define __sbnanalysis_ana_SBNOsc_Covariance__

/**
 * \file Covariance.h
 */

#include "json/json.h"
#include "core/Loader.hh"

#include <string>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cassert>

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>

class TTree;

namespace ana {
namespace SBNOsc {

class EventSample {
    
    public:
        /** Constructors. */
        EventSample();
        EventSample(TTree* _tree, float scaleFactor) : tree(_tree), fScaleFactor(scaleFactor) {}
        EventSample(std::vector<std::string> filenames, float fScaleFactor);
        EventSample(TFile* _file, float ScaleFactor, std::string Det, std::string Desc);
        
        TFile* file;            //!< File containing the tree
        TTree* tree;            //!< Event tree
        float fScaleFactor;     //!< Factor for POT (etc.) scaling
        std::string fDet;       //!< What detector it comes from
        std::string fDesc;      //!< (Very concise) Description of sample
    
};


class Covariance {
    
    public:
        Covariance(std::vector<EventSample> samples, char *configFileName);
        //SavePNGs(std::string directory);
        
        TH2D *covmat, *fcovmat, *corrmat;       // Covariance, fractional covariance and correlation matrices.
        
        TH1D *bkg_counts;                       // Counts of non-neutrino events. Won't be oscillated
        TH2D *nu_counts;                        // Counts of geniune neutrino events. For oscillation
        
        TH1D *CV_counts;                        // CV universe counts.
        std::vector <double> energies;          // Bin centres for bins in all hists
        
        std::vector <std::string> sample_order; // Description of samples, in order they were plotted.
        std::vector <int> sample_bins;          // Bin limits of samples.
    
    private:
    
        std::string fWeightKey;
        int fNumAltUnis;
        
        std::string fEnergyType;
        
        double fSelectionEfficiency, fRejectionEfficiency;
        
        std::map <std::string, std::vector <double> > fBins;
        
        std::map <std::string, float> fScaleTargets;
        
        std::string fOutputDirectory;
        
        
    
};

}   // namespace SBNOsc
}   // namespace ana

#endif// __sbnanalysis_ana_SBNOsc_Covariance__
