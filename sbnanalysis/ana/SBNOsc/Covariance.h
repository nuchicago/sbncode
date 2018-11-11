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
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

class TTree;

namespace ana {
namespace SBNOsc {

class EventSample {
    
    public:
    
        /** Constructors. */
        EventSample();
        EventSample(TTree* _tree, float scaleFactor) : tree(_tree), fScaleFactor(scaleFactor) {}
        EventSample(std::vector<std::string> filenames, float fScaleFactor);
        EventSample(TFile* _file, float ScaleFactor, std::string Det, std::string Desc, std::vector <double> bins, int scale_sample, std::string nutype);
        
        TFile* file;                //!< File containing the tree
        TTree* tree;                //!< Event tree
        float fScaleFactor;         //!< Factor for POT (etc.) scaling
        std::string fDet;           //!< What detector it comes from
        std::string fDesc;          //!< (Very concise) Description of sample
        std::vector <double> fBins; //!< Energy bin limits
        int fScaleSample;           //!< Scale to this sample (shape+rate chisq)?
        bool fNuType;               //!< Is this sample a neutrino sample
    
};


class Covariance {
    
    public:
        
        // Functions
        
        Covariance(std::vector <EventSample> samples, char *configFileName);
        
        void ScanEvents(), GetCovs(), GetCounts(), Write(std::string directory);
        
        // Output
        
        TH2D *cov, *fcov, *corr;                // Covariance, fractional covariance and correlation
                                                // matrices.
        
        std::vector <TH1D*> numu_counts, 
            numu_bkgs, nue_counts, nue_bkgs;    // For plotting pretty histograms
    
    private:
    
        // Configuration parameters
        
        std::string fWeightKey;
        int fNumAltUnis;
        
        std::string fEnergyType;
        
        double fSelectionEfficiency, fRejectionEfficiency;
        
        std::map <std::string, float> fScaleTargets;
        
        bool fSignalOnly;
        
        // Internal functions
        
        std::vector <double> GetUniWeights(std::map <std::string, std::vector <double> > weights, int n_unis);
        
        // Stored objects
        
        std::vector <int> sample_bins;
        int num_bins;
        
        std::vector <EventSample> ev_samples;
        
        std::vector <std::vector <double> > nu_counts;
        std::vector <double> bkg_counts;
        
        
    
};

}   // namespace SBNOsc
}   // namespace ana

#endif// __sbnanalysis_ana_SBNOsc_Covariance__
