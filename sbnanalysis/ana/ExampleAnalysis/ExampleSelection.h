#ifndef __sbnanalysis_ana_ExampleAnalysis_ExampleSelection__
#define __sbnanalysis_ana_ExampleAnalysis_ExampleSelection__

/* Includes */
#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include <vector>
#include <string>

/* Forward declarations */
class TH2D;
class TH1D;
class THStack;
class TAttFill;

/** All analysis code is defined in namespace "ana" */
namespace ana {
/** Code specific to the ExampleAnalysis. */
namespace ExampleAnalysis {

class ExampleSelection : public core::SelectionBase {
    
    public:
    
    ExampleSelection();
    void Initialize(Json::Value* config=NULL);
    void Finalize();
    bool ProcessEvent(gallery::Event& ev);

    protected:
    
    unsigned fEventCounter;         //!< Count processed events

    /** Configuration parameters */
    art::InputTag fTruthTag;        //!< art tag for MCTruth information
    art::InputTag fPartTag;   // for use with particle (non-nu) information
    art::InputTag fEventWgtTag;
    art::InputTag fFluxWgtTag;

    /** Custom data branches */
    //int fNuCount;  //!< Number of neutrino interactions in the event
    
    /** My stuff */
    
    // Detector
    std::string det;
    
    // Particle properties
    double ml, mn, mp, Eb;
    
    // Weights and universes
    int n_unis;
    
    // Reason to pass cut
    int nucount, nuindet, mu_passes, pi_passes, double_pass;
    
    // Histograms
    std::vector <std::string> unis1, unis2;
    //double min_nuE, max_nuE, nuE_interval, num_nuE_bins;
    
    THStack *len_stack, *count_stack;
    std::vector <TH1D*> len_hists;
    std::vector <std::vector <TH1D*> > nuE_hists;
    
    /* Note: output from the ProcessEvent stage will be (let n index neutrinos from 0 to N and u index universes from 0 to U):
    1. nuE: vector of length N containing energies of detected neutrinos
    2. why_pass: vector of length N containing info on why that detected neutrino passed the cut (0 if OK, 1 if because of long-tracked pion, i.e. mistake)
    3. weights: vector of length N of vectors of length U, element weights[n][u] containing the weight corresponding to neutrino n for universe u.
    */
    
    // TestCuts
    int numucount, numu_CC, numu_CC_lessthanmu, nuV_indet;
    int nu_passcut1, nu_passcut2, nu_passdouble;
    int mu_pass1, mu_pass2, pi_pass1, pi_pass2;
    int mu_contained, mu_escaped, pi_contained, pi_escaped;
    TH1D *mu_len_contained, *mu_len_escaped, *pi_len_contained, *pi_len_escaped;
    THStack *contained, *escaped;
};

}// namespace ExampleAnalysis
}// namespace ana

#endif// __sbnanalysis_ana_ExampleAnalysis_ExampleSelection__

