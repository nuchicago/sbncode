#include <iostream>
#include <vector>
#include <string>
#include <TMath.h>
#include <TLorentzVector.h>
#include <json/json.h>

#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "core/Event.hh"
#include "NumuSelection.h"
#include "Utilities.h"

namespace ana {
namespace SBNOsc {

bool in_SBND(std::vector <double> coords) {

    return (coords[0] > -189.15)*(coords[0] < -12.65)*(coords[1] > -190)*(coords[1] < 190)*(coords[2] > 10)*(coords[2] < 490) + (coords[0] > 12.65)*(coords[0] < 189.15)*(coords[1] > -190)*(coords[1] < 190)*(coords[2] > 10)*(coords[2] < 490);

}

bool in_MicroBooNE(std::vector <double> coords) {
    
    return (coords[0] > 8.45)*(coords[0] < 244.8)*(coords[1] > -105.53)*(coords[1] < 107.47)*(coords[2] > 10.1)*(coords[2] < 1026.9);
    
}

bool in_ICARUS(std::vector <double> coords) {
    
    return (coords[0] > -354.49)*(coords[0] < -226.29)*(coords[1] > -163.41)*(coords[1] < 133.41)*(coords[2] > -899.950652)*(coords[2] < 869.950652) + (coords[0] > -206.14)*(coords[0] < -77.94)*(coords[1] > -163.41)*(coords[1] < 133.41)*(coords[2] > -899.950652)*(coords[2] < 869.950652);
    
    // Still need to add other TPCs if the simulation starts including them...
    
}

bool in_detector(std::vector <double> coords, std::string det) {
    
    if (det == (std::string)"SBND") {
        return in_SBND(coords);
    } else if (det == (std::string)"MicroBooNE") {
        return in_MicroBooNE(coords);
    } else if (det == (std::string)"ICARUS") {
        return in_ICARUS(coords);
    } else {
        std::cout << "Not a valid detector!" << std::endl;
        return 0;
    }
    
}

std::vector <double> get_lens(auto mcpart, std::string det) {
    
    /* Loop over positions until one is outside of the detector. Add up the distance travelled as we go. 
       Check if the particle leaves detector. If so then retun in second position of return vector. Else, 
       return in first position. */
    
    int retind = 0;                     // Index of return vector to return length in
    double len = 0;
    std::vector <double> oldcoords;
    
    for (int i = 0; i < mcpart->NumberTrajectoryPoints(); i++) {
        
        TLorentzVector event = mcpart->Position(i);
        std::vector <double> newcoords = {event.X(), event.Y(), event.Z()};
        if (i == 0) { oldcoords = newcoords; }
        
        if (!in_detector(newcoords, det)) {
            retind = 1;
            break;
        }
        
        len += TMath::Sqrt(TMath::Power(newcoords[0] - oldcoords[0], 2) +
                           TMath::Power(newcoords[1] - oldcoords[1], 2) +
                           TMath::Power(newcoords[2] - oldcoords[2], 2));
        
        oldcoords = newcoords;
        
    }
    
    std::vector <double> ret_vec(2);
    ret_vec[retind] = len;
    
    return ret_vec;
    
}

NumuSelection::NumuSelection() : SelectionBase(), fEventCounter(0), fNuAll(0), fNuCount(0), fNuinFid(0) {}

void NumuSelection::Initialize(Json::Value* config) {
    
    // Load configuration parameters
    fTruthTag = { "generator" };
    fPartTag = { "largeant" };
    fDet = "";

    if (config) {
        fTruthTag = { (*config)["SBNOsc"].get("MCTruthTag", "generator").asString() };
        fPartTag = { (*config)["SBNOsc"].get("MCParticleTag", "largeant").asString() };
        fDet = (*config)["SBNOsc"].get("detector", "").asString();
    }
    
    std::cout << std::endl << std::endl << std::endl << "fDet is: " << fDet << "." << std::endl << std::endl << std::endl << std::endl;
    
    // For reconstructing energies
    mn = 0.93956541; mp = 0.93827208; ml = 0.10565837; Eb = 0.030;

}

bool NumuSelection::ProcessEvent(const gallery::Event& ev, std::vector<Event::Interaction>& reco) {
    
    if (fEventCounter % 100 == 0) {
        std::cout << std::endl << "NumuSelection: Processing event " << fEventCounter << " ("
        << fNuCount << " neutrinos selected of " << fNuinFid << " that reacted in fiducial vol.)."
        << std::endl << std::endl;
    }
    fEventCounter++;

    // Grab a data product from the event
    auto const& mctruths = *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);
    auto const& mctruth_handle = ev.getValidHandle<std::vector<simb::MCTruth>>(fTruthTag);
    
    /* Associations b/w MCTruths (neutrinos) and MCParticles */
    const art::FindManyP <simb::MCParticle, sim::GeneratedParticleInfo> find_many_mcparticle(mctruth_handle, ev, fPartTag);

    // Iterate through the neutrinos
    for (size_t i=0; i<mctruths.size(); i++) {
        
        auto const& mctruth = mctruths.at(i);
        const simb::MCNeutrino& nu = mctruth.GetNeutrino();
        
        /* Skip if vertex is outside fiducial volume */
        fNuAll++; fNuinFid++;
        std::vector <double> vertex = {nu.Nu().Vx(), nu.Nu().Vy(), nu.Nu().Vz()};
        if (!in_detector(vertex, fDet)) { fNuinFid--; continue; }
        
        /* Loop over associated MCParticles to see if this neutrino passes the cuts */
        std::vector <art::Ptr <simb::MCParticle> > const& mcparticle = find_many_mcparticle.at(i);
        
        std::vector <int> pass_cut(2); double nuE;
        
        for (auto const& mcpart : mcparticle) {
            
            // Look only at muons and pions
            if ((mcpart->PdgCode()*mcpart->PdgCode() == 13*13 && nu.CCNC() == simb::kCC) || 
                (mcpart->PdgCode()*mcpart->PdgCode() == 211*211 && nu.CCNC() == simb::kNC)) {
                
                // Cut off those whose total track length is less than 50cm
                if (mcpart->Trajectory().TotalLength() < 50) { continue; }
                
                // Muon (0) or pion (1)? Could create extra branch for this, potentially...
                int part_ind = 0;
                if (mcpart->PdgCode()*mcpart->PdgCode() == 211*211) { part_ind = 1; }
                
                // Store lens
                std::vector <double> lens_forcuts = get_lens(mcpart, fDet);
                
                // Pass cuts?
                if (lens_forcuts[0] > 50 || lens_forcuts[1] > 100) { 
                
                    pass_cut[part_ind] += 1;
                
                    double px = mcpart->Px(), py = mcpart->Py(), pz = mcpart->Pz(), E = mcpart->E(),
                        p = TMath::Sqrt(px*px + py*py + pz*pz), pxy = TMath::Sqrt(px*px + py*py),
                        theta = TMath::ATan(pxy/pz);
                    
                    nuE = 0.5 * (mp*mp - (mn - Eb)*(mn - Eb) - ml*ml + 2*(mn - Eb)*E) / ((mn - Eb) - E + p*TMath::Cos(theta));
                    
                }
                
            }
                
        } // End loop over MCParticles
        
        if (pass_cut[0] + pass_cut[1] == 1) {
            Event::Interaction interaction = TruthReco(mctruth);
            interaction.neutrino.energy = nuE;
            reco.push_back(interaction);
        }
        
    }

    bool selected = !reco.empty();

    if (selected) { fNuCount += reco.size(); }

    return selected;
}

void NumuSelection::Finalize() {
    
    std::cout << std::endl << std::endl << std::endl 
    << "In all, we had " << fNuAll << " neutrinos" << std::endl
    << "Of these, " << fNuinFid << " reacted in the fiducial volume (a proportion" << (double)fNuinFid/fNuAll << ")" << std::endl
    << "And we selected " << fNuCount << " of those (a proportion " << (double)fNuCount/fNuinFid << ")." 
    << std::endl << std::endl << std::endl;
    
    std::cout << fNuAll << " & " << fNuinFid << " (" << (double)fNuinFid/fNuAll << ") & " << fNuCount << " (" << (double)fNuCount/fNuinFid << ") " << std::endl << std::endl;
    
}


}// namespace SBNOsc
}// namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NumuSelection)