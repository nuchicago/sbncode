////////////////////////////////////////////////////////////////////////
// Class:       PandoraTesting
// Plugin Type: analyzer (art v2_08_04)
// File:        PandoraTesting_module.cc
//
// Generated at Mon Dec 18 08:00:14 2017 by Rhiannon Jones using cetskelgen
// from cetlib version v3_01_01.
////////////////////////////////////////////////////////////////////////

// Playing around with Pandoraw objects
// Taken from Rhiannon Jones github


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "TROOT.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"

namespace pndr {
  class PandoraTesting;
}


class pndr::PandoraTesting : public art::EDAnalyzer {
public:
  explicit PandoraTesting(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PandoraTesting(PandoraTesting const &) = delete;
  PandoraTesting(PandoraTesting &&) = delete;
  PandoraTesting & operator = (PandoraTesting const &) = delete;
  PandoraTesting & operator = (PandoraTesting &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p);

private:

  // Declare member data here.
  std::map< std::vector< int >, int > m_selection;
  float m_detectorHalfLengthX;
  float m_detectorHalfLengthY;
  float m_detectorHalfLengthZ;
  float m_coordinateOffsetX;
  float m_coordinateOffsetY;
  float m_coordinateOffsetZ;
  float m_selectedBorderX;
  float m_selectedBorderY;
  float m_selectedBorderZ;

  // counters 
  int events, fiducial, not_fiducial, topology;

  // nTuples
  TNtuple *fNt_vtx;
};


pndr::PandoraTesting::PandoraTesting(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  
  this->reconfigure(p);

}

void pndr::PandoraTesting::analyze(art::Event const & e)
{
 
  // Implementation of required member function here.
  // Get the MCTruth information 
  art::Handle< std::vector< simb::MCTruth > > mct_handle;
  e.getByLabel("generator", mct_handle );
  int mct_size = mct_handle->size();
  // Initialise to be true, and if any of the counters don't match, 
  // change to false
  bool contained_topology = true;

  typedef::std::map< std::vector< int >, int > topology_map;
  
  if(mct_handle.isValid() && mct_size) {
  
    // Loop over the truth info
    for(auto const& mct : (*mct_handle)) {
 
      // Check the neutrino came from the beam
      if(mct.Origin() != simb::kBeamNeutrino) continue;
 
      events++;

      //-------------------------------------------------------------------------
      //   Check the neutrino interaction vertex is within the fiducial volume
      //-------------------------------------------------------------------------
      float nuVtxX = mct.GetNeutrino().Lepton().Vx();
      float nuVtxY = mct.GetNeutrino().Lepton().Vy();
      float nuVtxZ = mct.GetNeutrino().Lepton().Vz();
 
 
      if (    (nuVtxX > (m_detectorHalfLengthX - m_coordinateOffsetX - m_selectedBorderX)) 
           || (nuVtxX < (-m_coordinateOffsetX + m_selectedBorderX)) 
           || (nuVtxY > (m_detectorHalfLengthY - m_coordinateOffsetY - m_selectedBorderY)) 
           || (nuVtxY < (-m_coordinateOffsetY + m_selectedBorderY)) 
           || (nuVtxZ > (m_detectorHalfLengthZ - m_coordinateOffsetZ - m_selectedBorderZ)) 
           || (nuVtxZ < (-m_coordinateOffsetZ + m_selectedBorderZ))){
 
        contained_topology = false;
   
        not_fiducial++;

      }
      else{
      
        fiducial++;

        
        //-------------------------------------------------------------------------
        //                           Check the interaction
        //-------------------------------------------------------------------------
   
        // Get the number of particles
        const int n_particles = mct.NParticles();
   
        // Loop over the map and get the elements
        for( topology_map::iterator it = m_selection.begin(); it != m_selection.end(); ++it ){

          std::vector< int > pdgVect = it->first;
          int                count   = it->second;

          // Initialise a counter for the number within the mct vector
          int particle_counter = 0;

          // Loop over the particles in the event
          for( int i = 0; i < n_particles; ++i ){

            // Find if the pdg code of the current particle is one of the ones in the map
            if( std::find( pdgVect.begin(), pdgVect.end(), mct.GetParticle(i).PdgCode() ) != pdgVect.end() ) ++particle_counter;

          }

          // If the counters don't match, return false
          if( particle_counter != count ){

            // Not topology event
            contained_topology =  false;

          }
        }
      }
    }
  }

  if( contained_topology ){

    topology++;
    
    // Get the PFParticle information 
    art::Handle< std::vector< recob::PFParticle > > pfp_handle;
    e.getByLabel("pandoraNu", pfp_handle );
    int pfp_size = pfp_handle->size();
    
    if( pfp_handle.isValid() && pfp_size ){
      
      // Loop over truth
      for( auto & mct : (*mct_handle) ){

        // True neutrino vertex ( primary vertex position)
        double nu_x, nu_y, nu_z;

        nu_x = mct.GetNeutrino().Lepton().Vx();
        nu_y = mct.GetNeutrino().Lepton().Vy();
        nu_z = mct.GetNeutrino().Lepton().Vz();
        
        // If the current particle is the primary, get the track/vertex
        // association and make a list of the distances from the true neutrino 
        // vertex
        art::FindMany< recob::Track  > ftrk( pfp_handle, e, "pandoraNu" );
        art::FindMany< recob::Vertex > fvtx( pfp_handle, e, "pandoraNu" );

        for( int i = 0; i < pfp_size; ++i ){
        
          art::Ptr< recob::PFParticle > pfp( pfp_handle, i );
        
          if( pfp->IsPrimary() ){
            
            // Define track and vertex handles
            std::vector<const recob::Track*>  trk_assn = ftrk.at(i);
            std::vector<const recob::Vertex*> vtx_assn = fvtx.at(i);
          
            // Get the distance from the associated reconstructed vertex and the
            // true neutrino vertex
            for( unsigned int j = 0; j < vtx_assn.size(); ++j ){
            
              double xyz[3];

              // Set array to be current vertex position
              vtx_assn[j]->XYZ(xyz);

              // x,y,z of current vertex
              double x, y, z, R;

              x = xyz[0];
              y = xyz[1];
              z = xyz[2];
              R =  sqrt( pow( ( x - nu_x ), 2 ) + pow( ( y - nu_y ), 2 ) + pow( ( z - nu_z ), 2 ) ); 
        
              fNt_vtx->Fill( x, y, z, nu_x, nu_y, nu_z, R );
            }
          }
        }
      }
    }
  }
}

void pndr::PandoraTesting::beginJob()
{
  
  // Implementation of optional member function here.
  // Counter initialisation
  events       = 0;
  fiducial     = 0;
  not_fiducial = 0;
  topology     = 0;

  // Initialise nTuples

  fNt_vtx = new TNtuple( "fNt_vtx", "Pandora vertex metrics", "x:y:z:nuX:nuY:nuZ:R" );
  fNt_vtx->SetDirectory(0);

}

void pndr::PandoraTesting::endJob()
{
  // Implementation of optional member function here.
  // Loop over everything and check it's working  
  // Show the chosen topology to filter on
  typedef::std::map< std::vector< int >, int > topology_map;

  std::cout << "=================================================================================" << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;

  std::cout << "Topology:" << std::endl;
  
  for ( topology_map::iterator it = m_selection.begin(); it != m_selection.end(); ++it ) {
    std::vector< int > pdgVect = it->first;
    int                count   = it->second;
    
    std::cout << "  " << count << " particles with PDG in [";
    for ( int pdg : pdgVect ) {
        std::cout << " " << pdg << " ";
    }
    std::cout << "]" << std::endl;

  }

  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  
  std::cout << " Number of events, total           : " << events       << std::endl;
  std::cout << " Number of events, fiducial        : " << fiducial     << std::endl;
  std::cout << " Number of events, not fiducial    : " << not_fiducial << std::endl;
  std::cout << " Number of events, fiducial, cc0pi : " << topology     << std::endl;

  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  
  std::cout << " Fraction of fiducial events       : " << fiducial / double( events )   << std::endl;
  std::cout << " Fraction of fiducial cc0pi events : " << topology / double( fiducial ) << std::endl;
  
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  std::cout << "=================================================================================" << std::endl;

  TFile *f = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/pandoratesting/pandoratesting/root/pandora_vertex_metrics.root", "RECREATE" );

  fNt_vtx->Write();

  f->Close();

  delete f;
  delete fNt_vtx;

}

void pndr::PandoraTesting::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.

   std::vector< int > blankVect;
   std::vector< std::vector< int > > input;
 
   std::vector< int > selection1 = p.get< std::vector< int > >("Selection1",        blankVect);
   if ( selection1.size() != 0 ) input.push_back(selection1);
 
   std::vector< int > selection2 = p.get< std::vector< int > >("Selection2",        blankVect);
   if ( selection2.size() != 0 ) input.push_back(selection2);
 
   std::vector< int > selection3 = p.get< std::vector< int > >("Selection3",        blankVect);
   if ( selection3.size() != 0 ) input.push_back(selection3);
 
 
   for ( auto & inputVect : input ) {
     if ( inputVect.size() < 2 ) {
       std::cerr << " Error: Selection vector must have at least 2 elements " <<    std::endl;
       std::cerr << "        First element:     Number of particles of PDG code(s)  specified " << std::endl;
       std::cerr << "        Remaining element: PDG codes to filter on " << std::   endl;
       exit(1);
     }
 
     int count = inputVect[0];
     inputVect.erase( inputVect.begin() );
 
     m_selection.insert( std::make_pair( inputVect, count ) );
   }
 
   m_detectorHalfLengthX = p.get<float>("DetectorHalfLengthX");
   m_detectorHalfLengthY = p.get<float>("DetectorHalfLengthY");
   m_detectorHalfLengthZ = p.get<float>("DetectorHalfLengthZ");
   m_coordinateOffsetX   = p.get<float>("CoordinateOffsetX");
   m_coordinateOffsetY   = p.get<float>("CoordinateOffsetY");
   m_coordinateOffsetZ   = p.get<float>("CoordinateOffsetZ");
   m_selectedBorderX     = p.get<float>("SelectedBorderX");
   m_selectedBorderY     = p.get<float>("SelectedBorderY");
   m_selectedBorderZ     = p.get<float>("SelectedBorderZ");
 
}

DEFINE_ART_MODULE(pndr::PandoraTesting)
