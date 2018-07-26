#ifndef __sbnanalysis_ana_NuMuSelection_SelectNu__
#define __sbnanalysis_ana_NuMuSelection_SelectNu__

// Includes
#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"

namespace ana {

  /** Code specific to the NuMuSelection. */
  namespace NuMuSelection {

/**
 * \class SelectNu
 */
class SelectNu : public core::SelectionBase {
public:
  /** Constructor. */
  SelectNu();

  void Initialize(Json::Value* config=NULL);

  void Finalize();

  bool ProcessEvent(gallery::Event& ev);

protected:
  TTree *_neutrinoTree;
  double _min_reconstructed_R;
  std::vector<double> *_min_reconstructed_Rs; 

  // config
  double _config_dist_cut;

};

  }  // namespace NuMuSelection
}  // namespace ana

#endif  // __sbnanalysis_ana_NuMuSelection_SelectNu__

