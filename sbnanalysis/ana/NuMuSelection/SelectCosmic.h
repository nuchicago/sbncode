#ifndef __sbnanalysis_ana_NuMuSelection_SelectCosmic__
#define __sbnanalysis_ana_NuMuSelection_SelectCosmic__

// Includes
#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"

namespace ana {

  /** Code specific to the NuMuSelection. */
  namespace NuMuSelection {

/**
 * \class SelectCosmic
 */
class SelectCosmic : public core::SelectionBase {
public:
  /** Constructor. */
  SelectCosmic();

  void Initialize(Json::Value* config=NULL);

  void Finalize();

  bool ProcessEvent(gallery::Event& ev);

protected:

};

  }  // namespace NuMuSelection
}  // namespace ana

#endif  // __sbnanalysis_ana_NuMuSelection_SelectCosmic__
