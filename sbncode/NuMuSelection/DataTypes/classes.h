#include "sbncode/NuMuSelection/DataTypes/AssocTruthInfo.h"
#include "canvas/Persistency/Common/Wrapper.h"
#include <vector>

namespace {
  struct dictionary {
    numuselection::AssocTruthInfo ato;
    std::vector<numuselection::AssocTruthInfo> v_ato;
    art::Wrapper<numuselection::AssocTruthInfo> w_ato;
    art::Wrapper<std::vector<numuselection::AssocTruthInfo>> w_v_ato;
  };
}
