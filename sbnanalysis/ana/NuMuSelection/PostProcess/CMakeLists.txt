# Define the subproject name
project(NuMuSelection)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# make dictionary for analysis classes
ROOT_GENERATE_DICTIONARY(NuMuSelection_PostProcess_dict OscChi2.h MakeHistos.h LINKDEF linkdef.h)
# and make it into a library
add_library(sbnanalysis_NuMuSelection_PostProcess SHARED NuMuSelection_PostProcess_dict.cxx OscChi2.cc MakeHistos.cc)

target_link_libraries(
  sbnanalysis_NuMuSelection_PostProcess
  NuMuSelection_NuMuUtil
  NuMuSelection_classes
  canvas
  cetlib
  gallery
  nusimdata_SimulationBase
  lardataobj_Simulation_dict
  larcoreobj_SummaryData
  lardataobj_RecoBase
  lardataobj_RecoBase_dict
  ${LARSIM_BASE_LIBRARY}
  ${LARSIM_BASE_DICT}
  ${ROOT_LIBRARIES}
)

# install
install(TARGETS sbnanalysis_NuMuSelection_PostProcess DESTINATION lib)

# and install the RecobInteraction dict
add_custom_command(TARGET NuMuSelection_PostProcess_dict POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        libNuMuSelection_PostProcess_dict_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)
