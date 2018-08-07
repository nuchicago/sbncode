# Define the subproject name
project(NuMuSelection)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(NuMuSelection_NuMuUtil SHARED NuMuUtil.cc)
target_link_libraries(
  NuMuSelection_NuMuUtil
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



install(TARGETS NuMuSelection_NuMuUtil DESTINATION lib)

add_subdirectory(DataTypes)
add_subdirectory(Selection)
