# Define the subproject name
project(SBNOsc)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories($ENV{LARSIM_INC})

find_library(LARSIM_BASE_LIBRARY larsim_EventWeight_Base HINTS $ENV{LARSIM_LIB})
find_library(LARSIM_BASE_DICT larsim_EventWeight_Base_dict HINTS $ENV{LARSIM_LIB})

# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_SBNOsc_Utilities SHARED Utilities.cxx)
add_library(sbnanalysis_SBNOsc_NueSelection SHARED NueSelection.cxx)
add_library(sbnanalysis_SBNOsc_NumuSelection SHARED NumuSelection.cxx)
add_library(sbnanalysis_SBNOsc_Covariance SHARED Covariance.cxx)
add_library(sbnanalysis_SBNOsc_Chi2Sensitivity SHARED Chi2Sensitivity.cxx)

# Define library linking
target_link_libraries(sbnanalysis_SBNOsc_NueSelection sbnanalysis_SBNOsc_Utilities)
target_link_libraries(sbnanalysis_SBNOsc_NumuSelection sbnanalysis_SBNOsc_Utilities)
target_link_libraries(sbnanalysis_SBNOsc_Chi2Sensitivity sbnanalysis_SBNOsc_Covariance)
target_link_libraries(${LARSIM_BASE_LIBRARY})
target_link_libraries(${LARSIM_BASE_DICT})

# Build executables
add_executable(sbn-cov sbn-cov.cxx)
target_link_libraries(
  jsoncpp
  sbn-cov
  sbnanalysis_SBNOsc_Covariance
  sbnanalysis_Event
  ${ROOT_LIBRARIES}
)

add_executable(sbn-sensitivity sbn-sensitivity.cxx)
target_link_libraries(
  jsoncpp
  sbn-sensitivity
  sbnanalysis_SBNOsc_Chi2Sensitivity
  sbnanalysis_SBNOsc_Covariance
  sbnanalysis_Event
  ${ROOT_LIBRARIES}
)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SBNOsc_Utilities DESTINATION lib)
install(TARGETS sbnanalysis_SBNOsc_NueSelection DESTINATION lib)
install(TARGETS sbnanalysis_SBNOsc_NumuSelection DESTINATION lib)
install(TARGETS sbnanalysis_SBNOsc_Covariance DESTINATION lib)
install(TARGETS sbnanalysis_SBNOsc_Chi2Sensitivity DESTINATION lib)

install(TARGETS sbn-cov DESTINATION bin/sbnosc)
install(TARGETS sbn-sensitivity DESTINATION bin/sbnosc)

