# Define the subproject name
project(NuMuSelection)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_NuMuSelection_SelectCosmic SHARED SelectCosmic.cxx)
add_library(sbnanalysis_NuMuSelection_SelectNu SHARED SelectNu.cxx)

# Define library linking
# target_link_libraries(sbnanalysis_ExampleAnalysis_ExampleSelection sbnanalysis_ExampleAnalysis_ExampleTools)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_NuMuSelection_SelectCosmic DESTINATION lib)
install(TARGETS sbnanalysis_NuMuSelection_SelectNu DESTINATION lib)

