# Define the subproject name
project(TruthSelection)

# Directories where includes are found
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Declare libraries to build
add_library(sbnanalysis_TruthSelection_TruthSelection SHARED TruthSelection.cxx)
add_library(sbnanalysis_TruthSelection_Selections SHARED Selections.cxx)

# Define library linking
target_link_libraries(sbnanalysis_TruthSelection_TruthSelection sbnanalysis_TruthSelection_Selections)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_TruthSelection_Selections DESTINATION lib)
install(TARGETS sbnanalysis_TruthSelection_TruthSelection DESTINATION lib)

