# Define the subproject name
project(NuMuSelection)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# make dictionary for analysis classes
ROOT_GENERATE_DICTIONARY(NuMuSelection_dict RecobInteraction.h LINKDEF linkdef.h)
# and make it into a library
add_library(NuMuSelection_classes NuMuSelection_dict.cxx)

# install
install(TARGETS NuMuSelection_classes DESTINATION lib)

# and install the RecobInteraction dict
add_custom_command(TARGET NuMuSelection_classes POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
	libNuMuSelection_dict_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)	


