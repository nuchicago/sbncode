# Define the subproject name
project(NuMuSelection)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(DataTypes)
add_subdirectory(Selection)
