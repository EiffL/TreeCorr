cmake_minimum_required(VERSION 2.8)


# Adding customized cmake module
list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include(ExternalProject)
include(FindPkgConfig)

project(treecorr)

# Find required packages
find_package(PythonInterp REQUIRED)
find_package(PythonLibsNew REQUIRED)

# Compile boost packages to match the specific system environment
include(BuildBoost)
include(BuildBoostNumpy)

# Define include and library directories
include_directories(
  ${Boost_INCLUDE_DIRS}
  ${BoostNumpy_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS}
)
link_directories(${BoostNumpy_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS})

# Compilation flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang, disabling OpenMP support
  set(CMAKE_CXX_FLAGS "-O3 -fomit-frame-pointer -fno-common -fPIC")
else()
  set(CMAKE_CXX_FLAGS "-O3 -fomit-frame-pointer -fno-common -fPIC -fopenmp")
endif()

# Build the treecorr module
include_directories(include)
add_library(treecorr SHARED src/BinnedCorr2.cpp src/BinnedCorr3.cpp src/Cell.cpp src/Field.cpp)
add_dependencies(treecorr BoostNumpy Boost)
target_link_libraries(treecorr ${BoostNumpy_LIBRARIES} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
set_target_properties(treecorr PROPERTIES SUFFIX .so)
set_target_properties(treecorr PROPERTIES PREFIX "")

# Installs the module inside the python package 
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})


install(TARGETS treecorr DESTINATION treecorr)
