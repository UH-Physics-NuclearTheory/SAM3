cmake_minimum_required(VERSION 3.16)

project(SAM3
        LANGUAGES CXX
        DESCRIPTION "Subensemble Acceptance Method 3.0 (SAM-3.0) reference implementation")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------------------------------------------------------------------------
# Eigen (the only dependency). Resolved in this order:
#   1. an explicit include directory:  -DEIGEN3_INCLUDE_DIR=/path/to/eigen
#   2. an installed Eigen3 CMake package (find_package)
#   3. downloaded automatically via FetchContent (disable with -DSAM3_FETCH_EIGEN=OFF)
# ---------------------------------------------------------------------------
option(SAM3_FETCH_EIGEN "Download Eigen via FetchContent if it is not found locally" ON)

if(EIGEN3_INCLUDE_DIR AND EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/Dense")
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
  message(STATUS "Eigen: using EIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR}")
else()
  find_package(Eigen3 QUIET)
  if(TARGET Eigen3::Eigen)
    message(STATUS "Eigen: found installed package (${Eigen3_DIR})")
  elseif(SAM3_FETCH_EIGEN)
    message(STATUS "Eigen: not found locally -- downloading 3.4.0 via FetchContent")
    if(POLICY CMP0135)
      cmake_policy(SET CMP0135 NEW)   # use extraction time for downloaded files
    endif()
    include(FetchContent)
    # Eigen is header-only: SOURCE_SUBDIR points at a directory with no
    # CMakeLists.txt, so MakeAvailable downloads/extracts but does NOT run
    # Eigen's own CMake project.
    FetchContent_Declare(eigen
      URL           https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
      URL_HASH      SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
      SOURCE_SUBDIR headers-only-no-build)
    FetchContent_MakeAvailable(eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${eigen_SOURCE_DIR}")
  else()
    message(FATAL_ERROR
      "Eigen not found. Install Eigen3, pass -DEIGEN3_INCLUDE_DIR=/path/to/eigen, "
      "or enable the download with -DSAM3_FETCH_EIGEN=ON.")
  endif()
endif()

# ---------------------------------------------------------------------------
# Header-only SAM-3.0 library (sam3.h).
# ---------------------------------------------------------------------------
add_library(sam3 INTERFACE)
target_include_directories(sam3 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(sam3 INTERFACE Eigen3::Eigen)
target_compile_features(sam3 INTERFACE cxx_std_17)

# ---------------------------------------------------------------------------
# Tests and cross-check tool.
# ---------------------------------------------------------------------------
# Runnable examples: minimal quick-start, and a worked acceptance model.
add_executable(example example.cpp)
target_link_libraries(example PRIVATE sam3)

add_executable(example_binomial example_binomial.cpp)
target_link_libraries(example_binomial PRIVATE sam3)

add_executable(test_sam3 test_sam3.cpp)
target_link_libraries(test_sam3 PRIVATE sam3)

# Emits a Mathematica data file for the C++ <-> notebook cross-check
# (see ../crosscheck.sh and ../Mathematica/crosscheck_cpp.wl).
add_executable(bridge_sam3 bridge_sam3.cpp)
target_link_libraries(bridge_sam3 PRIVATE sam3)

enable_testing()
add_test(NAME SAM3-tests COMMAND test_sam3)
