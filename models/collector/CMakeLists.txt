cmake_minimum_required(VERSION 2.8)
project(collector CXX)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" CXX11)
if(NOT CXX11)
  message(FATAL_ERROR "Your compiler does not support C++11")
endif()


# Find and use GSL
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/src/CMake)
FIND_PACKAGE(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

# Find and use Boost.
find_package(Boost 1.50 REQUIRED COMPONENTS filesystem program_options system)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# Search for the HDF5 C++ library.  The CMake FindHDF5 package does not work
# correctly for C++, because it will report success if it finds only the C
# library and not C++.  It also uses the rather unfortunate 'h5c++' compiler
# wrapper which sets a lot of unnecessary crap.
find_library(HDF5_CXX_LIB hdf5_cpp)
if("${HDF5_CXX_LIB}" STREQUAL "HDF5_CXX_LIB-NOTFOUND")
  message(STATUS "HDF5 C++ library not found")
else()
  message(STATUS "Found HDF5 C++ library: ${HDF5_CXX_LIB}")
  # Having found the C++ lib, also locate the header and C lib.
  find_path(HDF5_INCLUDE_DIRS H5Cpp.h)
  include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
  find_library(HDF5_C_LIB hdf5)
  # Also link to zlib for HDF5 compression.
  set(HDF5_LIBRARIES ${HDF5_CXX_LIB} ${HDF5_C_LIB} z)
  # Set preprocessor definition so HDF5 output is compiled.
  add_definitions(-DTRENTO_HDF5)
endif()

# Set default architecture flags to optimize for the current system,
# but don't override user settings.
if(NOT CMAKE_CXX_FLAGS MATCHES "-march|-mtune")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# warnings adapted from http://stackoverflow.com/a/9862800
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Wpedantic -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-conversion -Wsign-promo -Wstrict-overflow=2 -Wundef")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Werror")

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og")
endif()

set(LIBRARY_NAME "lib${PROJECT_NAME}")

add_subdirectory(src)
install(FILES calc_obs.py DESTINATION ${PROJECT_NAME})

