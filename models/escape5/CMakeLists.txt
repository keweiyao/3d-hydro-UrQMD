cmake_minimum_required(VERSION 2.8)
project(escape CXX Fortran)


# Search for the HDF5 C++ library.  The CMake FindHDF5 package does not work
# correctly for C++, because it will report success if it finds only the C
# library and not C++.  It also uses the rather unfortunate 'h5c++' compiler
# wrapper which sets a lot of unnecessary craps

set(CMAKE_Fortran_COMPILER "gfortran") 
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(ROOT)
if(ROOT_FOUND)
  message(STATUS "Found ROOT library: ${ROOT_LIBRARIES}")
  include_directories(${ROOT_INCLUDE_DIR})
  set(ROOT_LIBS ${ROOT_LIBRARIES})
else()
  message(FATAL_ERROR "root not found!")
endif()


# Set default architecture flags to optimize for the current system,
# but don't override user settings.
if(NOT CMAKE_CXX_FLAGS MATCHES "-march|-mtune")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# warnings adapted from http://stackoverflow.com/a/9862800
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Werror")

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_REQUIRED_FLAGS} --coverage")

add_subdirectory(src.visc)
install (DIRECTORY Tb DESTINATION ${PROJECT_NAME})
