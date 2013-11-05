# mpm CMakeLists.txt
# 
# Breannan Smith (smith@cs.columbia.edu) 
# Updated: 11/05/2013

cmake_minimum_required( VERSION 2.6 )

#project( MPM C CXX Fortran )
project( MPM C )

# If the user provided a CMAKE_BUILD_TYPE, create a cache variable
if( NOT CMAKE_BUILD_TYPE STREQUAL "" )
  set( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
       build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
       Release RelWithDebInfo MinSizeRel." FORCE )
# Otherwise default to debug
else( NOT CMAKE_BUILD_TYPE STREQUAL "" )
  set( CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of
       build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
       Release RelWithDebInfo MinSizeRel." FORCE )
endif( NOT CMAKE_BUILD_TYPE STREQUAL "" )

# Enforce the use of a valid build type
if( NOT CMAKE_BUILD_TYPE STREQUAL "None" AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
  if( NOT CMAKE_BUILD_TYPE STREQUAL "Release" AND NOT CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
    if( NOT CMAKE_BUILD_TYPE STREQUAL "MinSizeRel" )
      message (SEND_ERROR "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
    endif( NOT CMAKE_BUILD_TYPE STREQUAL "MinSizeRel" )
  endif( NOT CMAKE_BUILD_TYPE STREQUAL "Release" AND NOT CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
endif( NOT CMAKE_BUILD_TYPE STREQUAL "None" AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )

# Add a compiler define for diagnostic purposes
#add_definitions( -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} )

###############################################################################
### Set flags used across all builds (TODO: Check compiler types here!)

# Append an extra space if flags are already set
if( CMAKE_CXX_FLAGS )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} " )
endif( CMAKE_CXX_FLAGS )

# Enable aggressive warnings
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused")
# Enable c++0x support
#set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x" )
#set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c -std=gnu99 -O3 -Wall -Wstrict-prototypes -pedantic -g -pg -funroll-loops -I/usr/lib/openmpi/include/ -rdynamic" )

###############################################################################
### Set preprocessor defines used in all builds

add_definitions( -DCOLUMBIA )

# Disable Eigen's built in multithreading as we will use OpenMP around Eigen
#add_definitions( -DEIGEN_DONT_PARALLELIZE )
# Disable spurious warnings in Intel compilers
#add_definitions( -DEIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS )

# Optionally enable timers
#option( ENABLE_TIMERS "Enable the collection of timing stats for the simulation" OFF )
#if( ENABLE_TIMERS )
#  add_definitions( -DTIMERS_ENABLED )
#endif()

# Optionally enable resumable serialization of state
#option( ENABLE_SERIALIZATION "Enable resumable serialization of simulation state" OFF )

###############################################################################
### Set preprocessor defines used in debug builds

#if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
#  add_definitions( -DEIGEN_DONT_ALIGN ) # Also does -DEIGEN_DONT_VECTORIZE
#  add_definitions( -DEIGEN_INITIALIZE_MATRICES_BY_NAN )
#endif( CMAKE_BUILD_TYPE STREQUAL "Debug" )


###############################################################################
### Set paths, etc

# Create a default location for custom cmake modules
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake )

# Allow projects to reference source relative to the root of the project
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )


###############################################################################
### Allow the user to set a preferred BLAS and LAPACK implementation

#set( SUGGESTED_BLAS_LAPACK "" CACHE STRING "The desired LAPACK implementation (e.g. Generic, Apple, ACML, Intel10_64lp, Intel10_64lp_seq, ..." )


###############################################################################
### Allow OpenMP use throughout

#find_package( OpenMP )
#if( OPENMP_FOUND )
#  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
#endif()


###############################################################################
### Add libraries, plugins

# Allow targets to generate tests
#enable_testing()


# Optionally build in support for Qt4
#option( USE_QT4 "Build Qt4 problem drivers" OFF )

# Optionally build in support for Ipopt
#option( USE_IPOPT "Use the Ipopt solver" OFF )


# Core SDIC library
#add_subdirectory( SDIC )

# Core two dimensional ball simulation library
#add_subdirectory( TwoDBalls )

# Utilities shared by the command-line and Qt4 interface
#add_subdirectory( TwoDBallsUtils )

# Command line interface for two dimensional ball simulation
#add_subdirectory( TwoDBallsCLI )

# Qt4 interface for two dimensional ball simulation
#if( USE_QT4 )
#  add_subdirectory( TwoDBallsQt4 )
#endif()

# Unit tests for two dimensional ball simulation
#add_subdirectory( TwoDBallsTests )


# Core three dimensional rigid body simulation library
#add_subdirectory( ThreeDRigidBodies )

# Command line interface for three dimensional rigid body simulation
#add_subdirectory( ThreeDRigidBodiesCLI )

# Qt4 interface for the three dimensional rigid body simulation
#if( USE_QT4 )
#  add_subdirectory( ThreeDRigidBodiesQt4 )
#endif()

# Command line interface for computing inertia of 3D meshes
#add_subdirectory( ThreeDRigidBodiesInertiaComputer )

# Unit tests for two dimensional ball simulation
#add_subdirectory( ThreeDRigidBodiesTests )

add_subdirectory( mpm2d )
