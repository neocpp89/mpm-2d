# mpm2d CMakeLists.txt
# 
# Breannan Smith (smith@cs.columbia.edu) 
# Updated: 11/05/2013

include( CMakeSourceFiles.txt )

# libConfuse (http://www.nongnu.org/confuse/) is required
find_package( confuse REQUIRED )
include_directories( ${CONFUSE_INCLUDE_DIRS} )

# BLAS is required
find_package( BLAS REQUIRED )

# SuireSparse is required
find_package( SuiteSparse REQUIRED )
include_directories( ${CSPARSE_INCLUDE_DIR} )

# HDF5 is required
find_package( HDF5 REQUIRED )

add_executable( mpm_2d ${Headers} ${Sources} )

target_link_libraries( mpm_2d ${CONFUSE_LIBRARIES} ${BLAS_LIBRARIES} ${CSPARSE_LIBRARY} ${HDF5_C_LIBRARIES} )
