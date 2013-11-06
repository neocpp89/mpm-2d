# viz CMakeLists.txt
# 
# Breannan Smith (smith@cs.columbia.edu) 
# Updated: 11/05/2013

include( CMakeSourceFiles.txt )

# Freetype is required
find_package( Freetype REQUIRED )
include_directories( ${FREETYPE_INCLUDE_DIRS} )

# FTGL is required
find_package( FTGL REQUIRED )

# OpenGL is required
find_package( OpenGL REQUIRED )

# SDL is required
find_package( SDL REQUIRED )
include_directories( ${SDL_INCLUDE_DIR} )

# libpng is required
find_package( PNG REQUIRED )
include_directories( ${PNG_INCLUDE_DIRS} )

# libConfuse (http://www.nongnu.org/confuse/) is required
find_package( confuse REQUIRED )

# BLAS is required
#find_package( BLAS REQUIRED )

# SuireSparse is required
#find_package( SuiteSparse REQUIRED )
#include_directories( ${CSPARSE_INCLUDE_DIR} )

# HDF5 is required
#find_package( HDF5 REQUIRED )

add_executable( viz ${Headers} ${Sources} )

target_link_libraries( viz ${SDL_LIBRARY} ${CONFUSE_LIBRARIES} ${FREETYPE_LIBRARIES} ${FTGL_LIBRARY} ${OPENGL_LIBRARIES} )

#target_link_libraries( mpm_2d ${CONFUSE_LIBRARIES} ${BLAS_LIBRARIES} ${CSPARSE_LIBRARY} ${HDF5_C_LIBRARIES} )
