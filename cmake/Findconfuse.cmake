# Locate Confuse include paths and libraries
#
# $Id$
# Copyright (C) Sebastien Tricaud 2012
#
# CONFUSE_FOUND        - If libconfuse is found
# CONFUSE_INCLUDE_DIRS - Where include/confuse is found
# CONFUSE_LIBRARIES    - List of libraries when using libconfuse
# CONFUSE_DEFINITIONS  - List of definitions to be added when using libconfuse
#

set(CONFUSE_DEFINITIONS "")

if(CONFUSE_INLINE_BUILD)
	if(WIN32)
		set(CONFUSE_FOUND true)
		set(CONFUSE_INCLUDE_DIR "${confuse_SOURCE_DIR}\\src\\include")
		set(CONFUSE_LIBRARY "${confuse_SOURCE_DIR}\\src\\${CMAKE_BUILD_TYPE}\\confuse.lib")
	else(WIN32)
		set(CONFUSE_FOUND true)
		set(CONFUSE_INCLUDE_DIR "${confuse_SOURCE_DIR}/src/include")
		set(CONFUSE_LIBRARY "${confuse_SOURCE_DIR}/src/libconfuse.so")
	endif(WIN32)
else(CONFUSE_INLINE_BUILD)
	find_path(CONFUSE_INCLUDE_DIR confuse.h confuse/confuse.h
	          HINTS "../src/include" ${CONFUSE_INCLUDEDIR}
	          PATH_SUFFIXES confuse )

	find_library(CONFUSE_LIBRARY NAMES confuse
	             HINTS "../src/" ${CONFUSE_LIBDIR} 
		     PATH_SUFFIXES src/ )
endif(CONFUSE_INLINE_BUILD)

set(CONFUSE_LIBRARIES ${CONFUSE_LIBRARY})
set(CONFUSE_INCLUDE_DIRS ${CONFUSE_INCLUDE_DIR})

mark_as_advanced(CONFUSE_INCLUDE_DIRS CONFUSE_LIBRARIES)

