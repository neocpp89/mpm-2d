CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
PROJECT(mpm-2d)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -march=native -std=gnu99 -O3 -Wall -Wextra -Wstrict-prototypes")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -march=native -std=c++11 -O3 -Wall -Wextra")
# SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -Werror")

# for all binaries:
# SET(CMAKE_INSTALL_RPATH "$ORIGIN/../mpm/material/:$ORIGIN/../mpm/bc/:$ORIGIN/")

# default install is to a subdirectory of build dir
# override with `cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..`
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/installed" CACHE PATH "Default install path" FORCE)
endif()

# for main library
FIND_PACKAGE(CXSparse REQUIRED)

# for driver program
FIND_PACKAGE(confuse REQUIRED)
FIND_PACKAGE(DL REQUIRED)

# for visualization program
FIND_PACKAGE(FTGL)
FIND_PACKAGE(Freetype)
FIND_PACKAGE(OpenGL)
FIND_PACKAGE(png)
FIND_PACKAGE(SDL)

INCLUDE_DIRECTORIES(${CXSPARSE_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CONFUSE_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${DL_INCLUDE_DIR})
enable_testing()

if (CMAKE_CONFIGURATION_TYPES)
    add_custom_target(vtest COMMAND ${CMAKE_CTEST_COMMAND}
        --force-new-ctest-process --output-on-failure
        --build-config "$<CONFIGURATION>")
else ()
    add_custom_target(vtest COMMAND ${CMAKE_CTEST_COMMAND}
        --force-new-ctest-process --output-on-failure)
endif ()

INCLUDE_DIRECTORIES(
    ${PROJECT_SOURCE_DIR}/libmpm
    ${PROJECT_SOURCE_DIR}/src
)
ADD_SUBDIRECTORY(libmpm)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(materialsrc)
ADD_SUBDIRECTORY(bcsrc)
ADD_SUBDIRECTORY(tests)
if (FTGL_FOUND AND FREETYPE_FOUND
    AND OPENGL_FOUND
    AND PNG_FOUND
    AND SDL_FOUND)
    ADD_SUBDIRECTORY(vizsrc)
endif()

# create output directories for examples
file(MAKE_DIRECTORY examples/silo/output)
file(MAKE_DIRECTORY examples/silo/rycroft-output)
file(MAKE_DIRECTORY examples/heap/output)
file(MAKE_DIRECTORY examples/square-hopper/output)
file(MAKE_DIRECTORY examples/hopper/output)

# copy examples over
install(DIRECTORY examples DESTINATION mpm)

# generate documentation for libmpm
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

