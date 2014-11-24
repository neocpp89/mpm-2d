CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
PROJECT(mpm-2d)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
SET(CMAKE_C_FLAGS "-fPIC -march=native -std=gnu99 -g -O3 -Wall -Wextra")
SET(CMAKE_CXX_FLAGS "-fPIC -march=native -std=c++11 -g -O3 -Wall -Wextra")
# SET(CMAKE_C_FLAGS "-march=native -std=gnu99 -g -O3 -Wall -Wextra -Werror")

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

install(DIRECTORY examples DESTINATION mpm)
