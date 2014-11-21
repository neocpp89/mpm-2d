CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
PROJECT(mpm-2d)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
SET(CMAKE_C_FLAGS "-march=native -std=gnu99 -g -O3 -Wall -Wextra")
SET(CMAKE_CXX_FLAGS "-march=native -std=c++11 -g -O3 -Wall -Wextra")
# SET(CMAKE_C_FLAGS "-march=native -std=gnu99 -g -O3 -Wall -Wextra -Werror")

FIND_PACKAGE(CXSparse REQUIRED)
FIND_PACKAGE(confuse REQUIRED)
FIND_PACKAGE(DL REQUIRED)

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
    ${PROJECT_SOURCE_DIR}/src
)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(materialsrc)
ADD_SUBDIRECTORY(bcsrc)
if (FTGL_FOUND AND FREETYPE_FOUND
    AND OPENGL_FOUND
    AND PNG_FOUND
    AND SDL_FOUND)
    ADD_SUBDIRECTORY(vizsrc)
endif()
