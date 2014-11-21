ADD_EXECUTABLE(mpm_2d
    main.c
    element.c
    implicit.c
    interpolate.c
    loading.c
    map.c
    material.c
    particle.c
    process.c
    reader.c
    rtsafe.c
    tensor.c
    writer.c
)
target_link_libraries(mpm_2d ${CXSPARSE_LIBRARY})
target_link_libraries(mpm_2d m)
target_link_libraries(mpm_2d pthread)
target_link_libraries(mpm_2d ${CONFUSE_LIBRARY})
target_link_libraries(mpm_2d ${DL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(mpm_2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
