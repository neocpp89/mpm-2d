add_library(mpm
    element.c
    interpolate.c
    loading.c
    map.c
    material.c
    particle.c
    process_usl.c
    rtsafe.c
    tensor.c
)
target_include_directories(mpm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
