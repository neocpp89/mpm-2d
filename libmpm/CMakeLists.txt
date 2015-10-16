add_library(mpm
    element.c
    interpolate.c
    loading.c
    material.c
    particle.c
    process_usl.c
    tensor.c
)
target_include_directories(mpm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
