add_library(mpm
    element.c
    loading.c
    material.c
    process_usl.c
    tensor.c
    spmd.c
)
target_include_directories(mpm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
