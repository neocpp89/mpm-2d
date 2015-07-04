add_executable(particle_movement particle_movement.c)
target_link_libraries(particle_movement mpm)
add_test(test_particle_movement particle_movement)

add_executable(cgsolver cgsolver.c ../materialsrc/cs_cg.c)
include_directories(../materialsrc/)
target_link_libraries(cgsolver mpm)
target_link_libraries(cgsolver ${CXSPARSE_LIBRARY})
add_test(test_cgsolver cgsolver)

# A bit of hack since we need to sync with the material directory
set(LOCAL_MATERIALS
    dp_rd.c
    dp_ri.c
    dp_ri_mt.c
    g_local.c
    g_local_mu2.c
    g_local_mu2_ext.c
    g_local_mu2_plane_strain.c
    g_local_mu2_reverse.c
    g_local_mu2_true_2d.c
    isolin.c
)

set(NONLOCAL_MATERIALS
    g_nonlocal_avg.c
    g_nonlocal_bicgstab.c
    g_nonlocal.c
    g_nonlocal_mu2.c
    g_nonlocalxs71.c
)

foreach(matc ${LOCAL_MATERIALS})
    string(REGEX REPLACE "\\.c$" "" mat ${matc})
    add_executable(material_${mat} test_material.c ../materialsrc/${matc})
    target_link_libraries(material_${mat} pthread)
    target_link_libraries(material_${mat} m)
    target_link_libraries(material_${mat} mpm)
    configure_file(${mat}.cfg ${CMAKE_CURRENT_BINARY_DIR}/${mat}.cfg COPYONLY)
    add_test(test_material_${mat} material_${mat} ${mat}.cfg)
endforeach(matc)

# We need a different test for the nonlocal materials...
#foreach(matc ${NONLOCAL_MATERIALS})
#    string(REGEX REPLACE "\\.c$" "" mat ${matc})
#    add_executable(material_${mat} test_material.c ../materialsrc/${matc})
#    target_link_libraries(material_${mat} pthread)
#    target_link_libraries(material_${mat} m)
#    target_link_libraries(material_${mat} mpm)
#    target_link_libraries(material_${mat} ${CXSPARSE_LIBRARY})
#    add_test(test_material_${mat} material_${mat})
#endforeach(matc)
