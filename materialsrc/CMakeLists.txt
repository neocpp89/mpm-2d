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
    string(REGEX REPLACE "\\.c$" "" matso ${matc})
    add_library(${matso} SHARED ${matc})
    target_link_libraries(${matso} mpm)
    set_target_properties(${matso} PROPERTIES PREFIX "")
    install(TARGETS ${matso} LIBRARY DESTINATION mpm/material)
endforeach(matc)

foreach(matc ${NONLOCAL_MATERIALS})
    string(REGEX REPLACE "\\.c$" "" matso ${matc})
    add_library(${matso} SHARED ${matc})
    target_link_libraries(${matso} mpm)
    target_link_libraries(${matso} ${CXSPARSE_LIBRARY})
    set_target_properties(${matso} PROPERTIES PREFIX "")
    install(TARGETS ${matso} LIBRARY DESTINATION mpm/material)
endforeach(matc)
