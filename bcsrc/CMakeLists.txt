set(BOUNDARY_CONDITIONS
    box.c
    square_silo.c
)

foreach(bcc ${BOUNDARY_CONDITIONS})
    string(REGEX REPLACE "\\.c$" "" bcso ${bcc})
    add_library(${bcso} SHARED ${bcc})
    target_link_libraries(${bcso} mpm)
    set_target_properties(${bcso} PROPERTIES PREFIX "")
endforeach(bcc)

