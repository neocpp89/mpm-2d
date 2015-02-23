set(SOURCES
    main.c
    reader.c
    writer.c
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.c.in" "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.c" @ONLY)
list(APPEND SOURCES "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.c" GitSHA1.h)
ADD_EXECUTABLE(mpm_2d ${SOURCES})
target_link_libraries(mpm_2d mpm)
target_link_libraries(mpm_2d ${CXSPARSE_LIBRARY})
target_link_libraries(mpm_2d m)
target_link_libraries(mpm_2d rt) # for clock_gettime
target_link_libraries(mpm_2d pthread)
target_link_libraries(mpm_2d ${CONFUSE_LIBRARY})
target_link_libraries(mpm_2d ${DL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(mpm_2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(mpm_2d PROPERTIES INSTALL_RPATH "$ORIGIN/../mpm/material/:$ORIGIN/../mpm/bc/:$ORIGIN/../mpm/:$ORIGIN/")
install(TARGETS mpm_2d RUNTIME DESTINATION bin)
