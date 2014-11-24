ADD_EXECUTABLE(mpm_viz
    viz.cpp
    tokenizer.cpp
    viz_colormap.cpp
    viz_reader.cpp
)
target_include_directories(mpm_viz PUBLIC ${FTGL_INCLUDE_DIR})
target_include_directories(mpm_viz PUBLIC ${FREETYPE_INCLUDE_DIRS})
target_include_directories(mpm_viz PUBLIC ${OPENGL_INCLUDE_DIRS})
target_include_directories(mpm_viz PUBLIC ${PNG_INCLUDE_DIRS})
target_include_directories(mpm_viz PUBLIC ${SDL_INCLUDE_DIRS})
target_link_libraries(mpm_viz ${CONFUSE_LIBRARY})
target_link_libraries(mpm_viz ${FTGL_LIBRARIES})
target_link_libraries(mpm_viz ${FREETYPE_LIBRARIES})
target_link_libraries(mpm_viz ${OPENGL_LIBRARIES})
target_link_libraries(mpm_viz ${PNG_LIBRARIES})
target_link_libraries(mpm_viz ${SDL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(mpm_viz PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS mpm_viz RUNTIME DESTINATION bin)
