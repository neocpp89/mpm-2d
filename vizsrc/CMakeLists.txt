ADD_EXECUTABLE(viz
    viz.cpp
    tokenizer.cpp
    viz_colormap.cpp
    viz_reader.cpp
)
target_include_directories(viz PUBLIC ${FTGL_INCLUDE_DIR})
target_include_directories(viz PUBLIC ${FREETYPE_INCLUDE_DIRS})
target_include_directories(viz PUBLIC ${OPENGL_INCLUDE_DIRS})
target_include_directories(viz PUBLIC ${PNG_INCLUDE_DIRS})
target_include_directories(viz PUBLIC ${SDL_INCLUDE_DIRS})
target_link_libraries(viz ${CONFUSE_LIBRARY})
target_link_libraries(viz ${FTGL_LIBRARIES})
target_link_libraries(viz ${FREETYPE_LIBRARIES})
target_link_libraries(viz ${OPENGL_LIBRARIES})
target_link_libraries(viz ${PNG_LIBRARIES})
target_link_libraries(viz ${SDL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(viz PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
