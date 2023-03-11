add_library(${EXTENSION_NAME}
    ${CMAKE_CURRENT_SOURCE_DIR}/directed.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/undirected.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/general.cpp
)

target_include_directories(${EXTENSION_NAME} PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
set_target_properties(${EXTENSION_NAME} PROPERTIES 
    LINKER_LANGUAGE CXX
    POSITION_INDEPENDENT_CODE TRUE
)
target_link_libraries(${EXTENSION_NAME} BaseGraph::core)
