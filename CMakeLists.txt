set(EXTENSION_NAME metrics)
add_library(${EXTENSION_NAME} directed.cpp undirected.cpp general.cpp)
target_link_libraries(${EXTENSION_NAME} core)
add_extension(${EXTENSION_NAME})


if (BUILD_BINDINGS)
    pybind11_add_module(_${EXTENSION_NAME} "${CMAKE_CURRENT_SOURCE_DIR}/bind_metrics.cpp")
    target_link_libraries(_${EXTENSION_NAME} PRIVATE ${EXTENSION_NAME})
    add_pybind_extension(_${EXTENSION_NAME})
endif()

if (BUILD_TESTS)
    file(GLOB METRICS_TEST_SRC "${CMAKE_CURRENT_SOURCE_DIR}/test_metrics/*.cpp")

    add_executable(test_metrics ${METRICS_TEST_SRC})
    target_link_libraries(test_metrics ${EXTENSION_NAME})
    add_gtest(test_metrics)
endif()
