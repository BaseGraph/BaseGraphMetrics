find_package(GTest REQUIRED)

set(TEST_FILES_NAMES test_directed test_undirected)

foreach(TEST_NAME ${TEST_FILES_NAMES})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
    target_link_libraries(${TEST_NAME} ${EXTENSION_NAME} BaseGraph::core gtest gtest_main)
    add_test(${TEST_NAME} ${TEST_NAME})
endforeach()
