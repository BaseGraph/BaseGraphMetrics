find_package(pybind11 REQUIRED)

pybind11_add_module(_${EXTENSION_NAME} MODULE pybind.cpp)
target_link_libraries(_${EXTENSION_NAME} PRIVATE BaseGraph::core ${EXTENSION_NAME})

install(TARGETS _${EXTENSION_NAME} DESTINATION ./basegraph/${EXTENSION_NAME})
install(FILES ${EXTENSION_NAME}/__init__.py DESTINATION ./basegraph/${EXTENSION_NAME})
