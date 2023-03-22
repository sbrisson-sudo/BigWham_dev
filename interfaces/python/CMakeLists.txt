pybind11_add_module(bigwhamPybind bigwham_pybind.cc)
target_link_libraries(bigwhamPybind PRIVATE BigWhamShared)
set_target_properties(bigwhamPybind
  PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH ON
    INSTALL_RPATH "${CBLAS_LIBRARY_DIRS};${CMAKE_INSTALL_PREFIX}/lib"
    BUILD_RPATH "${CBLAS_LIBRARY_DIRS};${PROJECT_BINARY_DIR}"
)

set(SRC_PYBIND11 "${CMAKE_CURRENT_SOURCE_DIR}/bigwham4py.py"
                 "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py")
add_custom_target(copy_pybind11_files ALL COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_PYBIND11}
                                           ${CMAKE_BINARY_DIR}/interfaces/python)
