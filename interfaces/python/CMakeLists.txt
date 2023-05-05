pybind11_add_module(py_bigwham py_bigwham.cc)
target_link_libraries(py_bigwham PRIVATE BigWhamShared)
set_target_properties(py_bigwham
  PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH ON
    INSTALL_RPATH "${CBLAS_LIBRARY_DIRS};${CMAKE_INSTALL_PREFIX}/lib"
    BUILD_RPATH "${CBLAS_LIBRARY_DIRS};${PROJECT_BINARY_DIR}"
)

set(SRC_PYBIND11 "${CMAKE_CURRENT_SOURCE_DIR}/hmatrix.py"
                 "${CMAKE_CURRENT_SOURCE_DIR}/hmatrix_rectangular.py"
                 # "${CMAKE_CURRENT_SOURCE_DIR}/hmatrix_mpi.py"
                 "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py")

message(STATUS "Use the following command in terminal to use BigWham python interface:
export PYTHONPATH=\"${CMAKE_BINARY_DIR}/interfaces/python\"
                    OR
-- Use the following command in your python script
import sys
sys.path.append(\"${CMAKE_BINARY_DIR}/interfaces/python\")
")
add_custom_target(copy_pybind11_files ALL COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_PYBIND11}
                                           ${CMAKE_BINARY_DIR}/interfaces/python)
