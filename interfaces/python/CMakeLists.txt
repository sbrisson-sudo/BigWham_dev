pybind11_add_module(py_bigwham py_bigwham.cc)
target_link_libraries(py_bigwham PRIVATE BigWhamShared)
set_target_properties(py_bigwham 
  PROPERTIES 
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/interfaces/python/bigwham4py"
    INSTALL_RPATH_USE_LINK_PATH ON
    INSTALL_RPATH "${CBLAS_LIBRARY_DIRS};${CMAKE_INSTALL_PREFIX}/lib"
    BUILD_RPATH "${CBLAS_LIBRARY_DIRS};${PROJECT_BINARY_DIR}"
)

set(PYTHON_PACKAGE_DIR "${CMAKE_BINARY_DIR}/interfaces/python")
set(BIGWHAM4PY_DIR "${PYTHON_PACKAGE_DIR}/bigwham4py")

# Ensure the destination directories exist
file(MAKE_DIRECTORY ${PYTHON_PACKAGE_DIR})
file(MAKE_DIRECTORY ${BIGWHAM4PY_DIR})

# Define the source files
set(SRC_SETUP_PY "${CMAKE_CURRENT_SOURCE_DIR}/setup.py")
set(SRC_BIGWHAM4PY
    "${CMAKE_CURRENT_SOURCE_DIR}/bigwham4py/bigwham4py.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/bigwham4py/__init__.py")

# Copy setup.py to interfaces/python/
add_custom_target(copy_setup_py ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_SETUP_PY} ${PYTHON_PACKAGE_DIR}/setup.py
  DEPENDS ${SRC_SETUP_PY})

# Copy the bigwham4py package files to interfaces/python/bigwham4py/
add_custom_target(copy_bigwham4py_files ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_BIGWHAM4PY} ${BIGWHAM4PY_DIR}/
  DEPENDS ${SRC_BIGWHAM4PY})

# Ensure the compiled module is placed inside the `bigwham4py` package
set_target_properties(py_bigwham 
  PROPERTIES 
    LIBRARY_OUTPUT_DIRECTORY "${BIGWHAM4PY_DIR}"
    INSTALL_RPATH_USE_LINK_PATH ON
    INSTALL_RPATH "${CBLAS_LIBRARY_DIRS};${CMAKE_INSTALL_PREFIX}/lib"
    BUILD_RPATH "${CBLAS_LIBRARY_DIRS};${PROJECT_BINARY_DIR}"
)

message(STATUS "To use Bigwham python interface bigwham4py, you can :
1. add ${CMAKE_CURRENT_SOURCE_DIR}/bigwham4py to your python path (via the PYTHONPATH environment variable or via sys.path.append())
2. build and install it via pip, for that go in ${CMAKE_CURRENT_SOURCE_DIR} and execute :
python3 setup.py build
python3 -m pip install .
")