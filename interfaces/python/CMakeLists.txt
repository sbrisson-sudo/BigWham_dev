pybind11_add_module(bigwhamPybind bigwhamPybind.cpp)
target_link_libraries(bigwhamPybind PRIVATE BigWhamShared)

set(SRC_PYBIND11 "${CMAKE_CURRENT_SOURCE_DIR}/bigwham4py.py")
add_custom_target(copy_pybind11_files ALL COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_PYBIND11}
                                           ${CMAKE_BINARY_DIR}/interfaces/python)
