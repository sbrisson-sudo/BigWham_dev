find_package(ZLIB REQUIRED)

add_executable(rep_test
  rep_test.cpp
  cnpy.cpp)

target_include_directories(rep_test PRIVATE BigWhamShared ${ZLIB_INCLUDE_DIRS})
target_link_libraries(rep_test PRIVATE BigWhamShared ${ZLIB_LIBRARIES})

add_executable(rep_test_rowhmat
  rep_test_rowhmat.cpp
  cnpy.cpp)

target_include_directories(rep_test_rowhmat PRIVATE BigWhamShared ${ZLIB_INCLUDE_DIRS})
target_link_libraries(rep_test_rowhmat PRIVATE BigWhamShared ${ZLIB_LIBRARIES})

configure_file(generate_penny_mesh.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(mesh_conn.npy ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(mesh_coords.npy ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(rep_test.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(rep_test PRIVATE)
endif()
