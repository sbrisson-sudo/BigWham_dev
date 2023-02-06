add_executable(rep_test
  rep_test.cpp
  cnpy.cpp)

find_package(ZLIB REQUIRED)
target_include_directories(rep_test PRIVATE BigWhamShared ${ZLIB_INCLUDE_DIRS})
target_link_libraries(rep_test PRIVATE BigWhamShared ${ZLIB_LIBRARIES})

configure_file(generate_penny_mesh.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(rep_test PRIVATE -Wno-narrowing)
endif()
