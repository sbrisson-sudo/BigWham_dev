find_package(ZLIB REQUIRED)

option(BIGWHAM_OPENMP "OPENMP" ON) # MUST BE ON for openMP

# add_library(cnpy SHARED cnpy.cpp)
# target_include_directories(cnpy PRIVATE ${ZLIB_INCLUDE_DIRS})
# target_link_libraries(cnpy PRIVATE ${ZLIB_LIBRARIES})

add_executable(rep_test  rep_test.cpp)
target_link_libraries(rep_test PRIVATE BigWhamShared cnpy)

# add_executable(rep_test_rowhmat  rep_test_rowhmat.cpp)
# target_link_libraries(rep_test_rowhmat PRIVATE BigWhamShared cnpy)

# add_executable(build_hmat build_hmat.cc)
# target_link_libraries(build_hmat PRIVATE BigWhamShared cnpy)

add_executable(mat_vec mat_vec.cc)
target_link_libraries(mat_vec PRIVATE BigWhamShared cnpy)

if (USE_CUDA)
    add_executable(mat_vec_gpu mat_vec_gpu.cc)
    target_link_libraries(mat_vec_gpu PRIVATE BigWhamShared cnpy)
    target_link_libraries(mat_vec_gpu PRIVATE ${CUDA_LIBRARIES})

    add_executable(mat_vec_gpu_getMemoryRequirement mat_vec_gpu_getMemoryRequirement.cc)
    target_link_libraries(mat_vec_gpu_getMemoryRequirement PRIVATE BigWhamShared cnpy)
    target_link_libraries(mat_vec_gpu_getMemoryRequirement PRIVATE ${CUDA_LIBRARIES})
endif()

add_executable(mat_vec_selection mat_vec_selection.cc)
target_link_libraries(mat_vec_selection PRIVATE BigWhamShared cnpy)

add_executable(test_contiguous_array2d test_contiguous_array2d_vector.cpp)
target_include_directories(test_contiguous_array2d PUBLIC ../../il ../../src)


configure_file(generate_penny_tri_mesh.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(generate_line_segment_mesh.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(generate_2d_fracture_bulk_mesh.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# configure_file(mesh_conn.npy ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
# configure_file(mesh_coords.npy ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(rep_test.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

if (VTUNE_PROFILING)

    if (ITT_INCLUDE_DIR AND ITT_LIBRARY)
        message(STATUS "Found Intel ITT: ${ITT_LIBRARY}")
        add_library(ittnotify STATIC IMPORTED)
        set_target_properties(ittnotify PROPERTIES IMPORTED_LOCATION ${ITT_LIBRARY})
        target_include_directories(mat_vec PRIVATE ${ITT_INCLUDE_DIR})
        target_link_libraries(mat_vec PRIVATE ittnotify)
        target_compile_definitions(mat_vec PRIVATE USE_ITT=1)
    else()
        message(FATAL_ERROR "Intel ITT not found. Please install Intel VTune Profiler.")
    endif()

endif()

add_executable(dense_matrix_2DS0_H dense_matrix_2DS0_H.cc)
target_link_libraries(dense_matrix_2DS0_H PRIVATE InsideLoop cnpy BigWhamShared)
target_include_directories(dense_matrix_2DS0_H PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/il
    ${CMAKE_SOURCE_DIR}/external/cnpy
    ${CMAKE_SOURCE_DIR}/external/nlohmann
)

add_executable(dense_matrix_2DS0T0_V dense_matrix_2DS0T0_V.cc)
target_link_libraries(dense_matrix_2DS0T0_V PRIVATE InsideLoop cnpy BigWhamShared GSL::gsl)
target_include_directories(dense_matrix_2DS0T0_V PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/il
    ${CMAKE_SOURCE_DIR}/external/cnpy
    ${CMAKE_SOURCE_DIR}/external/nlohmann
)

add_executable(hr_pattern_2DS0T0 hr_pattern_2DS0T0.cc)
target_link_libraries(hr_pattern_2DS0T0 PRIVATE InsideLoop cnpy BigWhamShared GSL::gsl)
target_include_directories(hr_pattern_2DS0T0 PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/il
    ${CMAKE_SOURCE_DIR}/external/cnpy
    ${CMAKE_SOURCE_DIR}/external/nlohmann
)

add_executable(matvec_2DS0T0_V matvec_2DS0T0_V.cc)
target_link_libraries(matvec_2DS0T0_V PRIVATE InsideLoop cnpy BigWhamShared GSL::gsl)
target_include_directories(matvec_2DS0T0_V PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/il
    ${CMAKE_SOURCE_DIR}/external/cnpy
    ${CMAKE_SOURCE_DIR}/external/nlohmann
)