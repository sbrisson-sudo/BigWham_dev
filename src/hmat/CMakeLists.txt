set(TMP_SOURCES_FILES ${SOURCE_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/cluster/cluster.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/hmatrix/h_pattern.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/hmatrix/hmat.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/hmatrix/hmat_selection.cc
    # ${CMAKE_CURRENT_SOURCE_DIR}/square_matrix_generator.h
    )

if(USE_CUDA)
    list(APPEND TMP_SOURCES_FILES ${CMAKE_CURRENT_SOURCE_DIR}/hmatrix/hmat_cuda.cc)
    list(APPEND TMP_SOURCES_FILES ${CMAKE_CURRENT_SOURCE_DIR}/hmatrix/cuda_kernels_utilities.cu)
endif()

set(SOURCE_FILES ${TMP_SOURCES_FILES} PARENT_SCOPE)