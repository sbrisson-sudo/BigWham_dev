#####################################
# BIGWHAM                           #
# compilation configuration         #
#####################################
cmake_minimum_required(VERSION 3.5.1)

if(CMAKE_VERSION VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(BigWham VERSION 1.0
        DESCRIPTION "This is the library containing the wheels"
        LANGUAGES CXX)

# This is only for pybind11 to behave properly
set(CMAKE_CXX_STANDARD 17)


# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

message(STATUS "CMAKE_BUILD_TYPE :  ${CMAKE_BUILD_TYPE}")

#if(NOT CMAKE_BUILD_TYPE)
#    message(STATUS "Build Type not set, defaulting to Release..." )
#    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build Type")
#endif()

# set options: False = 0, True = 1  #
#####################################
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(BIGWHAM_OPENMP "OPENMP" ON) # MUST BE ON for openMP

option(BIGWHAM_PYTHON_INTERFACE "Build Pybind11 binding" ON)
option(BIGWHAM_JULIA_INTERFACE "Build Julia binding" OFF)
option(BIGWHAM_MATHEMATICA_INTERFACE "make mathematica interface" OFF)  # --> to add static option to CMAKE_CXX_FLAGS
option(BIGWHAM_TESTING "Build google tests" OFF)
option(BIGWHAM_HDF5 "Use HDF5 for checkpoint/restart" OFF)
option(BIGWHAM_DEV_EXAMPLE OFF)

#option(USE_STATIC "use static" OFF)  # --> to add static option to CMAKE_CXX_FLAGS

#####################################
# optimization options              #
#####################################
set(_optimization_flags "-O3 -DNDEBUG")
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(_optimization_flags "${_optimization_flags} -xHost")
else()
    if (BIGWHAM_MATHEMATICA_INTERFACE)
        set(_optimization_flags "${_optimization_flags}")
    else()
        set(_optimization_flags "${_optimization_flags} -march=native")
    endif()
endif()
set(CMAKE_CXX_FLAGS_RELEASE "${_optimization_flags}" CACHE STRING "Flags used by the compiler during release builds" FORCE)


#if(USE_STATIC)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
#endif()

################################################################################
# set compiler: clang, gcc, intel etc.
################################################################################
message(STATUS "Detecting the compiler:")
message(STATUS "  The compiler is:  ${CMAKE_CXX_COMPILER}")

if(BIGWHAM_MATHEMATICA_INTERFACE)
  set(IL_MATH_VENDOR "MKL_sequential" CACHE INTERNAL "" FORCE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(BIGWHAM_OPENMP OFF CACHE INTERNAL "" FORCE)
  endif()
endif()

add_subdirectory(il)
add_subdirectory(cnpy)

set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_SKIP_BUILD_RPATH OFF)
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${CBLAS_LIBRARY_DIRS}")
#set(CMAKE_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${CBLAS_LIBRARY_DIRS}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

################################################
# TO PROFILE WITH VTUNE
################################################
option(VTUNE_PROFILING "Enable Intel VTune Profiling" OFF)

if (VTUNE_PROFILING)
    find_path(ITT_INCLUDE_DIR ittnotify.h PATHS 
        $ENV{VTUNE_AMPLIFIER_DIR}/include
        /opt/intel/oneapi/vtune/latest/sdk/include
        /opt/intel/vtune_profiler/sdk/include
    )
    find_library(ITT_LIBRARY NAMES ittnotify PATHS 
        $ENV{VTUNE_AMPLIFIER_DIR}/lib64
        /opt/intel/oneapi/vtune/latest/sdk/lib64
        /opt/intel/vtune_profiler/sdk/lib64
    )
    message(STATUS "Found Intel ITT: ${ITT_INCLUDE_DIR}, ${ITT_LIBRARY}")
endif()

################################################
# TO LINK CUDA LIBRARIES
if (USE_CUDA)
    SET(CUDA_SEPARABLE_COMPILATION ON) # cf https://stackoverflow.com/questions/22115197/dynamic-parallelism-undefined-reference-to-cudaregisterlinkedbinary-linking

    find_package(CUDA REQUIRED)
    enable_language(CUDA)
    set(CMAKE_CUDA_ARCHITECTURES 70 75 80 86)

    message(STATUS "Building with CUDA support : ${CUDA_TOOLKIT_ROOT_DIR}")
    message(STATUS "CUDA_LIBRARIES: ${CUDA_LIBRARIES}")
    message(STATUS "CUDA_CUBLAS_LIBRARIES: ${CUDA_CUBLAS_LIBRARIES}")
    message(STATUS "CUDA_cusparse_LIBRARY: ${CUDA_cusparse_LIBRARY}")
    message(STATUS "CUDA_cudart_static_LIBRARY: ${CUDA_cudart_static_LIBRARY}")

    include_directories(
        ${CBLAS_INCLUDE_DIR}
        ${CUDA_INCLUDE_DIRS}
        )
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++14)

    # # CUDA NVTX : only for debugging purposes
    # set(NVTX_INCLUDE_DIR ${CUDA_TOOLKIT_ROOT_DIR}/nsight-compute-2024.1.1/host/target-linux-x64/nvtx/include)
    # message(STATUS "NVTX_INCLUDE_DIR: ${NVTX_INCLUDE_DIR}")

    # include_directories(${NVTX_INCLUDE_DIR})

    ################################################
    # TO LINK MAGMA LIBRARY
    # Allow the user to specify MAGMA_ROOT via command line
    option(MAGMA_ROOT "Path to MAGMA installation directory" "")

    # Check if MAGMA_ROOT was provided
    if(MAGMA_ROOT)
        # Set variables based on MAGMA_ROOT
        message(STATUS "MAGMA support enabled. Using MAGMA at: ${MAGMA_ROOT}")
        set(MAGMA_LIBRARY ${MAGMA_ROOT}/lib/libmagma.so)
        set(MAGMA_INCLUDE_DIR ${MAGMA_ROOT}/include)
        message(STATUS "MAGMA_LIBRARY : ${MAGMA_LIBRARY}")
        message(STATUS "MAGMA_INCLUDE_DIR : ${MAGMA_INCLUDE_DIR}")
    else()
        message(STATUS "MAGMA support disabled. To enable, specify -DMAGMA_ROOT=/path/to/magma")
    endif()
    ################################################

endif()
################################################


################################################################################
# file definitions
################################################################################

# Source files are added by src subdirectories
set(SOURCE_FILES "")
add_subdirectory(src)

# message(STATUS ${SOURCE_FILES})

add_library(BigWham OBJECT ${SOURCE_FILES})
target_compile_features(BigWham PUBLIC cxx_std_17)
target_compile_definitions(BigWham PUBLIC BIGWHAM_OPENMP)
target_link_libraries(BigWham PUBLIC InsideLoop cnpy)
set_target_properties(BigWham PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(BigWham PUBLIC src il cnpy)

if (USE_CUDA)
    target_link_libraries(BigWham PUBLIC
        dl pthread rt # https://stackoverflow.com/questions/62962539/linking-error-for-cuda-separate-compilation-mode-and-static-linking-of-cudart
        ${BLAS_LIBRARIES}
        ${CUDA_LIBRARIES}
        ${CUDA_CUBLAS_LIBRARIES}
        ${CUDA_cusparse_LIBRARY}
        ${CUDA_cudart_static_LIBRARY}
    )
    target_compile_definitions(BigWham PRIVATE USE_CUDA)

    # set_target_properties(BigWham PROPERTIES
    # CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET BigWham PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)

    if(MAGMA_ROOT)
        target_include_directories(BigWham PRIVATE ${MAGMA_INCLUDE_DIR})
        target_link_libraries(BigWham PRIVATE ${MAGMA_LIBRARY})
        target_compile_definitions(BigWham PRIVATE USE_MAGMA)
    endif()

endif()

add_library(BigWhamShared SHARED $<TARGET_OBJECTS:BigWham>)
add_library(BigWhamStatic STATIC $<TARGET_OBJECTS:BigWham>)
foreach(tgt BigWhamShared BigWhamStatic)
    target_link_libraries(${tgt} PUBLIC BigWham)
    set_target_properties(${tgt} PROPERTIES OUTPUT_NAME BigWham)
endforeach()

set(BIGWHAM_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src;${PROJECT_SOURCE_DIR}/il"
    CACHE PATH "BigWham include dirs")
mark_as_advanced(BIGWHAM_INCLUDE_DIRS)

add_executable(BigWham_exe src/main.cc)
target_link_libraries(BigWham_exe BigWhamShared)
include(BigWhamTools)

if(BIGWHAM_HDF5)
    enable_language(C)
    find_package(HDF5 REQUIRED COMPONENTS C)
    target_compile_definitions(BigWham PUBLIC BIGWHAM_HDF5)
    target_link_libraries(BigWham PUBLIC ${HDF5_LIBRARIES})
    target_include_directories(BigWham PUBLIC ${HDF5_INCLUDE_DIRS})


    if(HDF5_IS_PARALLEL)
        find_package(MPI REQUIRED COMPONENTS C)
        target_compile_options(BigWham PUBLIC ${MPI_C_COMPILE_OPTIONS})
        target_link_libraries(BigWham PUBLIC ${MPI_C_LIBRARIES})
        target_include_directories(BigWham PUBLIC ${MPI_C_INCLUDE_DIRS})
    endif()
endif()

################################################################
# Julia
################################################################
if(BIGWHAM_JULIA_INTERFACE)
    message(STATUS "Building Julia interface: ${BIGWHAM_JULIA_INTERFACE}")
    add_subdirectory(interfaces/julia)
endif()


################################################################
# Pybind --->use with PYTHON version >=3.7
################################################################
if(BIGWHAM_PYTHON_INTERFACE)
    message(STATUS "Building python interface: ${BIGWHAM_PYTHON_INTERFACE}")
    # find_package (Python3 COMPONENTS Interpreter Development)
    # message(STATUS "Python Executable: ${Python3_EXECUTABLE}")
    # message(STATUS "Python Version: ${Python3_VERSION}")
    # message(STATUS "Python Include Dir: ${Python3_INCLUDE_DIRS}")
    find_package (Python REQUIRED COMPONENTS Interpreter Development)
    message(STATUS "Python Executable: ${Python_EXECUTABLE}")
    message(STATUS "Python Version: ${Python_VERSION}")
    message(STATUS "Python Include Dir: ${Python_INCLUDE_DIRS}")
    set(PYBIND11_PYTHON_VERSION 3)
    add_external_package(pybind11 VERSION 2.13.0 CONFIG)
    add_subdirectory(interfaces/python)
endif()

################################################################################
# Mathematica interface
################################################################################
if(BIGWHAM_MATHEMATICA_INTERFACE)
    message(STATUS "Building the mathematica interface: ${BIGWHAM_MATHEMATICA_INTERFACE}")
    add_external_package(ltemplate VERSION 0.5.4)
    add_subdirectory(interfaces/mathematica)
endif()

################################################################################
# For Unit Tests
################################################################################
if(BIGWHAM_TESTING)
    set(TESTED_SOURCE_FILES ${SOURCE_FILES})
    set(UNIT_TEST_FILES
        il/gtest/src/gtest-all.cc
            src/elasticity/_test/penny_shaped_crack_analytical_sol.cc
            src/core/_test/FaceData_tests.cpp
            src/elasticity/fullspace_iso_3d_rectangle/_test/elastic3DR0_element_tests.cc
         # src/_test/elastic3DR0_element_benchmark.cpp
            src/elasticity/_test/penny_shaped_crack_analytical_sol_benchmark.cc
            src/elasticity/_test/elastic_square_matrix_generator_test.cc
            src/hmat/hmatrix/_test/Hmat_test_gaussian_mat.cpp
            src/core/_test/bemesh_test.cc
            src/core/_test/hierarchical_representation_test.cc
            src/elements/_test/triangle_test.cc
            src/elements/_test/rectangle_test.cc
            src/elements/_test/point_test.cc
            src/elements/_test/segment_test.cc
            src/elasticity/fullspace_iso_2d_segment/_test/plane_elasticity_segment_test.cc
            src/elasticity/fullspace_iso_2d_segment/_test/segment0_kernels_test.cc
            src/elasticity/fullspace_iso_sp3d_segment/_test/segment0_sp3d_test.cc
            src/elasticity/fullspace_iso_3d_rectangle/_test/rectangle0_kernels_test.cc
            src/elasticity/fullspace_iso_3d_rectangle/_test/rectangle0_mode1_kernels_test.cc
            src/elasticity/fullspace_iso_3d_triangle/_test/triangle0_kernels_test.cc
            src/elasticity/fullspace_iso_3d_triangle/_test/triangle2_kernels_test.cc
            src/io/_test/bigwham_io_helper_test.cc
            src/io/_test/bigwham_io_gen_test_2d.cc
            src/io/_test/bigwham_io_gen_test_3d.cc
            src/_test/bigwhamio_test_2d.cc
            src/elasticity/fullspace_iso_3d_triangle/_test/triangle0_tests.cc)

    add_executable(BigWhamElastUnitTest ${UNIT_TEST_FILES} src/_test/Gtests-main.cc)
    target_link_libraries(BigWhamElastUnitTest PUBLIC BigWhamShared)

    target_include_directories(BigWhamElastUnitTest
        PUBLIC ${CMAKE_SOURCE_DIR}/il/gtest
        )

    # For unit tests: The precondition of our fonctions are checked with assert
    # macros that terminate the program in debug mode. In order to test those macros
    # in our unit tests, the flag IL_UNIT_TEST turns them into exceptions that can
    # be caught and checked.
    target_compile_definitions(BigWhamElastUnitTest
        PRIVATE IL_UNIT_TEST
                IL_FOLDER='\"${IL_FOLDER}\"')
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(BigWhamElastUnitTest PRIVATE -Wno-terminate -Wno-return-type ) #-Wall -Wextra -pedantic
    endif()
endif()

################################################################################
# For compiling _dev_example, folder with examples for developers
################################################################################
if(BIGWHAM_DEV_EXAMPLES)
    add_subdirectory(src/_dev_example)
endif()


