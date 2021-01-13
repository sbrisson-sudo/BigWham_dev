cmake_minimum_required(VERSION 3.15)
project(BigWham CXX)
project(BigWham_exe CXX)
project(pyparty)
project(BigWhamElastUnitTest)
set(CMAKE_CXX_STANDARD 11)

#####################################
# set options: False = 0, True = 1  #
#####################################
set(CMAKE_VERBOSE_MAKEFILE 1)
set(IL_MKL 1)
set(IL_OPENBLAS 0)
set(IL_OPENMP 0)

set(USE_INTEL 0)  # compiler option
set(USE_STATIC 0) # to add static option to CMAKE_CXX_FLAGS

set(CMAKE_CXX_FLAGS_DEBUG     "${CMAKE_CXX_FLAGS_DEBUG}   -g -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE   "${CMAKE_CXX_FLAGS_RELEASE} -g -std=c++11 -O3 -mavx2")

if(USE_STATIC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
endif()


################################################################################
# choose compiler: clang, gcc, intel ???
################################################################################
#
MESSAGE("Detecting the compiler:")
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    MESSAGE("-- The compiler is: Clang")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    MESSAGE("-- The compiler is: GNU")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    MESSAGE("-- The compiler is: Intel")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    MESSAGE("-- The compiler is: MSVC")
endif()

if(APPLE)
    MESSAGE("-- located at /usr/bin/")
    set(CMAKE_C_COMPILER /usr/bin/clang)
    set(CMAKE_CXX_COMPILER /usr/bin/clang++)
elseif(UNIX)
    if(USE_INTEL)
        set(CMAKE_C_COMPILER /opt/intel/compilers_and_libraries/linux/bin/intel64/icc)
        set(CMAKE_CXX_COMPILER /opt/intel/compilers_and_libraries/linux/bin/intel64/icpc)
    else()
        set(CMAKE_C_COMPILER /usr/bin/gcc)
        set(CMAKE_CXX_COMPILER /usr/bin/g++)
    endif()
elseif(WIN32)
    if(MINGW)
        set(CMAKE_C_COMPILER C:/msys64/mingw64/bin/gcc)
        set(CMAKE_CXX_COMPILER C:/msys64/mingw64/bin/g++)
    endif()
endif()

################################################################################
# parallel framework
################################################################################

# For OpenMP
if(IL_OPENMP)

    if(APPLE)
        if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            MESSAGE("-- OPENMP ENABLED")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
            #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNUMBEROFTHREADS=4")
            #             link_directories(/opt/clang-3.9.0/lib)
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            MESSAGE("-- OPENMP ENABLED")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        endif()
    elseif(UNIX)
        if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            MESSAGE("-- OPENMP ENABLED")
            # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
            # link_directories(/opt/clang-3.9.0/lib)
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            MESSAGE("-- OPENMP ENABLED")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            MESSAGE("-- OPENMP ENABLED")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
        endif()
    elseif(WIN32)
        if(MINGW)
            MESSAGE("-- OPENMP ENABLED")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        endif()
    endif()
endif()

# For TBB
if(IL_TBB)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_TBB")
    # target_link_libraries(InsideLoop tbb)
endif()

################################################################################
# math framework
################################################################################

# For Intel MKL
if(IL_MKL)
    MESSAGE("Detecting MKL")
    if(APPLE)
        if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            MESSAGE("-- The MKL is found at: " $ENV{MKLROOT}/lib/)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DIL_MKL -DIL_BLAS")
            #           set(CMAKE_MKL_LIBRARIES mkl_intel_lp64 mkl_core mkl_sequential tbb dl pthread m) #iomp5
            set(CMAKE_MKL_LIBRARIES $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_core.a $ENV{MKLROOT}/lib/libmkl_sequential.a  tbb dl  pthread  m) #mkl_intel_thread mkl_sequential  iomp5
            #            link_directories($ENV{MKLROOT}/lib $ENV{MKLROOT}/../compiler/lib $ENV{TBBROOT}/lib)

        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DIL_MKL -DIL_BLAS")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
            set(CMAKE_MKL_LIBRARIES $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_core.a $ENV{MKLROOT}/lib/libmkl_sequential.a  tbb dl  pthread  m) #mkl_intel_thread mkl_sequential  iomp5
            link_directories($ENV{MKLROOT}/lib/intel64_lin)
            link_directories($ENV{MKLROOT}/lib $ENV{MKLROOT}/../compiler/lib )
            # link_directories(/opt/intel/mkl/lib/intel64/opt/intel/compilers_and_libraries/linux/lib/intel64)
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=sequential -DIL_MKL -DIL_BLAS")
            set(CMAKE_MKL_LIBRARIES $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_core.a $ENV{MKLROOT}/lib/libmkl_sequential.a  tbb dl  pthread  m) #mkl_intel_thread mkl_sequential  iomp5
        endif()
    elseif(UNIX)
        if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64  -DIL_MKL -DIL_BLAS")
            set(CMAKE_MKL_LIBRARIES mkl_intel_lp64 mkl_core mkl_intel_thread iomp5 dl pthread m)
            link_directories($ENV{MKLROOT}/lib $ENV{MKLROOT}/../compiler/lib )
            # link_directories(/opt/intel/mkl/lib/intel64 /opt/intel/compilers_and_libraries/linux/lib/intel64)
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -m64 -lmkl_sequential -DIL_MKL -DIL_BLAS") #-m64
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
            set(CMAKE_MKL_LIBRARIES $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a tbb dl pthread m)
            link_directories($ENV{MKLROOT}/lib/intel64 $ENV{MKLROOT}/../compilers_and_libraries/linux/lib/intel64 $ENV{TBBROOT}/lib/intel64/gcc4.8)
            #link_directories($ENV{MKLROOT}/lib $ENV{MKLROOT}/../compiler/lib )
            # link_directories(/opt/intel/mkl/lib/intel64 /opt/intel/compilers_and_libraries/linux/lib/intel64)
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=sequential -DIL_MKL -DIL_BLAS")# -mkl=sequential
            set(CMAKE_MKL_LIBRARIES $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a  tbb dl  pthread  m) #mkl_intel_thread mkl_sequential  iomp5
            link_directories($ENV{MKLROOT}/lib/intel64 $ENV{MKLROOT}/../compilers_and_libraries/linux/lib/intel64)
            #            link_directories(/opt/intel/mkl/lib/intel64 /opt/intel/compilers_and_libraries/linux/lib/intel64)
        endif()
    elseif(WIN32) # to fill up here ...
        if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            # to fill up here ...
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            # to fill up here ...
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            # to fill up here ...
        endif()
    endif()

    # For OpenBlas
elseif(IL_OPENBLAS)
    if(APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_MKL -DIL_OPENBLAS -DIL_BLAS")
        set(CMAKE_OPENBLAS_LIBRARIES "openblas")
    elseif(UNIX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_OPENBLAS -DIL_BLAS")
        set(CMAKE_OPENBLAS_LIBRARIES "openblas")
    elseif(WIN32)
        if(MINGW)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_OPENBLAS -DIL_BLAS")
            set(CMAKE_OPENBLAS_LIBRARIES "openblas")
        endif()
    endif()
endif()


set(SOURCE_FILES_DIR ${PROJECT_SOURCE_DIR}/src)

include_directories(il $ENV{TBBROOT}/include/tbb  src  Hmat-lib .)

MESSAGE("Using TBB:")
MESSAGE("-- The TBB is found at: " $ENV{TBBROOT})

set(SOURCE_FILES
        Hmat-lib/cluster/cluster.h
        Hmat-lib/cluster/cluster.cpp
        # 2D source code
        src/core/Mesh2D.cpp src/core/Mesh2D.h
        src/core/SegmentData.h
        src/elasticity/2d/ElasticS3DP0_element.cpp
        src/elasticity/2d/Elastic2DP1_element.cpp
        src/elasticity/2d/ElasticHMatrix2DP0.h
        src/elasticity/2d/ElasticHMatrix2DP1.h
        src/elasticity/jacobi_prec_assembly.cpp
        src/elasticity/PostProcessDDM.cpp
        src/elasticity/2d/FullMatrixAssembly2D.cpp
        src/elasticity/2d/FullMatrixAssembly2D.h
        # 3D source code
        Hmat-lib/linearAlgebra/factorization/lowRankApproximation.cpp
        src/core/Mesh3D.cpp src/core/Mesh3D.h
        src/core/FaceData.cpp
        #  - T6 element
        src/core/element_utilities.cpp
        src/elasticity/3d/Elastic3DT6_element.cpp
        src/elasticity/3d/h_potential.cpp
        src/elasticity/3d/tensor_utilities.cpp
        src/elasticity/3d/ElasticHMatrix3DT6.h
        #  - R0 element
#        src/elasticity/3d/Elastic3DR0_element.h
        src/elasticity/3d/Elastic3DR0_element.cpp
        src/elasticity/3d/ElasticHMatrix3DR0.h

        # Solver
        src/solvers/HIterativeSolverUtilities.h
        # BigWham Interface
        src/BigWham.h
        )

set(CMAKE_POSITION_INDEPENDENT_CODE True)
add_library(BigWham STATIC ${SOURCE_FILES} )
target_link_libraries(BigWham  ${CMAKE_MKL_LIBRARIES} ) #${CMAKE_OPENBLAS_LIBRARIES}

set_target_properties(BigWham PROPERTIES LINKER_LANGUAGE CXX)
add_executable(BigWham_exe  main.cpp )
include_directories(il $ENV{TBBROOT}/include/  )

target_link_libraries(BigWham_exe   BigWham ) #${CMAKE_OPENBLAS_LIBRARIES}

if(IL_MKL)
    if(APPLE)
        target_link_libraries(BigWham_exe   ${CMAKE_MKL_LIBRARIES}  -lpthread -lm -ldl)

        target_compile_definitions(BigWham PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
        target_include_directories(BigWham PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include $ENV{TBBROOT}/include)
    elseif(UNIX)
        target_link_libraries(BigWham_exe "-Wl,--start-group" ${CMAKE_MKL_LIBRARIES} "-Wl,--end-group" -lpthread -lm -ldl)

        target_compile_definitions(BigWham PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
        target_include_directories(BigWham PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include $ENV{TBBROOT}/include)
        # target_include_directories(HFPx2D PRIVATE ${CMAKE_SOURCE_DIR} /opt/intel/mkl/include)
    elseif(WIN32)
        target_compile_definitions(BigWham PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
        target_include_directories(BigWham PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)
    endif()
elseif(IL_OPENBLAS)
    target_compile_definitions(BigWham PRIVATE IL_OPENBLAS=1 IL_BLAS=1 IL_OPENMP=1)
    target_include_directories(BigWham PRIVATE ${CMAKE_SOURCE_DIR})
endif()

################################################################
# Pybind
################################################################
add_library(party SHARED
#        Hmat-lib/cluster/cluster.h
#        Hmat-lib/cluster/cluster.cpp
#        # 2D source code
#        src/core/Mesh2D.cpp src/core/Mesh2D.h
#        src/core/SegmentData.h
#        src/elasticity/2d/ElasticS3DP0_element.cpp
#        src/elasticity/2d/Elastic2DP1_element.cpp
#        src/elasticity/2d/ElasticHMatrix2DP0.h
#        src/elasticity/2d/ElasticHMatrix2DP1.h
#        src/elasticity/jacobi_prec_assembly.cpp
#        src/elasticity/PostProcessDDM.cpp
#        src/elasticity/2d/FullMatrixAssembly2D.cpp
#        src/elasticity/2d/FullMatrixAssembly2D.h
#        # 3D source code
#        Hmat-lib/linearAlgebra/factorization/lowRankApproximation.cpp
#        src/core/Mesh3D.cpp src/core/Mesh3D.h
#        src/core/FaceData.cpp
#        #  - T6 element
#        src/core/element_utilities.cpp
#        src/elasticity/3d/Elastic3DT6_element.cpp
#        src/elasticity/3d/h_potential.cpp
#        src/elasticity/3d/tensor_utilities.cpp
#        src/elasticity/3d/ElasticHMatrix3DT6.h
#        #  - R0 element
#        #        src/elasticity/3d/Elastic3DR0_element.h
#        src/elasticity/3d/Elastic3DR0_element.cpp
#        src/elasticity/3d/ElasticHMatrix3DR0.h
#
#        # Solver
#        src/solvers/HIterativeSolverUtilities.h
#        # BigWham Interface
#        src/BigWham.h
        src/test.cpp
        )
#target_link_libraries(party  ${CMAKE_MKL_LIBRARIES} ) #${CMAKE_OPENBLAS_LIBRARIES}
#if(IL_MKL)
#    if(APPLE)
#        target_compile_definitions(party PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
#        target_include_directories(party PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include $ENV{TBBROOT}/include)
#    elseif(UNIX)
#        target_compile_definitions(party PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
#        target_include_directories(party PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include $ENV{TBBROOT}/include)
#        # target_include_directories(HFPx2D PRIVATE ${CMAKE_SOURCE_DIR} /opt/intel/mkl/include)
#    elseif(WIN32)
#        target_compile_definitions(party PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
#        target_include_directories(party PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)
#    endif()
#elseif(IL_OPENBLAS)
#    target_compile_definitions(party PRIVATE IL_OPENBLAS=1 IL_BLAS=1 IL_OPENMP=1)
#    target_include_directories(party PRIVATE ${CMAKE_SOURCE_DIR})
#endif()

# use this
include_directories(${PROJECT_SOURCE_DIR}/pybind11-master/include/pybind11)
add_subdirectory(${PROJECT_SOURCE_DIR}/pybind11-master)
# or use
# find_package(pybind11 REQUIRED)

add_library(pyparty MODULE ${PROJECT_SOURCE_DIR}/src/pyparty.cpp)
target_link_libraries(pyparty PRIVATE pybind11::module party)

set_target_properties(pyparty PROPERTIES PREFIX "")
set_target_properties(pyparty PROPERTIES OUTPUT_NAME "pyparty")


add_custom_command(TARGET pyparty POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pyparty> /Users/carloperuzzo/PyFrac_git/src/pyparty.so
        COMMENT "Created  /Users/carloperuzzo/PyFrac_git/src/pyparty.so"
        )
add_custom_command(TARGET pyparty POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pyparty> /Users/carloperuzzo/bigwham/pybind_test/pyparty.so
        COMMENT "Created  /Users/carloperuzzo/bigwham/pybind_test/pyparty.so"
        )

################################################################################
# For Unit Tests
################################################################################
set(IL_FOLDER "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_FOLDER='\"${IL_FOLDER}\"'")
set(TESTED_SOURCE_FILES
        src/core/FaceData.cpp
        src/elasticity/3d/Elastic3DR0_element.cpp
        src/elasticity/3d/Elastic3DR0_element.h

        )
set(UNIT_TEST_FILES
        #src/elasticity/_test/elasticity_tests.cpp
        src/elasticity/_test/elastic3DR0_element_tests.cpp

        il/gtest/src/gtest-all.cc
        )

add_executable(BigWhamElastUnitTest ${TESTED_SOURCE_FILES} ${UNIT_TEST_FILES} test.cpp )

target_include_directories(BigWhamElastUnitTest PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/il/gtest)
target_link_libraries(BigWhamElastUnitTest ${CMAKE_MKL_LIBRARIES} pthread)

# For unit tests: The precondition of our fonctions are checked with assert
# macros that terminate the program in debug mode. In order to test those macros
# in our unit tests, the flag IL_UNIT_TEST turns them into exceptions that can
# be caught and checked.
target_compile_options(BigWhamElastUnitTest PRIVATE "-DIL_UNIT_TEST")