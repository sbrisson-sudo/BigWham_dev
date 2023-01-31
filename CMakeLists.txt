#####################################
# BIGWHAM                           #
# compilation configuration         #
#####################################

cmake_minimum_required(VERSION 3.7...3.19)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(BigWham_exe VERSION 1.0
        DESCRIPTION "This is the executable of BigWham"
        LANGUAGES CXX)

project(bigwhamPybind VERSION 1.0
        DESCRIPTION "This will create a shared library to be imported and used in Python"
        LANGUAGES CXX)

project(BigWhamElastUnitTest VERSION 1.0
        DESCRIPTION "This is for testing the wheels in BigWham"
        LANGUAGES CXX)

project(BigWham VERSION 1.0
        DESCRIPTION "This is the library containing the wheels"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

#####################################
# set options: False = 0, True = 1  #
#####################################
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
option(BUILD_PYTHON_BINDINGS "Build Pybind11 binding" ON)
option(BUILD_STATIC_LIB "Build STATIC .a library" ON)
option(BUILD_GOOGLE_TESTS "Build google tests" ON)

option(LLDB_ENABLE_PYTHON "LLDB python debugging" ON)
option(CMAKE_VERBOSE_MAKEFILE "makefile verbose" ON)
option(IL_MKL "MKL" ON)      # --> if IL_MKL=0 then it must be IL_OPENBLAS=1
option(IL_OPENBLAS "OPENBLAS" OFF)
option(IL_OPENMP "OPENMP" ON) # --> to activate openMP for bigwham (recommended but should have openmp lib installed)
option(IL_TBB "TBB" ON)      # --> this must be 1 at the moment

option(COPY_TO_SPECIFIC_LOCATION "bigwham copy to specific location" OFF) # --> update BIGWHAM_PYBIND_DEST if you want to copy the shared library to a specific destination

option(USE_INTEL "use intel" OFF)   # --> compiler option
option(USE_STATIC "use static" OFF)  # --> to add static option to CMAKE_CXX_FLAGS

set(MATHEMATICA_INTERFACE "make mathematica interface" ON)  # --> to add static option to CMAKE_CXX_FLAGS

#####################################
# optimization options              #
#####################################

set(CMAKE_CXX_FLAGS_DEBUG     "${CMAKE_CXX_FLAGS_DEBUG}   -std=c++17")
set(CMAKE_CXX_FLAGS_RELEASE   "${CMAKE_CXX_FLAGS_RELEASE}  -std=c++17 -O3  -mavx2")
#set(CMAKE_CXX_FLAGS_RELEASE   "${CMAKE_CXX_FLAGS_RELEASE}  -std=c++11 -Ofast -march=native -mtune=intel  -mavx2")

if(USE_STATIC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
endif()

################################################################################
# set compiler: clang, gcc, intel etc.
################################################################################

if(APPLE)
    #    MESSAGE("-- located at /usr/bin/")
    #    set(CMAKE_C_COMPILER /usr/bin/clang)
    #    set(CMAKE_CXX_COMPILER /usr/bin/clang++)
    MESSAGE("   -- located at /usr/local/opt/llvm/bin/clang")
    set(CMAKE_C_COMPILER /usr/local/opt/llvm/bin/clang)
    set(CMAKE_CXX_COMPILER /usr/local/opt/llvm/bin/clang++)
    #    MESSAGE("-- located at /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/")
    #    set(CMAKE_CXX_COMPILER /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c)
    #    set(CMAKE_CXX_COMPILER /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++)
elseif(UNIX)
    if(USE_INTEL)
        #set(CMAKE_C_COMPILER /opt/intel/compilers_and_libraries/linux/bin/intel64/icc)
        #set(CMAKE_CXX_COMPILER /opt/intel/compilers_and_libraries/linux/bin/intel64/icpc)
        #set(CMAKE_C_COMPILER /opt/intel/bin/icc)
        #set(CMAKE_CXX_COMPILER /opt/intel/bin/icpc)
    else()
         set(CMAKE_C_COMPILER /usr/bin/gcc)
         set(CMAKE_CXX_COMPILER /usr/bin/g++)
         #set(CMAKE_CXX_COMPILER /cm/local/apps/gcc/10.2.0/bin/c++)
         #set(CMAKE_C_COMPILER /cm/local/apps/gcc/10.2.0/bin/gcc)
    endif()
elseif(WIN32)
    if(MINGW)
        set(CMAKE_C_COMPILER C:/msys64/mingw64/bin/gcc)
        set(CMAKE_CXX_COMPILER C:/msys64/mingw64/bin/g++)
    endif()
endif()

MESSAGE("   Detecting the compiler:")
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    MESSAGE("   -- The compiler is:  ${CMAKE_CXX_COMPILER}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    MESSAGE("   -- The compiler is: ${CMAKE_CXX_COMPILER}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    MESSAGE("   -- The compiler is: ${CMAKE_CXX_COMPILER}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    MESSAGE("   -- The compiler is: ${CMAKE_CXX_COMPILER}")
endif()

################################################################################
# parallel framework
################################################################################

# For OpenMP
if(IL_OPENMP)
    if(APPLE)
        if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            MESSAGE("   -- OPENMP ENABLED")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp") #
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNUMBEROFTHREADS=8")
            set(OPENMP_LIBRARY /usr/local/opt/libomp/lib/)#libomp.dylib
            set(OPENMP_LIBRARIES /usr/local/opt/libomp/lib/libomp.dylib)
            link_directories(/usr/local/opt/libomp/lib/) #
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            MESSAGE("   -- OPENMP ENABLED")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        endif()
    elseif(UNIX)
        if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            MESSAGE("   -- OPENMP ENABLED")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
            link_directories(/opt/clang-3.9.0/lib)
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            MESSAGE("   -- OPENMP ENABLED")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fopenmp")
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            MESSAGE("   -- OPENMP ENABLED")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
        endif()
    elseif(WIN32)
        if(MINGW)
            MESSAGE("   -- OPENMP ENABLED")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        endif()
    endif()
endif()

# For TBB
if(IL_TBB)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_TBB")
endif()

################################################################################
# math framework
################################################################################

# For Intel MKL
if(IL_MKL)
    MESSAGE("   Detecting MKL ${IL_MKL} ${IL_OPENBLAS}")
    MESSAGE(${CMAKE_CXX_COMPILER_ID})
    set(MKLROOT $ENV{MKLROOT})
    if(APPLE)
        if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            MESSAGE("   -- The MKL is found at: " $ENV{MKLROOT}/lib/)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DIL_MKL") # -DIL_BLAS
            set(CMAKE_MKL_LIBRARIES $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_core.a $ENV{MKLROOT}/lib/libmkl_sequential.a  tbb dl  pthread  m)
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            MESSAGE("   -- The MKL is found at: " $ENV{MKLROOT}/lib/)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -m64 -lmkl_sequential -DIL_MKL ") #-m64 -DIL_BLAS
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
            set(CMAKE_MKL_LIBRARIES $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a tbb dl pthread m)
            link_directories($ENV{MKLROOT}/lib/intel64 $ENV{MKLROOT}/../compilers_and_libraries/linux/lib/intel64 $ENV{TBBROOT}/lib/intel64/gcc4.7)
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=sequential -DIL_MKL -DIL_BLAS")
            set(CMAKE_MKL_LIBRARIES $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_core.a $ENV{MKLROOT}/lib/libmkl_sequential.a  tbb dl  pthread  m) #mkl_intel_thread mkl_sequential  iomp5
        endif()
    elseif(UNIX)
        if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64  -DIL_MKL ") #-DIL_BLAS
            set(CMAKE_MKL_LIBRARIES mkl_intel_lp64 mkl_core mkl_sequential  tbb dl pthread m)
            link_directories($ENV{MKLROOT}/lib $ENV{MKLROOT}/../compiler/lib )
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -m64 -lmkl_sequential -DIL_MKL ") #-m64 -DIL_BLAS
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
            set(CMAKE_MKL_LIBRARIES $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a tbb dl pthread m)
            link_directories($ENV{MKLROOT}/lib/intel64 $ENV{TBBROOT}/lib/intel64/gcc4.8)
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=sequential -DIL_MKL -DIL_BLAS")# -mkl=sequential
            set(CMAKE_MKL_LIBRARIES $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a  tbb dl  pthread  m) #mkl_intel_thread mkl_sequential  iomp5
            link_directories($ENV{MKLROOT}/lib/intel64 $ENV{MKLROOT}/../compilers_and_libraries/linux/lib/intel64)
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
    MESSAGE("   -- The MKL is found at: " $ENV{MKLROOT})

    # For OpenBlas
elseif(IL_OPENBLAS)
    MESSAGE("   Detecting OPENBLAS ${IL_OPENBLAS}")
    if(APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_MKL -DIL_OPENBLAS -DIL_BLAS")
        set(CMAKE_OPENBLAS_LIBRARIES "openblas")
    elseif(UNIX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_OPENBLAS -DIL_BLAS")
        # set(CMAKE_OPENBLAS_LIBRARIES "openblas")
        set(CMAKE_OPENBLAS_LIBRARIES $ENV{OPENBLASROOT}/lib/libopenblas.a) 
        set(OPENBLASROOT $ENV{OPENBLASROOT}/include )
        message(${OPENBLASROOT})
    elseif(WIN32)
        if(MINGW)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_OPENBLAS -DIL_BLAS")
            set(CMAKE_OPENBLAS_LIBRARIES "openblas")
        endif()
    endif()
endif()

################################################################################
# file definitions
################################################################################

set(SOURCE_FILES_DIR ${PROJECT_SOURCE_DIR}/src)

if(IL_MKL)
    include_directories(il $ENV{TBBROOT}/include/tbb $ENV{MKLROOT}/include  src src/hmat .)
elseif(IL_OPENBLAS)
    include_directories(il $ENV{TBBROOT}/include/tbb ${OPENBLASROOT}  src src/hmat .)
endif()


set(SOURCE_FILES
        src/hmat/cluster/cluster.h
        src/hmat/cluster/cluster.cpp
        # 2D source code
        src/core/Mesh2D.cpp src/core/Mesh2D.h
        src/core/SegmentData.h
        src/core/SquareMatrixGenerator.h
        src/core/BEMesh.h
        src/elasticity/2d/ElasticS3DP0_element.cpp
        src/elasticity/2d/Elastic2DP1_element.cpp
        src/elasticity/PostProcessDDM_2d.cpp
        src/elasticity/PostProcessDDM_3d.cpp
        src/elasticity/2d/FullMatrixAssembly2D.cpp
        src/elasticity/2d/FullMatrixAssembly2D.h

        # 3D source code
        src/core/Mesh3D.cpp src/core/Mesh3D.h
        src/core/FaceData.cpp
        #  - T6 element
        src/core/element_utilities.cpp
        src/elasticity/3d/Elastic3DT6_element.cpp
        src/elasticity/3d/h_potential_3DT6.cpp
        src/elasticity/3d/tensor_utilities_3DT6.cpp
        src/elasticity/3d/ElasticHMatrix3DT6.h
        #  - T0 element
        src/elasticity/3d/Elastic3DT0_element.cpp
        src/elasticity/3d/ElasticHMatrix3DT0.h
        src/elasticity/3d/ElasticHMatrix3DT0_modes2and3.h
        #  - R0 element
        src/elasticity/3d/Elastic3DR0_common.cpp
        src/elasticity/3d/Elastic3DR0_element.cpp
        src/elasticity/3d/Elastic3DR0_mode1Cartesian_element.cpp
        src/elasticity/3d/ElasticHMatrix3DR0.h
        src/elasticity/3d/ElasticHMatrix3DR0displ.h
        src/elasticity/3d/ElasticHMatrix3DR0_mode1Cartesian.h
        src/elasticity/3d/ElasticHMatrix3DR0_modes2and3Cartesian.h
        # BigWham Interface
        src/BigWham.h
        #new
        src/core/BIE_Kernel.h
        src/core/BEMesh.h
        src/core/BoundaryElement.h
        src/elasticity/BIE_elastostatic.h

        src/core/hierarchical_representation.h
        src/hmat/hmatrix/toHPattern.cpp
        src/core/elements/Segment.cpp
        src/core/elements/Segment.h
        src/core/elements/Triangle.cpp
        src/core/elements/Triangle.h
        src/core/elements/Rectangle.cpp
        src/core/elements/Rectangle.h
        src/elasticity/2d/BIE_elastostatic_segment_0_impls.h
        src/elasticity/2d/BIE_elastostatic_segment_1_impls.h)


################################################################################
# buid static library 
################################################################################

if(BUILD_STATIC_LIB)
    set(CMAKE_POSITION_INDEPENDENT_CODE True)
    add_library(BigWham STATIC ${SOURCE_FILES} )
    
    if(IL_MKL)
        target_link_libraries(BigWham  ${CMAKE_MKL_LIBRARIES} )
    elseif(IL_OPENBLAS)
        target_link_libraries(BigWham  ${CMAKE_OPENBLAS_LIBRARIES} -lgfortran -lpthread -lm -ltbb)
    endif()
    
    set_target_properties(BigWham PROPERTIES LINKER_LANGUAGE CXX)

    add_executable(BigWham_exe  main.cpp )
   
    if(IL_MKL)
       include_directories(il $ENV{TBBROOT}/include/ $ENV{MKLROOT}/include/)
    elseif(IL_OPENBLAS)
       include_directories(il $ENV{TBBROOT}/include/tbb  ${OPENBLASROOT} )
    endif()
   
    target_link_libraries(BigWham_exe BigWham)

    if(IL_MKL)
        if(APPLE)
            target_link_directories(BigWham_exe PRIVATE ${OPENMP_LIBRARY})
            target_link_libraries(BigWham_exe  ${CMAKE_MKL_LIBRARIES} -lpthread -lm -ldl)
            target_compile_definitions(BigWham PRIVATE IL_MKL=1 IL_BLAS=0 IL_OPENMP=0)
            target_include_directories(BigWham PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include $ENV{TBBROOT}/include $ENV{TBBROOT}/lib )
        elseif(UNIX)
            target_link_libraries(BigWham_exe "-Wl,--start-group" ${CMAKE_MKL_LIBRARIES} "-Wl,--end-group" -lpthread -lm -ldl)
            target_compile_definitions(BigWham PRIVATE IL_MKL=1 IL_BLAS=0 IL_OPENMP=0)
            target_include_directories(BigWham PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include $ENV{TBBROOT}/include)
        elseif(WIN32) # not tested
            target_compile_definitions(BigWham PRIVATE IL_MKL=1 IL_BLAS=0 IL_OPENMP=0)
            target_include_directories(BigWham PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)
        endif()
    elseif(IL_OPENBLAS)
        #target_link_libraries(BigWham_exe BigWham  ${CMAKE_OPENBLAS_LIBRARIES} -lgfortran -lpthread -lm  -ltbb   )
        target_compile_definitions(BigWham PRIVATE IL_OPENBLAS=1 IL_BLAS=1 IL_OPENMP=0)
        target_include_directories(BigWham PRIVATE ${CMAKE_SOURCE_DIR}  ${OPENBLASROOT}   $ENV{TBBROOT}/include/tbb )
       
    endif()
endif()


################################################################
# Pybind --->use with PYTHON version >=3.7
################################################################

message("BUILD_PYTHON_BINDINGS ${BUILD_PYTHON_BINDINGS}")

if(BUILD_PYTHON_BINDINGS)
    ## Set root of python virtual env
    set(PYENV_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/bigwham_venv" CACHE PATH "root path of bigwham venv")
    find_package(Python COMPONENTS Interpreter Development)
    set(PYBIND11_PYTHON_VERSION 3.5)
    # message("The Python interpreter set by the user is : " ${Python_EXECUTABLE})
    # message("The Python Library set by the user is : " ${Python_LIBRARIES})
    set(pyenv_root_lib_path "${PYENV_ROOT}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/.pth")
    message("BigWham Python Library Path: ${pyenv_root_lib_path}")
    add_subdirectory(pybind11)

    # if (APPLE)
    #     set(PYTHON_EXECUTABLE "/Users/bricelecampion/opt/anaconda3/envs/python3.9/bin/python3.9")
    #     set(PYTHON_LIBRARY "/Users/bricelecampion/opt/anaconda3/envs/python3.9/lib/libpython3.9.dylib ")
    # elseif(LINUX)
    #     #set(PYTHON_EXECUTABLE "/home/peruzzo/anaconda3/envs/PyFrac/bin/python3.8")
    #     #set(PYTHON_LIBRARY "/home/peruzzo/anaconda3/envs/PyFrac/lib/libpython3.8.so")
    #     #set(PYTHON_LIBRARY "/home/peruzzo/anaconda3/envs/PyFrac/lib/libpython3.8.so")
    # elseif(WINDOWS)
    #     set(PYTHON_EXECUTABLE "/home/peruzzo/anaconda3/envs/PyFrac/bin/python3.8")
    #     set(PYTHON_LIBRARY "/home/peruzzo/anaconda3/envs/PyFrac/bin/python3.8/libpython3.9.dylib ")
    # endif()

    add_library(part SHARED  ${SOURCE_FILES})
    add_subdirectory(interfaces/python)

endif()

################################################################################
# Mathematica interface
################################################################################
message("THis is the mma interface flag", ${MATHEMATICA_INTERFACE})
if(MATHEMATICA_INTERFACE)
    add_subdirectory(interfaces/mathematica)
endif()

################################################################################
# For Unit Tests
################################################################################
if(BUILD_GOOGLE_TESTS)

   # set(IL_FOLDER "${CMAKE_SOURCE_DIR}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_FOLDER='\"${IL_FOLDER}\"'")
    set(TESTED_SOURCE_FILES ${SOURCE_FILES})
    set(UNIT_TEST_FILES
            il/gtest/src/gtest-all.cc
            # src/elasticity/_test/elasticity_tests.cpp
            src/_test/penny_shaped_crack_analytical_sol.h
            src/_test/penny_shaped_crack_analytical_sol.cpp
            src/core/_test/FaceData_tests.cpp
            src/elasticity/_test/elastic3DR0_element_tests.cpp
            #src/_test/elastic3DR0_element_benchmark.cpp
            src/_test/penny_shaped_crack_analytical_sol_benchmark.cpp
            src/core/_test/segment_test.cpp src/elasticity/BIE_elastostatic.h src/elasticity/_test/segment0_sp3d_tests.cpp src/core/_test/BEmesh_tests.cpp src/elasticity/_test/elastic_square_matrix_generator_tests.cpp src/core/_test/hierarchical_representation_tests.cpp src/hmat/hmatrix/_test/Hmat_test_gaussian_mat.cpp src/core/_test/triangle_tests.cpp src/core/_test/rectangle_tests.cpp src/_test/bigwhamio_tests_2d.cpp)

    add_executable(BigWhamElastUnitTest ${TESTED_SOURCE_FILES} ${UNIT_TEST_FILES} test.cpp )

    #target_include_directories(BigWhamElastUnitTest PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/il/gtest ${CMAKE_SOURCE_DIR}/il ${CMAKE_SOURCE_DIR}/Hmat-lib $ENV{MKLROOT}/include )
    target_include_directories(BigWhamElastUnitTest PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/il/gtest ${CMAKE_SOURCE_DIR}/il ${CMAKE_SOURCE_DIR}/Hmat-lib ${OPENBLASROOT})
    if(IL_MKL)
        if(APPLE)
            target_link_directories(BigWhamElastUnitTest PRIVATE ${OPENMP_LIBRARY})
            target_link_libraries(BigWhamElastUnitTest   ${CMAKE_MKL_LIBRARIES}  -lpthread -lm -ldl)
        elseif(UNIX)
            target_link_directories(BigWhamElastUnitTest PRIVATE ${OPENMP_LIBRARY})
            target_link_libraries(BigWhamElastUnitTest "-Wl,--start-group" ${CMAKE_MKL_LIBRARIES} "-Wl,--end-group" -lpthread -lm -ldl)

        elseif(WIN32) # not tested
            target_compile_definitions(BigWhamElastUnitTest PRIVATE IL_MKL=1 IL_BLAS=0 IL_OPENMP=0)
            target_include_directories(BigWhamElastUnitTest PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)
        endif()
    elseif(IL_OPENBLAS)
        #target_compile_definitions(BigWhamElastUnitTest PRIVATE IL_OPENBLAS=1 IL_BLAS=1 IL_OPENMP=0)
        target_link_directories(BigWhamElastUnitTest PRIVATE ${OPENMP_LIBRARY})
        target_link_libraries(BigWhamElastUnitTest PRIVATE  ${CMAKE_OPENBLAS_LIBRARIES} -lgfortran -lpthread -lm    -ltbb) 
        target_include_directories(BigWhamElastUnitTest PRIVATE ${CMAKE_SOURCE_DIR}  $ENV{TBBROOT}/include/tbb ${OPENBLASROOT}  )
    endif()

    # For unit tests: The precondition of our fonctions are checked with assert
    # macros that terminate the program in debug mode. In order to test those macros
    # in our unit tests, the flag IL_UNIT_TEST turns them into exceptions that can
    # be caught and checked.
    target_compile_options(BigWhamElastUnitTest PRIVATE "-DIL_UNIT_TEST")
endif()


################################################################################
# For running examples
################################################################################

# add_library(bigwham_shared SHARED  ${SOURCE_FILES})
# add_subdirectory(src/examples)
