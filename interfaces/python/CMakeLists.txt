if(IL_MKL)
    if(APPLE)
        target_link_directories(part PRIVATE ${OPENMP_LIBRARY})
        target_link_libraries(part ${CMAKE_MKL_LIBRARIES} -lpthread -lm -ldl)
        target_compile_definitions(part PRIVATE IL_MKL=1 IL_BLAS=0 IL_OPENMP=0)
        target_include_directories(part PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include $ENV{TBBROOT}/include)
    elseif(UNIX)
        target_link_directories(part PRIVATE ${OPENMP_LIBRARY})
        target_link_libraries(part "-Wl,--start-group" ${CMAKE_MKL_LIBRARIES} "-Wl,--end-group" -lpthread -lm -ldl)
        target_compile_definitions(part PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
        target_include_directories(part PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include $ENV{TBBROOT}/include)
    elseif(WIN32)
        target_compile_definitions(part PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
        target_include_directories(part PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)
    endif()
elseif(IL_OPENBLAS)
    target_link_libraries(part  ${CMAKE_OPENBLAS_LIBRARIES} -lgfortran -lpthread -lm    -ltbb   )
    target_compile_definitions(part PRIVATE IL_OPENBLAS=1 IL_BLAS=1 IL_OPENMP=1)
    target_include_directories(part PRIVATE ${CMAKE_SOURCE_DIR}  $ENV{TBBROOT}/include/tbb ${OPENBLASROOT})
endif()

pybind11_add_module(bigwhamPybind bigwhamPybind.cpp)
target_link_libraries(bigwhamPybind PRIVATE part)
set(SRC_PYBIND11 "${CMAKE_CURRENT_SOURCE_DIR}/bigwham4py.py")
add_custom_target(copy_pybind11_files ALL COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_PYBIND11}
                                           ${CMAKE_BINARY_DIR}/interfaces/python)

FILE(WRITE ${pyenv_root_lib_path} "${CMAKE_BINARY_DIR}/interfaces/python")
# target_compile_definitions(cmake_example
#                        PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
# add_subdirectory(${PROJECT_SOURCE_DIR}/interfaces/python/pybind11-master)
# include_directories(${PROJECT_SOURCE_DIR}/interfaces/python/pybind11-master/include/  )
# or use
# find_package(pybind11 REQUIRED)
# add_library(bigwhamPybind MODULE interfaces/python/src/bigwhamPybind.cpp)
# target_link_libraries(bigwhamPybind PRIVATE pybind11::module part)

# this command is needed to copy the library bigwhamPybind to a specific location
# add_custom_command(TARGET bigwhamPybind POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:bigwhamPybind> ${PROJECT_SOURCE_DIR}/interfaces/python/lib/bigwhamPybind.so
#         COMMENT "Created ${PROJECT_SOURCE_DIR}/interfaces/python/lib/bigwhamPybind.so")

# if(COPY_TO_SPECIFIC_LOCATION)
#       # this command is needed to copy the library bigwhamPybind to a specific location
#       add_custom_command(TARGET bigwhamPybind POST_BUILD
#            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:bigwhamPybind>  ${BIGWHAM_PYBIND_DEST}/bigwhamPybind.so
#            COMMENT "Created ${BIGWHAM_PYBIND_DEST}/bigwhamPybind.so ")
# endif()


