find_package(JlCxx)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

message(STATUS "Found JlCxx at ${JlCxx_location}")

include_directories(${J1Cxx_INCLUDE_DIRS})

add_library(jl_bigwham SHARED jl_bigwham.cc)
target_link_libraries(jl_bigwham J1Cxx::cxxwrap_julia)