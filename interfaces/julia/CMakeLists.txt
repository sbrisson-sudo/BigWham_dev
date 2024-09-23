set(JULIA_EXECUTABLE "julia")

set(JULIA_SCRIPT "import Pkg; Pkg.add(\"CxxWrap\", io=devnull);using CxxWrap;println(CxxWrap.prefix_path() * \"/lib/cmake/JlCxx\")")

execute_process(
    COMMAND ${JULIA_EXECUTABLE} -e "${JULIA_SCRIPT}"
    OUTPUT_VARIABLE JlCxx_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)


find_package(JlCxx REQUIRED)


get_target_property(JlCxx_INCLUDE_DIR JlCxx::cxxwrap_julia INTERFACE_INCLUDE_DIRECTORIES)

get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
message(STATUS "Found JlCxx Cmake Dir at ${JlCxx_DIR}")
message(STATUS "Found JlCxx at ${JlCxx_location}")

# message("Julia Library ${Julia_LIBRARY}")
# message("Julia Include ${Julia_INCLUDE_DIRS}")
# message("JlCxx CxxWrap Include ${JlCxx_INCLUDE_DIR}")

add_library(jl_bigwham SHARED jl_bigwham.cc)
target_link_libraries(jl_bigwham PRIVATE JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl BigWhamShared)

set(SRC_CXXJULIA "${CMAKE_CURRENT_SOURCE_DIR}/bigwham4jl.jl")

message(STATUS "Use the following command in terminal to use BigWham julia interface:
include(\"${CMAKE_BINARY_DIR}/interfaces/julia/bigwham4jl.jl\")
import .BigWham:BEMatrix
")
add_custom_target(copy_juliacxx_files ALL COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_CXXJULIA}
                                           ${CMAKE_BINARY_DIR}/interfaces/julia)
