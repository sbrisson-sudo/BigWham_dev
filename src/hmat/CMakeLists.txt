# add_subdirectory(elements)
set(SOURCE_FILES ${SOURCE_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/cluster/cluster.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/hmatrix/toHPattern.cpp
    PARENT_SCOPE
)
