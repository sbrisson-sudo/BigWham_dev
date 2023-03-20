add_subdirectory(core)
add_subdirectory(hmat)
add_subdirectory(elasticity)
## All subdirectories editing SOURCE_FILES should named before
set(SOURCE_FILES ${SOURCE_FILES}
  ${CMAKE_CURRENT_SOURCE_DIR}/bigwham_io.h
  PARENT_SCOPE)
