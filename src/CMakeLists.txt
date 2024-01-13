add_subdirectory(core)
add_subdirectory(elements)
add_subdirectory(elasticity)
add_subdirectory(hmat)
add_subdirectory(io)
## All subdirectories editing SOURCE_FILES should named before
set(SOURCE_FILES ${SOURCE_FILES}
  PARENT_SCOPE)
