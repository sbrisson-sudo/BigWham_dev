## All subdirectories editing SOURCE_FILES should named before
set(SOURCE_FILES ${SOURCE_FILES}
  ${CMAKE_CURRENT_SOURCE_DIR}/bigwham_io_gen.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bigwham_io_helper.cc
  PARENT_SCOPE)
