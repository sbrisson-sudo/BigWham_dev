## All subdirectories editing SOURCE_FILES should named before
set(SOURCE_FILES ${SOURCE_FILES}
  ## 3D Tri 0 Element
  ${CMAKE_CURRENT_SOURCE_DIR}/
  ${CMAKE_CURRENT_SOURCE_DIR}/3d/elastic3DT0_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/3d/bie_elastostatic_triangle_0_impls.cc
  PARENT_SCOPE)
