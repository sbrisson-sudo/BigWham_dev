## All subdirectories editing SOURCE_FILES should named before
set(SOURCE_FILES ${SOURCE_FILES}
  
  ## 3D Tri 0 Element
  ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_3d_triangle/elastic_3dT0_element.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_3d_triangle/elastic_3dT0_element_af.cc
        ## 3D Tri 0 Element
  ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_3d_triangle/elastic_3dT6_element.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_3d_triangle/element_utilities_3DT6.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_3d_triangle/h_potential_3DT6.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_3d_triangle/tensor_utilities_3DT6.cc

        #   ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_3d_triangle/bie_elastostatic_triangle_0_influence.h

  ## 3D Rect 0 Element
        ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_3d_rectangle/elastic_3dR0_common.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_3d_rectangle/elastic_3dR0_element.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_3d_rectangle/elastic_3dR0_mode1Cartesian_element.cc

  ## 2D Segment 0
     #   ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_2d_segment/bie_elastostatic_segment_0_influence.h
  ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_2d_segment/elastic_2dP0_segment.cc

  ## 2D Segment 1
      #  ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_2d_segment/bie_elastostatic_segment_1_influence.h
  ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_2d_segment/elastic_2dP1_segment.cc

  ## Simplified 3D  Segment 0
   ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_sp3d_segment/elastic_sp3dP0_element.cc
  
  ## AxisSymm Ring 
  ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_axisymmetry_flat_unidirectional/elliptic_integral.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fullspace_iso_axisymmetry_flat_unidirectional/elastic_axi3dP0_element.cc


        PARENT_SCOPE)
