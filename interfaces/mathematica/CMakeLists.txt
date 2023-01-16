## Root files
set(src_math_root_files "PacletInfo.m"
                        "BigWhamLink.nb"
                        "BigWhamLink-RunTests.nb")
add_custom_target(copy_math_root ALL COMMAND ${CMAKE_COMMAND} -E 
  copy_if_different ${src_math_root_files}
  ${CMAKE_BINARY_DIR}/interfaces/mathematica)

## library resource files
set(src_math_libresources_files "BigWhamLink/LibraryResources/Source/HMatExpr.h"
                                "BigWhamLink/LibraryResources/Source/Manager.h")
add_custom_target(copy_math_libresources ALL COMMAND ${CMAKE_COMMAND} -E 
  copy_if_different ${src_math_libresources_files}
  ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/LibraryResources/Source)

## copy ltemplate 
add_custom_target(copy_math_ltemplate ALL COMMAND ${CMAKE_COMMAND} -E 
  copy_if_different BigWhamLink/Kernel/LTemplate
  ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/Kernel)

## copy tests
add_custom_target(copy_math_tests ALL COMMAND ${CMAKE_COMMAND} -E 
  copy_if_different BigWhamLink/Tests
  ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/Tests)


## copy Examples
add_custom_target(copy_math_examples ALL COMMAND ${CMAKE_COMMAND} -E 
  copy_if_different BigWhamLink/Examples
  ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/Examples)


# 
configure_file(BigWhamLink/Kernel/BuildSettings.m
  ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/Kernel/BuildSettings.m)


configure_file(BigWhamLink/Kernel/BuildSettings.m
  ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/Kernel/BuildSettings.m)
