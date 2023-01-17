## Make directories
add_custom_target(make_math_root ALL COMMAND ${CMAKE_COMMAND} -E 
  make_directory
  ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink)
add_custom_target(make_math_libresources ALL COMMAND ${CMAKE_COMMAND} -E 
  make_directory
  ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/LibraryResources/Source)
add_custom_target(make_math_kernel ALL COMMAND ${CMAKE_COMMAND} -E 
  make_directory
  ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/Kernel)

## Root files
set(src_math_root_files "${CMAKE_CURRENT_SOURCE_DIR}/PacletInfo.m")
                        # "${CMAKE_CURRENT_SOURCE_DIR}/BigWhamLink.nb"
                        # "${CMAKE_CURRENT_SOURCE_DIR}/BigWhamLink-RunTests.nb")
add_custom_target(copy_math_root ALL COMMAND ${CMAKE_COMMAND} -E 
  copy_if_different ${src_math_root_files}
  ${CMAKE_BINARY_DIR}/interfaces/mathematica)

add_dependencies(copy_math_root BigWham)

## library resource files
set(src_math_libresources_files "${CMAKE_CURRENT_SOURCE_DIR}/BigWhamLink/LibraryResources/Source/HMatExpr.h"
                                "${CMAKE_CURRENT_SOURCE_DIR}/BigWhamLink/LibraryResources/Source/Manager.h")
add_custom_target(copy_math_libresources ALL COMMAND ${CMAKE_COMMAND} -E 
  copy_if_different ${src_math_libresources_files}
  ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/LibraryResources/Source)

# ## copy ltemplate 
# add_custom_target(copy_math_ltemplate ALL COMMAND ${CMAKE_COMMAND} -E 
#   copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/BigWhamLink/Kernel/LTemplate
#   ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/Kernel/LTemplate)

# ## copy tests
# add_custom_target(copy_math_tests ALL COMMAND ${CMAKE_COMMAND} -E 
#   copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/BigWhamLink/Tests
#   ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/Tests)


# ## copy Examples
# add_custom_target(copy_math_examples ALL COMMAND ${CMAKE_COMMAND} -E 
#   copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/BigWhamLink/Examples
#   ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/Examples)

# #
# set(BLA_STATIC ON)
# message(${BLA_STATIC})
# find_package(BLAS REQUIRED)

# 
set(LTemplateLibPath "${CMAKE_SOURCE_DIR}/LTemplate/LTemplate")
configure_file(BigWhamLink/Kernel/BuildSettings.m
  ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/Kernel/BuildSettings.m)


configure_file(BigWhamLink/Kernel/BigWhamLink.m
  ${CMAKE_BINARY_DIR}/interfaces/mathematica/BigWhamLink/Kernel/BigWhamLink.m)
