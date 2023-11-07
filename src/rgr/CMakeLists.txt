set(target_name matrix)
add_library(
        ${target_name} INTERFACE
  #      rgr/matrix/Matrix.hpp
)
include(CompileOptions)
target_include_directories(
        ${target_name}
        INTERFACE
        ${CMAKE_CURRENT_LIST_DIR}
)
