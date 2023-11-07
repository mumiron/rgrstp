set(target_name matrix_test)
add_executable(
        ${target_name}
        rgr/rgr/matrix.test.cpp
)

include(CompileOptions)
set_compile_options(${target_name})

target_link_libraries(
        ${target_name}
        PRIVATE
        matrix
        gtest_main
)

include(GoogleTest)
gtest_discover_tests(${target_name})
