cmake_minimum_required (VERSION 2.6)
project (AdvancedNumericalMethods)

file(GLOB_RECURSE sources cpp_src/*.cpp cpp_src/*.h)
file(GLOB_RECURSE test_sources cpp_tests/*.cpp)

add_executable(Matrix ${sources})
add_executable(Matrix_Tests ${test_sources})

target_compile_options(Matrix PUBLIC -fopenmp -lgomp -std=c++11 -Wall)
target_compile_options(Matrix_Tests PUBLIC -fopenmp -std=c++11 -Wall)

target_link_libraries(Matrix  PUBLIC -fopenmp)
target_link_libraries(Matrix_Tests PUBLIC -fopenmp)
