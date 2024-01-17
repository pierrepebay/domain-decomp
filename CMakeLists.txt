cmake_minimum_required(VERSION 3.27)
project(domain_decomp)

set(CMAKE_CXX_STANDARD 20)

add_library(jacobi_lib Jacobi.cpp Jacobi.h)
add_library(domain_lib Domain.cpp Domain.h)

add_executable(domain_decomp main.cpp)

target_link_libraries(domain_decomp jacobi_lib domain_lib)
