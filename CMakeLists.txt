cmake_minimum_required(VERSION 3.10)
project(fluidsim)
set(CMAKE_CXX_STANDARD 20)
add_executable(fluidsim
				src/main.cpp
				src/utils.hpp
				src/fluid.hpp
				src/solver.hpp
)

find_package(OpenMP REQUIRED)


set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

target_include_directories(fluidsim PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(fluidsim PRIVATE OpenMP::OpenMP_CXX)
