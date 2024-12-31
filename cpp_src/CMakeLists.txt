cmake_minimum_required(VERSION 3.5)
project(mwmrna LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
# If GCC, add -fconcepts-diagnostics-dep  th=2
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts-diagnostics-depth=2")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_subdirectory(mwmrna_lib)
add_subdirectory(exe)

