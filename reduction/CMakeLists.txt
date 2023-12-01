cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(reduction SHARED
        reduction.cpp
)

target_include_directories(reduction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
