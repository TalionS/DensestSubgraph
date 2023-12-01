cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(verification SHARED
        verification.cpp
)

target_include_directories(verification PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
