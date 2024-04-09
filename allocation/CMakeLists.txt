cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(allocation SHARED
        allocation.cpp
)

target_include_directories(allocation PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ../utility
)

target_link_libraries(allocation PRIVATE utility)