cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(extraction SHARED
        extraction.cpp
)

target_include_directories(extraction PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ../utility  # 这里假设 utility 库的头文件在同级的 utility 目录下
)

# 链接 utility 库
target_link_libraries(extraction PRIVATE utility)