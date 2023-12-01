cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(DensestSubgraph)

include_directories(
        ${PROJECT_SOURCE_DIR}/
)

add_executable(DensestSubgraph main.cpp)
target_link_libraries(DensestSubgraph utility reduction allocation extraction verification)

add_executable(debug test.cpp)
target_link_libraries(debug utility reduction allocation extraction verification)

add_subdirectory(utility)
add_subdirectory(reduction)
add_subdirectory(allocation)
add_subdirectory(extraction)
add_subdirectory(verification)
