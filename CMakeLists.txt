cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")

project(DensestSubgraph)
find_package(Boost 1.66.0)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

include_directories(
        ${PROJECT_SOURCE_DIR}/
)

add_executable(DensestSubgraph main.cpp)
target_link_libraries(DensestSubgraph utility reduction allocation extraction verification)

#verificationadd_executable(debug test.cpp)
#target_link_libraries(debug utility reduction allocation extraction verification)

add_subdirectory(utility)
add_subdirectory(reduction)
add_subdirectory(allocation)
add_subdirectory(extraction)
add_subdirectory(verification)
