cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--stack,104857600")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-stack_size,104857600")
endif()

project(DensestSubgraph)
find_package(Boost 1.66.0)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

#find_package(OpenMP REQUIRED)

include_directories(
        ${PROJECT_SOURCE_DIR}/
)

add_executable(DensestSubgraph main.cpp)
target_link_libraries(DensestSubgraph utility reduction allocation extraction verification)
if(OpenMP_CXX_FOUND)
    target_link_libraries(DensestSubgraph OpenMP::OpenMP_CXX)
endif()

#add_executable(debug test.cpp)
#target_link_libraries(debug utility reduction allocation extraction verification)
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(debug PUBLIC OpenMP::OpenMP_CXX)
#endif()

add_subdirectory(utility)
add_subdirectory(reduction)
add_subdirectory(allocation)
add_subdirectory(extraction)
add_subdirectory(verification)
