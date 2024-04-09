cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(utility SHARED
        graph.cpp
        flownetwork.cpp
        args.cpp
        ratioselection.cpp
        xycore.cpp
        lp.cpp
        wcore.cpp
        app.cpp
        heap.cpp
)

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(utility PUBLIC OpenMP::OpenMP_CXX)
endif()
target_include_directories(utility PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
