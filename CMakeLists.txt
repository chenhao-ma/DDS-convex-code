cmake_minimum_required(VERSION 3.10)
project(DDSapp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3")

add_executable(DDSapp main.cpp Graph.cpp Graph.h Args.cpp Args.h Edge.cpp Edge.h EdgeFN.cpp EdgeFN.h FlowNetwork.cpp FlowNetwork.h)