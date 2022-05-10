cmake_minimum_required(VERSION 3.20)
project(ParallelMatrix)

set(CMAKE_CXX_STANDARD 20)

add_executable(ParallelMatrix main.cpp Matrix.cpp Matrix.h utils.h)
