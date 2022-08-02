cmake_minimum_required(VERSION 3.21)
project(final C)

set(CMAKE_C_STANDARD 90)

add_executable(final main.c matrix_calculations.c spkmeansmodule.c spkmeans.h useful.h useful.c kmeans.h jacobi.c jacobi.h matrix_calculations.h)
