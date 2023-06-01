cmake_minimum_required(VERSION 3.25)
project(chopper_lib C)

set(CMAKE_C_STANDARD 11)

add_library(chopper_lib SHARED chopper-lib.c)
