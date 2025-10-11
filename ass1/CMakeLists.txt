cmake_minimum_required(VERSION 3.10)
project(ass1)

set(CMAKE_C_STANDARD 11)

include_directories(/opt/homebrew/Cellar/mpich/4.3.1/include)
link_directories(/opt/homebrew/Cellar/mpich/4.3.1/lib)

add_executable(main main.c
        koh.h
        koh.c)

target_link_libraries(main mpi)