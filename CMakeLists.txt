cmake_minimum_required(VERSION 3.20.0)

set (CMAKE_CXX_STANDARD 23)
project(kissing_searcher VERSION 0.1.0 LANGUAGES C CXX)


set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -g")


add_executable(kissing_searcher main.cpp)

target_include_directories(kissing_searcher PUBLIC .)