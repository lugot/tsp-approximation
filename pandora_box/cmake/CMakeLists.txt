cmake_minimum_required(VERSION 3.7)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(tsp_approximation)

find_package(CPLEX REQUIRED)
include_directories(include ${CPLEX_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include)

add_library(applib
    src/globals.c
    )

add_executable(tsp_approximation src/main.c)

target_link_libraries(tsp_approximation applib ${CPLEX_LIBRARIES})
