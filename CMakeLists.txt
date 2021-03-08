cmake_minimum_required(VERSION 3.9)
project(NGP C)

set(EXEC fortune)

file(GLOB_RECURSE sources "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")

include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/include
    )

add_library(applib
    src/inputs.c
    src/rbtree.c
    src/fortune.c
    src/point.c
    src/globals.c
    src/pqueue.c
    src/event.c
    src/inputs.c
    src/halfedge.c
    src/animation.c
    )

add_executable(${EXEC}
    src/main.c
    )

target_include_directories(${EXEC} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
set_target_properties(${EXEC} PROPERTIES
                      C_STANDARD 99
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

add_subdirectory(deps/BOV)
target_link_libraries(${EXEC} applib bov m)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${EXEC})
