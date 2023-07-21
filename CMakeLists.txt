cmake_minimum_required(VERSION 3.5)

project(tinyrender LANGUAGES CXX VERSION 0.0.1)
set(CMAKE_CXX_FLAGS "-std=c++17 -Wall")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_SOURCE_DIR}/src")
file(GLOB_RECURSE lib_SOURCES "src/*.cpp" "src/*/*.cpp")
add_library(tinyrender SHARED ${lib_SOURCES})
target_compile_definitions(tinyrender PRIVATE WIN_EXPORT)

file(GLOB_RECURSE test_SOURCES "test/*.cpp")
add_executable(line-test ${test_SOURCES})
target_include_directories(line-test PRIVATE "${PROJECT_SOURCE_DIR}/test")
target_link_libraries(line-test tinyrender)

add_custom_target(link_file ALL DEPENDS line-test)
add_custom_command (TARGET link_file
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink 
    "${PROJECT_SOURCE_DIR}/test/obj" "$<TARGET_FILE_DIR:line-test>/obj"
)

enable_testing()

add_test(
    NAME catch_test
    COMMAND $<TARGET_FILE:line-test> --success
)