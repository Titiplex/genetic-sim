cmake_minimum_required(VERSION 3.31)
project(Evo)

find_package(glfw3 REQUIRED)
find_package(glm CONFIG REQUIRED)

add_executable(evo main.cpp glad/src/glad.c)
target_include_directories(evo
        PRIVATE ${CMAKE_SOURCE_DIR}/include
        PRIVATE ${CMAKE_SOURCE_DIR}/glad/include
)
target_link_libraries(evo PRIVATE glfw opengl32 glm::glm)
