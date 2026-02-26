cmake_minimum_required(VERSION 3.16)
project(Evo)

find_package(glfw3 REQUIRED)
find_package(glm CONFIG REQUIRED)

add_executable(evo main.cpp glad/src/glad.c)
target_include_directories(evo
        PRIVATE ${CMAKE_SOURCE_DIR}/include
        PRIVATE ${CMAKE_SOURCE_DIR}/glad/include
)
if (WIN32)
    target_link_libraries(evo PRIVATE glfw opengl32 glm::glm)
else ()
    find_package(OpenGL REQUIRED)
    target_link_libraries(evo PRIVATE glfw OpenGL::GL glm::glm)
endif ()
