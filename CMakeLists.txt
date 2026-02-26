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


option(EVO_ENABLE_OPENMP "Enable OpenMP acceleration" ON)
if (EVO_ENABLE_OPENMP)
    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        target_link_libraries(evo PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(evo PRIVATE EVO_HAS_OPENMP=1)
    endif ()
endif ()
