cmake_minimum_required(VERSION 3.20)

project(genetic_sim_v3 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(EVO_ENABLE_WARNINGS "Enable strict warnings" ON)

# ---- Dependencies: GLFW + GLM via FetchContent ----
include(FetchContent)

FetchContent_Declare(
        glfw
        GIT_REPOSITORY https://github.com/glfw/glfw.git
        GIT_TAG 3.4
)
FetchContent_MakeAvailable(glfw)

FetchContent_Declare(
        glm
        GIT_REPOSITORY https://github.com/g-truc/glm.git
        GIT_TAG 1.0.1
)
FetchContent_MakeAvailable(glm)

include(FetchContent)

FetchContent_Declare(
        imgui
        GIT_REPOSITORY https://github.com/ocornut/imgui.git
        GIT_TAG docking
)
FetchContent_MakeAvailable(imgui)

# ---- GLAD (pre-generated) ----
add_library(glad STATIC
        third_party/glad/src/glad.c
)
target_include_directories(glad PUBLIC
        third_party/glad/include
)

# ---- ImGUI ----

add_library(imgui_lib
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_demo.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)

target_include_directories(imgui_lib PUBLIC
        ${imgui_SOURCE_DIR}
        ${imgui_SOURCE_DIR}/backends
)

target_link_libraries(imgui_lib PUBLIC glfw)

# ---- App ----
add_executable(genetic_sim
        src/main.cpp
        src/Environment.cpp
        src/Genetics.cpp
        src/Sim.cpp
        src/Renderer.cpp
)

target_include_directories(genetic_sim PRIVATE
        include
)

target_link_libraries(genetic_sim PRIVATE
        glad
        glfw
        glm
        imgui_lib
)

if (WIN32)
    target_link_libraries(genetic_sim PRIVATE opengl32)
elseif (APPLE)
    find_library(COCOA Cocoa REQUIRED)
    find_library(OpenGL OpenGL REQUIRED)
    find_library(IOKit IOKit REQUIRED)
    find_library(CoreVideo CoreVideo REQUIRED)
    target_link_libraries(genetic_sim PRIVATE ${COCOA} ${OpenGL} ${IOKit} ${CoreVideo})
else ()
    find_package(OpenGL REQUIRED)
    target_link_libraries(genetic_sim PRIVATE OpenGL::GL)
endif ()

if (EVO_ENABLE_WARNINGS)
    if (MSVC)
        target_compile_options(genetic_sim PRIVATE /W4 /permissive-)
    else ()
        target_compile_options(genetic_sim PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    endif ()
endif ()