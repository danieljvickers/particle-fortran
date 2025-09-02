cmake_minimum_required(VERSION 3.15)

# ---- Project setup ----
project(MyFortranProject
    VERSION 0.1
    LANGUAGES Fortran
)

# ---- Compiler settings ----
enable_language(Fortran)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

option(USE_GPU "Enable GPU offloading" OFF)
if(USE_GPU)
    add_compile_definitions(USE_GPU)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp -foffload=nvptx-none")
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
endif()

# ---- Executable target ----
add_subdirectory(src)
