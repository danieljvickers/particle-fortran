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

option(USE_GPU "Enable GPU offloading" ON)

if(USE_GPU)
    add_compile_definitions(USE_GPU)

    # Detect if we are using NVHPC
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
        # OpenMP GPU offload flags for NVHPC
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp=gpu -Minfo=accel")
    else()
        # Fallback for GNU or other compilers
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp -foffload=nvptx-none")
    endif()

else()
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp")
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
    endif()
endif()

# ---- Executable target ----
add_subdirectory(src)
