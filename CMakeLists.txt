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

# ---- Executable target ----
add_subdirectory(src)
