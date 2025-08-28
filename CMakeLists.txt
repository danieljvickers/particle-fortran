cmake_minimum_required(VERSION 3.15)

# ---- Project setup ----
project(MyFortranProject
    VERSION 0.1
    LANGUAGES Fortran
)

# ---- Compiler settings ----
enable_language(Fortran)

# Require MPI
find_package(MPI REQUIRED)

# Require OpenMP
find_package(OpenMP REQUIRED)

# ---- Executable target ----
add_executable(particle_fortran
    src/main.f90
)

# Link against MPI and OpenMP
target_link_libraries(particle_fortran
    PRIVATE
        MPI::MPI_Fortran
        OpenMP::OpenMP_Fortran
)

# Optional: stricter compiler flags (example for gfortran/ifort)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(particle_fortran PRIVATE -Wall -Wextra -O2)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(particle_fortran PRIVATE -warn all -O2)
endif()

# ---- Installation (optional) ----
install(TARGETS particle_fortran RUNTIME DESTINATION bin)
