# Add library for simulation modules
add_subdirectory(simulation)

# Create executable
add_executable(particle_fortran main.f90)

# Link to simulation library, MPI, OpenMP
target_link_libraries(particle_fortran
    PRIVATE
        simulation_lib
        MPI::MPI_Fortran
        OpenMP::OpenMP_Fortran
)

# Optional: compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(particle_fortran PRIVATE -Wall -Wextra -O2)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(particle_fortran PRIVATE -warn all -O2)
endif()

# Add compile definition for GPU mode
target_compile_definitions(particle_fortran PUBLIC USE_GPU)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")