# Library for all module code
add_library(simulation_lib
    m_constants.f90
    m_initialize.f90
    m_time_step.f90
    m_save.f90
)

target_link_libraries(simulation_lib PUBLIC OpenMP::OpenMP_Fortran)

set_target_properties(simulation_lib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
)

target_include_directories(simulation_lib PUBLIC
    ${CMAKE_BINARY_DIR}/mod
)

# Add compile definition for GPU mode
target_compile_definitions(simulation_lib PUBLIC USE_GPU)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")