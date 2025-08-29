# Library for all module code
add_library(simulation_lib
    m_constants.f90
    m_initialize.f90
    m_time_step.f90
    m_save.f90
)

set_target_properties(simulation_lib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
)

target_include_directories(simulation_lib PUBLIC
    ${CMAKE_BINARY_DIR}/mod
)