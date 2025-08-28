# Library for all module code
add_library(simulation_lib
    m_initialize.f90
)

set_target_properties(simulation_lib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
)

target_include_directories(simulation_lib PUBLIC
    ${CMAKE_BINARY_DIR}/mod
)