module constants
    implicit none

    ! mathematical constants
    real(8), parameter :: C_PI = 3.14159265358979323846

    ! physical constants
    real(8), parameter :: C_G = 6.674e-11
    real(8), parameter :: C_Density = 200  ! kg/m^3 density of the asteroids
    real(8), parameter :: C_M_s = 1.9891e30  ! mass of the central star, which is the mass of the sun for now in kg
    real(8), parameter :: C_R_s = 696e6  ! radius of the sun, in meters

end module constants