module initialize_data
    implicit none

    type :: Asteroid_t
        ! Positions
        real(8), allocatable :: x(:)
        real(8), allocatable :: y(:)

        ! Momentums
        real(8), allocatable :: px(:)
        real(8), allocatable :: py(:)

        ! Other variables
        real(8), allocatable :: m(:)  ! mass
        real(8), allocatable :: r(:)  ! radius, computed from mass
    end type Asteroid_t

    type(Asteroid_t) :: asteroids

    ! decleare the starting variables
    integer :: num_particles
    real(8) :: solar_mass, solar_radius
    real(8) :: mass_lower, mass_upper, radius_lower, radius_upper
    real(8) :: velocity_noise_bound

contains

    subroutine initialize_particles()
        integer :: i

        ! Variables used when initilizing the arrays
        num_particles = 1  ! number of particles in the simulation
        solar_mass = 1.9891e30  ! mass of the central star, which is the mass of the sun for now
        solar_radius = 696e6  ! radius of the sun, in meters
        mass_lower = 1e18  ! lower-bound mass of an astroid
        mass_upper = 1e19  ! upper-bound mass of an asteroid
        radius_lower = 1.082e11  ! lower-bound orbital radius of an asteroid, currently orbital radius of venus
        radius_upper = 1.5e11  ! upper-bound orbital radius of an asteroid, currently orbital radius of earth
        velocity_noise_bound = 0.1  ! the upper bound of the fraction of the velocity that will be perturbed from the ideal orbital velocity

        ! Allocate arrays
        allocate(asteroids%x(num_particles), asteroids%y(num_particles))  ! positions
        allocate(asteroids%px(num_particles), asteroids%py(num_particles))  ! momentums
        allocate(asteroids%m(num_particles), asteroids%r(num_particles))  ! other

        print *,  "Hello from initialize_data. The radius of the sun is ", solar_radius

    end subroutine initialize_particles

end module initialize_data
