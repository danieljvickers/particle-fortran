module initialize_data
    use constants
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

    subroutine get_random_number(val, lower_bound, upper_bound)
        real(8), intent(in)  :: lower_bound, upper_bound
        real(8), intent(out) :: val

        call random_number(val)
        val = lower_bound + (upper_bound - lower_bound) * val
    end subroutine get_random_number

    subroutine initialize_particles()
        integer :: i
        real(8) :: temp_val_1, temp_val_2  ! placeholders for random number generator

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

        call random_seed()
        do i = 1, num_particles
            ! randomly distribute the positions of the particles in the allowed anulus
            call get_random_number(temp_val_1, radius_lower**2, radius_upper**2)  ! temp radius
            temp_val_1 = sqrt(temp_val_1)
            call get_random_number(temp_val_2, dble(0.0), 2.0*C_PI)  ! temp angle
            asteroids%x(i) = temp_val_1 * cos(temp_val_2)
            asteroids%y(i) = temp_val_1 * sin(temp_val_2)

            call get_random_number(asteroids%m(i), mass_lower, mass_upper)  ! evenly distribute the mass
            asteroids%r(i) = (asteroids%m(i) / C_Density * 0.75 / C_PI)**(1.0/3.0)  ! compute the size of the object from the mass and density

            ! generate the momentums
            temp_val_2 = sqrt(C_G * solar_mass / temp_val_1) * asteroids%m(i)  ! compute the optimal orbital momentum magnitude from the mass and orbital radius
            call get_random_number(temp_val_2, temp_val_2 * (1-velocity_noise_bound), temp_val_2 * (1+velocity_noise_bound))  ! get a random momentum magnitude
            asteroids%px(i) = temp_val_2 * asteroids%y(i) / temp_val_1  ! assign x and y momentum to be in 
            asteroids%py(i) = -1.0 * temp_val_2 * asteroids%x(i) / temp_val_1

            print *, "x=", asteroids%x(i), "y=", asteroids%y(i)
            print *, "px=", asteroids%px(i), "py=", asteroids%py(i)
            print *, "r=", asteroids%r(i), "m=", asteroids%m(i)
        end do

    end subroutine initialize_particles

end module initialize_data
