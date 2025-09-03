module initialize_data
    use constants
    implicit none

    type :: particle_type
        ! Positions
        real(8), allocatable :: x(:)
        real(8), allocatable :: y(:)

        ! Momentums
        real(8), allocatable :: px(:)
        real(8), allocatable :: py(:)

        ! Accelerations
        real(8), allocatable :: ax(:)
        real(8), allocatable :: ay(:)

        ! Other variables
        real(8), allocatable :: m(:)  ! mass
        real(8), allocatable :: r(:)  ! radius, computed from mass
        logical, allocatable :: merged(:)
    end type particle_type

    type(particle_type) :: particles

    ! decleare the starting variables
    integer :: num_particles
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
        real(8) :: orbital_radius, start_angle, orbital_momentum  ! placeholders for random number generator

        ! Allocate arrays
        allocate(particles%x(num_particles), particles%y(num_particles))  ! positions
        allocate(particles%px(num_particles), particles%py(num_particles))  ! momentums
        allocate(particles%ax(num_particles), particles%ay(num_particles))  ! accelerations
        allocate(particles%m(num_particles), particles%r(num_particles))  ! other
        allocate(particles%merged(num_particles))  ! other

        call random_seed()
        do i = 1, num_particles
            ! randomly distribute the positions of the particles in the allowed anulus
            call get_random_number(orbital_radius, radius_lower**2, radius_upper**2)  ! temp radius
            orbital_radius = sqrt(orbital_radius)
            call get_random_number(start_angle, dble(0.0), 2.0*C_PI)  ! temp angle
            particles%x(i) = orbital_radius * cos(start_angle)
            particles%y(i) = orbital_radius * sin(start_angle)

            call get_random_number(particles%m(i), mass_lower, mass_upper)  ! evenly distribute the mass
            particles%r(i) = (particles%m(i) / C_Density * 0.75 / C_PI)**(1.0/3.0)  ! compute the size of the object from the mass and density

            ! generate the momentums
            orbital_momentum = sqrt(C_G * C_M_s / orbital_radius) * particles%m(i)  ! compute the optimal orbital momentum magnitude from the mass and orbital radius
            call get_random_number(orbital_momentum, orbital_momentum * (1-velocity_noise_bound), &
                orbital_momentum * (1+velocity_noise_bound))  ! get a random momentum magnitude
            particles%px(i) =  orbital_momentum * particles%y(i) / orbital_radius  ! assign x and y momentum to be in 
            particles%py(i) = -orbital_momentum * particles%x(i) / orbital_radius

            particles%merged(i) = .False.

        end do

    end subroutine initialize_particles

end module initialize_data
