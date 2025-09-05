module initialize_data
    use constants
    implicit none

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
    integer, allocatable :: merged(:)


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
        real(8) :: orbital_radius, start_angle, orbital_momentum_mag, random_momentum  ! placeholders for random number generator
        integer, allocatable :: seed(:)

        ! Allocate arrays
        allocate(x(num_particles), y(num_particles))  ! positions
        allocate(px(num_particles), py(num_particles))  ! momentums
        allocate(ax(num_particles), ay(num_particles))  ! accelerations
        allocate(m(num_particles), r(num_particles))  ! other
        allocate(merged(num_particles))  ! other
        allocate(seed(2))
        seed(1) = 6
        seed(2) = 9

        call random_seed(put=seed)
        do i = 1, num_particles
            ! randomly distribute the positions of the particles in the allowed anulus
            call get_random_number(orbital_radius, radius_lower**2, radius_upper**2)  ! temp radius
            orbital_radius = sqrt(orbital_radius)
            call get_random_number(start_angle, dble(0.0), 2.0*C_PI)  ! temp angle
            x(i) = orbital_radius * cos(start_angle)
            y(i) = orbital_radius * sin(start_angle)

            call get_random_number(m(i), mass_lower, mass_upper)  ! evenly distribute the mass
            r(i) = (m(i) / C_Density * 0.75 / C_PI)**(1.0/3.0)  ! compute the size of the object from the mass and density

            ! generate the momentums
            orbital_momentum_mag = sqrt(C_G * C_M_s / orbital_radius) * m(i)  ! compute the optimal orbital momentum magnitude from the mass and orbital radius
            px(i) =  orbital_momentum_mag * y(i) / orbital_radius  ! assign x and y momentum to be in 
            py(i) = -orbital_momentum_mag * x(i) / orbital_radius
            call get_random_number(random_momentum, -orbital_momentum_mag * velocity_noise_bound, &
                orbital_momentum_mag * velocity_noise_bound)  ! get a random momentum magnitude
            px(i) = px(i) + random_momentum
            call get_random_number(random_momentum, -orbital_momentum_mag * velocity_noise_bound, &
                orbital_momentum_mag * velocity_noise_bound)  ! get a random momentum magnitude
            py(i) = py(i) + random_momentum

            merged(i) = 0

        end do

    end subroutine initialize_particles

end module initialize_data
