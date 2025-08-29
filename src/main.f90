! program main
!     use mpi
!     use omp_lib
!     implicit none

!     integer :: ierr, rank, size
!     integer :: omp_rank, omp_size

!     ! Initialize MPI
!     call MPI_Init(ierr)
!     call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
!     call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

!     ! OpenMP parallel region
!     !$omp parallel private(omp_rank, omp_size)
!         omp_rank = omp_get_thread_num()
!         omp_size = omp_get_num_threads()

!         print *, 'Hello from MPI rank', rank, 'of', size, &
!                  'and OpenMP thread', omp_rank, 'of', omp_size
!     !$omp end parallel

!     ! Finalize MPI
!     call MPI_Finalize(ierr)

! end program main

program main
    use constants
    use initialize_data
    use time_step
    implicit none


    ! contains the variables that we will use to intialize
    real(8) :: escape_radius, dt
    integer :: num_time_steps, i, save_frequency
    
    ! Variables used when initilizing the arrays
    num_particles = 1  ! number of particles in the simulation
    mass_lower = 1e18  ! lower-bound mass of an astroid
    mass_upper = 1e19  ! upper-bound mass of an asteroid
    radius_lower = 1.082e11  ! lower-bound orbital radius of an asteroid, currently orbital radius of venus
    radius_upper = 1.5e11  ! upper-bound orbital radius of an asteroid, currently orbital radius of earth
    velocity_noise_bound = 0.1  ! the upper bound of the fraction of the velocity that will be perturbed from the ideal orbital velocity

    ! time step variables
    escape_radius = 40 * C_R_s
    num_time_steps = 1000
    dt = 60*60*6

    ! initialize particles
    call initialize_particles()

    ! take time steps
    call handle_collisions(particles, num_particles)   ! needs to be called once to handle all of the particles that may be overlapping
    do i = 1, num_time_steps
        call take_time_step(particles, num_particles)
        call handle_collisions(particles, num_particles)

        if (modulo(i, save_frequency) .eq. 0) then
            ! call the save function
            cycle
        end if
    end do

end program main