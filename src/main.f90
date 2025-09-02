! program hello_omp
!     use omp_lib
!     integer :: tid, nthreads
!     !$omp parallel private(tid)
!         tid = omp_get_thread_num()
!         nthreads = omp_get_num_threads()
!         print *, 'Thread ', tid, ' of ', nthreads
!     !$omp end parallel
! end program hello_omp


program main
    use omp_lib

    use constants
    use initialize_data
    use time_step
    use save_data
    implicit none


    ! contains the variables that we will use to intialize
    real(8) :: escape_radius, dt
    integer :: num_time_steps, i, save_frequency, time_frequency

    integer :: count_start, count_end, count_rate, total_count_start, total_count_end
    real(8) :: elapsed, total_elapsed

    integer :: merged_in_sun, flew_to_infinity, merged_together
    
    ! Variables used when initilizing the arrays
    num_particles = 8192  ! number of particles in the simulation
    mass_lower = 1e24  ! lower-bound mass of an astroid
    mass_upper = 1e25  ! upper-bound mass of an asteroid
    radius_lower = 1.082e11  ! lower-bound orbital radius of an asteroid, currently orbital radius of venus
    radius_upper = 1.5e11  ! upper-bound orbital radius of an asteroid, currently orbital radius of earth
    velocity_noise_bound = 0.1  ! the upper bound of the fraction of the velocity that will be perturbed from the ideal orbital velocity

    ! time step variables
    escape_radius = 20 * radius_upper
    num_time_steps = 100000
    dt = 60*60*6
    save_frequency = 10
    time_frequency = 100

    ! counting metrics
    merged_in_sun = 0
    flew_to_infinity = 0
    merged_together = 0

    ! initialize particles
    call initialize_particles()

    ! take time steps
    call system_clock(count_rate=count_rate)
    elapsed = dble(0.0)

    call system_clock(total_count_start)
    call handle_collisions(particles, num_particles, escape_radius, merged_in_sun, flew_to_infinity, merged_together)   ! needs to be called once to handle all of the particles that may be overlapping
    call save_all(particles, num_particles, "data", 0)
    do i = 1, num_time_steps
        call system_clock(count_start)

        call take_time_step(particles, num_particles, dt)
        call handle_collisions(particles, num_particles, escape_radius, merged_in_sun, flew_to_infinity, merged_together)

        call system_clock(count_end)
        elapsed = elapsed + real(count_end - count_start, 8) / real(count_rate, 8)

        if (modulo(i, save_frequency) .eq. 0) then
            call save_all(particles, num_particles, "data", i)
        end if
        if (modulo(i, time_frequency) .eq. 0) then
            print '(I5,A,I0,A,F8.4)', i, "    Elapsed time for previous " , time_frequency, " time steps is: ", elapsed
            elapsed = dble(0.0)
        end if
    end do
    call system_clock(total_count_end)

    total_elapsed = real(total_count_end - total_count_start, 8) / real(count_rate, 8)
    print *, "Average time per time step: ", total_elapsed / num_time_steps

    ! Count Metrics
    print *, "Particles Removed:    ", flew_to_infinity + merged_in_sun + merged_together
    print *, "Particles Fell into the Sun:    ", flew_to_infinity + merged_in_sun + merged_together
    print *, "Particles Flew Off to Infinity:    ", flew_to_infinity + merged_in_sun + merged_together
    print *, "Particles Merged Together:    ", flew_to_infinity + merged_in_sun + merged_together
    print *, "Particles remaining: ", num_particles - flew_to_infinity - merged_in_sun - merged_together

end program main
