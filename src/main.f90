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
    num_particles = 65536  ! number of particles in the simulation
    mass_lower = 1e26  ! lower-bound mass of an astroid
    mass_upper = 1e27  ! upper-bound mass of an asteroid
    radius_lower = 0.5e11  ! lower-bound orbital radius of an asteroid, currently orbital radius of venus
    radius_upper = 2.5e11  ! upper-bound orbital radius of an asteroid, currently orbital radius of earth
    velocity_noise_bound = 0.001  ! the upper bound of the fraction of the velocity that will be perturbed from the ideal orbital velocity

    ! time step variables
    escape_radius = 20 * radius_upper
    num_time_steps = 1000
    dt = 60*60*2
    save_frequency = 10
    time_frequency = 1000

    ! counting metrics
    merged_in_sun = 0
    flew_to_infinity = 0
    merged_together = 0

    ! initialize particles
    call initialize_particles()

#ifdef USE_GPU
    ! Move all particle arrays to the device once and keep them there
    !$omp target enter data map(to: x(1:num_particles), y(1:num_particles), &
    !$omp&                        px(1:num_particles), py(1:num_particles), &
    !$omp&                        ax(1:num_particles), ay(1:num_particles), &
    !$omp&                        m(1:num_particles), r(1:num_particles), &
    !$omp&                        merged(1:num_particles))
#endif

    ! take time steps
    call system_clock(count_rate=count_rate)
    elapsed = dble(0.0)

    call system_clock(total_count_start)
    ! call handle_collisions(particles, num_particles, escape_radius, merged_in_sun, flew_to_infinity, merged_together)   ! needs to be called once to handle all of the particles that may be overlapping
    call save_all("data", 0)
    do i = 1, num_time_steps
        call system_clock(count_start)

        call take_time_step(dt)
        ! call handle_collisions(particles, num_particles, escape_radius, merged_in_sun, flew_to_infinity, merged_together)

        call system_clock(count_end)
        elapsed = elapsed + real(count_end - count_start, 8) / real(count_rate, 8)

        if (modulo(i, save_frequency) .eq. 0) then
            call save_all("data", i)
        end if
        if (modulo(i, time_frequency) .eq. 0) then
            print '(I5,A,I0,A,F8.4)', i, "    Elapsed time for previous " , time_frequency, " time steps is: ", elapsed
            elapsed = dble(0.0)
        end if
    end do
    call system_clock(total_count_end)

#ifdef USE_GPU
    ! Bring back results
    !$omp target exit data map(from: x(1:num_particles), y(1:num_particles), &
    !$omp&                        px(1:num_particles), py(1:num_particles))
#endif

    total_elapsed = real(total_count_end - total_count_start, 8) / real(count_rate, 8)
    print *, "Average time per time step: ", total_elapsed / num_time_steps

    ! Count Metrics
    print *, "Particles Removed:              ", flew_to_infinity + merged_in_sun + merged_together
    print *, "Particles Fell into the Sun:    ", merged_in_sun
    print *, "Particles Flew Off to Infinity: ", flew_to_infinity
    print *, "Particles Merged Together:      ", merged_together
    print *, "Particles remaining:            ", num_particles - flew_to_infinity - merged_in_sun - merged_together

end program main
