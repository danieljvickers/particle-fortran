module time_step
    use omp_lib

    use constants
    use initialize_data
    implicit none

contains

#ifdef USE_GPU
    !$omp declare target
#endif
    subroutine get_distance(distance, x1, y1, x2, y2)

        real(8), intent(in) :: x1, y1, x2, y2
        real(8), intent(out) :: distance

        distance = sqrt((x1 - x2)**2 + (y1 - y2)**2)

    end subroutine get_distance
#ifndef USE_GPU
    !$omp end declare target
#endif


    subroutine handle_collisions( max_allowed_distance, merged_in_sun, flew_to_infinity, merged_together)

        integer :: i, j
        real(8) :: distance, x_com, y_com, combined_mass, merge_distance

        real(8), intent(in) :: max_allowed_distance
        integer, intent(inout) :: merged_in_sun, flew_to_infinity, merged_together

#ifdef USE_GPU
        
#else
        !$omp parallel do private(i, j, distance)
#endif
        do i = 1, num_particles
            if (merged(i)) cycle

            ! check for solar merges
            call get_distance(distance, x(i), y(i), dble(0.0), dble(0.0))
            if ((distance .le. C_R_s)) then
                ! TODO :: for now we assume the particle is so small that it has no mass comapred to sun
                merged_in_sun = merged_in_sun + 1
                merged(i) = .true.
            else if (distance .ge. max_allowed_distance) then
                flew_to_infinity = flew_to_infinity + 1
                merged(i) = .true.
            end if
        end do
#ifndef USE_GPU
        !$omp end parallel do
#endif

        ! check for particle-to-particle merging
        ! TODO :: Consider swaping to an iterative parallel search and then serial merge.
        !$omp parallel do private(i, j, distance, combined_mass)
        do i = 1, num_particles
            do j = i+1, num_particles
                if (merged(j)) then
                    cycle  ! skips if the second particles has collided
                end if

                call get_distance(distance, x(i), y(i), x(j), y(j))
                if ( distance .le. r(i) + r(j)) then  ! true if they should collide perfectly inelasically

                    !$omp critical
                    merged_together = merged_together + 1
                    ! add the momentum
                    px(i) = px(i) + px(j)
                    py(i) = py(i) + py(j)

                    ! compute the center of mass and move particles
                    combined_mass = m(i) + m(j)
                    x(i) = (x(i) * m(i) + x(j) * m(j)) / combined_mass
                    y(i) = (y(i) * m(i) + y(j) * m(j)) / combined_mass
                    
                    ! update mass and radius
                    m(i) = combined_mass
                    r(i) = (m(i) / C_Density * 0.75 / C_PI)**(1.0/3.0)

                    ! count the second particle as merged
                    merged(j) = .True.
                    !$omp end critical
                end if
            end do
        end do

    end subroutine handle_collisions


    subroutine take_time_step(dt)

        integer :: i, j, tid, n_threads
        real(8) :: acceleration_x, acceleration_y, velocity_x, velocity_y, distance

        real(8), intent(in) :: dt

#ifdef USE_GPU
        !$omp target teams distribute parallel do
#else
        !$omp parallel do default(shared) private(i, distance)
#endif
        do i = 1, num_particles
            if (merged(i)) cycle
            call get_distance(distance, x(i), y(i), dble(0.0), dble(0.0))
            ax(i) = -C_G * C_M_s * x(i) / (distance**3)
            ay(i) =  -C_G * C_M_s * y(i) / (distance**3)
        end do
#ifndef USE_GPU
        !$omp end parallel do
#endif

#ifdef USE_GPU
        !$omp target teams distribute parallel do collapse(2) default(shared) private(distance, acceleration_x, acceleration_y)
#else
        !$omp parallel do default(shared) private(i, j, acceleration_x, acceleration_y)
#endif
        do i = 1, num_particles
            do j = i+1, num_particles
                if (merged(j) .or. merged(i)) cycle

                call get_distance(distance, x(i), y(i), x(j), y(j))
                acceleration_x = -C_G * m(i) * m(j) * (x(i) - x(j)) / (distance**3)
                acceleration_y = -C_G * m(i) * m(j) * (y(i) - y(j)) / (distance**3)

                ! We update each acceleration, using attomics to prevent race conditions in parallel
                !$omp atomic update
                ax(i) = ax(i) + acceleration_x / m(i)
                !$omp atomic update
                ax(j) = ax(j) - acceleration_x / m(j)
                !$omp atomic update
                ay(i) = ay(i) + acceleration_y / m(i)
                !$omp atomic update
                ay(j) = ay(j) - acceleration_y / m(j)

            end do
        end do
#ifndef USE_GPU
        !$omp end parallel do
#endif

#ifdef USE_GPU
        !$omp target teams distribute parallel do default(shared) private(i, velocity_x, velocity_y)
#else
        !$omp parallel do default(shared) private(i, velocity_x, velocity_y)
#endif
        do i = 1, num_particles
            velocity_x = (px(i) / m(i)) + dt * ax(i)
            velocity_y = (py(i) / m(i)) + dt * ay(i)

            x(i)  = x(i)  + dt * velocity_x
            y(i)  = y(i)  + dt * velocity_y
            px(i) = m(i) * velocity_x
            py(i) = m(i) * velocity_y
        end do
#ifndef USE_GPU
        !$omp end parallel do
#endif

    end subroutine take_time_step

end module time_step
