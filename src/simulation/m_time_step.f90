module time_step
    use omp_lib

    use constants
    use initialize_data
    implicit none

contains

    subroutine get_distance(distance, x1, y1, x2, y2)

        real(8), intent(in) :: x1, y1, x2, y2
        real(8), intent(out) :: distance

        distance = sqrt((x1 - x2)**2 + (y1 - y2)**2)

    end subroutine get_distance


    subroutine handle_collisions(particles, num_particles, max_allowed_distance)

        integer :: i, j
        real(8) :: distance, x_com, y_com, combined_mass, merge_distance

        type(particle_t), intent(inout) :: particles
        integer, intent(in) :: num_particles
        real(8), intent(in) :: max_allowed_distance

        !$omp parallel do private(i, j, distance)
        do i = 1, num_particles
            if (particles%merged(i)) then
                cycle  ! skips if collided
            end if

            ! check for solar merges
            call get_distance(distance, particles%x(i), particles%y(i), dble(0.0), dble(0.0))
            if ((distance .le. C_R_s) .or. (distance .ge. max_allowed_distance)) then
                ! TODO :: for now we assume the particle is so small that it has no mass comapred to sun
                print *, "Merging into sun ", i
                particles%merged(i) = .true.
                cycle
            end if
        end do
        !$omp end parallel do

        ! check for particle-to-particle merging
        ! TODO :: Consider swaping to an iterative parallel search and then serial merge.
        !$omp parallel do private(i, j, distance, combined_mass)
        do i = 1, num_particles
            do j = i+1, num_particles
                if (particles%merged(j)) then
                    cycle  ! skips if the second particles has collided
                end if

                call get_distance(distance, particles%x(i), particles%y(i), particles%x(j), particles%y(j))
                if ( distance .le. particles%r(i) + particles%r(j)) then  ! true if they should collide perfectly inelasically

                    !$omp critical
                    ! add the momentum
                    particles%px(i) = particles%px(i) + particles%px(j)
                    particles%py(i) = particles%py(i) + particles%py(j)

                    ! compute the center of mass and move particles
                    particles%x(i) = (particles%x(i) * particles%m(i) + particles%x(j) * particles%m(j)) / combined_mass
                    particles%y(i) = (particles%y(i) * particles%m(i) + particles%y(j) * particles%m(j)) / combined_mass
                    
                    ! update mass and radius
                    combined_mass = particles%m(i) + particles%m(j)
                    particles%m(i) = combined_mass
                    particles%r(i) = (particles%m(i) / C_Density * 0.75 / C_PI)**(1.0/3.0)

                    ! count the second particle as merged
                    particles%merged(j) = .True.
                    !$omp end critical
                end if
            end do
        end do

    end subroutine handle_collisions


    subroutine take_time_step(particles, num_particles, dt)

        integer :: i, j, tid, n_threads
        real(8), allocatable :: ax_local(:,:), ay_local(:,:)
        real(8) :: acceleration_x, acceleration_y, velocity_x, velocity_y, distance

        type(particle_t), intent(inout) :: particles
        integer, intent(in) :: num_particles
        real(8), intent(in) :: dt

        !$omp parallel
        n_threads = omp_get_num_threads()
        !$omp end parallel

        ! reset the memory on acceleration
        allocate(ax_local(num_particles, n_threads))
        allocate(ay_local(num_particles, n_threads))
        ax_local = 0.0d0
        ay_local = 0.0d0

        !$omp parallel default(shared) private(i,j,distance,acceleration_x,acceleration_y,tid)
        tid = omp_get_thread_num()

        !$omp do schedule(dynamic)
        do i = 1, num_particles

            call get_distance(distance, particles%x(i), particles%y(i), dble(0.0), dble(0.0))
            ax_local(i, tid+1) = ax_local(i, tid+1) - (C_G * C_M_s * particles%x(i) / (distance**3))
            ay_local(i, tid+1) = ay_local(i, tid+1) - (C_G * C_M_s * particles%y(i) / (distance**3))

            do j = i+1, num_particles
                if (particles%merged(j)) cycle

                call get_distance(distance, particles%x(i), particles%y(i), particles%x(j), particles%y(j))
                acceleration_x = -C_G * particles%m(i) * particles%m(j) * (particles%x(i) - particles%x(j)) / (distance**3)
                acceleration_y = -C_G * particles%m(i) * particles%m(j) * (particles%y(i) - particles%y(j)) / (distance**3)

                ! We update each acceleration, using attomics to prevent race conditions in parallel
                ax_local(i, tid+1) = ax_local(i, tid+1) + acceleration_x / particles%m(i)
                ax_local(j, tid+1) = ax_local(j, tid+1) - acceleration_x / particles%m(j)
                ay_local(i, tid+1) = ay_local(i, tid+1) + acceleration_y / particles%m(i)
                ay_local(j, tid+1) = ay_local(j, tid+1) - acceleration_y / particles%m(j)

            end do
        end do
        !$omp end do
        !$omp end parallel

        ! Reduce private arrays into the shared particles%ax, %ay
        particles%ax = 0.0d0
        particles%ay = 0.0d0
        do tid = 1, n_threads
            particles%ax = particles%ax + ax_local(:, tid)
            particles%ay = particles%ay + ay_local(:, tid)
        end do

        deallocate(ax_local, ay_local)

        ! use those accelerations to compute the updated positions and momentums using eulers method
        !$omp parallel do default(shared) private(i, velocity_x, velocity_y)
        do i = 1, num_particles

            velocity_x = (particles%px(i) / particles%m(i)) + dt * particles%ax(i)
            velocity_y = (particles%py(i) / particles%m(i)) + dt * particles%ay(i)

            particles%x(i) = particles%x(i) + dt * velocity_x
            particles%y(i) = particles%y(i) + dt * velocity_y
            particles%px(i) = particles%m(i) * velocity_x
            particles%py(i) = particles%m(i) * velocity_y

        end do
        !$omp end parallel do

    end subroutine take_time_step

end module time_step
