module time_step
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
        real(8) :: distance, x_com, y_com, combined_mass

        type(particle_t), intent(inout) :: particles
        integer, intent(in) :: num_particles
        real(8), intent(in) :: max_allowed_distance

        do i = 1, num_particles - 1
            if (particles%merged(i)) then
                cycle  ! skips if the first particles has collided
            end if

            ! check for solar merges
            call get_distance(distance, particles%x(i), particles%y(i), dble(0.0), dble(0.0))
            if ((distance .le. C_R_s) .or. (distance .ge. max_allowed_distance)) then
                ! TODO :: for now we assume the particle is so small that it has no mass comapred to sun
                particles%merged(i) = .true.
                cycle
            end if

            ! check for particle-to-particle merging
            do j = i+1, num_particles
                if (particles%merged(j)) then
                    cycle  ! skips if the second particles has collided
                end if

                call get_distance(distance, particles%x(i), particles%y(i), particles%x(j), particles%y(j))
                if (distance .le. particles%r(i) + particles%r(j)) then  ! true if they should collide perfectly inelasically
                    ! add the momentum
                    particles%px(i) = particles%px(i) + particles%px(j)
                    particles%py(i) = particles%py(i) + particles%py(j)

                    ! compute the center of mass and move particles
                    combined_mass = particles%m(i) + particles%m(j)
                    particles%x(i) = (particles%x(i) * particles%m(i) + particles%x(j) * particles%m(j)) / combined_mass
                    particles%y(i) = (particles%y(i) * particles%m(i) + particles%y(j) * particles%m(j)) / combined_mass
                    
                    ! update mass and radius
                    particles%m(i) = combined_mass
                    particles%r(i) = (particles%m(i) / C_Density * 0.75 / C_PI)**(1.0/3.0)

                    ! count the second particle as merged
                    particles%merged(j) = .True.
                end if
            end do
        end do

    end subroutine handle_collisions


    subroutine take_time_step(particles, num_particles, dt)

        integer :: i, j
        real(8) :: acceleration_x, acceleration_y, velocity_x, velocity_y, distance

        type(particle_t), intent(inout) :: particles
        integer, intent(in) :: num_particles
        real(8), intent(in) :: dt

        ! reset the memory on acceleration
        do i = 1, num_particles
            particles%ax(i) = dble(0.0)
            particles%ay(i) = dble(0.0)
        end do

        do i = 1, num_particles

            call get_distance(distance, particles%x(i), particles%y(i), dble(0.0), dble(0.0))
            particles%ax(i) = particles%ax(i) - (C_G * C_M_s * particles%x(i) / (distance**3))
            particles%ay(i) = particles%ay(i) - (C_G * C_M_s * particles%y(i) / (distance**3))

            do j = i+1, num_particles
                if (particles%merged(j)) then
                    cycle  ! skips if the second particles has collided
                end if

                acceleration_x = -C_G * particles%m(i) * particles%m(j) * (particles%x(i) - particles%x(j)) / (distance**3)
                acceleration_y = -C_G * particles%m(i) * particles%m(j) * (particles%y(i) - particles%y(j)) / (distance**3)

                particles%ax(i) = particles%ax(i) + acceleration_x / particles%m(i)
                particles%ax(j) = particles%ax(j) - acceleration_x / particles%m(j)
                particles%ay(i) = particles%ay(i) + acceleration_y / particles%m(i)
                particles%ay(j) = particles%ay(j) - acceleration_y / particles%m(j)
            end do

        end do

        ! use those accelerations to compute the updated positions and momentums using eulers method
        do i = 1, num_particles

            velocity_x = (particles%px(i) / particles%m(i)) + dt * particles%ax(i)
            velocity_y = (particles%py(i) / particles%m(i)) + dt * particles%ay(i)

            particles%x(i) = particles%x(i) + dt * velocity_x
            particles%y(i) = particles%y(i) + dt * velocity_y
            particles%px(i) = particles%m(i) * velocity_x
            particles%py(i) = particles%m(i) * velocity_y

        end do

    end subroutine take_time_step

end module time_step
