module time_step
    use constants
    use initialize_data
    implicit none

contains

    subroutine handle_collisions(particles, num_particles)

        integer :: i, j
        real(8) :: distance, x_com, y_com, combined_mass

        type(particle_t), intent(inout) :: particles
        integer, intent(in) :: num_particles

        do i = 1, num_particles - 1
            if (particles%merged(i)) then
                cycle  ! skips if the first particles has collided
            end if

            ! check for solar merges
            distance = sqrt(particles%x(i)**2 + particles%y(i)**2)
            if (distance .le. C_R_s) then
                ! TODO :: for now we assume the particle is so small that it has no mass comapred to sun
                particles%merged(i) = .true.
                cycle
            end if

            ! check for particle-to-particle merging
            do j = i+1, num_particles
                if (particles%merged(j)) then
                    cycle  ! skips if the second particles has collided
                end if

                distance = sqrt((particles%x(i) - particles%x(j))**2 + (particles%y(i) - particles%y(j))**2)
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

    
    subroutine take_time_step(particles, num_particles)

    end subroutine take_time_step

end module time_step
