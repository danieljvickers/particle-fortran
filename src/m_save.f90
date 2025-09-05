module save_data
    use initialize_data
    implicit none

contains

subroutine save_all(directory, index)
    
    integer, intent(in) :: index
    character(len=*), intent(in) :: directory

    character(len=256) :: filename
    integer :: unit, ios, i

    ! update data on GPU
#ifdef USE_GPU
    !$omp target update from(x, y, px, py, r, merged)
#endif
    
    write(filename, '(A,"/",I0,".bin")') trim(directory), index

    ! Open the file as binary (unformatted stream)
    open(newunit=unit, file=filename, form="unformatted", access="stream", &
         action="write", status="replace", iostat=ios)
    if (ios /= 0) then
        print *, "Error opening file: ", trim(filename)
        return
    end if

    ! Write interleaved doubles: x1,y1,x2,y2,...
    do i = 1, num_particles
        if (merged(i) .eq. 0) then
            write(unit) x(i), y(i), px(i), py(i), r(i)
        end if
    end do

    ! Close file
    close(unit)

end subroutine save_all


subroutine coallese_particles()

    integer :: num_shifted, i

    num_shifted = 0

#ifdef USE_GPU
    !$omp target update from(x, y, px, py, r, m, merged)
#endif

    do i=1, num_particles
        if (merged(i) .eq. 1) then
            num_shifted = num_shifted + 1
            merged(i) = 0
        else
            x(i-num_shifted) = x(i)
            y(i-num_shifted) = y(i)
            px(i-num_shifted) = px(i)
            py(i-num_shifted) = py(i)
            r(i-num_shifted) = r(i)
            m(i-num_shifted) = m(i)
        end if
    end do

    num_particles = num_particles - num_shifted

#ifdef USE_GPU
    !$omp target update to(x, y, px, py, r, m, merged)
#endif

end subroutine coallese_particles

end module save_data