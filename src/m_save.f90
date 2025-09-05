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

end module save_data