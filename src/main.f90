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
    implicit none


    ! contains the variables that we will use to intialize
    real(8) :: escape_factor
    
    escape_factor = 40  ! the fractional distance from the star, relative to radius_upper, that the asteroid is considered "escaped" 

    ! initialize particles
    call initialize_particles()

    ! take time steps

end program main