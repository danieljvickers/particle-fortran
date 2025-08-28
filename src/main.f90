program main
    use mpi
    use omp_lib
    implicit none

    integer :: ierr, rank, size
    integer :: omp_rank, omp_size

    ! Initialize MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    ! OpenMP parallel region
    !$omp parallel private(omp_rank, omp_size)
        omp_rank = omp_get_thread_num()
        omp_size = omp_get_num_threads()

        print *, 'Hello from MPI rank', rank, 'of', size, &
                 'and OpenMP thread', omp_rank, 'of', omp_size
    !$omp end parallel

    ! Finalize MPI
    call MPI_Finalize(ierr)

end program main
