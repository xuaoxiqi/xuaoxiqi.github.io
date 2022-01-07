
    module commondata
        use mpi
        implicit none
        integer :: myID, nProc, ierr, idest
        integer :: istatus(MPI_STATUS_SIZE)
        real(8), parameter :: Pi=4.0d0*datan(1.0d0)
        integer, parameter :: nx=50000001
        real(8), parameter :: dx=1.0d0/dble(nx-1)
        real(8) :: X(nx)
    end module commondata
    
    
    program main
    use commondata
    use mpi
    use omp_lib
    implicit none
    integer :: namelen
    character(len=20) :: proc_name
    integer :: myMaxThreads, threadID
    
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nProc, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myID, ierr)
    
    if (myID.EQ.0) then
        write(*,*) "Max Running processes:", nProc
        write(*,*) " "
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)    

    call MPI_GET_PROCESSOR_NAME(proc_name, namelen, ierr)
    write(*,*) "From process #", myID," on ", proc_name
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if(myID.EQ.0) then
        write(*,*) " "
    endif

#ifdef _OPENMP
    call OMP_set_num_threads(8)
    myMaxThreads = OMP_get_max_threads()
    if(myID.EQ.0) then
        write(*,*) "Max Running threads:",myMaxThreads
        write(*,*) " "
    endif    

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    !$omp parallel default(none) shared(myID, proc_name) private(threadID)
    threadID = omp_get_thread_num()
    write(*,*) "From process #", myID, "; Thread #", threadID
    !$omp end parallel
#endif

    call trapezoidalRule()

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)    
    if(myID.EQ.0) then
        write(*,*) "Successfully: simulation completed!"
    endif
    
    call MPI_FINALIZE(ierr)
    
    stop
    end program main
    
    
    subroutine trapezoidalRule()
    use commondata
    use mpi
    implicit none
    integer :: i
    real(8) :: myPi, myPiTemp
    
    myPiTemp = 0.0d0
    !$omp parallel do default(none) shared(X,nProc) private(i,myID) reduction(+:myPiTemp)
    do i=myID+1,nx,nProc
        X(i) = (dble(i)-0.5d0)*dx
        myPiTemp = myPiTemp+4.0d0/(X(i)*X(i)+1.0d0)*dx
    enddo
    !$omp end parallel do
    
    call MPI_ALLREDUCE(myPiTemp, myPi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    if(myID.EQ.0) then
        write(*,*) "Pi (calculation) =", myPi
        write(*,*) "Pi (machine)     =", Pi
    endif
    
    return
    end subroutine trapezoidalRule
null
