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
        write(*,100) "=====Max Running processes:", nProc
    endif
    
    call MPI_GET_PROCESSOR_NAME(proc_name, namelen, ierr)
    write(*,*) "Greetings from process #", myID," on ", proc_name 
    
#ifdef _OPENMP
    write(*,*) "Starting OpenMP >>>>>>"
    call OMP_set_num_threads(8)
    myMaxThreads = OMP_get_max_threads()
    write(*,100) "-----Max Running threads:",myMaxThreads
    
    !$omp parallel default(none) private(threadID)
    threadID = omp_get_thread_num()
    write(*,100) "Greetings from thread #", threadID
    !$omp end parallel
    write(*,*) "    "
#endif
100 format(1X,A,I3)

    call trapezoidalRule()
    
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