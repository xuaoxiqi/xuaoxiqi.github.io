
    module commondata
        use mpi
        implicit none
        real(kind=8), parameter :: Pi=4.0d0*datan(1.0d0)
        integer(kind=4), parameter :: nx=6
        real(kind=8), parameter :: dx=1.0d0/dble(nx-1)
        real(kind=8), allocatable :: Xcenter(:)
        
        !----for MPI----!
        integer(kind=4) :: myID, nProc, ierr, iRoot
        integer(kind=4), allocatable :: start1d(:), end1d(:), count1d(:)
        integer(kind=4) :: istatus(MPI_STATUS_SIZE)
        integer(kind=4) :: nxLocal
        integer(kind=4) :: iStart, iEnd
    end module commondata
    
    
    program main
    use commondata
    use mpi
    use omp_lib
    implicit none
    integer(kind=4) :: myMaxThreads, threadID
    
    !----------------
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nProc, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myID, ierr)
    
    if (myID.EQ.0) then
        write(*,*) "Grids: From 1 to nx =", nx
        write(*,*) "Integration range: [0, 1]"
        write(*,*) "    "
    
        write(*,*) "Max Running processes:", nProc
        write(*,*) " "
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)    

#ifdef _OPENMP
    call OMP_set_num_threads(8)
    myMaxThreads = OMP_get_max_threads()
    if(myID.EQ.0) then
        write(*,*) "Max Running threads:",myMaxThreads
        write(*,*) " "
    endif    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    !$omp parallel default(none) shared(myID) private(threadID)
    threadID = omp_get_thread_num()
    write(*,*) "From process #", myID, "; Thread #", threadID
    !$omp end parallel
#endif
    
    allocate (start1d(0:nProc-1))
    allocate (end1d(0:nProc-1))
    allocate (count1d(0:nProc-1))
    
    nxLocal = nx/nProc+1
    allocate (Xcenter(1:nxLocal))

    !--------------------------------------------------------------------    
    call StartEnd(1, nx-1)

    iStart = 1
    iEnd = count1d(myID)
    !--------------------------------------------------------------------    
    
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
    integer(kind=4) :: i
    real(kind=8) :: myPi, myPiTemp
    
    myPiTemp = 0.0d0
    !$omp parallel do default(none) shared(Xcenter, myID, start1d, iStart, iEnd) private(i) reduction(+:myPiTemp)
    do i = iStart, iEnd
        Xcenter(i) = (dble(start1d(myID)-1+i)-0.5d0)*dx
        !~ write(*, *) "Xcenter(i):", real(Xcenter(i))
        myPiTemp = myPiTemp+4.0d0/(Xcenter(i)*Xcenter(i)+1.0d0)*dx
    enddo
    !$omp end parallel do
    write(*,*) "From process #", myID
    write(*,*) "    "
    
    iRoot = 0
    call MPI_REDUCE(myPiTemp, myPi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, iRoot, MPI_COMM_WORLD, ierr)
    
    if(myID.EQ.0) then
        write(*,*) "Pi (calculation) =", myPi
        write(*,*) "Pi (machine)     =", Pi
    endif
    
    return
    end subroutine trapezoidalRule


    subroutine StartEnd(iS1, iS2)
    use commondata
    implicit none
    integer(kind=4) :: leng, iBlock
    integer(kind=4) :: ir
    integer(kind=4) :: iS1, iS2
    integer(kind=4) :: i
    
    leng = iS2-iS1+1
    iBlock = leng/nProc
    ir= leng-iBlock*nProc
    
    do i = 0, nProc-1
        
        if(i.LT.ir) then  
            count1d(i) = iBlock+1
            start1d(i) = iS1+i*(iBlock+1) 
            end1d(i) = start1d(i)+count1d(i)-1

        else
            count1d(i) = iBlock   
            start1d(i) = iS1+i*iBlock+ir   
            end1d(i) = start1d(i)+count1d(i)-1  
        endif

    enddo
    
    return
    end subroutine StartEnd