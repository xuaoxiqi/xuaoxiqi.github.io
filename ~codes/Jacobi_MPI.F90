
    module commondata
        use mpi
        implicit none
        integer, parameter :: nx=800, ny=800
        real(8) :: A(nx,ny), Anew(nx,ny)
        real(8) :: errorMax, maxGlobal
        integer :: itc
        real(8), parameter :: eps=1e-6
        integer, parameter :: itc_max=20000
        
        integer :: myID, nProc
        integer :: iStart, iEnd
        integer :: iStartMinus1, iEndPlus1
        integer :: iStartLocal, iEndLocal
        integer, parameter :: maxProc=48
        integer :: gStart(0:maxProc-1), gEnd(0:maxProc-1), gCount(0:maxProc-1)
        integer :: ierr, iRoot
        integer :: iStatus(MPI_STATUS_SIZE)
        integer :: l_nbr, r_nbr
    end module commondata


    program main
    use commondata
    use mpi
    implicit none
    integer :: i, j
    real(8) :: xp(nx), yp(ny)
    real(8) :: start, finish
    real(8) :: start2, finish2
    real(8) :: clock
    integer :: myCount
    integer :: kount, iDest
    integer :: iStart1, kount1
    integer :: iSrc

    !----------------
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nProc, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myID, ierr)

    if(myID.EQ.0) then
        write(*,*) 'nProc =',  nProc
    endif
    !----------------
    
    call startEnd(1, ny) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    iStart = gStart(myID)
    iEnd = gEnd(myID)
    myCount = gCount(myID)
    write(*,*) "myID, iStart, iEnd=", myID, iStart, iEnd
    !----------------

    iEndLocal = iEnd
    if(myID.EQ.nProc-1) iEndLocal = iEnd-1
    
    iStartLocal = iStart
    if(myID.EQ.0) iStartLocal = 2
    
    iStartMinus1 = iStart-1
    iEndPlus1 = iEnd+1
    
    l_nbr = myID-1
    r_nbr = myID+1
    if(myID.EQ.0) l_nbr = MPI_PROC_NULL
    if(myID.EQ.nProc-1) r_nbr = MPI_PROC_NULL
    
    if(myID.EQ.0) then
        do i=1,nx
            Xp(i) = (i-1)
        enddo
        do j=1,ny
            Yp(j) = (j-1)
        enddo

        A = 0.0d0
        Anew = 0.0d0

        do i=1,nx
            A(i,1) = 0.0d0
            A(i,ny) = 100.0d0
        enddo
        do j=1,ny
            A(1,j) = 0.0d0
            A(nx,j) = 0.0d0
        enddo

    endif
    
    iRoot = 0
    call MPI_BCAST(A, nx*ny, MPI_DOUBLE_PRECISION, iRoot, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(Anew, nx*ny, MPI_DOUBLE_PRECISION, iRoot, MPI_COMM_WORLD, ierr)

    errorMax = 10.0d0
    
    itc = 0
        
    do while( (itc.LE.itc_max) )

        call calA()

        if(MOD(itc,2000).EQ.0) then
            call errorA()
        endif

        call updateA()

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                    
        itc = itc+1
        
    enddo

    !collect data
    if(myID.NE.0) then
        kount = gCount(myID)
        iDest = 0
        call MPI_SEND(A(1,iStart), kount*nx, MPI_DOUBLE_PRECISION, iDest, 110, MPI_COMM_WORLD, ierr)
    else
        do iSrc=1,nProc-1
            iStart1 = gStart(iSrc)
            kount1 = gCount(iSrc)
            call MPI_RECV(A(1,iStart1), kount1*nx, MPI_DOUBLE_PRECISION, iSrc, 110, MPI_COMM_WORLD, iStatus, ierr)
        enddo
    endif
    
    if(myID.EQ.0) then
        open(unit=02,file='Jacobi.dat',status='unknown')
        write(02,*) 'TITLE="Laplace Solver"'
        write(02,*) 'VARIABLES="X" "Y" "Numerical" '
        write(02,101) nx, ny
        do j=1,ny
            do i=1,nx
                write(02,100) xp(i), yp(j), A(i,j)
            enddo
        enddo
    100 format(1x,2(e11.4,' '),10(e13.6,' '))
    101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
        close(02)
    endif
    
    write(*,*) "myID =",myID

    call MPI_FINALIZE(ierr)

    end program main


    subroutine calA()
    use commondata
    implicit none
    integer :: i, j

    call MPI_SENDRECV(A(1,iEnd), nx, MPI_DOUBLE_PRECISION, r_nbr, 20, &
                                A(1,iStartMinus1), nx, MPI_DOUBLE_PRECISION, l_nbr, 20, &
                                MPI_COMM_WORLD, iStatus, ierr)
    call MPI_SENDRECV(A(1,iStart), nx, MPI_DOUBLE_PRECISION, l_nbr, 20, &
                                A(1,iEndPlus1), nx, MPI_DOUBLE_PRECISION, r_nbr, 20, &
                                MPI_COMM_WORLD, iStatus, ierr)
                        
    !!do j=2,ny-1
    do j=iStartLocal,iEndLocal
        do i=2,nx-1
            Anew(i,j) = 0.25d0*( A(i+1,j)+A(i-1,j)+A(i,j+1)+A(i,j-1) )
        enddo
    enddo

    return
    end subroutine calA


    subroutine updateA()
    use commondata
    implicit none
    integer :: i, j

    !!do j=2,ny-1
    do j=iStartLocal,iEndLocal
        do i=2,nx-1
            A(i,j) = Anew(i,j)
        enddo
    enddo

    return
    end subroutine updateA


    subroutine errorA
    use commondata
    implicit none
    integer :: i, j

    errorMax = 0
    !!do j=2,ny-1
    do j=iStartLocal,iEndLocal
        do i=2,nx-1
            errorMax = max( errorMax, abs(Anew(i,j)-A(i,j)) )
        enddo
    enddo

    call MPI_ALLREDUCE(errorMax, maxGlobal, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    
    if(myID.EQ.0) then
        write(*,*) itc, maxGlobal
    endif
    
    return
    end subroutine errorA
    
    
    subroutine StartEnd(iS1, iS2)
    use commondata
    implicit none
    integer :: leng, iBlock
    integer :: ir
    integer :: id, iS1, iS2
    integer :: i
    
    leng = iS2-iS1+1
    iBlock = leng/nProc
    ir= leng-iBlock*nProc
    do i=0,nProc-1
        if(i.LT.ir) then
            gStart(i) = iS1+i*iBlock+i
            gEnd(i) = gStart(i)+iBlock
        else
            gStart(i) = iS1+i*iBlock+ir
            gEnd(i) = gStart(i)+iBlock-1
        endif
        
        if(leng.LT.1) then
            gStart(i) = 1
            gEnd(i) = 0
            write(*,*) "Error!! Robust test......"
            stop
        endif
        
        gCount(i) = gEnd(i)-gStart(i)+1
    enddo
    
    return
    end subroutine StartEnd
