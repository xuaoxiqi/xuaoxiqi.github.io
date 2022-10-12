
    module commondata
        use mpi
        implicit none
        integer(kind=4), parameter :: nx=401, ny=401
        real(kind=8), allocatable :: Xp(:), Yp(:)
        real(kind=8), allocatable :: A(:,:), Anew(:,:)
        real(kind=8), allocatable :: AGlobal(:,:)
        real(kind=8) :: errorGloabl
        integer(kind=4) :: itc
        real(kind=8), parameter :: eps=1e-6
        integer(kind=4), parameter :: itc_max=200000
        
        !----for MPI----!
        integer(kind=4) :: myID, nProc, ierr, iRoot
        integer(kind=4), allocatable :: start1d(:), end1d(:), count1d(:), displ1d(:)
        integer(kind=4), allocatable :: start2d(:), end2d(:), count2d(:), displ2d(:)
        integer(kind=4) :: iStatus(MPI_STATUS_SIZE)
        integer(kind=4) :: iStart, iEnd
        integer(kind=4) :: iStartMinus1, iEndPlus1
        integer(kind=4) :: iStartLocal, iEndLocal
        integer(kind=4) :: leftNeighbor, rightNeighbor
        integer(kind=4) :: nyLocal
    end module commondata


    program main
    use commondata
    use mpi
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: timeStart, timeEnd
    character(len=24) :: ctime, string
    integer(kind=4) :: time
    integer(kind=4) :: iSrc, iDest

    !----------------
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nProc, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myID, ierr)
    
    !----------------
    if (myID.EQ.0) then
        string = ctime( time() )
        write(*,*) 'Start: ', string
        write(*,*) "Max Running processes:", nProc
        write(*,*) " "
    endif
    !----------------
    
    allocate (start1d(0:nProc-1))
    allocate (end1d(0:nProc-1))
    allocate (count1d(0:nProc-1))
    allocate (displ1d(0:nProc-1))
    
    allocate (start2d(0:nProc-1))
    allocate (end2d(0:nProc-1))
    allocate (count2d(0:nProc-1))
    allocate (displ2d(0:nProc-1))
    
    nyLocal = ny/nProc+1
    allocate (Xp(nx))
    allocate (Yp(ny))
    
    allocate (A(nx, 0:nyLocal+1))
    allocate (Anew(nx, nyLocal))
    allocate (AGlobal(nx, ny))
    !--------------------------------------------------------------------    
    call StartEnd(1, ny)

    iStart = 1
    iEnd = count1d(myID)
    !----------------
    
    iStartLocal = iStart
    if(myID.EQ.0) iStartLocal = 2
    
    iEndLocal = iEnd
    if(myID.EQ.nProc-1) iEndLocal = iEnd-1
    
    !================================================
    iStartMinus1 = iStart-1
    iEndPlus1 = iEnd+1
    
    leftNeighbor = myID-1
    rightNeighbor = myID+1
    if(myID.EQ.0) leftNeighbor = MPI_PROC_NULL
    if(myID.EQ.nProc-1) rightNeighbor = MPI_PROC_NULL
    !================================================
    
    if(myID.EQ.0) then
        do i = 1, nx
            Xp(i) = dble(i-1)/dble(nx-1)
        enddo
        do j = 1, ny
            Yp(j) = dble(j-1)/dble(ny-1)
        enddo

        AGlobal = 0.0d0
        do i=1,nx
            AGlobal(i,1) = 0.0d0
            AGlobal(i,ny) = 1.0d0
        enddo
        do j=1,ny
            AGlobal(1,j) = 0.0d0
            AGlobal(nx,j) = 0.0d0
        enddo
    endif
    
    errorGloabl = 10.0d0
    itc = 0
    
    A = 0.0d0
    Anew = 0.0d0
        
    call MPI_SCATTERV(AGlobal, count2d, displ2d, MPI_DOUBLE_PRECISION,  &
                                    A(1,1), count2d(myID), MPI_DOUBLE_PRECISION, iRoot, MPI_COMM_WORLD, ierr)

    !----------------------------------------------------------
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    timeStart = MPI_WTIME()
    !----------------------------------------------------------
    do while( (errorGloabl.GT.eps).AND.(itc.LE.itc_max) )
    
        itc = itc+1
        
        call calA()

        if(MOD(itc,2000).EQ.0) then
            call errorA()
        endif

        call updateA()

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    enddo
    !----------------------------------------------------------
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    timeEnd = MPI_WTIME()
    !----------------------------------------------------------
    
    iDest = 0
    call MPI_GATHERV(A(1,1), count2d(myID), MPI_DOUBLE_PRECISION, &
                           AGlobal, count2d, displ2d, MPI_DOUBLE_PRECISION, iDest, MPI_COMM_WORLD, ierr)
     
    if(myID.EQ.0) then
    
        write(*,*) "Successfully: DNS completed!"
        string = ctime( time() )
        write(*,*) 'End:   ', string
    
        write(*,*) "Time (CPU) =", real(timeEnd-timeStart),"seconds"
    
        call outputA_datFile()

        call outputA_pltFile()
    endif
    
    deallocate (AGlobal)
    deallocate (A)
    deallocate (Anew)
    deallocate (Xp)
    deallocate (Yp)

    call MPI_FINALIZE(ierr)

    stop
    end program main


    subroutine calA()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    call MPI_SENDRECV(A(1,iEnd), nx, MPI_DOUBLE_PRECISION, rightNeighbor, 20, &
                                A(1,iStartMinus1), nx, MPI_DOUBLE_PRECISION, leftNeighbor, 20, &
                                MPI_COMM_WORLD, iStatus, ierr)
    call MPI_SENDRECV(A(1,iStart), nx, MPI_DOUBLE_PRECISION, leftNeighbor, 20, &
                                A(1,iEndPlus1), nx, MPI_DOUBLE_PRECISION, rightNeighbor, 20, &
                                MPI_COMM_WORLD, iStatus, ierr)
                        
    do j=iStartLocal,iEndLocal
        do i = 2, nx-1
            Anew(i,j) = 0.25d0*( A(i+1,j)+A(i-1,j)+A(i,j+1)+A(i,j-1) )
        enddo
    enddo

    return
    end subroutine calA


    subroutine updateA()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    do j = iStartLocal, iEndLocal
        do i = 2, nx-1
            A(i,j) = Anew(i,j)
        enddo
    enddo

    return
    end subroutine updateA


    subroutine errorA
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: errorLocal

    errorLocal = 0.0d0
    do j=iStartLocal,iEndLocal
        do i=2,nx-1
            errorLocal = errorLocal+dabs(Anew(i,j)-A(i,j))
        enddo
    enddo

    call MPI_ALLREDUCE(errorLocal, errorGloabl, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    
    if(myID.EQ.0) then
        write(*,*) itc, errorGloabl
    endif
    
    return
    end subroutine errorA
    
    
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
            !-----------------------------------------------------------
            count2d(i) = (iBlock+1)*nx
            start2d(i) = iS1+i*(iBlock+1)*nx
            end2d(i) = start2d(i)+count2d(i)-1

        else
            count1d(i) = iBlock   
            start1d(i) = iS1+i*iBlock+ir   
            end1d(i) = start1d(i)+count1d(i)-1  
            !-----------------------------------------------------------
            count2d(i) = iBlock*nx 
            start2d(i) = iS1+i*iBlock*nx+ir*nx  
            end2d(i) = start2d(i)+count2d(i)-1  
        endif

        displ1d(i) = start1d(i)-iS1  
        displ2d(i) = start2d(i)-iS1  
        
    enddo
    
    return
    end subroutine StartEnd
    
    
    subroutine outputA_datFile()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    
    open(unit=02,file='Jacobi.dat',status='unknown')
    write(02,*) 'TITLE="Laplace Solver"'
    write(02,*) 'VARIABLES="X" "Y" "Numerical" '
    write(02,101) nx, ny
    do j = 1, ny
        do i = 1, nx
            write(02,100) Xp(i), Yp(j), AGlobal(i,j)
        enddo
    enddo
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(02)
    
    return
    end subroutine outputA_datFile


    subroutine outputA_pltFile()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    character(len=9) :: B2
    REAL(kind=4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3
    integer(kind=4), parameter :: kmax=1
    character(len=40) :: zoneName
    
    open(unit=41,file='Jacobi.plt', access='stream', form='unformatted')

    !---------------------------------------------
    zoneMarker= 299.0
    eohMarker = 357.0

    !I. HEAD SECTION--------------------------------------
    !c--Magic number, Version number
    write(41) '#!TDV101'

    !c--Integer value of 1
    write(41) 1

    Title='MyFirst'
    call dumpstring(title)

    !c-- Number of variables in this data file
    write(41) 3

    !c-- Variable names.
    V1='X'
    call dumpstring(V1)
    V2='Y'
    call dumpstring(V2)
    V3='Numerical'
    call dumpstring(V3)

    !c-----Zones-----------------------------

    !c--------Zone marker. Value = 299.0
    write(41) zoneMarker

    !--------Zone name.
    zoneName='ZONE 001'
    call dumpstring(zoneName)

    !---------Zone Color
    write(41) -1

    !---------ZoneType
    write(41) 0

    !---------DataPacking 0=Block, 1=Point
    write(41) 1

    !---------Specify Var Location. 0 = Do not specify, all data
    !---------is located at the nodes. 1 = Specify
    write(41) 0

    !---------Number of user defined face neighbor connections
    ! (value >= 0)
    write(41) 0

    !---------IMax,JMax,KMax
    write(41) nx
    write(41) ny
    write(41) kmax

    !-----------1=Auxiliary name/value pair to follow
    !-----------0=No more Auxiliar name/value pairs.
    write(41) 0
    write(41) eohMarker

    !----zone ------------------------------------------------------------
    write(41) zoneMarker

    !--------variable data format, 1=Float, 2=Double, 3=LongInt,4=ShortInt, 5=Byte, 6=Bit
    write(41) 1
    write(41) 1
    write(41) 1

    !--------Has variable sharing 0 = no, 1 = yes.
    write(41) 0

    !----------Zone number to share connectivity list with (-1 = no
    ! sharing).
    write(41) -1

    !---------------------------------------------------------------------
    do k=1,kmax
        do j=1,ny
            do i=1,nx
                write(41) real(Xp(i))
                write(41) real(Yp(j))
                write(41) real(AGlobal(i,j))
            end do
        end do
    enddo
    close(41)
    !---------------------------------------------------------------------
    
    return
    end subroutine outputA_pltFile
    
    subroutine dumpstring(instring)
    implicit none
    character(len=40) instring
    integer(kind=4) :: stringLength
    integer(kind=4) :: ii
    integer(kind=4) :: I

    stringLength=LEN_TRIM(instring)
    do ii=1,stringLength
        I=ICHAR(instring(ii:ii))
        write(41) I
    end do
    write(41) 0

    return
    end subroutine dumpstring