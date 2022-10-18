
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
        integer(kind=4) :: myID, nProc, ierr, iRoot, iSrc, iDest
        integer(kind=4) :: iStatus(MPI_STATUS_SIZE)
        integer(kind=4) :: iStartX, iEndX, iStartY, iEndY
        integer(kind=4) :: nxLocal, nyLocal
        integer(kind=4), allocatable :: countX(:), countY(:)
        integer(kind=4), allocatable :: i_start_global(:), j_start_global(:)
        
        integer(kind=4) :: dims(1:2), coords(1:2)
        integer(kind=4) :: MPI_COMM_2D, myID_2D, column_y
        integer(Kind=4) :: nbr_left, nbr_right, nbr_top, nbr_bottom
        logical :: periods(2)
        data periods/2*.FALSE./
    end module commondata


    program main
    use commondata
    use mpi
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: timeStart, timeEnd
    character(len=24) :: ctime, string
    integer(kind=4) :: time
    real(kind=8), allocatable :: tempA(:,:)

    !----------------
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nProc, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myID, ierr)
    
    !----------------
    if (myID.EQ.0) then
        string = ctime( time() )
        write(*,*) 'Start: ', string
        write(*,*) "Max Running processes:", nProc
    endif
    !----------------
    
    dims = 0
    call MPI_DIMS_CREATE (nProc, 2, dims, ierr)
    call MPI_CART_CREATE (MPI_COMM_WORLD, 2, dims, periods, .TRUE., MPI_COMM_2D, ierr)
    if (myID.EQ.0) then
        write(*,*) "CPU in each dimension is: ", dims(1), " &", dims(2)
    endif
        
    ! get myID in 2D 
    call MPI_COMM_RANK (MPI_COMM_2D, myID_2D, ierr)
    
    ! get the neighbors
    call MPI_CART_SHIFT (MPI_COMM_2D, 0, 1, nbr_left, nbr_right, ierr)
    call MPI_CART_SHIFT (MPI_COMM_2D, 1, 1, nbr_bottom, nbr_top, ierr)
    
    ! get the subdomain information
    call MPI_CART_GET (MPI_COMM_2D, 2, dims, periods, coords, ierr)
    
    allocate (countX(0: nProc-1))
    allocate (countY(0: nProc-1))
    countX = 0
    countY = 0
    
    allocate (i_start_global(0: nProc-1))
    allocate (j_start_global(0: nProc-1))
    i_start_global = 0
    j_start_global = 0
    
    ! determine sub-domain size in each dimension    
    !------------------------------------------------------------------------------------------------
    call StartEnd(1, nx, 1, ny, dims, coords)
    !------------------------------------------------------------------------------------------------
    nxLocal = countX(myID_2D)
    nyLocal = countY(myID_2D)
    
    iStartX = 1
    if(coords(1).EQ.0) iStartX = 2
    iEndX = nxLocal
    if(coords(1).EQ.dims(1)-1) iEndX = nxLocal-1
    
    iStartY = 1
    if(coords(2).EQ.0) iStartY = 2
    iEndY = nyLocal
    if(coords(2).EQ.dims(2)-1) iEndY = nyLocal-1

    allocate (Xp(1:nx))
    allocate (Yp(1:ny))
    
    ! allocate array with ghost layers outside
    allocate (A(0:nxLocal+1, 0:nyLocal+1))
    allocate (Anew(0:nxLocal+1, 0:nyLocal+1))
    
    ! construct the datatype for exchange in the y-direction (non-contiguous)
    call MPI_TYPE_VECTOR(nyLocal, 1, nxLocal+2, MPI_DOUBLE_PRECISION, column_y, ierr)
    call MPI_TYPE_COMMIT(column_y, ierr)
    !--------------------------------------------------------------------    

    do i = 1, nx
        Xp(i) = dble(i-1)/dble(nx-1)
    enddo
    do j = 1, ny
        Yp(j) = dble(j-1)/dble(ny-1)
    enddo
    
    errorGloabl = 10.0d0
    itc = 0
    
    A = 0.0d0
    Anew = 0.0d0
    ! top boundary
    if (coords(2).EQ. dims(2)-1) then
        do i=0, nxLocal+1
            A(i, nyLocal) = 1.0d0
        enddo
    endif
    
    !----------------------------------------------------------
    call MPI_BARRIER(MPI_COMM_2D, ierr)
    timeStart = MPI_WTIME()
    !----------------------------------------------------------
    do while( (errorGloabl.GT.eps).AND.(itc.LE.itc_max) )
    
        itc = itc+1
        
        call calA()

        if(MOD(itc,2000).EQ.0) then
            call errorA()
        endif

        call updateA()

        call MPI_BARRIER(MPI_COMM_2D, ierr)

    enddo
    !----------------------------------------------------------
    call MPI_BARRIER(MPI_COMM_2D, ierr)
    timeEnd = MPI_WTIME()
    !----------------------------------------------------------
    
    if(myID.EQ.0) then
        allocate (AGlobal(nx, ny))
        
        AGlobal = 0.0d0
        do j = 1, nyLocal
            do i = 1, nxLocal 
                AGlobal(i+i_start_global(0)-1, j+j_start_global(0)-1) = A(i, j)
            enddo
        enddo
        
        do iSrc = 1, nProc-1

            allocate (tempA(0:countX(iSrc)+1, 0:countY(iSrc)+1))
            tempA = 0.0d0
            
            call MPI_RECV(tempA, (countX(iSrc)+2)*(countY(iSrc)+2), MPI_DOUBLE_PRECISION, iSrc, 110, MPI_COMM_2D, iStatus, ierr)

            do j = 1, countY(iSrc)
                do i = 1, countX(iSrc)
                    AGlobal(i+i_start_global(iSrc)-1, j+j_start_global(iSrc)-1) = tempA(i, j)
                enddo
            enddo
            
            deallocate(tempA)
        enddo
        
    else
        iDest = 0
        call MPI_SEND(A, (countX(myID_2D)+2)*(countY(myID_2D)+2), MPI_DOUBLE_PRECISION, iDest, 110, MPI_COMM_2D, ierr)
    endif
    
    if(myID.EQ.0) then
    
        write(*,*) "Successfully: DNS completed!"
        string = ctime( time() )
        write(*,*) 'End:   ', string
    
        write(*,*) "Time (CPU) =", real(timeEnd-timeStart),"seconds"
    
        call outputA_datFile()

        call outputA_pltFile()

        deallocate (AGlobal)
    endif
        
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
                                
    call MPI_SENDRECV(A(1, nyLocal), nxLocal, MPI_DOUBLE_PRECISION, nbr_top, 0, &
                                    A(1, 0), nxLocal, MPI_DOUBLE_PRECISION, nbr_bottom, 0, &
                                    MPI_COMM_2D, iStatus, ierr)
    call MPI_SENDRECV(A(1, 1), nxLocal, MPI_DOUBLE_PRECISION, nbr_bottom, 1, &
                                    A(1, nyLocal+1), nxLocal, MPI_DOUBLE_PRECISION, nbr_top, 1, &
                                    MPI_COMM_2D, iStatus, ierr)
    call MPI_SENDRECV(A(nxLocal, 1), 1, column_y, nbr_right, 0, &
                                    A(0, 1), 1, column_y, nbr_left, 0, &
                                    MPI_COMM_2D, iStatus, ierr)
    call MPI_SENDRECV(A(1, 1), 1, column_y, nbr_left, 1, &
                                    A(nxLocal+1, 1), 1, column_y, nbr_right, 1, &
                                    MPI_COMM_2D, iStatus, ierr)
                        
    do j = iStartY, iEndY
        do i = iStartX, iEndX
            Anew(i,j) = 0.25d0*( A(i+1,j)+A(i-1,j)+A(i,j+1)+A(i,j-1) )
        enddo
    enddo

    return
    end subroutine calA


    subroutine updateA()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    do j = iStartY, iEndY
        do i = iStartX, iEndX
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
    do j = iStartY, iEndY
        do i = iStartX, iEndX
            errorLocal = errorLocal+dabs(Anew(i,j)-A(i,j))
        enddo
    enddo

    call MPI_ALLREDUCE(errorLocal, errorGloabl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_2D, ierr)
    
    if(myID.EQ.0) then
        write(*,*) itc, errorGloabl
    endif
    
    return
    end subroutine errorA


    subroutine StartEnd(iSx1, iSx2, iSy1, iSy2, num_process, rank)
    use commondata
    implicit none
    integer(kind=4) :: num_process(1:2), rank(1:2), iBlock(1:2)
    integer(kind=4) :: leng(1:2), ir(1:2)
    integer(kind=4) :: iSx1, iSx2, iSy1, iSy2
    integer(kind=4) :: startX, startY, localCountX, localCountY
    integer(kind=4) :: i, j
    integer(kind=4) :: startPoint(1:2), endPoint(1:2)
    integer(kind=4) :: dir
    
    startPoint(1) = iSx1
    endPoint(1) = iSx2
    
    startPoint(2) = iSy1
    endPoint(2) = iSy2
    
    do dir = 1, 2
        leng(dir) = endPoint(dir)-startPoint(dir)+1
        iBlock(dir) = leng(dir)/dims(dir)
        ir(dir) = leng(dir)-iBlock(dir)*dims(dir)
    enddo
    
    if(rank(1).LT.ir(1)) then
        localCountX = iBlock(1)+1
        startX = startPoint(1)+rank(1)*(iBlock(1)+1) 
    elseif(rank(1).GE.ir(1)) then
        localCountX = iBlock(1)
        startX = startPoint(1)+rank(1)*iBlock(1)+ir(1)
    endif

    if(rank(2).LT.ir(2)) then
        localCountY = iBlock(2)+1
        startY = startPoint(2)+rank(2)*(iBlock(2)+1) 
    elseif(rank(2).GE.ir(2)) then
        localCountY = iBlock(2)
        startY = startPoint(2)+rank(2)*iBlock(2)+ir(2)
    endif
        
    call MPI_ALLGATHER(localCountX, 1, MPI_INTEGER, countX, 1, MPI_INTEGER, MPI_COMM_2D, ierr)
    call MPI_ALLGATHER(localCountY, 1, MPI_INTEGER, countY, 1, MPI_INTEGER, MPI_COMM_2D, ierr)
    
    call MPI_ALLGATHER(startX, 1, MPI_INTEGER, i_start_global, 1, MPI_INTEGER, MPI_COMM_2D, ierr)
    call MPI_ALLGATHER(startY, 1, MPI_INTEGER, j_start_global, 1, MPI_INTEGER, MPI_COMM_2D, ierr)

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