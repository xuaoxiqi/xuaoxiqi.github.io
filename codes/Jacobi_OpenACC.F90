
    module commondata
        implicit none
        integer(kind=4), parameter :: nx=401, ny=401
        real(kind=8), allocatable :: Xp(:), Yp(:)
        real(kind=8), allocatable :: A(:,:), Anew(:,:)
        real(kind=8) :: error
        integer(kind=4) :: itc
        real(kind=8), parameter :: eps=1e-6
        integer(kind=4), parameter :: itc_max=200000
    end module commondata


    program main
    use commondata
    implicit none
    integer :: i, j
    real(kind=8) :: timeStart, timeEnd
    INTEGER(kind=4) :: time
    character(len=24) :: ctime, string
    
    allocate (Xp(nx))
    allocate (Yp(ny))
    allocate (A(nx, ny))
    allocate (Anew(nx, ny))

    do i = 1, nx
        Xp(i) = dble(i-1)/dble(nx-1)
    enddo
    do j = 1, ny
        Yp(j) = dble(j-1)/dble(ny-1)
    enddo

    A = 0.0d0
    Anew = 0.0d0

    do i=1,nx
        A(i,1) = 0.0d0
        A(i,ny) = 1.0d0
    enddo
    do j=1,ny
        A(1,j) = 0.0d0
        A(nx,j) = 0.0d0
    enddo

    error = 10.0d0
    itc = 0
    
    string = ctime( time() )
    write(*,*) 'Start: ', string
    
    call CPU_TIME(timeStart)
    
    !$acc data copy(A) create(Anew)
    do while( (error.GT.eps).AND.(itc.LE.itc_max) )

        itc = itc+1
        
        call calA()

        if(MOD(itc,2000).EQ.0) then
            call errorA()
        endif

        call updateA()

    enddo
    !$acc end data

    call CPU_TIME(timeEnd)
    
    write(00,*) "Time = ", real(timeEnd-timeStart), "s"

    string = ctime( time() )
    write(00,*) 'End:   ', string
    
    call outputA_datFile()

    call outputA_pltFile()

    write(*,*) "Successfully: simulation completed!"
    
    stop
    end program main
    

    subroutine calA()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$acc parallel loop present(A, Anew)
    do j = 2, ny-1
        !$acc loop
        do i = 2, nx-1
            Anew(i,j) = 0.25d0*( A(i+1,j)+A(i-1,j)+A(i,j+1)+A(i,j-1) )
        enddo
    enddo
    !$acc end parallel
    
    return
    end subroutine calA


    subroutine updateA()
    use commondata
    implicit none
    integer(kind=4) :: i, j 

    !$acc parallel loop present(A, Anew)
    do j = 2, ny-1
        !$acc loop
        do i = 2, nx-1
            A(i,j) = Anew(i,j)
        enddo
    enddo
    !$acc end parallel

    return
    end subroutine updateA


    subroutine errorA()
    use commondata
    implicit none
    integer(kind=4) :: i, j 

    error = 0
    !$acc parallel loop reduction(+:error) present(A, Anew)
    do j = 2, ny-1
        !$acc loop reduction(+:error)
        do i = 2, nx-1
            error = error+abs(Anew(i,j)-A(i,j))
        enddo
    enddo
    !$acc end parallel

    write(*,*) itc, error

    return
    end subroutine errorA


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
            write(02,100) Xp(i), Yp(j), A(i,j)
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
                write(41) real(A(i,j))
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