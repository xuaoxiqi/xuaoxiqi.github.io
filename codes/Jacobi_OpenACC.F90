
    module commondata
        implicit none
        integer, parameter :: nx=2000, ny=2000
        real(8) :: A(nx,ny), Anew(nx,ny)
        real(8) :: error
        integer :: itc
        real(8), parameter :: eps=1e-6
        integer, parameter :: itc_max=20000
    end module commondata


    program main
    use commondata
    implicit none
    integer :: i, j
    real(8) :: xp(nx), yp(ny)
    real(8) :: start, finish

    do i=1,nx
        Xp(i) = (i-1)
    enddo
    do j=1,ny
        Yp(j) = (j-1)
    enddo

    A = 0.0d0
    Anew = 0

    do i=1,nx
        A(i,1) = 0.0d0
        A(i,ny) = 100.0d0
    enddo
    do j=1,ny
        A(1,j) = 0.0d0
        A(nx,j) = 0.0d0
    enddo

    error = 10.0d0
    itc = 0
    
    call CPU_TIME(start)
    
    !$acc data copy(A) create(Anew)
    do while( (error.GT.eps).AND.(itc.LE.itc_max) )

        call calA()

        if(MOD(itc,2000).EQ.0) then
            call errorA()
        endif

        call updateA()

        itc = itc+1

    enddo
    !$acc end data

    call CPU_TIME(finish)
    
    write(*,*) "Time (CPU) = ",finish-start, "s"

    open(unit=02,file='Jacobi.dat',status='unknown')
    write(02,*) 'TITLE="Lid Driven Cavity(MRT)"'
    write(02,*) 'VARIABLES="X" "Y" "A" '
    write(02,101) nx, ny
    do j=1,ny
        do i=1,nx
            write(02,100) xp(i), yp(j), A(i,j)
        enddo
    enddo
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(02)

    end program main


    subroutine calA()
    use commondata
    implicit none
    integer :: i, j

    !$acc parallel loop reduction(max:error)
    do j=2,ny-1
        !$acc loop reduction(max:error)
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

    !$acc parallel loop
    do j=2,ny-1
        !$acc loop
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

    error = 0
    !$acc parallel loop reduction(max:error)
    do j=2,ny-1
        !$acc loop reduction(max:error)
        do i=2,nx-1
            error = max( error, abs(Anew(i,j)-A(i,j)) )
        enddo
    enddo

    write(*,*) itc, error

    return
    end subroutine errorA
