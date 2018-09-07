
    module commondata
        implicit none
        integer, parameter :: nx=800, ny=800
        real(8) :: A(nx,ny), Anew(nx,ny)
        real(8) :: error
        integer :: itc
        real(8), parameter :: eps=1e-6
        integer, parameter :: itc_max=20000
    end module commondata


    program main
    use commondata
    use omp_lib
    implicit none
    integer :: i, j
    real(8) :: xp(nx), yp(ny)
    integer :: myMaxThreads
    real(8) :: start, finish
    real(8) :: start2, finish2

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

#ifdef _OPENMP
    call OMP_set_num_threads(24)
    myMaxThreads = OMP_get_max_threads()
    write(*,*) "Max Threads=",myMaxThreads
#endif

    call CPU_TIME(start)
#ifdef _OPENMP
    start2 = OMP_get_wtime()
#endif

    do while( (error.GT.eps).AND.(itc.LE.itc_max) )

        call calA()

        if(MOD(itc,2000).EQ.0) then
            call errorA()
        endif

        call updateA()

        itc = itc+1

    enddo
    
    call CPU_TIME(finish)
#ifdef _OPENMP
    finish2 = OMP_get_wtime()
#endif
    
    write(*,*) "Time (CPU) = ",finish-start, "s"
#ifdef _OPENMP
    write(*,*) "Time (OMP) = ", finish2-start2, "s"
#endif


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

    end program main


    subroutine calA()
    use commondata
    implicit none
    integer :: i, j

    !$omp parallel do default(none) shared(A,Anew) private(i,j)
    do j=2,ny-1
        do i=2,nx-1
            Anew(i,j) = 0.25d0*( A(i+1,j)+A(i-1,j)+A(i,j+1)+A(i,j-1) )
        enddo
    enddo
    !$omp end parallel do

    return
    end subroutine calA


    subroutine updateA()
    use commondata
    implicit none
    integer :: i, j

    !$omp parallel do default(none) shared(A,Anew) private(i,j)
    do j=2,ny-1
        do i=2,nx-1
            A(i,j) = Anew(i,j)
        enddo
    enddo
    !$omp end parallel do

    return
    end subroutine updateA


    subroutine errorA
    use commondata
    implicit none
    integer :: i, j

    error = 0
    !$omp parallel do shared(A,Anew) private(i,j) reduction(max:error)
    do j=2,ny-1
        do i=2,nx-1
            error = max( error, abs(Anew(i,j)-A(i,j)) )
        enddo
    enddo
    !$omp end parallel do

    write(*,*) itc, error

    return
    end subroutine errorA
