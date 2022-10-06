
    module commondata
        implicit none
        real(kind=8), parameter :: Pi=4.0d0*datan(1.0d0)
        integer(kind=4), parameter :: nx=6
        real(kind=8), parameter :: dx=1.0d0/dble(nx-1)
        real(kind=8) :: Xcenter(1:nx-1)
    end module commondata
    
    
    program main
    use commondata
    use omp_lib
    implicit none
    integer(kind=4) :: myMaxThreads
    integer(kind=4) :: threadID
    
    write(*,*) "Grids: From 1 to nx =", nx
    write(*,*) "Integration range: [0, 1]"
    write(*,*) "    "
    
#ifdef _OPENMP
    write(*,*) "Starting OpenMP >>>>>>"
    call OMP_set_num_threads(8)
    myMaxThreads = OMP_get_max_threads()
    write(*,100) "-----Max Running threads:", myMaxThreads
    
    !$omp parallel default(none) private(threadID)
    threadID = omp_get_thread_num()
    write(*,100) "Greetings from thread #", threadID
    !$omp end parallel
    write(*,*) "    "
#endif
100 format(1X,A,I3)

    call trapezoidalRule()
    
    write(*,*) "Successfully: simulation completed!"
    
    stop
    end program main
    
    
    subroutine trapezoidalRule()
    use commondata
    implicit none
    real(kind=8) :: f 
    integer(kind=4) :: i
    real(kind=8) :: myPi
    
    myPi = 0.0d0
    !$omp parallel do default(none) shared(Xcenter) private(i) reduction(+:myPi)
    do i = 1, nx-1
        Xcenter(i) = (dble(i)-0.5d0)*dx
        myPi = myPi+f(Xcenter(i))*dx
    enddo
    !$omp end parallel do
    
    write(*,*) "Pi (calculation) =", myPi
    write(*,*) "Pi (machine)     =", Pi
    
    return
    end subroutine trapezoidalRule
    
    
    function f(x)
    implicit none
    real(kind=8) :: f
    real(kind=8) :: x
    
    f = 4.0d0/(x*x+1.0d0)
    
    return
    end function f