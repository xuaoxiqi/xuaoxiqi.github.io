    module commondata
        implicit none
        real(8), parameter :: Pi=4.0d0*datan(1.0d0)
        integer, parameter :: nx=50000001
        real(8), parameter :: dx=1.0d0/dble(nx-1)
        real(8) :: X(nx)
    end module commondata
    
    
    program main
    use commondata
    use omp_lib
    implicit none
    integer :: myMaxThreads
    integer :: threadID
    
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
    
    write(*,*) "Successfully: simulation completed!"
    
    stop
    end program main
    
    
    subroutine trapezoidalRule()
    use commondata
    implicit none
    real(8) :: f 
    integer :: i
    real(8) :: myPi
    
    myPi = 0.0d0
    !$omp parallel do default(none) shared(X) private(i) reduction(+:myPi)
    do i=1,nx
        X(i) = (dble(i)-0.5d0)*dx
        myPi = myPi+f(X(i))*dx
    enddo
    !$omp end parallel do
    
    write(*,*) "Pi (calculation) =", myPi
    write(*,*) "Pi (machine)     =", Pi
    
    return
    end subroutine trapezoidalRule
    
    function f(x)
    implicit none
    real(8) :: f
    real(8) :: x
    
    f = 4.0d0/(x*x+1.0d0)
    
    return
    end function f