    module commondata
        implicit none
        real(8), parameter :: Pi=4.0d0*datan(1.0d0)
        integer, parameter :: nx=50000001
        real(8), parameter :: dx=1.0d0/dble(nx-1)
        real(8) :: X(nx)
    end module commondata
    
    
    program main
    use commondata
    implicit none
    
    !$acc data create(X)
    call trapezoidalRule()
    !$acc end data
    
    write(*,*) "Successfully: simulation completed!"
    
    stop
    end program main
    
    
    subroutine trapezoidalRule()
    use commondata
    implicit none
    integer :: i
    real(8) :: myPi
    
    myPi = 0.0d0
    !$acc parallel loop reduction(+:myPi) present(X)
    do i=1,nx
        X(i) = (dble(i)-0.5d0)*dx
        myPi = myPi+4.0d0/(X(i)*X(i)+1.0d0)*dx
    enddo
    !$acc end parallel
    
    write(*,*) "Pi (calculation) =", myPi
    write(*,*) "Pi (machine)     =", Pi
    
    return
    end subroutine trapezoidalRule
