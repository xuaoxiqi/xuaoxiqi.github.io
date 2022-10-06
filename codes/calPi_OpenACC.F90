    module commondata
        implicit none
        real(kind=8), parameter :: Pi=4.0d0*datan(1.0d0)
        integer(kind=4), parameter :: nx=6
        real(kind=8), parameter :: dx=1.0d0/dble(nx-1)
        real(kind=8) :: Xcenter(1:nx-1)
    end module commondata
    
    
    program main
    use commondata
    implicit none
    
    write(*,*) "Grids: From 1 to nx =", nx
    write(*,*) "Integration range: [0, 1]"
    write(*,*) "    "
    
    !$acc data create(Xcenter)
    call trapezoidalRule()
    !$acc end data
    
    write(*,*) "Successfully: simulation completed!"
    
    stop
    end program main
    
    
    subroutine trapezoidalRule()
    use commondata
    implicit none
    integer(kind=4) :: i
    real(kind=8) :: myPi
    
    myPi = 0.0d0
    !$acc parallel loop reduction(+:myPi) present(Xcenter)
    do i = 1, nx-1
        Xcenter(i) = (dble(i)-0.5d0)*dx
        myPi = myPi+4.0d0/(Xcenter(i)*Xcenter(i)+1.0d0)*dx
    enddo
    !$acc end parallel
    
    write(*,*) "Pi (calculation) =", myPi
    write(*,*) "Pi (machine)     =", Pi
    
    return
    end subroutine trapezoidalRule
