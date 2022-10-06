    module commondata
        integer(kind=4), parameter :: nProc=3
        integer(kind=4), allocatable :: start1d(:), end1d(:), count1d(:)
    end module commondata


    program main
    use commondata
    implicit none
    integer(kind=4), parameter :: n=8
    integer(kind=4) :: local_n
    integer(kind=4) :: i, i_start_global

    write(*,*) "Solution  1: >>>"
    call StartEnd(1, n)
    
    write(*,*) "    "
    write(*,*) "Solution  2: >>>"
    do i = 0, nProc-1
        call decompose_1d(n, local_n, i, nProc, i_start_global)
        write(*,*) "start1d =", i_start_global
        write(*,*) "    "
    enddo

    stop
    end program main
    

!*!*Begin subroutine StartEnd    
    subroutine StartEnd(iS1, iS2)
    use commondata
    implicit none
    integer(kind=4) :: leng, iBlock
    integer(kind=4) :: ir
    integer(kind=4) :: iS1, iS2
    integer(kind=4) :: i
    
    allocate(start1d(0:nProc-1))
    allocate(end1d(0:nProc-1))
    allocate(count1d(0:nProc-1))
    
    leng = iS2-iS1+1
    iBlock = leng/nProc
    ir= leng-iBlock*nProc
    
    write(*,*) "****************************"
    write(*,*) "nProc  =", nProc
    write(*,*) "leng   =", leng
    write(*,*) "iBlock =", iBlock
    write(*,*) "ir =", ir
    write(*,*) "****************************"
    
    do i = 0, nProc-1
        
        if(i.LT.ir) then  
            count1d(i) = iBlock+1
            start1d(i) = iS1+i*(iBlock+1) 
            end1d(i) = start1d(i)+count1d(i)-1

        else
            count1d(i) = iBlock   
            start1d(i) = iS1+i*iBlock+ir   
            end1d(i) = start1d(i)+count1d(i)-1  
        endif
        
        write(*,*) "    "
        write(*,*) "From process #", i
        write(*,100) "start1d =", start1d(i)
        write(*,100) "end1d   =", end1d(i)
        write(*,100) "count1d =", count1d(i)

    enddo
    
100 format(1A, I3)
    
    return
    end subroutine StartEnd
!*!*End subroutine StartEnd



    subroutine decompose_1d(total_n, local_n, rank, num_process, i_start_global)
        implicit none
        integer, intent(in) :: total_n, rank, num_process
        integer, intent(out) :: local_n, i_start_global

        local_n = total_n / num_process
        if (rank < MOD(total_n, num_process)) then
            local_n = local_n + 1
        endif
        write(*,*) "local_n =", local_n
        
        if (local_n > total_n / num_process) then ! --- 5 5 '5' 4 4 4
            i_start_global = local_n * rank
        else                    ! --- 5 5 5 4 '4' 4
            i_start_global = local_n * rank + mod(total_n, num_process)
        endif

    end subroutine decompose_1d