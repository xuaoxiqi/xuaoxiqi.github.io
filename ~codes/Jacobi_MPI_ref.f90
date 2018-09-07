    module commondata
        real(8), parameter :: a=1.0d0, b=1.0d0
        integer(4), parameter :: nx=2000, ny=2000
        integer(4), parameter :: nit=20000
    end module commondata

program jacobi_parallel_poison_f90
    !
    ! solve the discrete 2-D Poisson equation
    ! with Parallel Jacobi iteration method
    !
    use mpi
    use commondata
    implicit none

    integer(4) :: npx, npy
    integer(4) :: nlx, nly, np, myid, myidx, myidy
    integer(4) :: i, j, err1, ierr, namelen
    integer(4),parameter :: nsize=1024*1024*64/8
    real(8) :: hx, hy, d, dx, dy, t0, t1
    real(8) :: u(nsize), ut(nsize), f(nsize)
    character(LEN=MPI_MAX_PROCESSOR_NAME) :: hostname

! ----------------------------------------------------------------------
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

    if (myid == 0) then
        npx = int(sqrt(real(np)))
        npy = int(sqrt(real(np)))
        write(*,*) "a=",real(a),"   ,b=",real(b)
        write(*,*) "nx=",nx,"   ,ny=",ny
        write(*,*) "Max. interation=",nit
        write(*,*) "npx=",npx," ,npy=",npy
        if ( npx*npy /= np) then
            write(*,'(1x,A)') 'Error: npx*npy not equal to np!'
            call MPI_ABORT(MPI_COMM_WORLD, err1, ierr)
        endif
    endif

    call MPI_BCAST(npx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(npy,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    ! compute 2-D processor coordinates: myidx, myidy
    !
    myidx = mod(myid, npx)
    myidy = myid / npx

    ! compute number of local grid points: nlx, nly
    nlx = (nx-1) / npx + 1
    if ( myidx < mod(nx-1, npx))  nlx = nlx + 1
    nly = (ny-1) / npy + 1
    if ( myidy < mod(ny-1, npy))  nly = nly + 1

    ! check buffer
    if ((nlx+1)*(nly+1) > nsize ) then
        write(*,'(1x,A)') 'Error: insufficient buffer!(enlarge nsize)'
        call MPI_ABORT(MPI_COMM_WORLD, err1, ierr)
    endif

    hx = a / dble(nx)
    hy = b / dble(ny)

    ! inital the data
    call init(nlx,nly,npx,npy,myidx,myidy,hx,hy,u,ut,f)

    t0 = MPI_WTIME()

    ! Jacobi iteration
    call jacobi(nlx,nly,npx,npy,myidx,myidy,hx,hy,u,f)
    
    t1 = MPI_WTIME()

    if (myid == 0 ) write(*,'(1x,A,f12.4)') , 'Elapsed time: ', t1-t0

    ! check the solution
    call checksol(nlx, nly, myidx, myidy, u, ut)

    call MPI_FINALIZE(ierr)

!----------------------------------------------------------------------
contains
    function g(x,y)
        real(8) :: g, x, y
        g = -(x*x + y*y)/4d0
    end function g


    subroutine init(nlx,nly,npx,npy,myidx,myidy,hx,hy,u,ut,f)
    use mpi
    use commondata
    implicit none
    integer(4) :: nlx, nly, npx, npy, myidx, myidy
    integer(4) :: x0, y0, i, j, ierr
    real(8) :: u(0:nlx, 0:nly), ut(0:nlx, 0:nly), f(nlx-1, nly-1)
    real(8) :: hx, hy, d, dx, dy, xi, yj

    ! compute the coordinates of (0, 0) of the subdomain
    x0 = myidx * (nx-1)/npx + min(myidx, mod(nx-1, npx))
    y0 = myidy * (ny-1)/npy + min(myidy, mod(ny-1, npy))

    d  = 1d0 / (2/(hx*hx) + 2/(hy*hy))
    dx = d/(hx*hx)
    dy = d/(hy*hy)

    ! true solution: ut(x,y) = - (x^2 + y^2)/4
    do j = 0, nlx
    do i = 0, nly
        xi = hx*dble(i+x0)
        yj = hy*dble(j+y0)
        ut(i,j) = -0.25d0 * (xi*xi + yj*yj)
    enddo
    enddo

    ! initial guess: u(x,y) = 0
    ! right hand side : f(x,y) = 1 * d
    do j = 1, nly - 1
    do i = 1, nlx - 1
        u(i,j) = 0d0
        f(i,j) = 1d0 * d
    enddo
    enddo

    ! boundary condition: g(x,y) = -(x^2 + y^2)/4
    if ( myidx .eq. 0) then
        do j = 0, nly
            u(0, j) = ut(0, j)
        enddo
    endif
    if ( myidx .eq. npx-1) then
        do j = 0, nly
            u(nlx, j) = ut(nlx, j)
        enddo
    endif
    if ( myidy .eq. 0) then
        do i = 0, nlx
            u(i, 0) = ut(i, 0)
        enddo
    endif
    if ( myidy .eq. npy-1) then
        do i = 0, nlx
            u(i, nly) = ut(i, nly)
        enddo
    endif

    return
    end subroutine init

 ! ----------------------------------------------------------------------
    subroutine jacobi(nlx,nly,npx,npy,myidx,myidy,hx,hy,u,f)
    use mpi
    use commondata
    implicit none
    integer(4) :: nlx, nly, npx, npy, myidx, myidy
    integer(4) :: left, right, above, below, type1, type2, i, j, iter
    integer(4) :: status(MPI_STATUS_SIZE), err1, ierr
    real(8)    :: hx, hy, d, dx, dy, res0, res
    real(8)    :: u(0:nlx, 0:nly), unew(0:nlx, 0:nly), f(nlx-1, nly-1)

    d  = 1d0 / (2/(hx*hx) + 2/(hy*hy))
    dx = d/(hx*hx)
    dy = d/(hy*hy)

    left = myidx - 1
    if ( left >= 0) then
        left = myidy * npx + left
    else
        left = MPI_PROC_NULL
    endif

    right = myidx + 1
    if ( right < npx) then
        right = myidy * npx + right
    else
        right = MPI_PROC_NULL
    endif

    above = myidy + 1
    if ( above < npy) then
        above = above * npx + myidx
    else
        above = MPI_PROC_NULL
    endif

    below = myidy -1
    if ( below >= 0) then
        below = below * npx + myidx
    else
        below = MPI_PROC_NULL
    endif

    ! create two new mpi datatype
    call MPI_TYPE_CONTIGUOUS(nlx-1, MPI_DOUBLE_PRECISION, type1, ierr)
    call MPI_TYPE_VECTOR(nly-1, 1, nlx+1, MPI_DOUBLE_PRECISION, type2, ierr)
    call MPI_TYPE_COMMIT(type1, ierr)
    call MPI_TYPE_COMMIT(type2, ierr)

    ! begin iteration loop
    do iter = 1, nit
        if(MOD(iter,100).EQ.0) then
            if(myid.EQ.0) then
                write(*,*) "iter=",iter
            endif
        endif
        
        do j = 1, nly - 1
        do i = 1, nlx - 1
            unew(i,j) = f(i,j) + dx * (u(i+1,j) + u(i-1,j)) + dy * (u(i,j+1) + u(i,j-1))
        enddo
        enddo

        ! update approximation
        do j = 1, nly - 1
        do i = 1, nlx - 1
            u(i,j) = unew(i,j)
        enddo
        enddo

        ! update inner boundary values
        ! u(1:nlx-1, 1) -> below, u(1:nlx-1, 0) <- below
        call MPI_SENDRECV(u(1,1), 1, type1, below, 111, u(1,0), 1, type1, below, 111, MPI_COMM_WORLD, status, ierr)

        ! u(1:nlx-1,nly-1) -> above, u(1:nlx-1,nly) <- above
        call MPI_SENDRECV(u(1,nly-1), 1, type1, above, 111, u(1,nly), 1, type1, above, 111, MPI_COMM_WORLD, status, ierr)

        ! u(1,1:nly-1) -> left, u(0,1:nly-1) <- left
        call MPI_SENDRECV(u(1,1), 1, type2, left,  111, u(0,1), 1, type2, left,  111, MPI_COMM_WORLD, status, ierr)

        ! u(nlx-1,1:nly-1) -> right, u(nlx,1:nly-1) <- right
        call MPI_SENDRECV(u(nlx-1,1), 1, type2, right, 111, u(nlx,1), 1, type2, right, 111, MPI_COMM_WORLD, status, ierr)

        ! check residual
        if ( .true. ) then
            res0 = 0d0
            do j = 1, nly - 1
            do i = 1, nlx - 1
                res = f(i,j) - u(i,j) + dx * (u(i+1,j) + u(i-1,j)) + dy * (u(i,j+1) + u(i,j-1))
                res0 = max(res0, abs(res/d))
            enddo
            enddo
            call MPI_REDUCE(res0, res, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
            ! if ( myidx+myidy .eq. 0 .and. mod(iter,10) .eq. 0)
            ! & write(*,'(A5,i5,A12,f20.12)') 'iter=', iter, ',  residual=', res
        endif
    enddo
    call  MPI_TYPE_FREE(type1, ierr)
    call  MPI_TYPE_FREE(type2, ierr)

    return
    end subroutine jacobi

! ----------------------------------------------------------------------
    subroutine checksol(nlx, nly, myidx, myidy, u, ut)
    use mpi
    implicit none
    integer(4) :: nlx, nly, myidx, myidy, i, j, ierr
    real(8)    :: u(0:nlx, 0:nly), ut(0:nlx, 0:nly), err0, err

    err0 = 0d0
    do j = 1, nly - 1
    do i = 1, nlx - 1
        err0 = max(err0, abs(ut(i,j) - u(i,j)))
    enddo
    enddo

    call MPI_REDUCE(err0, err, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    if ( myidx+myidy == 0) write(*,'(A, f20.16)'), 'Error=', err

    return
    end subroutine checksol

end program

