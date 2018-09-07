
        program main
        use omp_lib
        implicit none
        integer, parameter :: n=10000
        integer :: X(n)
        integer :: A(n), B(n)
        integer :: C(n), D(n)
        integer :: i, j
        integer :: k        
        
#ifdef _OPENMP
        call OMP_set_num_threads(12)
        write(*,*) "I am OpenMP"
#endif

        write(*,*) "Initialize A and B..."
        write(*,*) "Initialize C and D..."
        !$omp parallel do default(none) shared(A,C) private(i) 
        do i=1,n
            A(i) = i
            C(i) = i
        enddo
        !$omp end parallel do
        
        do i=1,n
            B(i) = i
            D(i) = i
            X(i) = 2*i-1
        enddo

        write(*,*) "Update A and B..."
        !!$omp parallel do default(none) shared(X,A) private(i,j) 
        !$omp parallel do default(none) shared(X) reduction(+:A) private(i,j)
        do i=1,n
            do j=i+1,n
                !!$omp atomic
                A(i) = A(i)+X(i)
                !!$omp atomic
                A(j) = A(j)-X(i)
            enddo
        enddo
        !$omp end parallel do
        
        do i=1,n
            do j=i+1,n
                B(i) = B(i)+X(i)
                B(j) = B(j)-X(i)
            enddo
        enddo
        
        write(*,*) "Update C and D..."
        !$omp parallel do default(none) shared(C) private(i)
        do i=1,n
            C(i) = C(i)+2
        enddo
        !$omp end parallel do
        do i=1,n
            D(i) = D(i)+2
        enddo

        write(*,*) "Compare A and B..."
        do i=1,n
            if(A(i).NE.B(i)) write(*,*) i, A(i), B(i)
        enddo
        
        write(*,*) "Compare C and D..."
        do i=1,n
            if(C(i).NE.D(i)) write(*,*) i, C(i), D(i)
        enddo
        
        stop
        end program main
