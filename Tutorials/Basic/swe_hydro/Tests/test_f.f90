! test some of the fortran functions

subroutine test_cons_to_prim(passed) bind(C, name="test_cons_to_prim")

    ! test cons_to_prim function
    use compute_flux_module, only : cons_to_prim

    logical, intent(out) :: passed

    integer, parameter :: lo(3) = (/ 1, 1, 1 /), hi(3) = (/ 100, 1, 1 /), Ncomp = 5
    double precision :: U_cons(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    double precision :: U_prim(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    double precision :: U_prim_test(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    double precision :: p(lo(1):hi(1))
    double precision :: h(lo(1):hi(1))
    double precision, parameter :: gamma = 5.0d0/3.0d0, M = 1.0d0, R = 100.0d0
    double precision :: dx(3) = (/ 1.0d-3, 1.0d-3, 1.0d-3 /), gamma_z, v2(lo(1):hi(1)), W(lo(1):hi(1)), alpha0, rand
    integer i

    alpha0 = sqrt(1.0d0 - 2.0d0 * M / R)
    gamma_z = (alpha0 + M * dx(3) * lo(3) / (R * alpha0))**2

    do i = 1, hi(1)
        call random_number(rand)
        U_prim(i, 1, 1, 1) = 10.d0 * rand
        U_prim(i, 1, 1, 2) = 0.8d0 * rand - 0.4d0
        U_prim(i, 1, 1, 3) = rand - 0.5d0
        U_prim(i, 1, 1, 4) = rand - 0.5d0
        U_prim(i, 1, 1, 5) = 1.0d-2 * rand
    end do

    ! calculate conservative variables
    v2 = U_prim(:,1,1,2)**2 + U_prim(:,1,1,3)**2 + &
         gamma_z * U_prim(:,1,1,4)**2

    W = sqrt(1.0d0 / (1.0d0 - v2))

    h = 1.0d0 + gamma * U_prim(:,1,1,5)
    p = (gamma - 1.0d0) * U_prim(:,1,1,1) * U_prim(:,1,1,5)

    U_cons(:,1,1,1) = U_prim(:,1,1,1) * W
    do i = 2, 4
        U_cons(:,1,1,i) = U_prim(:,1,1,1) * h * W**2 * U_prim(:,1,1,i)
    end do
    U_cons(:,1,1,5) = U_prim(:,1,1,1) * W * (h * W - 1.0d0) - p

    call cons_to_prim(U_cons, lo, hi, U_prim_test, lo, hi, p, lo, hi, lo, hi, Ncomp, gamma, alpha0, M, R, dx)

    if (any(abs(U_prim_test - U_prim) > 1.0d-10)) then
        write(*,*) "cons_to_prim failed :'("
        do i = lo(1), hi(1)
            !write(*,*) "U_prim: ", U_prim(i,:,:,:)
            !write(*,*) "U_prim_test: ", U_prim_test(i,:,:,:)
            !write(*,*) "p: ", p(i)
            write(*,*) "delta U_prim: ", abs(U_prim_test(i,:,:,:) - U_prim(i,:,:,:))
        end do
        passed = .false.
    else
        write(*,*) "cons_to_prim passed :D"
        do i = lo(1), hi(1)
            !write(*,*) "delta U_prim: ", abs(U_prim_test(i,:,:,:) - U_prim(i,:,:,:))
        end do
        passed = .true.
    end if

end subroutine test_cons_to_prim
