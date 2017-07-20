! test some of the fortran functions

subroutine test_cons_to_prim(passed) bind(C, name="test_cons_to_prim")
        ! test cons_to_prim function
    use compute_flux_module, only : cons_to_prim

    implicit none

    logical, intent(inout) :: passed

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
    end if

end subroutine test_cons_to_prim

subroutine test_rhoh_from_p(passed) bind(C, name="test_rhoh_from_p")
    implicit none
    logical, intent(inout) :: passed

    integer, parameter :: lo(3) = (/ 1,1,1 /), hi(3) = (/ 5, 1, 1 /)
    double precision :: rhoh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision :: rhoh_test(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision :: p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision, parameter :: gamma = 5.0d0 / 3.0d0

    rho(:, 1, 1) = (/ 1.0d-3, 1.0d-3, 1.0d3, 1.0d3, 1.124d0 /)
    p(:, 1,1) = (/ 1.0d-3, 1.0d3, 1.0d-3, 1.0d3, 13.12d0 /)
    rhoh(:,1,1) = (/ 3.5d-3, 2500.001d0, 1000.0025d0, 3500d0, 33.924d0 /)

    call rhoh_from_p(rhoh_test, p, rho, gamma, lo, hi, lo, hi)

    if (any(abs(rhoh_test - rhoh) > 1.d-5)) then
        write(*,*) "rhoh_from_p failed :("
        passed = .false.
    else
        write(*,*) "rhoh_from_p passed :D"
    end if

end subroutine test_rhoh_from_p

subroutine test_p_from_rhoh(passed) bind(C, name="test_p_from_rhoh")
    implicit none
    logical, intent(inout) :: passed

    integer, parameter :: lo(3) = (/ 1,1,1 /), hi(3) = (/ 5, 1, 1 /)
    double precision :: rhoh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision :: p_test(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision :: p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision, parameter :: gamma = 5.0d0 / 3.0d0

    rho(:, 1, 1) = (/ 1.0d-3, 1.0d-3, 1.0d3, 1.0d3, 1.124d0 /)
    p(:, 1,1) = (/ 1.0d-3, 1.0d3, 1.0d-3, 1.0d3, 13.12d0 /)
    rhoh(:,1,1) = (/ 3.5d-3, 2500.001d0, 1000.0025d0, 3500d0, 33.924d0 /)

    call p_from_rhoh(rhoh, p_test, rho, gamma, lo, hi)

    if (any(abs(p_test - p) > 1.d-5)) then
        write(*,*) "p_from_rhoh failed :("
        passed = .false.
    else
        write(*,*) "p_from_rhoh passed :D"
    end if

end subroutine test_p_from_rhoh

subroutine test_p_from_rho_eps(passed) bind(C, name="test_p_from_rho_eps")
    implicit none
    logical, intent(inout) :: passed

    integer, parameter :: lo(3) = (/ 1,1,1 /), hi(3) = (/ 5, 1, 1 /)
    double precision :: eps(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision :: p_test(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision :: p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision, parameter :: gamma = 5.0d0 / 3.0d0

    rho(:, 1, 1) = (/ 1.0d-3, 0.1d0, 1.0d3, 1.0d3, 1.124d0 /)
    p(:,1,1) = (/ 6.6666666667d-5, 6.6666666667d-5, 666666.66666667d0, 666.6666667d0, 25.420384d0 /)
    eps(:,1,1) = (/ 0.1d0, 1.d-3, 1.0d3, 1.0d0, 33.924d0 /)

    call p_from_rho_eps(rho, eps, p_test, gamma, lo, hi)

    if (any(abs(p_test - p) > 1.d-5)) then
        write(*,*) "p_from_rho_eps failed :("
        !write(*,*) "delta p: ", abs(p_test - p)
        passed = .false.
    else
        write(*,*) "p_from_rho_eps passed :D"
    end if

end subroutine test_p_from_rho_eps

subroutine test_calc_gamma_up(passed) bind(C, name="test_calc_gamma_up")
    implicit none
    logical, intent(inout) :: passed

end subroutine test_calc_gamma_up

subroutine test_calc_gamma_down(passed) bind(C, name="test_calc_gamma_down")
    implicit none
    logical, intent(inout) :: passed

end subroutine test_calc_gamma_down

subroutine test_gr_sources(passed) bind(C, name="test_gr_sources")
    implicit none
    logical, intent(inout) :: passed

end subroutine test_gr_sources

subroutine test_W_swe(passed) bind(C, name="test_W_swe")
    implicit none
    logical, intent(inout) :: passed

    integer, parameter :: lo(3) = (/ 1, 1, 1 /), hi(3) = (/ 100, 1, 1 /), Ncomp = 3
    double precision :: U_cons(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    double precision :: U_prim(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    double precision :: W(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision :: W_test(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision, parameter :: gamma = 5.0d0/3.0d0, M = 1.0d0, R = 100.0d0
    double precision :: dx(3) = (/ 1.0d-3, 1.0d-3, 1.0d-3 /), gamma_up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 9), v2(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), alpha0, rand
    integer i

    alpha0 = sqrt(1.0d0 - 2.0d0 * M / R)
    gamma_up(:,:,:,:) = 0.0d0
    gamma_up(:,:,:,1) = 1.0d0
    gamma_up(:,:,:,5) = 1.0d0
    gamma_up(:,:,:,9) = (alpha0 + M * dx(3) * lo(3) / (R * alpha0))**2

    do i = 1, hi(1)
        call random_number(rand)
        U_prim(i, 1, 1, 1) = 0.5d0 * rand + 1.1d0
        U_prim(i, 1, 1, 2) = 1.2d0 * rand - 0.6d0
        U_prim(i, 1, 1, 3) = 1.2d0 * rand - 0.6d0
    end do

    ! calculate conservative variables
    v2(:,1,1) = U_prim(:,1,1,2)**2 + U_prim(:,1,1,3)**2

    W = sqrt(1.0d0 / (1.0d0 - v2))

    U_cons(:,1,1,1) = U_prim(:,1,1,1) * W(:,1,1)
    do i = 2, 3
        U_cons(:,1,1,i) = U_prim(:,1,1,1) * W(:,1,1)**2 * U_prim(:,1,1,i)
    end do

    call W_swe(U_cons, lo, hi, lo, hi, Ncomp, gamma_up, lo, hi, W_test)

    if (any(abs(W_test - W) > 1.d-5)) then
        write(*,*) "W_swe failed :("
        !write(*,*) "delta p: ", abs(p_test - p)
        passed = .false.
    else
        write(*,*) "W_swe passed :D"
    end if

end subroutine test_W_swe

subroutine test_swe_from_comp(passed) bind(C, name="test_swe_from_comp")
    implicit none
    logical, intent(inout) :: passed

end subroutine test_swe_from_comp

subroutine test_calc_gamma_swe(passed) bind(C, name="test_calc_gamma_swe")
    implicit none
    logical, intent(inout) :: passed

end subroutine test_calc_gamma_swe

subroutine test_comp_from_swe(passed) bind(C, name="test_comp_from_swe")
    implicit none
    logical, intent(inout) :: passed

end subroutine test_comp_from_swe

subroutine test_gr_swe_flux(passed) bind(C, name="test_gr_swe_flux")
    use compute_flux_module, only : gr_swe_flux
    implicit none
    logical, intent(inout) :: passed

    integer, parameter :: lo(3) = (/ 1, 1, 1 /), hi(3) = (/ 6, 1, 1 /), Ncomp = 3
    double precision :: U_cons(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, Ncomp)
    double precision :: f(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, Ncomp)
    double precision :: f_test(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    double precision :: alpha(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1)
    integer :: i, dirs(lo(1):hi(1)), lims(3)

    U_cons(:,:,:,:) = 0.0d0

    U_cons(1,1,1,:) = (/ 1.d0,  0.d0,  0.d0 /)
    U_cons(2,1,1,:) = U_cons(1,1,1,:)
    U_cons(3,1,1,:) = (/ 1.0d-3,0.5d0,0.0d0/)
    U_cons(4,1,1,:) = U_cons(3,1,1,:)
    U_cons(5,1,1,:) = (/ 1.d3,0.5d0,0.5d0 /)
    U_cons(6,1,1,:) = U_cons(5,1,1,:)
    dirs(:) = (/ 0,1,0,1,0,1 /)
    f(1,1,1,:) = (/ -0.01111111d0,0.005d0, -0.0d0 /)
    f(2,1,1,:) = (/ 0.02222222d0,  0.0d0        ,  0.005d0 /)
    f(3,1,1,:) = (/ 0.00078889d0,  0.39444295d0,  0.0d0 /)
    f(4,1,1,:) = (/ 2.22222222d-04,   1.11111111d-01,   0.0d0 /)
    f(5,1,1,:) = (/ -1.10706112d2,   4.99999742d5,  -5.53530559d-02 /)
    f(6,1,1,:) = (/ 2.22627221d2,   1.11313611d-01,   4.99999909d5 /)

    do i = lo(1), hi(1)
        lims = (/ i,1,1 /)
        call gr_swe_flux(U_cons, f_test, lo-1, hi+1, lims, lims, Ncomp, dirs(i), alpha)
    end do

    if (any(abs(f_test - f(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),:)) > 1.d-5)) then
        write(*,*) "gr_swe_flux failed :("
        write(*,*) "delta f: ", abs(f_test - f(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),:))
        passed = .false.
    else
        write(*,*) "gr_swe_flux passed :D"
    end if

end subroutine test_gr_swe_flux

subroutine test_f_of_p(passed) bind(C, name="test_f_of_p")
    use compute_flux_module, only : f_of_p
    implicit none
    logical, intent(inout) :: passed

    integer, parameter :: lo(3) = (/ 1, 1, 1 /), hi(3) = (/ 4, 1, 1 /), Ncomp = 5
    double precision :: U_cons(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    double precision :: f(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision :: p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision :: f_test(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision, parameter :: gamma = 5.0d0/3.0d0
    double precision :: gamma_up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 9)
    integer i

    do i = lo(1), hi(1)
        gamma_up(i,1,1,:) = (/ 0.80999862d0,  0.0d0,  0.0d0,  0.0d0,  0.80999862d0, 0.0d0,  0.0d0,  0.0d0,  0.80999862d0 /)
    end do

    U_cons(1,1,1,:) = (/ 1.0d-3, 0.0d0, 0.0d0, 0.0d0, 3.0d0 /)
    U_cons(2,1,1,:) = (/ 1.0d-3, 0.4d0, -0.4d0, 0.4d0, 1.0d3 /)
    U_cons(3,1,1,:) = (/ 1.0d3, 0.0d0, 0.0d0, 0.0d0, 1.0d-3 /)
    U_cons(4,1,1,:) = (/ 5.0d0, 0.3d0, 0.1d0, 0.4d0, 1.0d0 /)
    p(:,1,1) = (/ 2.0d0, 50.0d0, 20.0d0, 1.0d0 /)
    f(:,1,1) = (/ 0.d0, 616.66641981029716d0, -19.99933333333335d0, -0.34621947550289045d0 /)

    do i = lo(1), hi(1)
        call f_of_p(f_test(i,1,1), p(i,1,1), U_cons(i,1,1,:), Ncomp, gamma, gamma_up(i,1,1,:))
    end do

    if (any(abs(f_test - f) > 1.d-5)) then
        write(*,*) "f_of_p failed :("
        write(*,*) "delta f: ", abs(f_test - f)
        passed = .false.
    else
        write(*,*) "f_of_p passed :D"
    end if

end subroutine test_f_of_p

subroutine test_zbrent(passed) bind(C, name="test_zbrent")
    use compute_flux_module, only : zbrent
    implicit none
    logical, intent(inout) :: passed

end subroutine test_zbrent

subroutine test_gr_comp_flux(passed) bind(C, name="test_gr_comp_flux")
    use compute_flux_module, only : gr_comp_flux
    implicit none
    logical, intent(inout) :: passed

    integer, parameter :: lo(3) = (/ 1, 1, 1 /), hi(3) = (/ 8, 1, 1 /), Ncomp = 5
    double precision :: U_cons(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, Ncomp)
    double precision :: f(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, Ncomp)
    double precision :: f_test(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    double precision, parameter :: gamma = 5.0d0/3.0d0, alpha0 = 0.9, R = 100.d0, M = 1.0d0
    double precision :: gamma_up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 9), dx(3) = (/ 1.0d-1, 1.0d-1, 1.0d-1 /), alpha(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1)
    integer :: i, dirs(lo(1):hi(1)), lims(3)

    gamma_up(:,:,:,:) = 0.0d0
    gamma_up(:,:,:,1) = 1.0d0
    gamma_up(:,:,:,5) = 1.0d0
    gamma_up(:,:,:,9) = (alpha0 + M * dx(3) * lo(3) / (R * alpha0))**2

    U_cons(:,:,:,:) = 0.0d0

    U_cons(1,1,1,:) = (/ 1.d0,  0.d0,  0.d0,  0.d0,  1.d0 /)
    U_cons(2,1,1,:) = U_cons(1,1,1,:)
    U_cons(3,1,1,:) = U_cons(1,1,1,:)
    U_cons(4,1,1,:) = (/ 1.13133438d0,  1.02393398d0,  1.02393398d0,  1.02393398d0,  1.61511222d0 /)
    U_cons(5,1,1,:) = U_cons(4,1,1,:)
    U_cons(6,1,1,:) = U_cons(4,1,1,:)
    U_cons(7,1,1,:) = U_cons(4,1,1,:)
    U_cons(8,1,1,:) = (/ 0.01012376d0,  0.00104199d0, -0.00104199d0,  0.00104199d0,  0.00022944d0 /)
    dirs(:) = (/ 0,1,2,0,1,0,1,2 /)
    f(1,1,1,:) = (/ -0.11111111d0,  0.66666667d0, -0.0d0, -0.0d0, -0.11111111d0 /)
    f(2,1,1,:) = (/ 0.22222222d0,  0.0d0,  0.66666667d0,  0.0d0,  0.22222222d0 /)
    f(3,1,1,:) = (/ -0.33333333d0, -0.0d0,  0.0d0,  0.66666667d0, -0.33333333d0 /)
    f(4,1,1,:) = (/ 0.14920997d0,  0.80171176d0,  0.13504509d0,  0.13504509d0,  0.41301469d0 /)
    f(5,1,1,:) = (/ 0.52632143d0,  0.47635642d0,  1.14302308d0,  0.47635642d0,  0.95138543d0 /)
    f(6,1,1,:) = (/ 0.00014921d0,  0.00080171d0, -0.00013505d0,  0.00013505d0,  0.00041301d0 /)
    f(7,1,1,:) = (/ -2.35061459d-05,  -2.12746488d-05,   6.87941315d-04, -2.12746488d-05,  -2.33557774d-04 /)
    f(8,1,1,:) = (/ -2.55456349d-03,  -2.62928173d-04,   2.62928173d-04, -1.96261506d-04,  -5.12293429d-05 /)

    do i = lo(1), hi(1)
        lims = (/ i,1,1 /)
        call gr_comp_flux(U_cons, f_test, lims, lims, Ncomp, dirs(i), gamma, lo-1, hi+1, alpha, dx, alpha0, M, R)
    end do

    if (any(abs(f_test - f(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),:)) > 1.d-5)) then
        write(*,*) "gr_comp_flux failed :("
        write(*,*) "delta f: ", abs(f_test - f(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),:))
        passed = .false.
    else
        write(*,*) "gr_comp_flux passed :D"
    end if

end subroutine test_gr_comp_flux
