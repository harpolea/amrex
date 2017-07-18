subroutine W_swe(q, slo, shi, lo, hi, Ncomp, gamma_up, glo, ghi, W)
    ! Calculate lorentz factor
    implicit none

    integer, intent(in) :: slo(3), shi(3), lo(3), hi(3), Ncomp, glo(3), ghi(3)
    double precision, intent(in) :: q(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), Ncomp)
    double precision, intent(in) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), slo(3):shi(3), 9)
    double precision, intent(out) :: W(lo(1):hi(1), lo(2):hi(2), slo(3):shi(3))

    integer i,j,l

    do l = slo(3), shi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                if (q(i,j,l,1) < 1.d-20) then
                    W(i,j,l) = 1.0d0
                else
                    W(i,j,l) = sqrt((q(i,j,l,2)**2 * gamma_up(i,j,l,1)+&
                          2.0d0 * q(i,j,l,2) * q(i,j,l,3) * &
                          gamma_up(i,j,l,2) + q(i,j,l,3)**2 * &
                          gamma_up(i,j,l,5)) / q(i,j,l,1)**2 + 1.0d0)
                end if
                ! nan check
                if (W(i,j,l) /= W(i,j,l)) then
                    W(i,j,l) = 1.0d0
                end if
          end do
      end do
    end do

end subroutine W_swe

subroutine swe_from_comp(U_prim, prlo, prhi, U_swe, slo, shi, p_comp, &
     pclo, pchi, p_swe, lo, hi, n_cons_comp, n_swe_comp, &
     alpha0, M, R, dx) bind(C, name="swe_from_comp")
    ! Assume nlayers = 1 as 2d
    implicit none

    integer, intent(in) :: n_cons_comp, n_swe_comp
    integer, intent(in) :: prlo(3), prhi(3), slo(3), shi(3), pclo(3), pchi(3), lo(3), hi(3)
    double precision, intent(inout)  :: U_swe(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), n_swe_comp)
    double precision, intent(in) :: U_prim(prlo(1):prhi(1), prlo(2):prhi(2), prlo(3):prhi(3), n_cons_comp)
    double precision, intent(in) :: p_comp(pclo(1):pchi(1), pclo(2):pchi(2), pclo(3):pchi(3))
    double precision, intent(in) :: p_swe(slo(3):shi(3))
    double precision, intent(in) :: alpha0, M, R, dx(3)

    double precision gamma_up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 9)
    double precision h_comp(lo(3):hi(3)), ssq
    integer neighbour, minl(1)
    double precision zfrac, W
    integer i, j, k
    double precision dz, h

    call calc_gamma_up(gamma_up, lo, hi, lo, hi, alpha0, M, R, dx)

    dz = dx(3)

    do k = lo(3), hi(3)
        h_comp(k) = (hi(3) - lo(3) - k) * dz
    end do

    !write(*,*) "p_swe = ", p_swe(slo(3))

    do k = slo(3), shi(3)
        ! neighbour is the comp layer above
        ! check if this layer is above or below
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                ! find nearest layer
                minl = minloc(abs(p_comp(i,j,lo(3):hi(3)) - p_swe(k)))
                neighbour = minl(1) + lo(3) - 1

                if (p_comp(i,j,neighbour) > p_swe(k) .and. neighbour > lo(3)) then
                    neighbour = neighbour - 1
                end if
                if (neighbour >= hi(3)) then
                    neighbour = neighbour - 1
                end if
            zfrac = 1.0d0 - (p_swe(k) - p_comp(i,j,neighbour)) / (p_comp(i,j,neighbour+1) - p_comp(i,j,neighbour))

            if (abs(zfrac) > 1.e9 .or. zfrac /= zfrac) then
                zfrac = 1.0d0
            end if

            !write(*,*) "zfrac, neighbour", zfrac, neighbour, h_comp(neighbour)

            ! now interpolate and stick primitive compressible variables in U_comp
            ! TODO: slope limit here?
            U_swe(i,j,k,1) = h_comp(neighbour) * zfrac + &
                h_comp(neighbour+1) * (1.d0 - zfrac)
            U_swe(i,j,k,2:3) = U_prim(i,j,neighbour,2:3) * zfrac + &
                U_prim(i,j,neighbour+1,2:3) * (1.d0 - zfrac)
            h = h_comp(neighbour) * zfrac + &
                h_comp(neighbour+1) * (1.d0 - zfrac)
            if (n_swe_comp == 4) then
                U_swe(i,j,k,4) = h
            end if

            ! interpolate W
            ! NOTE: do I interpolate the primitive velocities then calculate W
            ! or do as I've done here and calculate W then interpolate??
            ssq = U_prim(i,j,neighbour,2)**2 * gamma_up(i,j,neighbour,1) + &
                2.0d0 * U_prim(i,j,neighbour,2) * U_prim(i,j,neighbour,3) * &
                    gamma_up(i,j,neighbour,2) + &
                2.0d0 * U_prim(i,j,neighbour,2) * U_prim(i,j,neighbour,4) * &
                    gamma_up(i,j,neighbour,3) + &
                U_prim(i,j,neighbour,3)**2 * gamma_up(i,j,neighbour,5) + &
                2.0d0 * U_prim(i,j,neighbour,3) * U_prim(i,j,neighbour,4) * &
                    gamma_up(i,j,neighbour,6) + &
                U_prim(i,j,neighbour,4)**2 * gamma_up(i,j,neighbour,9)

            ssq = 1.0d0 / sqrt(1.0d0 - ssq)

            W = U_prim(i,j,neighbour+1,2)**2 * gamma_up(i,j,neighbour+1,1) + &
                2.0d0 * U_prim(i,j,neighbour+1,2) * U_prim(i,j,neighbour+1,3) *&
                    gamma_up(i,j,neighbour+1,2) + &
                2.0d0 * U_prim(i,j,neighbour+1,2) * U_prim(i,j,neighbour+1,4) *&
                    gamma_up(i,j,neighbour+1,3) + &
                U_prim(i,j,neighbour+1,3)**2 * gamma_up(i,j,neighbour+1,5) + &
                2.0d0 * U_prim(i,j,neighbour+1,3) * U_prim(i,j,neighbour+1,4) *&
                    gamma_up(i,j,neighbour+1,6) + &
                U_prim(i,j,neighbour+1,4)**2 * gamma_up(i,j,neighbour+1,9)

            W = ssq * zfrac + (1.0d0 - zfrac) / sqrt(1.0d0 - W)
            U_swe(i,j,k,1) = -log(alpha0 + M * h / (0.04*R * alpha0)) * W
            U_swe(i,j,k,2) = U_swe(i,j,k,1) * W * U_swe(i,j,k,2)
            U_swe(i,j,k,3) = U_swe(i,j,k,1) * W * U_swe(i,j,k,3)
            end do
        end do
    end do

    write(*,*) "U_swe", U_swe(lo(1), lo(2), slo(3), :)
    write(*,*) "U_prim", U_prim(lo(1), lo(2), lo(3), :)

    ! calculate conserved variables
    !gamma_up = 0.0d0
    !gamma_up(:,:,:,1) = 1.0d0
    !gamma_up(:,:,:,5) = 1.0d0
    ! don't need this component
    !gamma_down(:,:,:,9) = 1.0d0 / (alpha0 + M * U_swe(:,:,:,4) / (R**2 * alpha0))**2

    !W(:,:,:) = U_swe(:,:,:,2)**2 * gamma_down(:,:,:,1) + &
    !    2.0d0 * U_swe(:,:,:,2) * U_swe(:,:,:,3) * gamma_down(:,:,:,2) + &
    !    U_swe(:,:,:,3)**2 * gamma_down(:,:,:,5)
    !W = 1.0d0 / sqrt (1.0d0 - W)

end subroutine swe_from_comp

subroutine calc_gamma_up_swe(U, slo, shi, lo, hi, n_comp, gamma_up)
    implicit none

    integer, intent(in) :: n_comp
    integer, intent(in) :: slo(3), shi(3), lo(3), hi(3)
    double precision, intent(in) :: U(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), n_comp)
    double precision, intent(out) :: gamma_up(lo(1):hi(1), lo(2):hi(2), slo(3):shi(3), 9)

    gamma_up = 0.0d0

    gamma_up(:,:,:,1) = 1.0d0
    gamma_up(:,:,:,5) = 1.0d0
    gamma_up(:,:,:,9) = exp(-2.0d0 * U(lo(1):hi(1), lo(2):hi(2), :,1))

end subroutine calc_gamma_up_swe

subroutine comp_from_swe(U_comp, clo, chi, U_swe, slo, shi, p, rho, lo, hi, n_cons_comp, n_swe_comp, gamma, dx, alpha0, M, R, nghost) bind(C, name="comp_from_swe")
    ! TODO: what do I do about vertical velocity component????
    use slope_module, only: slopez
    implicit none

    integer, intent(in) :: n_cons_comp, n_swe_comp, nghost
    integer, intent(in) :: clo(3), chi(3), slo(3), shi(3), lo(3), hi(3)
    double precision, intent(in)  :: U_swe(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), n_swe_comp)
    double precision, intent(out) :: U_comp(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3), n_cons_comp)
    double precision, intent(in) :: p(slo(3)+1:shi(3)-1)
    double precision, intent(in) :: rho(slo(3):shi(3))
    double precision, intent(in)  :: gamma, dx(3), alpha0, M, R

    double precision h_swe(lo(1)-nghost:hi(1)+nghost, lo(2)-nghost:hi(2)+nghost, slo(3):shi(3))
    double precision v_swe(lo(1)-nghost:hi(1)+nghost, lo(2)-nghost:hi(2)+nghost, lo(3)-nghost:hi(3)+nghost, 2)
    double precision h_comp(lo(3)-nghost:hi(3)+nghost)
    double precision zfrac
    integer neighbour, minl(1)
    integer i, j, k, nlo(3), nhi(3)
    double precision W(lo(1)-nghost:hi(1)+nghost, lo(2)-nghost:hi(2)+nghost, min(slo(3),lo(3))-nghost:max(shi(3),hi(3))+nghost)
    double precision rhoh(lo(1)-nghost:hi(1)+nghost, lo(2)-nghost:hi(2)+nghost, lo(3)-nghost:hi(3)+nghost)

    double precision gamma_up_swe(lo(1)-nghost:hi(1)+nghost, lo(2)-nghost:hi(2)+nghost, slo(3):shi(3), 9)
    double precision gamma_up(lo(1)-nghost:hi(1)+nghost, lo(2)-nghost:hi(2)+nghost, lo(3)-nghost:hi(3)+nghost, 9)
    double precision dz

    dz = dx(3)
    nlo = lo - nghost
    nhi = hi + nghost

    write(*,*) "comp from swe"
    write(*,*) "U_swe: ", U_swe(lo(1)+3, lo(2)+3, slo(3)+3, :)

    call calc_gamma_up_swe(U_swe, slo, shi, nlo, nhi, n_swe_comp, gamma_up_swe)
    call calc_gamma_up(gamma_up, nlo, nhi, nlo, nhi, alpha0, M, R, dx)
    call W_swe(U_swe, slo, shi, nlo, nhi, n_swe_comp, gamma_up_swe, nlo, nhi, W(:,:,slo(3):shi(3)))

    !write(*,*) "gamma_up", gamma_up_swe(lo(1), lo(2), lo(3), 7:9)
    !write(*,*) "W_swe", W(lo(1), lo(2), slo(3))
    !write(*,*) "p", p

    do k = nlo(3), nhi(3)
        h_comp(k) = (hi(3) - lo(3) - k) * dz
    end do

    !write(*,*) "hcomp: ", h_comp

    if (n_swe_comp > 3) then
        h_swe = U_swe(nlo(1):nhi(1),nlo(2):nhi(2),slo(3):shi(3),4)
    else
        ! HACK
        h_swe = alpha0 * 0.04*R / M * (exp(-U_swe(nlo(1):nhi(1), nlo(2):nhi(2),slo(3):shi(3),1)) - alpha0)
    end if

    do k = slo(3), shi(3)
        do j = nlo(2), nhi(2)
            do i = nlo(1), nhi(1)
                if (U_swe(i,j,k,1) < 1.d-20) then
                    v_swe(i,j,k,1) = U_swe(i,j,k,2)
                    v_swe(i,j,k,2) = U_swe(i,j,k,3)
                else
                    v_swe(i,j,k,1) = U_swe(i,j,k,2) / &
                        (W(i,j,k) * U_swe(i,j,k,1))
                    v_swe(i,j,k,2) = U_swe(i,j,k,3) / &
                        (W(i,j,k) * U_swe(i,j,k,1))
                end if
            end do
        end do
    end do

    U_comp(:,:,:,:) = 0.0d0

    ! calculate layer fracs and interpolate
    do k = nlo(3), nhi(3)
        ! neighbour is the swe layer above
        ! check if this layer is above or below
        do j = nlo(2), nhi(2)
            do i = nlo(1), nhi(1)
                ! find nearest layer
                minl = minloc(abs(h_swe(i,j,slo(3)+1:shi(3)-1) - h_comp(k)))
                neighbour = minl(1) + slo(3)
                if (h_swe(i,j,neighbour) < h_comp(k) .and. neighbour > lo(3)) then
                    neighbour = neighbour - 1
                end if
                if (neighbour >= shi(3)-1) then
                    neighbour = neighbour - 1
                end if
                zfrac = 1.0d0 - (h_swe(i,j,neighbour) - h_comp(k)) / dz
                zfrac = max(0.0d0, zfrac)

                !write(*,*) "zfrac, neighbour, hswe", zfrac, neighbour, h_swe(i,j,slo(3)+1:shi(3)-1)

                ! now interpolate and stick primitive compressible variables in U_comp
                ! TODO: slope limit here?
                U_comp(i,j,k,1) = rho(neighbour) * zfrac + &
                    rho(neighbour+1) * (1.0d0 - zfrac)
                U_comp(i,j,k,2:3) = v_swe(i,j,neighbour,:) * zfrac + &
                    v_swe(i,j,neighbour+1,:) * (1.0d0 - zfrac)

                U_comp(i,j,k,5) = p(neighbour) * zfrac + &
                    p(neighbour+1) * (1.0d0 - zfrac)

                W(i,j,k) = U_comp(i,j,k,2)**2*gamma_up(i,j,k,1) + &
                    2.0d0 * U_comp(i,j,k,2) * U_comp(i,j,k,3) * &
                        gamma_up(i,j,k,2) + &
                    2.0d0 * U_comp(i,j,k,2) * U_comp(i,j,k,4) * &
                        gamma_up(i,j,k,3) + &
                    U_comp(i,j,k,3)**2 * gamma_up(i,j,k,5) + &
                    2.0d0 * U_comp(i,j,k,3) * U_comp(i,j,k,4) * &
                        gamma_up(i,j,k,6) + &
                    U_comp(i,j,k,4)**2 * gamma_up(i,j,k,9)
                !write(*,*) "U_comp = ", U_comp(i,j,k,:)
                W(i,j,k) = 1.0d0 / sqrt(1.0d0 - W(i,j,k))
            end do
        end do
    end do

    !write(*,*) "U_comp", U_comp(lo(1), lo(2), lo(3), :)
    !write(*,*) "W", W(lo(1), lo(2), lo(3))

    call rhoh_from_p(rhoh, U_comp(:,:,:,5), U_comp(:,:,:,1), gamma, clo, chi, nlo, nhi)

    do k = nlo(3), nhi(3)
        do j = nlo(2), nhi(2)
            do i = nlo(1), nhi(1)
                U_comp(i,j,k,1) = U_comp(i,j,k,1) * W(i,j,k)
                U_comp(i,j,k,2) = rhoh(i,j,k) * W(i,j,k)**2 * U_comp(i,j,k,2)
                U_comp(i,j,k,3) = rhoh(i,j,k) * W(i,j,k)**2 * U_comp(i,j,k,3)
                U_comp(i,j,k,4) = rhoh(i,j,k) * W(i,j,k)**2 * U_comp(i,j,k,4)
                U_comp(i,j,k,5) = rhoh(i,j,k) * W(i,j,k)**2 - &
                                    U_comp(i,j,k,5) - U_comp(i,j,k,1)
            end do
        end do
    end do

    write(*,*) "U_comp", U_comp(lo(1), lo(2), lo(3), :)

end subroutine comp_from_swe

subroutine rhoh_from_p(rhoh, p, rho, gamma, clo, chi, lo, hi)
    implicit none

    integer, intent(in) :: clo(3), chi(3), lo(3), hi(3)
    double precision, intent(in)  :: p(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    double precision, intent(in)  :: rho(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    double precision, intent(in)  :: gamma
    double precision, intent(out)  :: rhoh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

    rhoh = rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) + gamma * p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) / (gamma - 1.0d0)

end subroutine rhoh_from_p


subroutine p_from_rhoh(rhoh, p, rho, gamma, lo, hi)
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    double precision, intent(out)  :: p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision, intent(in)  :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision, intent(in)  :: gamma
    double precision, intent(in)  :: rhoh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

    p = (rhoh - rho) * (gamma - 1.0d0) / gamma

end subroutine p_from_rhoh

subroutine p_from_rho_eps(rho, eps, p, gamma, lo, hi)
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    double precision, intent(out)  :: p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision, intent(in)  :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision, intent(in)  :: gamma
    double precision, intent(in)  :: eps(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

    p = (gamma - 1.0d0) * rho * eps

end subroutine p_from_rho_eps

subroutine calc_gamma_up(gamma_up, glo, ghi, lo, hi, alpha0, M, R, dx)
    implicit none

    integer, intent(in) :: glo(3), ghi(3), lo(3), hi(3)
    double precision, intent(out)  :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9)
    double precision, intent(in)  :: alpha0, M, R
    double precision, intent(in)  :: dx(3)

    integer k

    gamma_up(:,:,:,:) = 0.0d0
    gamma_up(:,:,:,1) = 1.0d0
    gamma_up(:,:,:,5) = 1.0d0

    do k = lo(3), hi(3)
        gamma_up(:,:,k,9) = (alpha0 + M * k * dx(3) / (0.04*R * alpha0))**2
    end do
end subroutine calc_gamma_up

subroutine gr_sources(S, slo, shi, U, ulo, uhi, p, plo, phi, alpha, alo, ahi, gamma_up, glo, ghi, M, R, gamma, Ncomp, lo, hi, dx)
    implicit none

    integer, intent(inout) :: slo(3), shi(3), ulo(3), uhi(3), plo(3), phi(3), alo(3), ahi(3), glo(3), ghi(3), lo(3), hi(3), Ncomp
    double precision, intent(out)  :: S(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3))
    double precision, intent(in)  :: U(ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3), Ncomp)
    double precision, intent(in)  :: p(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
    double precision, intent(in)  :: alpha(alo(1):ahi(1), alo(2):ahi(2), alo(3):ahi(3))
    double precision, intent(in)  :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9)
    double precision, intent(in)  :: M, R, gamma, dx(3)

    double precision Ssq(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision h(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision W2(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    double precision S_temp
    integer i,j,k

    Ssq(:,:,:) = U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 2)**2 * &
        gamma_up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),1) + &
        2.0d0 * U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 2) * &
        U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 3) * &
        gamma_up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),2) + &
        2.0d0 * U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 2) * &
        U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 4) * &
        gamma_up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),3) + &
        U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 3)**2 * &
        gamma_up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),5) + &
        2.0d0 * U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 3) * &
        U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 4) * &
        gamma_up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),6) + &
        U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 4)**2 * &
        gamma_up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),9)

    h = 1.0d0 + gamma * &
        (sqrt((U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 5) + &
        p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) + &
        U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1))**2 - Ssq) - &
        p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) * &
        (U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 5) + &
        p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) + &
        U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)) / &
        sqrt((U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 5) + &
        p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) + &
        U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1))**2 - Ssq) - &
        U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)) / &
        U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)

    W2 = 1.0d0 + Ssq / &
        (U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1) * h)**2

    do i = lo(1), hi(1)
        do j = lo(2), hi(2)
            do k = lo(3), hi(3)
                S_temp = - M / R**2 * &
                    (alpha(i,j,k) * U(i,j,k, 4)**2 / W2(i,j,k) + &
                    (U(i,j,k, 5) + p(i,j,k) + U(i,j,k, 1)) / &
                    alpha(i,j,k))
                if (S_temp == S_temp) then ! not nan
                    S(i,j,k) = S(i,j,k) + dx(3) * alpha(i,j,k) * S_temp
                end if

            end do
        end do
    end do

end subroutine gr_sources
