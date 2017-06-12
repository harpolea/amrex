subroutine phi(r, ph, size)
    ! superbee limiter
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: size
    real(amrex_real), intent(in) :: r(size)
    real(amrex_real), intent(out) :: ph(size)

    integer i

    ph = 0.0d0

    do  i = 1, size
        ! Van Leer MC
        if (r(i) > 0.d0) then
            ph(i) = min(2.0d0 * r(i) / (1.0d0 + r(i)), 2.0d0 / (1.0d0 + r(i)))
        else
            ph(i) = 0.0d0
        end if

        ! superbee
        !if (r(i) >= 1.0d0) then
        !    ph(i) = min(2.0d0, min(r(i), 2.0d0 / (1.0d0 + r(i))))
        !elseif ( r(i) >= 0.5d0 ) then
        !    ph(i) = 1.0d0
        !elseif ( r(i) > 0.0d0 ) then
        !    ph(i) = 2.0d0 * r(i)
        !end if
    end do

end subroutine phi

subroutine swe_from_comp(U_prim, U_swe, p_comp, p_swe, lo, hi, n_cons_comp, n_swe_comp, nlayers, gamma_up, alpha0, M, R, dz)
    ! Assume nlayers = 1 as 2d
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: n_cons_comp, n_swe_comp, nlayers
    integer, intent(in) :: lo(3), hi(3)
    real(amrex_real), intent(inout)  :: U_swe(lo(1):hi(1), lo(2):hi(2), nlayers,  n_swe_comp)
    real(amrex_real), intent(in) :: U_prim(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), n_cons_comp)
    real(amrex_real), intent(in) :: p_comp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real), intent(in) :: p_swe(nlayers)
    real(amrex_real), intent(in) :: alpha0, M, R, dz
    real(amrex_real), intent(in) :: gamma_up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 9)

    real(amrex_real) h_comp(lo(3):hi(3)), ssq
    real(amrex_real) h_swe(lo(1):hi(1), lo(2):hi(2), nlayers)
    integer neighbour, minl(1)
    real(amrex_real) zfrac
    real(amrex_real) W(lo(1):hi(1), lo(2):hi(2), nlayers)
    integer i, j, k

    do k = lo(3), hi(3)
        h_comp(k) = (hi(3) - lo(3) - k) * dz
    end do

    do k = 1, nlayers
        ! neighbour is the comp layer above
        ! check if this layer is above or below
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                ! find nearest layer
                minl = minloc(abs(p_comp(i,j,:) - p_swe(k)))
                neighbour = minl(1)
                if (p_comp(i,j,neighbour) > p_swe(k)) then
                    neighbour = neighbour - 1
                end if
            zfrac = 1.0d0 - (p_swe(k) - p_comp(i,j,neighbour)) / (p_comp(i,j,neighbour+1) - p_comp(i,j,neighbour))

            ! now interpolate and stick primitive compressible variables in U_comp
            ! TODO: slope limit here?
            U_swe(i,j,k,1) = h_comp(neighbour) * zfrac + &
                h_comp(neighbour+1) * (1.d0 - zfrac)
            U_swe(i,j,k,2:3) = U_prim(i,j,neighbour,2:3) * zfrac + &
                U_prim(i,j,neighbour+1,2:3) * (1.d0 - zfrac)
            U_swe(i,j,k,4) = h_comp(neighbour) * zfrac + &
                h_comp(neighbour+1) * (1.d0 - zfrac)

            ! interpolate W
            ! NOTE: do I interpolate the primitive velocities then calculate W
            ! or do as I've done here and calculate W then interpolate??
            ssq = U_prim(i,j,neighbour,2)**2 * gamma_up(i,j,neighbour,1) + &
                2.0d0 * U_swe(i,j,neighbour,2) * U_prim(i,j,neighbour,3) * &
                    gamma_up(i,j,neighbour,2) + &
                2.0d0 * U_swe(i,j,neighbour,2) * U_prim(i,j,neighbour,4) * &
                    gamma_up(i,j,neighbour,3) + &
                U_prim(i,j,neighbour,3)**2 * gamma_up(i,j,neighbour,5) + &
                2.0d0 * U_swe(i,j,neighbour,3) * U_prim(i,j,neighbour,4) * &
                    gamma_up(i,j,neighbour,6) + &
                U_prim(i,j,neighbour,4)**2 * gamma_up(i,j,neighbour,9)

            ssq = 1.0d0 / sqrt(1.0d0 - ssq)

            W(i,j,k) = &
                U_prim(i,j,neighbour+1,2)**2 * gamma_up(i,j,neighbour+1,1) + &
                2.0d0 * U_swe(i,j,neighbour+1,2) * U_prim(i,j,neighbour+1,3) *&
                    gamma_up(i,j,neighbour+1,2) + &
                2.0d0 * U_swe(i,j,neighbour+1,2) * U_prim(i,j,neighbour+1,4) *&
                    gamma_up(i,j,neighbour+1,3) + &
                U_prim(i,j,neighbour+1,3)**2 * gamma_up(i,j,neighbour+1,5) + &
                2.0d0 * U_swe(i,j,neighbour+1,3) * U_prim(i,j,neighbour+1,4) *&
                    gamma_up(i,j,neighbour+1,6) + &
                U_prim(i,j,neighbour+1,4)**2 * gamma_up(i,j,neighbour+1,9)

            W(i,j,k) = ssq * zfrac + (1.0d0 - zfrac) / sqrt(W(i,j,k))
            end do
        end do
    end do

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

    U_swe(:,:,:,1) = -log(alpha0 + M * U_swe(:,:,:,4) / (R**2 * alpha0)) * W
    U_swe(:,:,:,2) = U_swe(:,:,:,1) * W * U_swe(:,:,:,2)
    U_swe(:,:,:,3) = U_swe(:,:,:,1) * W * U_swe(:,:,:,3)

end subroutine swe_from_comp

subroutine calc_gamma_up_swe(U, lo, hi, nlayers, n_comp, gamma_up)
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: n_comp, nlayers
    integer, intent(in) :: lo(3), hi(3)
    real(amrex_real), intent(in) :: U(lo(1):hi(1), lo(2):hi(2), nlayers, n_comp)
    real(amrex_real), intent(out) :: gamma_up(lo(1):hi(1), lo(2):hi(2), nlayers, 9)

    gamma_up = 0.0d0

    gamma_up(:,:,:,1) = 1.0d0
    gamma_up(:,:,:,5) = 1.0d0
    gamma_up(:,:,:,9) = exp(-2.0d0 * U(:,:,:,1))

end subroutine calc_gamma_up_swe

subroutine comp_from_swe(U_comp, U_swe, p, rho, lo, hi, n_cons_comp, n_swe_comp, nlayers, gamma, gamma_up, glo, ghi, dz)
    ! TODO: what do I do about vertical velocity component????
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: n_cons_comp, n_swe_comp, nlayers
    integer, intent(in) :: lo(3), hi(3), glo(3), ghi(3)
    real(amrex_real), intent(in)  :: U_swe(lo(1):hi(1), lo(2):hi(2), nlayers, n_swe_comp)
    real(amrex_real), intent(out) :: U_comp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), n_cons_comp)
    real(amrex_real), intent(in) :: p(nlayers)
    real(amrex_real), intent(in) :: rho(nlayers)
    real(amrex_real), intent(in)  :: gamma, dz
    real(amrex_real), intent(in)  :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9)

    real(amrex_real) h_swe(lo(1):hi(1), lo(2):hi(2), nlayers)
    real(amrex_real) v_swe(lo(1):hi(1), lo(2):hi(2), nlayers, 2)
    real(amrex_real) h_comp(lo(3):hi(3))
    real(amrex_real) zfrac
    integer neighbour, minl(1)
    integer i, j, k
    real(amrex_real) W(lo(1):hi(1), lo(2):hi(2), max(nlayers, hi(3)-lo(3)))
    real(amrex_real) rhoh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real) gamma_up_swe(lo(1):hi(1), lo(2):hi(2), nlayers, 9)

    call calc_gamma_up_swe(U_swe, lo, hi, nlayers, n_swe_comp, gamma_up_swe)

    call W_swe(U_swe, lo, hi, nlayers, n_swe_comp, gamma_up_swe, lo, hi, W(:,:,1:nlayers))

    do k = lo(3), hi(3)
        h_comp(k) = (hi(3) - lo(3) - k) * dz
    end do

    h_swe = U_swe(:,:,:,4)
    v_swe(:,:,:,1) = U_swe(:,:,:,2) / (W(:,:,1:nlayers) * U_swe(:,:,:,1))
    v_swe(:,:,:,2) = U_swe(:,:,:,3) / (W(:,:,1:nlayers) * U_swe(:,:,:,1))

    ! calculate layer fracs and interpolate
    do k = lo(3), hi(3)
        ! neighbour is the swe layer above
        ! check if this layer is above or below
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                ! find nearest layer
                minl = minloc(abs(h_swe(i,j,:) - h_comp(k)))
                neighbour = minl(1)
                if (h_swe(i,j,neighbour) < h_comp(k)) then
                    neighbour = neighbour - 1
                end if
            zfrac = 1.0d0 - (h_swe(i,j,neighbour) - h_comp(k)) / dz

            ! now interpolate and stick primitive compressible variables in U_comp
            ! TODO: slope limit here?
            U_comp(i,j,k,1) = rho(neighbour) * zfrac + &
                rho(neighbour+1) * (1.0d0 - zfrac)
            U_comp(i,j,k,2:3) = v_swe(i,j,neighbour,:) * zfrac + &
                v_swe(i,j,neighbour+1,:) * (1.0d0 - zfrac)
            U_comp(i,j,k,5) = p(neighbour) * zfrac + &
                p(neighbour+1) * (1.0d0 - zfrac)
            end do
        end do

        W(:,:,k) = U_comp(:,:,k,2)**2*gamma_up(lo(1):hi(1),lo(2):hi(2),k,1) + &
            2.0d0 * U_comp(:,:,k,2) * U_comp(:,:,k,3) * &
                gamma_up(lo(1):hi(1),lo(2):hi(2),k,2) + &
            2.0d0 * U_comp(:,:,k,2) * U_comp(:,:,k,4) * &
                gamma_up(lo(1):hi(1),lo(2):hi(2),k,3) + &
            U_comp(:,:,k,3)**2 * gamma_up(lo(1):hi(1),lo(2):hi(2),k,5) + &
            2.0d0 * U_comp(:,:,k,3) * U_comp(:,:,k,4) * &
                gamma_up(lo(1):hi(1),lo(2):hi(2),k,6) + &
            U_comp(:,:,k,4)**2 * gamma_up(lo(1):hi(1),lo(2):hi(2),k,9)
        W(:,:,k) = 1.0d0 / sqrt(1.0d0 - W(:,:,k))
    end do

    call rhoh_from_p(rhoh, U_comp(:,:,:,5), U_comp(:,:,:,1), gamma, lo, hi)

    U_comp(:,:,:,1) = U_comp(:,:,:,1) * W(:,:,lo(3):hi(3))
    U_comp(:,:,:,2) = rhoh * W(:,:,lo(3):hi(3))**2 * U_comp(:,:,:,2)
    U_comp(:,:,:,3) = rhoh * W(:,:,lo(3):hi(3))**2 * U_comp(:,:,:,3)
    U_comp(:,:,:,4) = rhoh * W(:,:,lo(3):hi(3))**2 * U_comp(:,:,:,4)
    U_comp(:,:,:,5) = rhoh * W(:,:,lo(3):hi(3))**2 - U_comp(:,:,:,5) - &
                      U_comp(:,:,:,1)

end subroutine comp_from_swe

subroutine rhoh_from_p(rhoh, p, rho, gamma, lo, hi)
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    real(amrex_real), intent(in)  :: p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real), intent(in)  :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real), intent(in)  :: gamma
    real(amrex_real), intent(out)  :: rhoh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

    rhoh = rho + gamma * p / (gamma - 1.0d0)

end subroutine rhoh_from_p


subroutine p_from_rhoh(rhoh, p, rho, gamma, lo, hi)
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    real(amrex_real), intent(out)  :: p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real), intent(in)  :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real), intent(in)  :: gamma
    real(amrex_real), intent(in)  :: rhoh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

    p = (rhoh - rho) * (gamma - 1.0d0) / gamma

end subroutine p_from_rhoh

subroutine p_from_rho_eps(rho, eps, p, gamma, lo, hi)
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    real(amrex_real), intent(out)  :: p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real), intent(in)  :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real), intent(in)  :: gamma
    real(amrex_real), intent(in)  :: eps(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

    p = (gamma - 1.0d0) * rho * eps

end subroutine p_from_rho_eps
