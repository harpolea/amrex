subroutine f_of_p(f, p, U, Ncomp, gamma, gamma_up)
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: Ncomp
    real(amrex_real), intent(in)  :: U(Ncomp), p, gamma, gamma_up(9)
    real(amrex_real), intent(out) :: f

    real(amrex_real) sq

    sq = sqrt((U(5) + p + U(1))**2 - U(2)**2*gamma_up(1)-&
        2.0d0 * U(2) * U(3) * gamma_up(2) -&
        2.0d0 * U(2) * U(4) * gamma_up(3) -&
        U(3)**2 * gamma_up(5) -&
        2.0 * U(3) * U(4) * gamma_up(6) -&
        U(4)**2 * gamma_up(9))

    f = (gamma - 1.0d0) * sq / (U(5) + p + U(1)) * &
        (sq - p * (U(5) + p + U(1)) / sq - U(1)) - p

end subroutine f_of_p

subroutine zbrent(p, x1, b, U, Ncomp, gamma, gamma_up)
    ! route finder using brent's method
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: Ncomp
    real(amrex_real), intent(out) :: p
    real(amrex_real), intent(in)  :: U(Ncomp), gamma, gamma_up(9), x1
    real(amrex_real), intent(inout) :: b

    real(amrex_real), parameter :: TOL = 1.0d-12
    integer, parameter :: ITMAX = 100

    real(amrex_real) a, c, d, fa, fb, fc, fs, s
    logical mflag, con1, con2, con3, con4, con5
    integer i

    a = x1
    c = 0.0d0
    d = 0.0d0
    call f_of_p(fa, a, U, Ncomp, gamma, gamma_up)
    call f_of_p(fb, b, U, Ncomp, gamma, gamma_up)
    fc = 0.0d0

    if (fa * fb >= 0.0d0) then
        p = b
        return
    end if

    if (abs(fa) < abs(fb)) then
        d = a
        a = b
        b = d

        d = fa
        fa = fb
        fb = d
    end if

    c = a
    fc = fa

    mflag = .true.

    do i = 1, ITMAX
        if (fa /= fc .and. fb /= fc) then
            s = a*fb*fc / ((fa-fb) * (fa-fc)) + b*fa*fc / ((fb-fa)*(fb-fc)) +&
                c*fa*fb / ((fc-fa)*(fc-fb))
        else
            s = b - fb * (b-a) / (fb-fa)
        end if

        con1 = .false.

        if (0.25d0 * (3.0d0 * a + b) < b) then
            if ( s < 0.25d0 * (3.0d0 * a + b) .or. s > b) then
                con1 = .true.
            end if
        else if (s < b .or. s > 0.25d0  * (3.0d0 * a + b)) then
            con1 = .true.
        end if

        con2 = mflag .and. abs(s - b) >= 0.5d0 * abs(b-c)

        con3 = (.not. mflag) .and. abs(s-b) >= 0.5d0 * abs(c-d)

        con4 = mflag .and. abs(b-c) < TOL

        con5 = (.not. mflag) .and. abs(c-d) < TOL

        if (con1 .or. con2 .or. con3 .or. con4 .or. con5) then
            s = 0.5d0 * (a + b)
            mflag = .true.
        else
            mflag = .false.
        end if

        call f_of_p(fs, s, U, Ncomp, gamma, gamma_up)

        if (abs(fa) < abs(fb)) then
            d = a
            a = b
            b = d

            d = fa
            fa = fb
            fb = d
        end if

        d = c
        c = b
        fc = fb

        if (fa * fs < 0.0d0) then
            b = s
            fb = fs
        else
            a = s
            fa = fs
        end if

        if (fb == 0.0d0 .or. fs == 0.0d0 .or. abs(b-a) < TOL) then
            p = b
            return
        end if

    end do

    p = x1

end subroutine zbrent

subroutine cons_to_prim(U, U_prim, p, lo, hi, Ncomp, gamma, gamma_up, glo, ghi)
    ! convert from conserved variables (D, Sx, Sy, tau) to primitive variables (rho, v^x, v^y, eps). Also outputs the pressure
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: Ncomp
    integer, intent(in) :: lo(3), hi(3), glo(3), ghi(3)
    real(amrex_real), intent(in)  :: U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    real(amrex_real), intent(out) :: U_prim(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    real(amrex_real), intent(out) :: p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real), intent(in)  :: gamma
    real(amrex_real), intent(in)  :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9)

    real(amrex_real) pmin, pmax, ssq, q(Ncomp), fmin, fmax, sq, h, W2
    integer i, j, k

    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                q = U(i, j, k, :)
                ssq = q(2)**2 * gamma_up(i,j,k,1) + &
                    2.0d0 * q(2) * q(3) * gamma_up(i,j,k,2) + &
                    2.0d0 * q(2) * q(4) * gamma_up(i,j,k,3) + &
                    q(3)**2 * gamma_up(i,j,k,5) + &
                    2.0d0 * q(3) * q(4) * gamma_up(i,j,k,6) + &
                    q(4)**2 * gamma_up(i,j,k,9)

                pmin = (1.0d0 - ssq)**2 * q(5) * (gamma - 1.0d0)
                pmax = (gamma - 1.0d0) * (q(5) + q(1)) / (2.0d0 - gamma)

                if (pmin < 0.0d0) then
                    pmin = 0.d0
                end if

                if (pmax < 0.d0 .or. pmax < pmin) then
                    pmax = 1.0d0
                end if

                call f_of_p(fmin, pmin, q, Ncomp, gamma, gamma_up(i,j,k,:))
                call f_of_p(fmax, pmax, q, Ncomp, gamma, gamma_up(i,j,k,:))

                if (fmin * fmax > 0.0d0) then
                    pmin = 0.d0
                end if

                call f_of_p(fmin, pmin, q, Ncomp, gamma, gamma_up(i,j,k,:))

                if (fmin * fmax > 0.0d0) then
                    pmax = pmax * 10.d0
                end if

                call zbrent(p(i,j,k), pmin, pmax, q, Ncomp, gamma, gamma_up(i,j,k,:))

                if (p(i,j,k) /= p(i,j,k) .or. p(i,j,k) < 0.0d0 .or. p(i,j,k) > 1.0d0) then
                    p(i,j,k) = abs((gamma - 1.0d0) * (q(5) + q(1)) / (2.0d0 - gamma))

                    if (p(i,j,k) > 1.0d0) then
                        p(i,j,k) = 1.0d0
                    end if
                end if

                sq = sqrt((q(5) + p(i,j,k) + q(1))**2 - ssq)

                if (sq /= sq) then
                    sq = q(5) + p(i,j,k) + q(1)
                end if

                h = 1.0d0 + gamma * (sq - p(i,j,k) * (q(5) + p(i,j,k) + q(1)) / sq - q(1)) / q(1)
                W2 = 1.0d0 + ssq / (q(1) * h)**2

                U_prim(i,j,k,1) = q(1) * sq / (q(5) + p(i,j,k) + q(1))
                U_prim(i,j,k,2) = (gamma_up(i,j,k,1)*q(2) + gamma_up(i,j,k,2) * q(3) + gamma_up(i,j,k,3) * q(4)) /&
                    (W2 * h * U_prim(i,j,k,1))
                U_prim(i,j,k,3) = (gamma_up(i,j,k,4)*q(2) + gamma_up(i,j,k,5) * q(3) + gamma_up(i,j,k,6) * q(4)) /&
                    (W2 * h * U_prim(i,j,k,1))
                U_prim(i,j,k,4) = (gamma_up(i,j,k,3)*q(2) + gamma_up(i,j,k,6) * q(3) + gamma_up(i,j,k,9) * q(4)) /&
                    (W2 * h * U_prim(i,j,k,1))
                U_prim(i,j,k,5) = (h - 1.0) / gamma

            end do
        end do
    end do

end subroutine cons_to_prim

subroutine comp_flux(U, f, lo, hi, Ncomp, dir, gamma, gamma_up, glo, ghi, beta, alpha)
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: Ncomp, dir
    integer, intent(in) :: lo(3), hi(3), glo(3), ghi(3)
    real(amrex_real), intent(in)  :: U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    real(amrex_real), intent(out) :: f(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    real(amrex_real), intent(in)  :: gamma
    real(amrex_real), intent(in)  :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9)
    real(amrex_real), intent(in)  :: beta(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 3)
    real(amrex_real), intent(in)  :: alpha(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3))

    integer i, j, k
    real(amrex_real), parameter :: half = 0.5d0
    real(amrex_real) p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real) U_prim(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)

    call cons_to_prim(U, U_prim, p, lo, hi, Ncomp, gamma, gamma_up, glo, ghi)

    if (dir == 0) then
        f(:,:,:,1) = U(:,:,:,1) * (U_prim(:,:,:,2) -&
             beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        f(:,:,:,2) = U(:,:,:,2) * (U_prim(:,:,:,2) -&
             beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))&
              + p
        f(:,:,:,3) = U(:,:,:,3) * (U_prim(:,:,:,2) -&
             beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        f(:,:,:,4) = U(:,:,:,4) * (U_prim(:,:,:,2) -&
             beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        f(:,:,:,5) = U(:,:,:,5) * (U_prim(:,:,:,2) -&
             beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))&
              + p * U_prim(:,:,:,2)
    else if (dir == 1) then
        !write(*,*) "y, alpha", alpha(lo(1),lo(2),lo(3))
        f(:,:,:,1) = U(:,:,:,1) * (U_prim(:,:,:,3) -&
             beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        f(:,:,:,2) = U(:,:,:,2) * (U_prim(:,:,:,3) - &
            beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        f(:,:,:,3) = U(:,:,:,3) * (U_prim(:,:,:,3) -&
             beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))&
              + p
        f(:,:,:,4) = U(:,:,:,4) * (U_prim(:,:,:,3) - &
             beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        f(:,:,:,5) = U(:,:,:,5) * (U_prim(:,:,:,3) -&
             beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))&
              + p * U_prim(:,:,:,3)
    else
        !write(*,*) "z, alpha", alpha(lo(1),lo(2),lo(3))
        f(:,:,:,1) = U(:,:,:,1) * (U_prim(:,:,:,4) -&
             beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        f(:,:,:,2) = U(:,:,:,2) * (U_prim(:,:,:,4) - &
            beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        f(:,:,:,3) = U(:,:,:,3) * (U_prim(:,:,:,4) -&
             beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        f(:,:,:,4) = U(:,:,:,4) * (U_prim(:,:,:,4) -&
             beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))&
            + p
        f(:,:,:,5) = U(:,:,:,5) * (U_prim(:,:,:,4) -&
             beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) / alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))&
              + p * U_prim(:,:,:,4)
    end if

end subroutine comp_flux

subroutine compute_flux (U, datalo, datahi, lo, hi, Ncomp,&
     fluxx, fluxy, fluxz, flo, fhi, dx, dt, gamma, gamma_up, glo, ghi, beta, alpha)

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(3), hi(3), Ncomp, datalo(3), datahi(3), flo(3), fhi(3), glo(3), ghi(3)
  real(amrex_real), intent(in)    :: U(datalo(1):datahi(1), datalo(2):datahi(2), datalo(3):datahi(3), Ncomp)
  real(amrex_real), intent(inout) :: fluxx( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3), Ncomp)
  real(amrex_real), intent(inout) :: fluxy( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3), Ncomp)
  real(amrex_real), intent(inout) :: fluxz( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3), Ncomp)
  real(amrex_real), intent(in)    :: dt, dx(3), gamma
  real(amrex_real), intent(in)  :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9)
  real(amrex_real), intent(in)  :: beta(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 3)
  real(amrex_real), intent(in)  :: alpha(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3))

  ! local variables
  integer i,j,k
  real(amrex_real) S_upwind(Ncomp), S_downwind(Ncomp), S(Ncomp), r(Ncomp), ph(Ncomp), f_p(Ncomp), f_m(Ncomp)
  real(amrex_real), parameter :: half = 0.5d0, alph = 0.9d0

  real(amrex_real) Up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, Ncomp)
  real(amrex_real) Um(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, Ncomp)
  real(amrex_real) fp(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, Ncomp)
  real(amrex_real) fm(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, Ncomp)

      !write(*,*) "alpha", alpha(lo(1),lo(2),lo(3))
  ! x fluxes
  do        k = lo(3)-1, hi(3)+1
      do    j = lo(2)-1, hi(2)+1
         do i = lo(1)-1, hi(1)+1
             S_upwind = U(i+1,j,k,:) - U(i,j,k,:)
             S_downwind = U(i,j,k,:) - U(i-1,j,k,:)
             S = half * (S_upwind + S_downwind)

             r = S_upwind / S_downwind

             call phi(r, ph, Ncomp)

             S = S * ph

             Up(i,j,k,:) = U(i,j,k,:) + half * S
             Um(i,j,k,:) = U(i,j,k,:) - half * S
         end do
     end do
 end do

 !write(*,*) g, dt
 !return

 call comp_flux(Up, fp, lo-1, hi+1, Ncomp, 0, gamma, gamma_up, glo, ghi, beta, alpha)
 call comp_flux(Um, fm, lo-1, hi+1, Ncomp, 0, gamma, gamma_up, glo, ghi, beta, alpha)

  do        k = lo(3), hi(3)
      do    j = lo(2), hi(2)
         do i = lo(1), hi(1)
            f_p = half * (fp(i,j,k,:) + fm(i+1,j,k,:) + alph * (Up(i,j,k,:) - Um(i+1,j,k,:)))
            f_m = half * (fp(i-1,j,k,:) + fm(i,j,k,:) + alph * (Up(i-1,j,k,:) - Um(i,j,k,:)))

            fluxx(i,j,k,:) = - alpha(i,j,k) * (f_p - f_m)
         end do
      end do
  end do

  ! y fluxes
  do        k = lo(3)-1, hi(3)+1
      do    j = lo(2)-1, hi(2)+1
         do i = lo(1)-1, hi(1)+1
             S_upwind = U(i,j+1,k,:) - U(i,j,k,:)
             S_downwind = U(i,j,k,:) - U(i,j-1,k,:)
             S = half * (S_upwind + S_downwind)

             r = S_upwind / S_downwind

             call phi(r, ph, Ncomp)

             S = S * ph

             Up(i,j,k,:) = U(i,j,k,:) + half * S
             Um(i,j,k,:) = U(i,j,k,:) - half * S
         end do
     end do
 end do

 call comp_flux(Up, fp, lo-1, hi+1, Ncomp, 1, gamma, gamma_up, glo, ghi, beta, alpha)
 call comp_flux(Um, fm, lo-1, hi+1, Ncomp, 1, gamma, gamma_up, glo, ghi, beta, alpha)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             f_p = half * (fp(i,j,k,:) + fm(i,j+1,k,:) + alph * (Up(i,j,k,:) - Um(i,j+1,k,:)))
             f_m = half * (fp(i,j-1,k,:) + fm(i,j,k,:) + alph * (Up(i,j-1,k,:) - Um(i,j,k,:)))

             fluxy(i,j,k,:) = -alpha(i,j,k) * (f_p - f_m)
          end do
       end do
   end do

   ! z fluxes
   do        k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
              S_upwind = U(i,j,k+1,:) - U(i,j,k,:)
              S_downwind = U(i,j,k,:) - U(i,j,k-1,:)
              S = half * (S_upwind + S_downwind)

              r = S_upwind / S_downwind

              call phi(r, ph, Ncomp)

              S = S * ph

              Up(i,j,k,:) = U(i,j,k,:) + half * S
              Um(i,j,k,:) = U(i,j,k,:) - half * S
          end do
      end do
  end do

  call comp_flux(Up, fp, lo-1, hi+1, Ncomp, 2, gamma, gamma_up, glo, ghi, beta, alpha)
  call comp_flux(Um, fm, lo-1, hi+1, Ncomp, 2, gamma, gamma_up, glo, ghi, beta, alpha)

     do       k = lo(3), hi(3)
        do    j = lo(2), hi(2)
           do i = lo(1), hi(1)
              f_p = half * (fp(i,j,k,:) + fm(i,j,k+1,:) + alph * (Up(i,j,k,:) - Um(i,j,k+1,:)))
              f_m = half * (fp(i,j,k-1,:) + fm(i,j,k,:) + alph * (Up(i,j,k-1,:) - Um(i,j,k,:)))

              fluxz(i,j,k,:) = -alpha(i,j,k) * (f_p - f_m)
           end do
        end do
    end do

  !write(*,*) fluxx(5,5,1)

end subroutine compute_flux

subroutine grav_sources(U, datalo, datahi, lo, hi, Ncomp, f, flo, fhi, dx,&
     gamma, gamma_up, glo, ghi, alpha, M, R)

    use amrex_fort_module, only : amrex_real
    implicit none

    integer datalo(3), datahi(3), lo(3), hi(3), flo(3), fhi(3), Ncomp, glo(3), ghi(3)
    real(amrex_real), intent(in)    :: U(datalo(1):datahi(1), datalo(2):datahi(2), datalo(3):datahi(3), Ncomp)
    real(amrex_real), intent(inout) :: f(flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3), Ncomp)
    real(amrex_real), intent(in)    :: dx(3), gamma
    real(amrex_real), intent(in)  :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3),9)
    real(amrex_real), intent(in) :: alpha(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3))
    real(amrex_real), intent(in) :: M, R

    real(amrex_real) U_prim(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    real(amrex_real) p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real) ssq(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real) sq(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real) rhoh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real) W2(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

    call cons_to_prim(U, U_prim, p, lo, hi, Ncomp, gamma, gamma_up, glo, ghi)

    ssq = U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2)**2 *&
        gamma_up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) + &
        2.0d0 * U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) *&
         U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) *&
          gamma_up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) + &
        2.0d0 * U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) *&
         U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),4) *&
          gamma_up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) + &
        U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3)**2 *&
         gamma_up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5) + &
        2.0d0 * U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) *&
         U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),4) *&
          gamma_up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),6) + &
        U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),4)**2 *&
         gamma_up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),9)

    sq = sqrt((U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5) + p + U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))**2 - ssq)

    rhoh = U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) + &
        gamma * (sq - p * (U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5) + p +&
         U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)) / sq -&
          U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    W2 = 1.0d0 + ssq / rhoh**2

    f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),4) = - M / R**2 *&
        (alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) *&
         U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),4)**2 / W2 +&
         (U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5) + p +&
          U(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)) /&
           alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

end subroutine grav_sources


subroutine update_data (lo, hi, Ncomp, dataold, polo, pohi, &
     datanew, pnlo, pnhi, dx, dt, gamma, gamma_up, glo, ghi, beta, blo, bhi, alpha, alo, ahi, alpha0, M, R) bind(C, name="update_data")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(3), hi(3), Ncomp, polo(3), pohi(3), pnlo(3), pnhi(3), glo(3), ghi(3), blo(3), bhi(3), alo(3), ahi(3)
  real(amrex_real), intent(in)    :: dataold(polo(1):pohi(1),polo(2):pohi(2),polo(3):pohi(3), Ncomp)
  real(amrex_real), intent(inout) :: datanew(pnlo(1):pnhi(1),pnlo(2):pnhi(2),pnlo(3):pnhi(3), Ncomp)
  real(amrex_real), intent(in)    :: dx(3), gamma
  real(amrex_real), value         :: dt
  real(amrex_real), intent(in)    :: gamma_up(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),9)
  real(amrex_real), intent(in)    :: beta(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3)
  real(amrex_real), intent(in)    :: alpha(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
  real(amrex_real), value         :: alpha0, M, R

  ! local variables
  integer i,j,k
  real(amrex_real) :: dtdx(3)

  real(amrex_real) fluxx(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, lo(3)-4:hi(3)+4, Ncomp)
  real(amrex_real) fluxy(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, lo(3)-4:hi(3)+4, Ncomp)
  real(amrex_real) fluxz(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, lo(3)-4:hi(3)+4, Ncomp)

  dtdx = dt/dx

  !write(*,*) "dataold", dataold(lo(1), lo(2), lo(3),:)
  !write(*,*) "alpha", alpha(lo(1)-4,lo(2)-4,lo(3)-4)

  call compute_flux (dataold, polo, pohi, lo-4, hi+4, Ncomp, &
       fluxx, fluxy, fluxz, lo-4, hi+4, dx, dt, gamma, gamma_up, glo, ghi, beta, alpha)

    do      k = lo(3)-4, hi(3)+4
       do   j = lo(2)-4, hi(2)+4
         do i = lo(1)-4, hi(1)+4
            datanew(i,j,k,:) = dataold(i,j,k,:) &
                 + fluxx(i,j,k,:) * dtdx(1) + fluxy(i,j,k,:) * dtdx(2) &
                 + fluxz(i,j,k,:) * dtdx(3)
         end do
      end do
  end do

  !write(*,*) "datanew", datanew(lo(1), lo(2), lo(3),:)
  !write(*,*) "fluxx*dtdx", fluxx(lo(1), lo(2), lo(3),:)*dtdx(1)
  !write(*,*) "fluxy", fluxy(lo(1), lo(2), lo(3),:)*dtdx(2)
  !write(*,*) "fluxz", fluxz(lo(1), lo(2), lo(3),:)*dtdx(3)


  call grav_sources(datanew, pnlo, pnhi, lo-4, hi+4, Ncomp, fluxx, lo-4, hi+4,&
   dx, gamma, gamma_up, glo, ghi, alpha, M, R)

  do      k = lo(3)-4, hi(3)+4
     do   j = lo(2)-4, hi(2)+4
       do i = lo(1)-4, hi(1)+4
          datanew(i,j,k,4) = datanew(i,j,k,4) + fluxx(i,j,k,4) * dt
       end do
    end do
  end do

  !write(*,*) "datanew, fluxx*dt", datanew(lo(1), lo(2), lo(3),:), fluxx(lo(1), lo(2), lo(3),4)*dt

  call compute_flux (datanew, pnlo, pnhi, lo-2, hi+2, Ncomp, &
       fluxx, fluxy, fluxz, lo-4, hi+4, dx, dt, gamma, gamma_up, glo, ghi, beta, alpha)

   do       k = lo(3)-2, hi(3)+2
       do   j = lo(2)-2, hi(2)+2
         do i = lo(1)-2, hi(1)+2
            datanew(i,j,k,:) = 0.25d0 * (3.d0 * dataold(i,j,k,:) + &
                 datanew(i,j,k,:) + fluxx(i,j,k,:) * dtdx(1) + &
                 fluxy(i,j,k,:) * dtdx(2) + fluxz(i,j,k,:) * dtdx(3))
         end do
      end do
  end do

  call grav_sources(datanew, pnlo, pnhi, lo-2, hi+2, Ncomp, fluxx, lo-4, hi+4,&
   dx, gamma, gamma_up, glo, ghi, alpha, M, R)

  do      k = lo(3)-2, hi(3)+2
     do   j = lo(2)-2, hi(2)+2
       do i = lo(1)-2, hi(1)+2
          datanew(i,j,k,4) = datanew(i,j,k,4) + 0.25d0 * fluxx(i,j,k,4) * dt
       end do
    end do
  end do

  call compute_flux (datanew, pnlo, pnhi, lo, hi, Ncomp, &
       fluxx, fluxy, fluxz, lo-4, hi+4, dx, dt, gamma, gamma_up, glo, ghi, beta, alpha)

  do        k = lo(3), hi(3)
      do    j = lo(2), hi(2)
         do i = lo(1), hi(1)
            datanew(i,j,k,:) = 1.0d0 / 3.0d0 * (dataold(i,j,k,:) + &
                 2.0d0 * (datanew(i,j,k,:) + fluxx(i,j,k,:) * dtdx(1) + &
                 fluxy(i,j,k,:) * dtdx(2) + fluxz(i,j,k,:) * dtdx(3)))
         end do
      end do
  end do

  call grav_sources(datanew, pnlo, pnhi, lo, hi, Ncomp, fluxx, lo-4, hi+4, dx,&
   gamma, gamma_up, glo, ghi, alpha, M, R)

  do      k = lo(3), hi(3)
     do   j = lo(2), hi(2)
       do i = lo(1), hi(1)
          datanew(i,j,k,4) = datanew(i,j,k,4) + &
                             2.0d0 / 3.0d0 * fluxx(i,j,k,4) * dt
       end do
    end do
  end do

end subroutine update_data
