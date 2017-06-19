module compute_flux_module

  implicit none

  private

  public :: compute_flux_2d

contains

  subroutine compute_flux_2d(lo, hi, dt, dx, &
                             phi, ph_lo, ph_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi, &
                             phi_p, phi_m, fp, fm, slope, glo, ghi, Ncomp, gr)

    use slope_module, only: slopex, slopey

    integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2), Ncomp
    double precision, intent(in) :: dt, dx(2)
    integer, intent(in) :: ph_lo(2), ph_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),Ncomp)
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),Ncomp)
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),Ncomp)
    logical, intent(in) :: gr
    double precision, dimension(glo(1):ghi(1),glo(2):ghi(2),Ncomp) :: &
         phi_p, phi_m, fp, fm, slope

    integer :: i, j
    double precision :: dxdt(2), f_p(Ncomp), f_m(Ncomp)
    double precision, parameter :: g = 1.0d-3, gamma = 5.0d0/3.0d0, alph = 0.5d0
    double precision :: alpha(glo(1):ghi(1), glo(2):ghi(2))

    dxdt = dx/dt

    call slopex(lo-1, hi+1, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi, Ncomp)

    ! compute phi on x faces
    do    j = lo(2), hi(2)
       do i = lo(1)-1  , hi(1)+1
           phi_p(i,j,:) = phi(i,j,:) + 0.5d0 * slope(i,j,:)
           phi_m(i,j,:) = phi(i,j,:) - 0.5d0 * slope(i,j,:)
       end do
    end do

    if (Ncomp < 4) then
        if (gr) then
            call gr_swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, .true., alpha)
            call gr_swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, .true., alpha)
        else
            call swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, g, .true.)
            call swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, g, .true.)
        end if
    else
        if (gr) then
            call gr_comp_flux(phi_p, fp, lo, hi, Ncomp, 0, gamma, glo, ghi, alpha)
            call gr_comp_flux(phi_m, fm, lo, hi, Ncomp, 0, gamma, glo, ghi, alpha)
        else
            call comp_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, gamma, .true.)
            call comp_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, gamma, .true.)
        end if
    end if

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          f_p = 0.5d0 * (fp(i,j,:) + fm(i+1,j,:) + alph * dxdt(1) * (phi_p(i,j,:) - phi_m(i+1,j,:)))
          f_m = 0.5d0 * (fp(i-1,j,:) + fm(i,j,:) + alph * dxdt(1) * (phi_p(i-1,j,:) - phi_m(i,j,:)))

          flxx(i,j,:) = -(f_p - f_m)

          if (gr) then
              flxx(i,j,:) = flxx(i,j,:) * alpha(i,j)
          end if
       end do
    end do

    call slopey(lo-1, hi+1, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi, Ncomp)

    ! compute phi on y faces
    do    j = lo(2)-1  , hi(2)+1
       do i = lo(1), hi(1)
           phi_p(i,j,:) = phi(i,j,:) + 0.5d0 * slope(i,j,:)
           phi_m(i,j,:) = phi(i,j,:) - 0.5d0 * slope(i,j,:)
       end do
    end do

    if (Ncomp < 4) then
        if (gr) then
            call gr_swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, .false., alpha)
            call gr_swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, .false., alpha)
        else
            call swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, g, .false.)
            call swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, g, .false.)
        end if
    else
        if (gr) then
            call gr_comp_flux(phi_p, fp, lo, hi, Ncomp, 1, gamma, glo, ghi, alpha)
            call gr_comp_flux(phi_m, fm, lo, hi, Ncomp, 1, gamma, glo, ghi, alpha)
        else
            call comp_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, gamma, .false.)
            call comp_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, gamma, .false.)
        end if
    end if

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          f_p = 0.5d0 * (fp(i,j,:) + fm(i,j+1,:) + alph * dxdt(2) * (phi_p(i,j,:) - phi_m(i,j+1,:)))
          f_m = 0.5d0 * (fp(i,j-1,:) + fm(i,j,:) + alph * dxdt(2) * (phi_p(i,j-1,:) - phi_m(i,j,:)))

          flxy(i,j,:) = -(f_p - f_m)
          if (gr) then
              flxy(i,j,:) = flxy(i,j,:) * alpha(i,j)
          end if
       end do
    end do

  end subroutine compute_flux_2d

  subroutine swe_flux(U, f, glo, ghi, lo, hi, Ncomp, g, x_dir)
      implicit none

      integer, intent(in) :: Ncomp
      double precision, intent(in) :: g
      integer, intent(in) :: glo(2), ghi(2), lo(2), hi(2)
      double precision, intent(in)  :: U(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      logical, intent(in) :: x_dir

      integer i, j

      if (x_dir) then
          do    j = lo(2), hi(2)
             do i = lo(1)-1, hi(1)+1
                f(i,j,1) =  U(i,j,2)
                f(i,j,2) = U(i,j,2)**2/U(i,j,1) + 0.5d0 * g * U(i,j,1)**2
                f(i,j,3) =  U(i,j,2) * U(i,j,3) / U(i,j,1)
             end do
          end do
      else
          do    j = lo(2)-1, hi(2)+1
             do i = lo(1), hi(1)
                f(i,j,1) =  U(i,j,3)
                f(i,j,2) = U(i,j,2) * U(i,j,3) / U(i,j,1)
                f(i,j,3) =  U(i,j,3)**2/U(i,j,1) + 0.5d0 * g * U(i,j,1)**2
             end do
          end do
      end if

  end subroutine swe_flux

  subroutine W_swe(q, lo, hi, Ncomp, gamma_up, glo, ghi, W)
      ! Calculate lorentz factor
      implicit none

      integer, intent(in) :: lo(2), hi(2), Ncomp, glo(2), ghi(2)
      double precision, intent(in) :: q(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(in) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), 9)
      double precision, intent(out) :: W(glo(1):ghi(1), glo(2):ghi(2))

      integer i,j

      do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              W(i,j) = sqrt((q(i,j,2)**2 * gamma_up(i,j,1)+ &
                    2.0d0 * q(i,j,2) * q(i,j,3) * gamma_up(i,j,2) + &
                    q(i,j,3)**2 * gamma_up(i,j,5)) / q(i,j,1)**2 + 1.0d0)
              ! nan check
              if (W(i,j) /= W(i,j)) then
                  W(i,j) = 1.0d0
              end if
        end do
    end do

  end subroutine W_swe

  subroutine gr_swe_flux(U, f, glo, ghi, lo, hi, Ncomp, x_dir, alpha)
      implicit none

      integer, intent(in) :: Ncomp
      integer, intent(in) :: glo(2), ghi(2), lo(2), hi(2)
      double precision, intent(in)  :: U(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      logical, intent(in) :: x_dir
      double precision, intent(out) :: alpha(glo(1):ghi(1),glo(2):ghi(2))

      integer :: i, j
      double precision :: v(2), W(glo(1):ghi(1),glo(2):ghi(2))
      double precision :: beta(glo(1):ghi(1),glo(2):ghi(2),3)
      double precision :: gamma_up(glo(1):ghi(1),glo(2):ghi(2),9)

      alpha = exp(-U(:,:,1))
      beta = 0.0d0

      gamma_up = 0.0d0
      gamma_up(:,:,1) = 1.0d0
      gamma_up(:,:,5) = 1.0d0

      call W_swe(U, lo-1, hi+1, Ncomp, gamma_up, glo, ghi, W)

      if (x_dir) then
          do    j = lo(2), hi(2)
             do i = lo(1)-1, hi(1)+1
                v(1) = U(i,j,2) / (U(i,j,1) * W(i,j))
                v(2) = U(i,j,3) / (U(i,j,1) * W(i,j))

                f(i,j,1) =  U(i,j,1) * &
                      (v(1)*gamma_up(i,j,1) + v(2)*gamma_up(i,j,2) -&
                       beta(i,j,1) / alpha(i,j))
                f(i,j,2) = U(i,j,2) * &
                      (v(1)*gamma_up(i,j,1) + v(2)*gamma_up(i,j,2) -&
                       beta(i,j,1) / alpha(i,j)) + 0.5d0 * U(i,j,1)**2 / W(i,j)**2
                f(i,j,3) =  U(i,j,3) * &
                      (v(1)*gamma_up(i,j,1) + v(2)*gamma_up(i,j,2) -&
                       beta(i,j,1) / alpha(i,j))

                f(i,j,:) = f(i,j,:)
             end do
          end do
      else
          do    j = lo(2)-1, hi(2)+1
             do i = lo(1), hi(1)
                v(1) = U(i,j,2) / (U(i,j,1) * W(i,j))
                v(2) = U(i,j,3) / (U(i,j,1) * W(i,j))

                f(i,j,1) = U(i,j,1) * &
                      (v(1)*gamma_up(i,j,2) + v(2)*gamma_up(i,j,5) -&
                       beta(i,j,2) / alpha(i,j))
                f(i,j,2) = U(i,j,2) * &
                      (v(1)*gamma_up(i,j,2) + v(2)*gamma_up(i,j,5) -&
                       beta(i,j,2) / alpha(i,j))
                f(i,j,3) =  U(i,j,3) * &
                      (v(1)*gamma_up(i,j,2) + v(2)*gamma_up(i,j,5) -&
                       beta(i,j,2) / alpha(i,j)) +  0.5d0 * U(i,j,1)**2 / W(i,j)**2

                f(i,j,:) = f(i,j,:)
             end do
          end do
      end if

  end subroutine gr_swe_flux

  subroutine comp_flux(U, f, glo, ghi, lo, hi, Ncomp, gamma, x_dir)
      implicit none

      integer, intent(in) :: Ncomp
      double precision, intent(in) :: gamma
      integer, intent(in) :: glo(2), ghi(2), lo(2), hi(2)
      double precision, intent(in)  :: U(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      logical, intent(in) :: x_dir

      integer i, j
      double precision :: p(glo(1):ghi(1), glo(2):ghi(2))

      p = (gamma - 1.d0) * (U(:,:,4) - 0.5d0 * (U(:,:,2)**2 + U(:,:,3)**2) / U(:,:,1))

      if (x_dir) then
          do    j = lo(2), hi(2)
             do i = lo(1)-1, hi(1)+1
                f(i,j,1) =  U(i,j,2)
                f(i,j,2) = U(i,j,2)**2 / U(i,j,1) + p(i,j)
                f(i,j,3) =  U(i,j,2) * U(i,j,3) / U(i,j,1)
                f(i,j,4) = (U(i,j,4) + p(i,j)) * U(i,j,2) / U(i,j,1)
             end do
          end do
      else
          do    j = lo(2)-1, hi(2)+1
             do i = lo(1), hi(1)
                f(i,j,1) =  U(i,j,3)
                f(i,j,2) = U(i,j,2) * U(i,j,3) / U(i,j,1)
                f(i,j,3) =  U(i,j,3)**2 / U(i,j,1) + p(i,j)
                f(i,j,4) = (U(i,j,4) + p(i,j)) * U(i,j,3) / U(i,j,1)
             end do
          end do
      end if

  end subroutine comp_flux

  subroutine f_of_p(f, p, U, Ncomp, gamma, gamma_up)
      implicit none

      integer, intent(in) :: Ncomp
      double precision, intent(in)  :: U(Ncomp), p, gamma, gamma_up(9)
      double precision, intent(out) :: f

      double precision :: sq

      sq = sqrt((U(4) + p + U(1))**2 - U(2)**2*gamma_up(1)-&
          2.0d0 * U(2) * U(3) * gamma_up(2) -&
          U(3)**2 * gamma_up(4))

      f = (gamma - 1.0d0) * sq / (U(4) + p + U(1)) * &
          (sq - p * (U(4) + p + U(1)) / sq - U(1)) - p

  end subroutine f_of_p

  subroutine zbrent(p, x1, b, U, Ncomp, gamma, gamma_up)
      ! route finder using brent's method
      implicit none

      integer, intent(in) :: Ncomp
      double precision, intent(out) :: p
      double precision, intent(in)  :: U(Ncomp), gamma, gamma_up(9), x1
      double precision, intent(inout) :: b

      double precision, parameter :: TOL = 1.0d-12
      integer, parameter :: ITMAX = 100

      double precision a, c, d, fa, fb, fc, fs, s
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
      implicit none

      integer, intent(in) :: Ncomp
      integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2)
      double precision, intent(in)  :: U(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(out) :: U_prim(lo(1):hi(1), lo(2):hi(2), Ncomp)
      double precision, intent(out) :: p(lo(1):hi(1), lo(2):hi(2))
      double precision, intent(in)  :: gamma
      double precision, intent(in)  :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), 9)

      double precision :: pmin, pmax, ssq, q(Ncomp), fmin, fmax, sq, h, W2
      integer :: i, j

      do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              q = U(i,j,:)
              ssq = q(2)**2 * gamma_up(i,j,1) + &
                  2.0d0 * q(2) * q(3) * gamma_up(i,j,2) + &
                  q(3)**2 * gamma_up(i,j,5)

              pmin = (1.0d0 - ssq)**2 * q(4) * (gamma - 1.0d0)
              pmax = (gamma - 1.0d0) * (q(4) + q(1)) / (2.0d0 - gamma)

              if (pmin < 0.0d0) then
                  pmin = 0.d0
              end if

              if (pmax < 0.d0 .or. pmax < pmin) then
                  pmax = 1.0d0
              end if

              call f_of_p(fmin, pmin, q, Ncomp, gamma, gamma_up(i,j,:))
              call f_of_p(fmax, pmax, q, Ncomp, gamma, gamma_up(i,j,:))

              if (fmin * fmax > 0.0d0) then
                  pmin = 0.d0
              end if

              call f_of_p(fmin, pmin, q, Ncomp, gamma, gamma_up(i,j,:))

              if (fmin * fmax > 0.0d0) then
                  pmax = pmax * 10.d0
              end if

              call zbrent(p(i,j), pmin, pmax, q, Ncomp, gamma, gamma_up(i,j,:))

              if (p(i,j) /= p(i,j) .or. p(i,j) < 0.0d0 .or. p(i,j) > 1.0d0) then
                  p(i,j) = abs((gamma - 1.0d0) * (q(4) + q(1)) / (2.0d0 - gamma))

                  if (p(i,j) > 1.0d0) then
                      p(i,j) = 1.0d0
                  end if
              end if

              sq = sqrt((q(4) + p(i,j) + q(1))**2 - ssq)

              if (sq /= sq) then
                  sq = q(4) + p(i,j) + q(1)
              end if

              h = 1.0d0 + gamma * (sq - p(i,j) * (q(4) + p(i,j) + q(1)) / sq - q(1)) / q(1)
              W2 = 1.0d0 + ssq / (q(1) * h)**2

              U_prim(i,j,1) = q(1) * sq / (q(4) + p(i,j) + q(1))
              U_prim(i,j,2) = (gamma_up(i,j,1)*q(2) + gamma_up(i,j,2) * q(3)) /&
                  (W2 * h * U_prim(i,j,1))
              U_prim(i,j,3) = (gamma_up(i,j,4)*q(2) + gamma_up(i,j,5) * q(3)) /&
                  (W2 * h * U_prim(i,j,1))
              U_prim(i,j,4) = (h - 1.0d0) / gamma

          end do
      end do

  end subroutine cons_to_prim

  subroutine gr_comp_flux(U, f, lo, hi, Ncomp, dir, gamma, glo, ghi, alpha)
      implicit none

      integer, intent(in) :: Ncomp, dir
      integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2)
      double precision, intent(in)  :: U(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(in)  :: gamma
      double precision, intent(out) :: alpha(glo(1):ghi(1),glo(2):ghi(2))

      integer :: i,j
      double precision :: p(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1)
      double precision :: U_prim(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
      double precision :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), 9)
      double precision :: beta(glo(1):ghi(1), glo(2):ghi(2), 3)
      double precision, parameter :: M = 1.0d0, R = 100.0d0

      alpha = sqrt(1.0d0 - 2.0d0 * M / (R+1.0))
      beta = 0.0d0
      gamma_up = 0.0d0
      gamma_up(:,:,1) = 1.0d0
      gamma_up(:,:,5) = 1.0d0

      call cons_to_prim(U, U_prim, p, lo-1, hi+1, Ncomp, gamma, gamma_up, glo, ghi)

      if (dir == 0) then
          do j = lo(2), hi(2)
              do i = lo(1)-1, hi(1)+1
                  f(i,j,1) = U(i,j,1) * (U_prim(i,j,2) - &
                       beta(i,j,1) / alpha(i,j))
                  f(i,j,2) = U(i,j,2) * (U_prim(i,j,2) - &
                       beta(i,j,1) / alpha(i,j)) + p(i,j)
                  f(i,j,3) = U(i,j,3) * (U_prim(i,j,2) - &
                       beta(i,j,1) / alpha(i,j))
                  f(i,j,4) = U(i,j,4) * (U_prim(i,j,2) - &
                       beta(i,j,1) / alpha(i,j)) + p(i,j) * U_prim(i,j,2)

                  f(i,j,:) = f(i,j,:)
               end do
           end do
      else
          do j = lo(2)-1, hi(2)+1
              do i = lo(1), hi(1)
                  !write(*,*) "y, alpha", alpha(lo(1),lo(2),lo(3))
                  f(i,j,1) = U(i,j,1) * (U_prim(i,j,3) - &
                       beta(i,j,2) / alpha(i,j))
                  f(i,j,2) = U(i,j,2) * (U_prim(i,j,3) - &
                      beta(i,j,2) / alpha(i,j))
                  f(i,j,3) = U(i,j,3) * (U_prim(i,j,3) - &
                       beta(i,j,2) / alpha(i,j)) + p(i,j)
                  f(i,j,4) = U(i,j,4) * (U_prim(i,j,3) - &
                       beta(i,j,2) / alpha(i,j)) + p(i,j) * U_prim(i,j,3)

                  f(i,j,:) = f(i,j,:)
               end do
           end do
      end if
  end subroutine gr_comp_flux

end module compute_flux_module
