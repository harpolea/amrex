module compute_flux_module

  implicit none

  private

  public :: compute_flux_2d

contains

  subroutine compute_flux_2d(lo, hi, dt, dx, &
                             phi, ph_lo, ph_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi, &
                             phi_p, phi_m, fp, fm, slope, glo, ghi, Ncomp)

    use slope_module, only: slopex, slopey

    integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2), Ncomp
    double precision, intent(in) :: dt, dx(2)
    integer, intent(in) :: ph_lo(2), ph_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),Ncomp)
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),Ncomp)
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),Ncomp)
    double precision, dimension(glo(1):ghi(1),glo(2):ghi(2),Ncomp) :: &
         phi_p, phi_m, fp, fm, slope

    integer :: i, j
    double precision :: dxdt(2), f_p(Ncomp), f_m(Ncomp)
    double precision, parameter :: g = 1.0d-3, gamma = 5.0d0/3.0d0, alpha = 0.5d0
    logical :: gr = .true.

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
            call gr_swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, .true.)
            call gr_swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, .true.)
        else
            call swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, g, .true.)
            call swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, g, .true.)
        end if
    else
        call comp_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, gamma, .true.)
        call comp_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, gamma, .true.)
    end if

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          f_p = 0.5d0 * (fp(i,j,:) + fm(i+1,j,:) + alpha * dxdt(1) * (phi_p(i,j,:) - phi_m(i+1,j,:)))
          f_m = 0.5d0 * (fp(i-1,j,:) + fm(i,j,:) + alpha * dxdt(1) * (phi_p(i-1,j,:) - phi_m(i,j,:)))

          flxx(i,j,:) = -(f_p - f_m)
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
            call gr_swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, .false.)
            call gr_swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, .false.)
        else
            call swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, g, .false.)
            call swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, g, .false.)
        end if
    else
        call comp_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, gamma, .false.)
        call comp_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, gamma, .false.)
    end if

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          f_p = 0.5d0 * (fp(i,j,:) + fm(i,j+1,:) + alpha * dxdt(2) * (phi_p(i,j,:) - phi_m(i,j+1,:)))
          f_m = 0.5d0 * (fp(i,j-1,:) + fm(i,j,:) + alpha * dxdt(2) * (phi_p(i,j-1,:) - phi_m(i,j,:)))

          flxy(i,j,:) = -(f_p - f_m)
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

  subroutine gr_swe_flux(U, f, glo, ghi, lo, hi, Ncomp, x_dir)
      implicit none

      integer, intent(in) :: Ncomp
      integer, intent(in) :: glo(2), ghi(2), lo(2), hi(2)
      double precision, intent(in)  :: U(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      logical, intent(in) :: x_dir

      integer :: i, j
      double precision :: v(2), W(glo(1):ghi(1),glo(2):ghi(2))
      double precision :: alpha(glo(1):ghi(1),glo(2):ghi(2))
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

end module compute_flux_module
