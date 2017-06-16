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
    double precision, parameter :: g = 1.0d-3, alpha = 0.5d0

    dxdt = dx/dt

    call slopex(lo-1, hi+1, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi, Ncomp)
    !call sslope(phi, ph_lo, ph_hi, slope, glo, ghi, lo-1, hi+1, Ncomp, .true.)

    ! compute phi on x faces
    do    j = lo(2), hi(2)
       do i = lo(1)-1  , hi(1)+1
           phi_p(i,j,:) = phi(i,j,:) + 0.5d0 * slope(i,j,:)
           phi_m(i,j,:) = phi(i,j,:) - 0.5d0 * slope(i,j,:)
       end do
    end do

    call swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, g, .true.)
    call swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, g, .true.)

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
    !call sslope(phi, ph_lo, ph_hi, slope, glo, ghi, lo-1, hi+1, Ncomp, .false.)

    ! compute phi on y faces
    do    j = lo(2)-1  , hi(2)+1
       do i = lo(1), hi(1)
           phi_p(i,j,:) = phi(i,j,:) + 0.5d0 * slope(i,j,:)
           phi_m(i,j,:) = phi(i,j,:) - 0.5d0 * slope(i,j,:)
       end do
    end do


    call swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, g, .false.)
    call swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, g, .false.)

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
      double precision, parameter :: half = 0.5d0

      if (x_dir) then
          do    j = lo(2), hi(2)
             do i = lo(1)-1, hi(1)+1
                f(i,j,1) =  U(i,j,2)
                f(i,j,2) = U(i,j,2)**2/U(i,j,1) + half * g * U(i,j,1)**2
                f(i,j,3) =  U(i,j,2) * U(i,j,3) / U(i,j,1)
             end do
          end do
      else
          do    j = lo(2)-1, hi(2)+1
             do i = lo(1), hi(1)
                f(i,j,1) =  U(i,j,3)
                f(i,j,2) = U(i,j,2) * U(i,j,3) / U(i,j,1)
                f(i,j,3) =  U(i,j,3)**2/U(i,j,1) + half * g * U(i,j,1)**2
             end do
          end do
      end if

  end subroutine swe_flux

  subroutine phi(r, ph, size)
      implicit none

      integer, intent(in) :: size
      double precision, intent(in) :: r(size)
      double precision, intent(out) :: ph(size)

      integer i

      ph = 0.0d0

      do  i = 1, size
          ! Van Leer MC
          if (r(i) > 0.d0) then
              ph(i) = min(2.0d0 * r(i) / (1.0d0 + r(i)), 2.0d0 / (1.0d0 + r(i)))
          else
              ph(i) = 0.0d0
          end if
      end do

  end subroutine phi

  subroutine sslope(q, qlo, qhi, s, slo, shi, lo, hi, Ncomp, x_dir)
      implicit none

      integer, intent(in) :: qlo(2), qhi(2), slo(2), shi(2), lo(2), hi(2), Ncomp
      double precision, intent(in) :: q(qlo(1):qhi(1), qlo(2):qhi(2), Ncomp)
      double precision, intent(out) :: s(slo(1):shi(1), slo(2):shi(2), Ncomp)
      logical, intent(in) :: x_dir

      integer i, j
      double precision S_upwind(Ncomp), S_downwind(Ncomp), r(Ncomp), ph(Ncomp)

      do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              if (x_dir) then
                  S_upwind = q(i+1,j,:) - q(i,j,:)
                  S_downwind = q(i,j,:) - q(i-1,j,:)
              else
                  S_upwind = q(i,j+1,:) - q(i,j,:)
                  S_downwind = q(i,j,:) - q(i,j-1,:)
              end if

              s(i,j,:) = 0.5d0 * (S_upwind + S_downwind)

              r = S_upwind / S_downwind

              call phi(r, ph, Ncomp)

              s(i,j,:) = s(i,j,:) * ph
          end do
      end do


  end subroutine sslope

end module compute_flux_module
