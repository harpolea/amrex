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
    double precision :: dtdx(2), dxdt(2), f_p(Ncomp), f_m(Ncomp)
    double precision, parameter :: g = 1.0d0, alpha = 0.5d0

    dtdx = dt/dx
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
      use amrex_fort_module, only : amrex_real
      implicit none

      integer, intent(in) :: Ncomp
      real(amrex_real), intent(in) :: g
      integer, intent(in) :: glo(2), ghi(2), lo(2), hi(2)
      real(amrex_real), intent(in)  :: U(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      real(amrex_real), intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      logical, intent(in) :: x_dir

      integer i, j
      real(amrex_real), parameter :: half = 0.5d0

      if (x_dir) then
          do    j = lo(2), hi(2)
             do i = lo(1)-1, hi(1)+1
                f(i,j,1) =  U(i,j,2)
                f(i,j,2) = U(i,j,2)**2/U(i,j,1) + half * g * U(i,j,1)**2
                f(i,j,3) =  U(i,j,2) * U(i,j,3)/U(i,j,1)
             end do
          end do
      else
          do    j = lo(2)-1, hi(2)+1
             do i = lo(1), hi(1)
                f(i,j,1) =  U(i,j,3)
                f(i,j,2) = U(i,j,2) * U(i,j,3)/U(i,j,1)
                f(i,j,3) =  U(i,j,3)**2/U(i,j,1) + half * g * U(i,j,1)**2
             end do
          end do
      end if

  end subroutine swe_flux

end module compute_flux_module
