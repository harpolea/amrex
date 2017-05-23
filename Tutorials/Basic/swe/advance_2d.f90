subroutine phi(r, ph)
    use amrex_fort_module, only : amrex_real
    implicit none

    real(amrex_real), intent(in) :: r
    real(amrex_real), intent(out) :: ph

    ph = 0.0d0

    if (r >= 1.0d0) then
        ph = min(2.0d0, min(r, 2.0d0 / (1.0d0 + r)))
    elseif ( r >= 0.5d0 ) then
        ph = 1.0d0
    elseif ( r > 0.0d0 ) then
        ph = 2.0d0 * r
    end if

end subroutine phi

subroutine swe_flux(U, datalo, datahi, f, flo, fhi, lo, hi, Ncomp, g, x_dir)
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: Ncomp
    real(amrex_real), intent(in) :: g
    integer, intent(in) :: datalo(2), datahi(2), flo(2), fhi(2), lo(2), hi(2)
    real(amrex_real), intent(in)    :: U(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
    real(amrex_real), intent(out) :: f(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
    logical, intent(in) :: x_dir

    integer i, j

    if (x_dir) then
        do    j = lo(2), hi(2)
           do i = lo(1)-1, hi(1)+1
              f(i,j,1) =  U(i,j,2)
              f(i,j,2) = U(i,j,2)**2/U(i,j,1) + 0.5d0 * g * U(i,j,1)**2
              f(i,j,3) =  U(i,j,2) * U(i,j,3)
           end do
        end do
    else
        do    j = lo(2)-1, hi(2)+1
           do i = lo(1), hi(1)
              f(i,j,1) =  U(i,j,3)
              f(i,j,2) = U(i,j,2) * U(i,j,3)
              f(i,j,3) =  U(i,j,3)**2/U(i,j,1) + 0.5d0 * g * U(i,j,1)**2
           end do
        end do
    end if

end subroutine swe_flux

subroutine compute_flux (U, datalo, datahi, lo, hi, Ncomp, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx, g, dt) bind(C, name="compute_flux")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), Ncomp, datalo(2), datahi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: U(datalo(1):datahi(1), datalo(2):datahi(2), Ncomp)
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2), Ncomp)
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2), Ncomp)
  real(amrex_real), intent(in)    :: dt, dx(2)
  real(amrex_real), intent(inout) :: g

  ! local variables
  integer i,j,k
  real(amrex_real) S_upwind, S_downwind, S, r, ph, f_p, f_m
  real(amrex_real), parameter :: half = 0.5d0

  real(amrex_real) Up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
  real(amrex_real) Um(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
  real(amrex_real) fp(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
  real(amrex_real) fm(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
  real(amrex_real) dxdt(2)

  dxdt = dx/dt

  !write(*,*) lo, hi, datalo, datahi, fxlo, fxhi
  !return

  ! x fluxes
  do k = 1, Ncomp
      do    j = lo(2), hi(2)
         do i = lo(1)-1, hi(1)+1
             S_upwind = U(i+1,j,k) - U(i,j,k)
             S_downwind = U(i,j,k) - U(i-1,j,k)
             S = half * (S_upwind + S_downwind)

             if (abs(S_downwind) > 1.0d-10) then
                 r = S_upwind / S_downwind
             else
                 r = 1.0d6
             end if

             call phi(r, ph)

             S = S * ph

             Up(i,j,k) = U(i,j,k) + half * S
             Um(i,j,k) = U(i,j,k) - half * S
         end do
     end do
 end do

 !write(*,*) g, dt
 !return

 call swe_flux(Up, datalo, datahi, fp, fxlo, fxhi, lo, hi, Ncomp, g, .true.)
 call swe_flux(Um, datalo, datahi, fm, fxlo, fxhi, lo, hi, Ncomp, g, .true.)

 do k = 1, Ncomp
      do    j = lo(2), hi(2)
         do i = lo(1), hi(1)
            f_p = half * (fp(i,j,k) + fm(i+1,j,k) + dxdt(1) * (Up(i,j,k) - Um(i+1,j,k)))
            f_m = half * (fp(i-1,j,k) + fm(i,j,k) + dxdt(1) * (Up(i-1,j,k) - Um(i,j,k)))

            fluxx(i,j,k) = -(f_p - f_m) / dx(1)
         end do
      end do
  end do

  do k = 1, Ncomp
      do    j = lo(2)-1, hi(2)+1
         do i = lo(1), hi(1)
             S_upwind = U(i,j+1,k) - U(i,j,k)
             S_downwind = U(i,j,k) - U(i,j-1,k)
             S = half * (S_upwind + S_downwind)

             if (abs(S_downwind) > 1.0d-10) then
                 r = S_upwind / S_downwind
             else
                 r = 1.0d6
             end if

             call phi(r, ph)

             S = S * ph

             Up(i,j,k) = U(i,j,k) + half * S
             Um(i,j,k) = U(i,j,k) - half * S
         end do
     end do
 end do

 call swe_flux(Up, datalo, datahi, fp, fylo, fyhi, lo, hi, Ncomp, g, .false.)
 call swe_flux(Um, datalo, datahi, fm, fylo, fyhi, lo, hi, Ncomp, g, .false.)

  ! y-fluxes
  do k = 1, Ncomp
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             f_p = half * (fp(i,j,k) + fm(i,j+1,k) + dxdt(2) * (Up(i,j,k) - Um(i,j+1,k)))
             f_m = half * (fp(i,j-1,k) + fm(i,j,k) + dxdt(2) * (Up(i,j-1,k) - Um(i,j,k)))

             fluxy(i,j,k) = -(f_p - f_m) / dx(2)
          end do
       end do
   end do

  !write(*,*) fluxx(5,5,1)

end subroutine compute_flux


subroutine update_data (lo, hi, Ncomp, dataold, polo, pohi, datanew, pnlo, pnhi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx, dt) bind(C, name="update_data")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), Ncomp, polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: dataold(polo(1):pohi(1),polo(2):pohi(2), Ncomp)
  real(amrex_real), intent(inout) :: datanew(pnlo(1):pnhi(1),pnlo(2):pnhi(2), Ncomp)
  real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2), Ncomp)
  real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2), Ncomp)
  real(amrex_real), intent(in)    :: dx(2)
  real(amrex_real), value         :: dt

  ! local variables
  integer i,j,k
  real(amrex_real) :: dtdx(2)

  dtdx = dt/dx
  do        k = 1, Ncomp
      do    j = lo(2), hi(2)
         do i = lo(1), hi(1)

            datanew(i,j,k) = dataold(i,j,k) &
                 + dt * (fluxx(i,j  ,k) + fluxy(i  ,j,k))

         end do
      end do
  end do

end subroutine update_data
