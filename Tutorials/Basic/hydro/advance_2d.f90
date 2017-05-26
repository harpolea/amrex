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
        if (r(i) >= 1.0d0) then
            ph(i) = min(2.0d0, min(r(i), 2.0d0 / (1.0d0 + r(i))))
        elseif ( r(i) >= 0.5d0 ) then
            ph(i) = 1.0d0
        elseif ( r(i) > 0.0d0 ) then
            ph(i) = 2.0d0 * r(i)
        end if
    end do

end subroutine phi


subroutine comp_flux(U, f, lo, hi, Ncomp, x_dir, gamma)
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: Ncomp
    integer, intent(in) :: lo(2), hi(2)
    real(amrex_real), intent(in)  :: U(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
    real(amrex_real), intent(out) :: f(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
    logical, intent(in) :: x_dir
    real(amrex_real), intent(in)  :: gamma

    integer i, j
    real(amrex_real), parameter :: half = 0.5d0
    real(amrex_real) p(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1)

    p = (gamma - 1.d0) * (U(:,:,4) - half * (U(:,:,2)**2 + U(:,:,3)**2) / U(:,:,1))

    if (x_dir) then
        f(:,:,1) = U(:,:,2)
        f(:,:,2) = U(:,:,2)**2 / U(:,:,1) + p
        f(:,:,3) = U(:,:,2) * U(:,:,3) / U(:,:,1)
        f(:,:,4) = (U(:,:,4) + p) * U(:,:,2) / U(:,:,1)
    else
        f(:,:,1) = U(:,:,3)
        f(:,:,2) = U(:,:,2) * U(:,:,3) / U(:,:,1)
        f(:,:,3) = U(:,:,3)**2 / U(:,:,1) + p
        f(:,:,4) = (U(:,:,4) + p) * U(:,:,3) / U(:,:,1)
    end if

end subroutine comp_flux

subroutine compute_flux (U, datalo, datahi, lo, hi, Ncomp,&
     fluxx, fluxy, flo, fhi, dx, dt, gamma)! bind(C, name="compute_flux")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), Ncomp, datalo(2), datahi(2), flo(2), fhi(2)
  real(amrex_real), intent(in)    :: U(datalo(1):datahi(1), datalo(2):datahi(2), Ncomp)
  real(amrex_real), intent(inout) :: fluxx( flo(1): fhi(1), flo(2): fhi(2), Ncomp)
  real(amrex_real), intent(inout) :: fluxy( flo(1): fhi(1), flo(2): fhi(2), Ncomp)
  real(amrex_real), intent(in)    :: dt, dx(2), gamma

  ! local variables
  integer i,j,k
  real(amrex_real) S_upwind(Ncomp), S_downwind(Ncomp), S(Ncomp), r(Ncomp), ph(Ncomp), f_p(Ncomp), f_m(Ncomp)
  real(amrex_real), parameter :: half = 0.5d0, alph = 0.5d0

  real(amrex_real) Up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
  real(amrex_real) Um(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
  real(amrex_real) fp(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
  real(amrex_real) fm(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
  real(amrex_real) dxdt(2)

  dxdt = dx/dt

  ! x fluxes
  do    j = lo(2), hi(2)
     do i = lo(1)-1, hi(1)+1
         S_upwind = U(i+1,j,:) - U(i,j,:)
         S_downwind = U(i,j,:) - U(i-1,j,:)
         S = half * (S_upwind + S_downwind)

         r = S_upwind / S_downwind

         call phi(r, ph, Ncomp)

         S = S * ph

         Up(i,j,:) = U(i,j,:) + half * S
         Um(i,j,:) = U(i,j,:) - half * S
     end do
 end do

 !write(*,*) g, dt
 !return

 call comp_flux(Up, fp, lo, hi, Ncomp, .true., gamma)
 call comp_flux(Um, fm, lo, hi, Ncomp, .true., gamma)

  do    j = lo(2), hi(2)
     do i = lo(1), hi(1)
        f_p = half * (fp(i,j,:) + fm(i+1,j,:) + alph * dxdt(1) * (Up(i,j,:) - Um(i+1,j,:)))
        f_m = half * (fp(i-1,j,:) + fm(i,j,:) + alph * dxdt(1) * (Up(i-1,j,:) - Um(i,j,:)))

        fluxx(i,j,:) = -(f_p - f_m)
     end do
  end do

  ! y fluxes
  do    j = lo(2)-1, hi(2)+1
     do i = lo(1), hi(1)
         S_upwind = U(i,j+1,:) - U(i,j,:)
         S_downwind = U(i,j,:) - U(i,j-1,:)
         S = half * (S_upwind + S_downwind)

         r = S_upwind / S_downwind

         call phi(r, ph, Ncomp)

         S = S * ph

         Up(i,j,:) = U(i,j,:) + half * S
         Um(i,j,:) = U(i,j,:) - half * S
     end do
 end do

 call comp_flux(Up, fp, lo, hi, Ncomp, .false., gamma)
 call comp_flux(Um, fm, lo, hi, Ncomp, .false., gamma)

   do    j = lo(2), hi(2)
      do i = lo(1), hi(1)
         f_p = half * (fp(i,j,:) + fm(i,j+1,:) + alph * dxdt(2) * (Up(i,j,:) - Um(i,j+1,:)))
         f_m = half * (fp(i,j-1,:) + fm(i,j,:) + alph * dxdt(2) * (Up(i,j-1,:) - Um(i,j,:)))

         fluxy(i,j,:) = -(f_p - f_m)
      end do
   end do

  !write(*,*) fluxx(5,5,1)

end subroutine compute_flux


subroutine update_data (lo, hi, Ncomp,dataold, polo, pohi, &
     datanew, pnlo, pnhi, dx, dt, gamma) bind(C, name="update_data")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), Ncomp, polo(2), pohi(2), pnlo(2), pnhi(2)
  real(amrex_real), intent(in)    :: dataold(polo(1):pohi(1),polo(2):pohi(2), Ncomp)
  real(amrex_real), intent(inout) :: datanew(pnlo(1):pnhi(1),pnlo(2):pnhi(2), Ncomp)
  real(amrex_real), intent(in)    :: dx(2), gamma
  real(amrex_real), value         :: dt

  ! local variables
  integer i,j
  real(amrex_real) :: dtdx(2)

  real(amrex_real) fluxx(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, Ncomp)
  real(amrex_real) fluxy(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, Ncomp)

  dtdx = dt/dx

  call compute_flux (dataold, polo, pohi, lo-4, hi+4, Ncomp, &
       fluxx, fluxy, lo-4, hi+4, dx, dt, gamma)

   do   j = lo(2)-4, hi(2)+4
     do i = lo(1)-4, hi(1)+4
        datanew(i,j,:) = dataold(i,j,:) &
             + fluxx(i,j,:) * dtdx(1) + fluxy(i,j,:) * dtdx(2)
     end do
  end do


  call compute_flux (datanew, pnlo, pnhi, lo-2, hi+2, Ncomp, &
       fluxx, fluxy, lo-4, hi+4, dx, dt, gamma)

   do   j = lo(2)-2, hi(2)+2
     do i = lo(1)-2, hi(1)+2
        datanew(i,j,:) = 0.25d0 * (3.d0 * dataold(i,j,:) + &
             datanew(i,j,:) + fluxx(i,j,:) * dtdx(1) + fluxy(i,j,:) * dtdx(2))
     end do
  end do

  call compute_flux (datanew, pnlo, pnhi, lo, hi, Ncomp, &
       fluxx, fluxy, lo-4, hi+4, dx, dt, gamma)

  do    j = lo(2), hi(2)
     do i = lo(1), hi(1)
        datanew(i,j,:) = 1.0d0 / 3.0d0 * (dataold(i,j,:) + &
             2.0d0 * (datanew(i,j,:) + fluxx(i,j,:) * dtdx(1) + fluxy(i,j,:) * dtdx(2)))
     end do
  end do

end subroutine update_data
