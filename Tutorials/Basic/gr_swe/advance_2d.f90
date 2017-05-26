subroutine phi(r, ph, size)
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

subroutine swe_flux(U, f, lo, hi, Ncomp, nlayers, g, x_dir)
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: Ncomp, nlayers
    real(amrex_real), intent(in) :: g
    integer, intent(in) :: lo(2), hi(2)
    real(amrex_real), intent(in)  :: U(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp*nlayers)
    real(amrex_real), intent(out) :: f(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp*nlayers)
    logical, intent(in) :: x_dir

    integer i, j, l, m
    real(amrex_real), parameter :: half = 0.5d0

    do l = 0, nlayers-1
        if (x_dir) then
            do    j = lo(2), hi(2)
               do i = lo(1)-1, hi(1)+1
                  f(i,j,l*Ncomp+1) =  U(i,j,l*Ncomp+2)
                  f(i,j,l*Ncomp+2) = U(i,j,l*Ncomp+2)**2/U(i,j,l*Ncomp+1) + half * g * U(i,j,l*Ncomp+1)**2
                  f(i,j,l*Ncomp+3) =  U(i,j,l*Ncomp+2) * U(i,j,l*Ncomp+3)
               end do
            end do
        else
            do    j = lo(2)-1, hi(2)+1
               do i = lo(1), hi(1)
                  f(i,j,l*Ncomp+1) =  U(i,j,l*Ncomp+3)
                  f(i,j,l*Ncomp+2) = U(i,j,l*Ncomp+2) * U(i,j,l*Ncomp+3)
                  f(i,j,l*Ncomp+3) =  U(i,j,l*Ncomp+3)**2/U(i,j,l*Ncomp+1) + half * g * U(i,j,l*Ncomp+1)**2
               end do
            end do
        end if
    end do

end subroutine swe_flux

subroutine compute_flux (U, datalo, datahi, lo, hi, Ncomp, nlayers,&
     fluxx, fluxy, flo, fhi, dx, g, dt)! bind(C, name="compute_flux")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), Ncomp, nlayers, datalo(2), datahi(2), flo(2), fhi(2)
  real(amrex_real), intent(in)    :: U(datalo(1):datahi(1), datalo(2):datahi(2), Ncomp*nlayers)
  real(amrex_real), intent(inout) :: fluxx( flo(1): fhi(1), flo(2): fhi(2), Ncomp*nlayers)
  real(amrex_real), intent(inout) :: fluxy( flo(1): fhi(1), flo(2): fhi(2), Ncomp*nlayers)
  real(amrex_real), intent(in)    :: dt, dx(2)
  real(amrex_real), intent(inout) :: g

  ! local variables
  integer i,j,k
  real(amrex_real) S_upwind(Ncomp*nlayers), S_downwind(Ncomp*nlayers), S(Ncomp*nlayers), r(Ncomp*nlayers), ph(Ncomp*nlayers), f_p(Ncomp*nlayers), f_m(Ncomp*nlayers)
  real(amrex_real), parameter :: half = 0.5d0, alpha = 0.5d0

  real(amrex_real) Up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp*nlayers)
  real(amrex_real) Um(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp*nlayers)
  real(amrex_real) fp(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp*nlayers)
  real(amrex_real) fm(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp*nlayers)
  real(amrex_real) dxdt(2)

  dxdt = dx/dt

  ! x fluxes
  do    j = lo(2), hi(2)
     do i = lo(1)-1, hi(1)+1
         S_upwind = U(i+1,j,:) - U(i,j,:)
         S_downwind = U(i,j,:) - U(i-1,j,:)
         S = half * (S_upwind + S_downwind)

         r = S_upwind / S_downwind

         call phi(r, ph, Ncomp*nlayers)

         S = S * ph

         Up(i,j,:) = U(i,j,:) + half * S
         Um(i,j,:) = U(i,j,:) - half * S
     end do
 end do

 !write(*,*) g, dt
 !return

 call swe_flux(Up, fp, lo, hi, Ncomp, nlayers, g, .true.)
 call swe_flux(Um, fm, lo, hi, Ncomp, nlayers, g, .true.)

  do    j = lo(2), hi(2)
     do i = lo(1), hi(1)
        f_p = half * (fp(i,j,:) + fm(i+1,j,:) + alpha * dxdt(1) * (Up(i,j,:) - Um(i+1,j,:)))
        f_m = half * (fp(i-1,j,:) + fm(i,j,:) + alpha * dxdt(1) * (Up(i-1,j,:) - Um(i,j,:)))

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

         call phi(r, ph, Ncomp*nlayers)

         S = S * ph

         Up(i,j,:) = U(i,j,:) + half * S
         Um(i,j,:) = U(i,j,:) - half * S
     end do
 end do

 call swe_flux(Up, fp, lo, hi, Ncomp, nlayers, g, .false.)
 call swe_flux(Um, fm, lo, hi, Ncomp, nlayers, g, .false.)

   do    j = lo(2), hi(2)
      do i = lo(1), hi(1)
         f_p = half * (fp(i,j,:) + fm(i,j+1,:) + alpha * dxdt(2) * (Up(i,j,:) - Um(i,j+1,:)))
         f_m = half * (fp(i,j-1,:) + fm(i,j,:) + alpha * dxdt(2) * (Up(i,j-1,:) - Um(i,j,:)))

         fluxy(i,j,:) = -(f_p - f_m)
      end do
   end do

  !write(*,*) fluxx(5,5,1)

end subroutine compute_flux

subroutine buoyancy(U, datalo, datahi, lo, hi, Ncomp, nlayers,&
     buoy, flo, fhi, dx, g)

     use amrex_fort_module, only : amrex_real
     implicit none

     integer lo(2), hi(2), Ncomp, nlayers, datalo(2), datahi(2), flo(2), fhi(2)
     real(amrex_real), intent(in)    :: U(datalo(1):datahi(1), datalo(2):datahi(2), Ncomp*nlayers)
     real(amrex_real), intent(inout) :: buoy( flo(1): fhi(1), flo(2): fhi(2), Ncomp*nlayers)
     real(amrex_real), intent(in)    :: dx(2)
     real(amrex_real), intent(inout) :: g

     integer i,j,l,m

     buoy(:,:,:) = 0.0d0

     do l = 0, nlayers-1
         do    j = lo(2), hi(2)
            do i = lo(1), hi(1)
                do m = l+1, nlayers-1
                   buoy(i,j,l*Ncomp+2) = buoy(i,j,l*Ncomp+2) + g * U(i,j,l*Ncomp+1) * (U(i+1,j,m*Ncomp+1) - U(i-1,j,m*Ncomp+1)) / (2.0d0 * dx(1))

                   buoy(i,j,l*Ncomp+3) = buoy(i,j,l*Ncomp+3) + g * U(i,j,l*Ncomp+1) * (U(i,j+1,m*Ncomp+1) - U(i,j-1,m*Ncomp+1)) / (2.0d0 * dx(2))
                end do

                ! TODO: need to define densities somewhere
                do m = 0, l-1
                   buoy(i,j,l*Ncomp+2) = buoy(i,j,l*Ncomp+2) + g * U(i,j,l*Ncomp+1) * (U(i+1,j,m*Ncomp+1) - U(i-1,j,m*Ncomp+1)) / (2.0d0 * dx(1))

                   buoy(i,j,l*Ncomp+3) = buoy(i,j,l*Ncomp+3) + g * U(i,j,l*Ncomp+1) * (U(i,j+1,m*Ncomp+1) - U(i,j-1,m*Ncomp+1)) / (2.0d0 * dx(2))
                end do
            end do
         end do
     end do

end subroutine buoyancy


subroutine update_data (lo, hi, Ncomp, nlayers, dataold, polo, pohi, &
     datanew, pnlo, pnhi, dx, g, dt) bind(C, name="update_data")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), Ncomp, nlayers, polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: dataold(polo(1):pohi(1),polo(2):pohi(2), Ncomp*nlayers)
  real(amrex_real), intent(inout) :: datanew(pnlo(1):pnhi(1),pnlo(2):pnhi(2), Ncomp*nlayers)
  real(amrex_real), intent(in)    :: dx(2)
  real(amrex_real), intent(inout) :: g
  real(amrex_real), value         :: dt

  ! local variables
  integer i,j
  real(amrex_real) :: dtdx(2)

  real(amrex_real) fluxx(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, Ncomp*nlayers)
  real(amrex_real) fluxy(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, Ncomp*nlayers)

  dtdx = dt/dx

  call compute_flux (dataold, polo, pohi, lo-4, hi+4, Ncomp, nlayers, &
       fluxx, fluxy, lo-4, hi+4, dx, g, dt)

   do   j = lo(2)-4, hi(2)+4
     do i = lo(1)-4, hi(1)+4
        datanew(i,j,:) = dataold(i,j,:) &
             + fluxx(i,j,:) * dtdx(1) + fluxy(i,j,:) * dtdx(2)
     end do
  end do

  call buoyancy(datanew, polo, pohi, lo-4, hi+4, Ncomp, nlayers, &
       fluxx, lo-4, hi+4, dx, g)

   do   j = lo(2)-4, hi(2)+4
     do i = lo(1)-4, hi(1)+4
        datanew(i,j,:) = datanew(i,j,:) + dt * fluxx(i,j,:)
     end do
  end do

  call compute_flux (datanew, pnlo, pnhi, lo-2, hi+2, Ncomp, nlayers, &
       fluxx, fluxy, lo-4, hi+4, dx, g, dt)

   do   j = lo(2)-2, hi(2)+2
     do i = lo(1)-2, hi(1)+2
        datanew(i,j,:) = 0.25d0 * (3.d0 * dataold(i,j,:) + &
             datanew(i,j,:) + fluxx(i,j,:) * dtdx(1) + fluxy(i,j,:) * dtdx(2))
     end do
  end do

  call buoyancy(datanew, polo, pohi, lo-2, hi+2, Ncomp, nlayers, &
       fluxx, lo-4, hi+4, dx, g)

   do   j = lo(2)-2, hi(2)+2
     do i = lo(1)-2, hi(1)+2
        datanew(i,j,:) = datanew(i,j,:) + 0.25d0 * dt * fluxx(i,j,:)
     end do
  end do

  call compute_flux (datanew, pnlo, pnhi, lo, hi, Ncomp, nlayers, &
       fluxx, fluxy, lo-4, hi+4, dx, g, dt)

  do    j = lo(2), hi(2)
     do i = lo(1), hi(1)
        datanew(i,j,:) = 1.0d0 / 3.0d0 * (dataold(i,j,:) + &
             2.0d0 * (datanew(i,j,:) + fluxx(i,j,:) * dtdx(1) + fluxy(i,j,:) * dtdx(2)))
     end do
  end do

  call buoyancy(datanew, polo, pohi, lo, hi, Ncomp, nlayers, &
       fluxx, lo-4, hi+4, dx, g)

   do   j = lo(2), hi(2)
     do i = lo(1), hi(1)
        datanew(i,j,:) = datanew(i,j,:) + dt * fluxx(i,j,:) * 2.0d0/3.0d0
     end do
  end do

end subroutine update_data
