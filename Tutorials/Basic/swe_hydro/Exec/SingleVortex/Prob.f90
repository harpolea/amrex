
subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  double precision, intent(in) :: problo(*), probhi(*)

  ! nothing needs to be done here

end subroutine amrex_probinit


subroutine initdata(level, time, lo, hi, &
     phi, phi_lo, phi_hi, &
     dx, prob_lo, Ncomp) bind(C, name="initdata")

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), phi_lo(3), phi_hi(3), Ncomp
  double precision, intent(in) :: time
  double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1), &
                                        phi_lo(2):phi_hi(2), &
                                        phi_lo(3):phi_hi(3), Ncomp)
  double precision, intent(in) :: dx(3), prob_lo(3)

  integer          :: dm
  integer          :: i,j,k
  double precision :: x,y,z,r2,alpha
  double precision, parameter ::  alpha0 = 0.9899494937d0, M = 1.0d0, R = 100.0d0

  if (phi_lo(3) .eq. 0 .and. phi_hi(3) .eq. 0) then
     dm = 2
  else
     dm = 3
  end if

  phi(:,:,:,:) = 0.0d0

  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
        do i=lo(1),hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

           if ( dm.eq. 2) then
              r2 = ((x-0.5d0)**2 + (y-0.75d0)**2) / 0.01d0
           else
              r2 = ((x-0.5d0)**2 + (y-0.75d0)**2 + (z-0.5d0)**2) / 0.01d0
           end if
           alpha = sqrt(1.0d0 - 2.0d0 * M / (1.0d0 + exp(-r2) + R)) !alpha0 + M * (1.0d0 + exp(-r2)) / (R**2 * alpha0)
           !phi(i,j,k,1) = 1.d0 + exp(-r2)
           if (Ncomp == 4) then ! compressible
               phi(i,j,k,1) = 1.d0 + exp(-r2)
               phi(i,j,k,4) = phi(i,j,k,1)
           else
               phi(i,j,k,1) = -log(alpha)
           end if
        end do
     end do
  end do
  !$omp end parallel do

end subroutine initdata