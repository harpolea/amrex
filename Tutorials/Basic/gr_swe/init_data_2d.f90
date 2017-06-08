subroutine init_data(U, philo, phihi, lo, hi, Ncomp, nlayers, dx, prob_lo, prob_hi, gamma_up, glo, ghi, beta, blo, bhi, alpha, alo, ahi, alpha0, M, R) bind(C, name="init_data")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2), glo(2), ghi(2), blo(2), bhi(2), alo(2), ahi(2), Ncomp, nlayers
  real(amrex_real), intent(inout) :: U(philo(1):phihi(1),philo(2):phihi(2), Ncomp*nlayers)
  real(amrex_real), intent(in   ) :: dx(2), alpha0, M, R
  real(amrex_real), intent(in   ) :: prob_lo(2)
  real(amrex_real), intent(in   ) :: prob_hi(2)
  real(amrex_real), intent(inout) :: beta(blo(1):bhi(1),blo(2):bhi(2), 3*nlayers)
  real(amrex_real), intent(inout) :: gamma_up(glo(1):ghi(1),glo(2):ghi(2), 9*nlayers)
  real(amrex_real), intent(inout) :: alpha(alo(1):ahi(1),alo(2):ahi(2), nlayers)

  integer          :: i,j,k,l
  double precision :: x,y,r2,h

  do j = philo(2), phihi(2)
     y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
     do i = philo(1), phihi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

        !r2 = (x**2 + y**2) / 0.01d0

        !U(i,j,1) = 1.d0 + 0.1d0*exp(-r2)

        r2 = (x**2 + y**2)

        do l = 0, nlayers-1
            h = (nlayers - 1.0d0 - l) * 1.5d0 + 1.0d0

            if (r2 < 0.2**2) then
                h = h + 1.0d0
            end if

            alpha(i,j,l+1) = alpha0 + M * h / (R**2 * alpha0)
            U(i,j,l*Ncomp+1) = -log(alpha(i,j,l+1))
            U(i,j,l*Ncomp+4) = h
        end do

     end do
  end do

  ! do velocities
  do l = 0, nlayers-1
      do k = 2, Ncomp
          U(:,:,l*Ncomp+k) = 0.d0
      end do
  end do

  do l = 0, nlayers-1
      ! make sure alpha is defined in ghost cells
      h = (nlayers - 1.0d0 - l) * 0.5d0 + 1.0d0
      alpha(alo(1):lo(1)-1, alo(2):lo(2)-1, l+1) = &
            alpha0 + M * h / (R**2 * alpha0)
      alpha(hi(1)+1:ahi(1), hi(2)+1:ahi(2), l+1) = &
            alpha0 + M * h / (R**2 * alpha0)

      gamma_up(:,:,l*9+1) = 1.0d0
      gamma_up(:,:,l*9+5) = 1.0d0
      gamma_up(:,:,l*9+9) = alpha(:,:,l+1)**2
  end do

end subroutine init_data
