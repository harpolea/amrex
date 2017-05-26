subroutine init_data(U, philo, phihi, lo, hi, Ncomp, nlayers, dx, prob_lo, prob_hi, gamma_up, glo, ghi, beta, blo, bhi, alpha, alo, ahi) bind(C, name="init_data")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2), glo(2), ghi(2), blo(2), bhi(2), alo(2), ahi(2), Ncomp, nlayers
  real(amrex_real), intent(inout) :: U(philo(1):phihi(1),philo(2):phihi(2), Ncomp*nlayers)
  real(amrex_real), intent(in   ) :: dx(2)
  real(amrex_real), intent(in   ) :: prob_lo(2)
  real(amrex_real), intent(in   ) :: prob_hi(2)
  real(amrex_real), intent(inout) :: beta(blo(1):bhi(1),blo(2):bhi(2), 3*nlayers)
  real(amrex_real), intent(inout) :: gamma_up(glo(1):ghi(1),glo(2):ghi(2), 9*nlayers)
  real(amrex_real), intent(inout) :: alpha(alo(1):ahi(1),alo(2):ahi(2), nlayers)

  integer          :: i,j,k,l
  double precision :: x,y,r2, R

  do j = lo(2), hi(2)
     y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
     do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

        !r2 = (x**2 + y**2) / 0.01d0

        !U(i,j,1) = 1.d0 + 0.1d0*exp(-r2)

        r2 = (x**2 + y**2)

        do l = 0, nlayers-1
            if (r2 < 0.2**2) then
                U(i,j,l*Ncomp+1) = 4.0d0 - l*1.5d0
            else
                U(i,j,l*Ncomp+1) = 3.0d0 - l*1.5d0
            end if
        end do

        ! do velocities
        do l = 0, nlayers-1
            do k = 2, Ncomp
                U(i,j,l*Ncomp+k) = 0.d0
            end do
        end do

     end do
  end do

  do l = 0, nlayers-1
      alpha(:,:, l+1) = exp(-U(i,j,l*Ncomp+1))
      gamma_up(:,:,l*Ncomp+1) = 1.0d0
      gamma_up(:,:,l*Ncomp+5) = 1.0d0
      gamma_up(:,:,l*Ncomp+9) = exp(-2.0d0 * U(i,j,l*Ncomp+1))
  end do

end subroutine init_data
