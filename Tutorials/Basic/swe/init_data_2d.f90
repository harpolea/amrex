subroutine init_data(U, philo, phihi, lo, hi, Ncomp, dx, prob_lo, prob_hi) bind(C, name="init_data")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2), Ncomp
  real(amrex_real), intent(inout) :: U(philo(1):phihi(1),philo(2):phihi(2), Ncomp)
  real(amrex_real), intent(in   ) :: dx(2)
  real(amrex_real), intent(in   ) :: prob_lo(2)
  real(amrex_real), intent(in   ) :: prob_hi(2)

  integer          :: i,j,k
  double precision :: x,y,r2

  do j = lo(2), hi(2)
     y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
     do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

        !r2 = (x**2 + y**2) / 0.01d0

        !U(i,j,1) = 1.d0 + 0.1d0*exp(-r2)

        r2 = (x**2 + y**2)

        if (r2 < 0.2**2) then
            U(i,j,1) = 2.0d0
        else
            U(i,j,1) = 1.0d0
        end if

        ! do velocities
        do k = 2, Ncomp
            U(i,j,k) = 0.d0
        end do

     end do
  end do

end subroutine init_data
