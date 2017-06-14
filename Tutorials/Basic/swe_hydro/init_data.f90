subroutine init_swe_data(U, philo, phihi, lo, hi, Ncomp, dx, prob_lo, prob_hi, gamma_up, glo, ghi, beta, blo, bhi, alpha, alo, ahi, alpha0, M, R) bind(C, name="init_swe_data")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3), philo(3), phihi(3), Ncomp,  glo(3), ghi(3), blo(3), bhi(3), alo(3), ahi(3)
  real(amrex_real), intent(inout) :: U(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3), Ncomp)
  real(amrex_real), intent(in   ) :: dx(3)
  real(amrex_real), intent(in   ) :: prob_lo(3)
  real(amrex_real), intent(in   ) :: prob_hi(3)
  real(amrex_real), intent(inout) :: beta(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3), 3)
  real(amrex_real), intent(inout) :: gamma_up(glo(1):ghi(1),glo(2):ghi(2), glo(3):ghi(3),9)
  real(amrex_real), intent(inout) :: alpha(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
  real(amrex_real), value         :: alpha0, M, R

  integer          :: i,j,k,l
  double precision :: x,y,z,r2

  do k = lo(3)-6, hi(3)+6
      z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
      do j = lo(2)-6, hi(2)+6
         y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
         do i = lo(1)-6, hi(1)+6
            x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

            !r2 = (x**2 + y**2) / 0.01d0

            !U(i,j,1) = 1.d0 + 0.1d0*exp(-r2)

            r2 = (x**2 + y**2 + z**2)

            if (r2 < 0.2**2) then
                U(i,j,k,1) = 1.5d0
                U(i,j,k,5) = 1.5d0
            else
                U(i,j,k,1) = 1.0d0
                U(i,j,k,5) = 1.0d0
            end if

            ! do velocities
            do l = 2, Ncomp-1
                U(i,j,k,l) = 0.d0
            end do

         end do
      end do
      alpha(:,:,k) = alpha0 + M * z / (R**2 * alpha0)
  end do


      !write(*,*) "alpha", alpha(lo(1),lo(2),lo(3))

  gamma_up(:,:,:,1) = 1.0d0
  gamma_up(:,:,:,5) = 1.0d0

  gamma_up(:,:,:,9) = alpha(:,:,:)**2


end subroutine init_swe_data

subroutine init_comp_data(U, philo, phihi, lo, hi, Ncomp, dx, prob_lo, prob_hi, gamma_up, glo, ghi, beta, blo, bhi, alpha, alo, ahi, alpha0, M, R) bind(C, name="init_comp_data")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3), philo(3), phihi(3), Ncomp,  glo(3), ghi(3), blo(3), bhi(3), alo(3), ahi(3)
  real(amrex_real), intent(inout) :: U(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3), Ncomp)
  real(amrex_real), intent(in   ) :: dx(3)
  real(amrex_real), intent(in   ) :: prob_lo(3)
  real(amrex_real), intent(in   ) :: prob_hi(3)
  real(amrex_real), intent(inout) :: beta(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3), 3)
  real(amrex_real), intent(inout) :: gamma_up(glo(1):ghi(1),glo(2):ghi(2), glo(3):ghi(3),9)
  real(amrex_real), intent(inout) :: alpha(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
  real(amrex_real), value         :: alpha0, M, R

  integer          :: i,j,k,l
  double precision :: x,y,z,r2

  do k = lo(3)-6, hi(3)+6
      z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
      do j = lo(2)-6, hi(2)+6
         y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
         do i = lo(1)-6, hi(1)+6
            x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

            !r2 = (x**2 + y**2) / 0.01d0

            !U(i,j,1) = 1.d0 + 0.1d0*exp(-r2)

            r2 = (x**2 + y**2 + z**2)

            if (r2 < 0.2**2) then
                U(i,j,k,1) = 1.5d0
                U(i,j,k,5) = 1.5d0
            else
                U(i,j,k,1) = 1.0d0
                U(i,j,k,5) = 1.0d0
            end if

            ! do velocities
            do l = 2, Ncomp-1
                U(i,j,k,l) = 0.d0
            end do

         end do
      end do
      alpha(:,:,k) = alpha0 + M * z / (R**2 * alpha0)
  end do


      !write(*,*) "alpha", alpha(lo(1),lo(2),lo(3))

  gamma_up(:,:,:,1) = 1.0d0
  gamma_up(:,:,:,5) = 1.0d0

  gamma_up(:,:,:,9) = alpha(:,:,:)**2


end subroutine init_comp_data
