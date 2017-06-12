subroutine W_swe(q, lo, hi, nlayers, Ncomp, gamma_up, glo, ghi, W)
    ! Calculate lorentz factor
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: lo(2), hi(2), nlayers, Ncomp, glo(2), ghi(2)
    real(amrex_real), intent(in) :: q(lo(1):hi(1), lo(2):hi(2), nlayers*Ncomp)
    real(amrex_real), intent(in) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), 9*nlayers)
    real(amrex_real), intent(out) :: W(lo(1):hi(1), lo(2):hi(2), nlayers)

    integer i,j,l

    do l = 0, nlayers - 1
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                W(i,j,l+1) = sqrt((q(i,j,l*Ncomp+2)**2 * gamma_up(i,j,l*9+1)+&
                      2.0d0 * q(i,j,l*Ncomp+2) * q(i,j,l*Ncomp+3) * &
                      gamma_up(i,j,l*9+2) + q(i,j,l*Ncomp+3)**2 * &
                      gamma_up(i,j,l*9+5)) / q(i,j,l*Ncomp+1)**2 + 1.0d0)
                ! nan check
                if (W(i,j,l+1) /= W(i,j,l+1)) then
                    W(i,j,l+1) = 1.0d0
                end if
          end do
      end do
    end do

end subroutine W_swe

subroutine swe_flux(U, f, lo, hi, Ncomp, nlayers, x_dir,&
                    gamma_up, glo, ghi, beta, alpha)
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: Ncomp, nlayers
    integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2)
    real(amrex_real), intent(in)  :: U(lo(1):hi(1), lo(2):hi(2), Ncomp*nlayers)
    real(amrex_real), intent(out) :: f(lo(1):hi(1), lo(2):hi(2), Ncomp*nlayers)
    logical, intent(in) :: x_dir
    real(amrex_real), intent(in) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), 3*3*nlayers)
    real(amrex_real), intent(in) :: beta(glo(1):ghi(1), glo(2):ghi(2), 3*nlayers)
    real(amrex_real), intent(in) :: alpha(glo(1):ghi(1), glo(2):ghi(2), nlayers)

    integer i, j, l
    real(amrex_real) v(2)
    real(amrex_real) W(lo(1):hi(1), lo(2):hi(2), nlayers)
    real(amrex_real), parameter :: half = 0.5d0


    call W_swe(U, lo, hi, nlayers, Ncomp, gamma_up, glo, ghi, W)

    do l = 0, nlayers-1
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)

                v(1) = U(i,j,l*Ncomp+2) / (U(i,j,l*Ncomp+1) * W(i,j,l+1))
                v(2) = U(i,j,l*Ncomp+3) / (U(i,j,l*Ncomp+1) * W(i,j,l+1))

                if (x_dir) then
                  f(i,j,l*Ncomp+1) = U(i,j,l*Ncomp+1) * &
                        (v(1)*gamma_up(i,j,l*9+1) + v(2)*gamma_up(i,j,l*9+2) -&
                         beta(i,j,l*3+1) / alpha(i,j,l+1))
                  f(i,j,l*Ncomp+2) = U(i,j,l*Ncomp+2) * &
                        (v(1)*gamma_up(i,j,l*9+1) + v(2)*gamma_up(i,j,l*9+2) -&
                         beta(i,j,l*3+1) / alpha(i,j,l+1)) + half * U(i,j,l*Ncomp+1)**2 / W(i,j,l+1)**2
                  f(i,j,l*Ncomp+3) = U(i,j,l*Ncomp+3) * &
                        (v(1)*gamma_up(i,j,l*9+1) + v(2)*gamma_up(i,j,l*9+2) -&
                         beta(i,j,l*3+1) / alpha(i,j,l+1))
                else
                  f(i,j,l*Ncomp+1) = U(i,j,l*Ncomp+1) * &
                        (v(1)*gamma_up(i,j,l*9+2) + v(2)*gamma_up(i,j,l*9+5) -&
                         beta(i,j,l*3+2) / alpha(i,j,l+1))
                  f(i,j,l*Ncomp+2) = U(i,j,l*Ncomp+2) * &
                        (v(1)*gamma_up(i,j,l*9+2) + v(2)*gamma_up(i,j,l*9+5) -&
                         beta(i,j,l*3+2) / alpha(i,j,l+1))
                  f(i,j,l*Ncomp+3) = U(i,j,l*Ncomp+3) * &
                        (v(1)*gamma_up(i,j,l*9+2) + v(2)*gamma_up(i,j,l*9+5) -&
                         beta(i,j,l*3+2) / alpha(i,j,l+1)) +  half * U(i,j,l*Ncomp+1)**2 / W(i,j,l+1)**2
                end if
            end do
        end do
    end do

end subroutine swe_flux

subroutine compute_swe_flux (U, datalo, datahi, lo, hi, Ncomp, nlayers,&
     fluxx, fluxy, flo, fhi, dx, dt, gamma_up, glo, ghi, beta, alpha)! bind(C, name="compute_flux")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), Ncomp, nlayers, datalo(2), datahi(2), flo(2), fhi(2), glo(2), ghi(2)
  real(amrex_real), intent(in)    :: U(datalo(1):datahi(1), datalo(2):datahi(2), Ncomp*nlayers)
  real(amrex_real), intent(inout) :: fluxx( flo(1): fhi(1), flo(2): fhi(2), Ncomp*nlayers)
  real(amrex_real), intent(inout) :: fluxy( flo(1): fhi(1), flo(2): fhi(2), Ncomp*nlayers)
  real(amrex_real), intent(in)    :: dt, dx(2)

  ! local variables
  integer i,j,k
  real(amrex_real) S_upwind(Ncomp*nlayers), S_downwind(Ncomp*nlayers), S(Ncomp*nlayers), r(Ncomp*nlayers), ph(Ncomp*nlayers), f_p(Ncomp*nlayers), f_m(Ncomp*nlayers)
  real(amrex_real), parameter :: half = 0.5d0, alph = 0.5d0
  real(amrex_real), intent(in) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2),3*3*nlayers)
  real(amrex_real), intent(in) :: beta(glo(1):ghi(1), glo(2):ghi(2),3*nlayers)
  real(amrex_real), intent(in) :: alpha(glo(1):ghi(1), glo(2):ghi(2), nlayers)

  real(amrex_real) Up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp*nlayers)
  real(amrex_real) Um(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp*nlayers)
  real(amrex_real) fp(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp*nlayers)
  real(amrex_real) fm(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp*nlayers)
  real(amrex_real) dxdt(2)

  dxdt = dx/dt

  ! x fluxes
  do    j = lo(2)-1, hi(2)+1
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

 call swe_flux(Up, fp, lo-1, hi+1, Ncomp, nlayers, .true., gamma_up, glo, ghi, beta, alpha)
 call swe_flux(Um, fm, lo-1, hi+1, Ncomp, nlayers, .true., gamma_up, glo, ghi, beta, alpha)

  do    j = lo(2), hi(2)
     do i = lo(1), hi(1)
        f_p = half * (fp(i,j,:) + fm(i+1,j,:) + &
            alph * dxdt(1) * (Up(i,j,:) - Um(i+1,j,:)))
        f_m = half * (fp(i-1,j,:) + fm(i,j,:) + &
            alph * dxdt(1) * (Up(i-1,j,:) - Um(i,j,:)))

        fluxx(i,j,:) = -(f_p - f_m)
     end do
  end do

  ! y fluxes
  do    j = lo(2)-1, hi(2)+1
     do i = lo(1)-1, hi(1)+1
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

 call swe_flux(Up, fp, lo-1, hi+1, Ncomp, nlayers, .false., gamma_up, glo, ghi, beta, alpha)
 call swe_flux(Um, fm, lo-1, hi+1, Ncomp, nlayers, .false., gamma_up, glo, ghi, beta, alpha)

   do    j = lo(2), hi(2)
      do i = lo(1), hi(1)
         f_p = half * (fp(i,j,:) + fm(i,j+1,:) + &
            alph * dxdt(2) * (Up(i,j,:) - Um(i,j+1,:)))
         f_m = half * (fp(i,j-1,:) + fm(i,j,:) + &
            alph * dxdt(2) * (Up(i,j-1,:) - Um(i,j,:)))

         fluxy(i,j,:) = -(f_p - f_m)
      end do
   end do

  !write(*,*) fluxx(5,5,1)

end subroutine compute_swe_flux

subroutine buoyancy(U, datalo, datahi, lo, hi, Ncomp, nlayers,&
     buoy, flo, fhi, dx)

     use amrex_fort_module, only : amrex_real
     implicit none

     integer lo(2), hi(2), Ncomp, nlayers, datalo(2), datahi(2), flo(2), fhi(2)
     real(amrex_real), intent(in)    :: U(datalo(1):datahi(1), datalo(2):datahi(2), Ncomp*nlayers)
     real(amrex_real), intent(inout) :: buoy( flo(1): fhi(1), flo(2): fhi(2), Ncomp*nlayers)
     real(amrex_real), intent(in)    :: dx(2)

     integer i,j,l,m

     buoy(:,:,:) = 0.0d0

     do l = 0, nlayers-1
         do    j = lo(2), hi(2)
            do i = lo(1), hi(1)
                do m = l+1, nlayers-1
                   buoy(i,j,l*Ncomp+2) = buoy(i,j,l*Ncomp+2) + &
                        U(i,j,l*Ncomp+1) * (U(i+1,j,m*Ncomp+1) - &
                        U(i-1,j,m*Ncomp+1)) / (2.0d0 * dx(1))

                   buoy(i,j,l*Ncomp+3) = buoy(i,j,l*Ncomp+3) + &
                        U(i,j,l*Ncomp+1) * (U(i,j+1,m*Ncomp+1) -&
                        U(i,j-1,m*Ncomp+1)) / (2.0d0 * dx(2))
                end do

                ! TODO: need to define densities somewhere
                do m = 0, l-1
                   buoy(i,j,l*Ncomp+2) = buoy(i,j,l*Ncomp+2) + &
                        U(i,j,l*Ncomp+1) * (U(i+1,j,m*Ncomp+1) - &
                        U(i-1,j,m*Ncomp+1)) / (2.0d0 * dx(1))

                   buoy(i,j,l*Ncomp+3) = buoy(i,j,l*Ncomp+3) + &
                        U(i,j,l*Ncomp+1) * (U(i,j+1,m*Ncomp+1) - &
                        U(i,j-1,m*Ncomp+1)) / (2.0d0 * dx(2))
                end do
            end do
         end do
     end do

end subroutine buoyancy

subroutine height(alpha, U, alpha0, M, R, lo, hi, nlayers, Ncomp)

    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) ::  lo(2), hi(2), nlayers, Ncomp
    real(amrex_real), intent(in)    :: alpha(lo(1):hi(1),lo(2):hi(2), nlayers)
    real(amrex_real), intent(inout)    :: U(lo(1):hi(1),lo(2):hi(2), Ncomp*nlayers)
    real(amrex_real), intent(in) :: alpha0, M, R

    integer l

    do l = 0, nlayers-1
        U(:,:,l*Ncomp+4) = alpha0 * R**2 / M * (alpha(:,:,l+1) - alpha0)
    end do


end subroutine height

subroutine update_swe_data (lo, hi, Ncomp, nlayers, dataold, polo, pohi, &
     datanew, pnlo, pnhi, dx, dt, gamma_up, glo, ghi, beta, blo, bhi, alpha, alo, ahi, alpha0, M, R) bind(C, name="update_swe_data")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), Ncomp, nlayers, polo(2), pohi(2), pnlo(2), pnhi(2), glo(2), ghi(2), blo(2), bhi(2), alo(2), ahi(2)
  real(amrex_real), intent(in)    :: dataold(polo(1):pohi(1),polo(2):pohi(2), Ncomp*nlayers)
  real(amrex_real), intent(inout) :: datanew(pnlo(1):pnhi(1),pnlo(2):pnhi(2), Ncomp*nlayers)
  real(amrex_real), intent(in)    :: dx(2), alpha0, M, R
  real(amrex_real), value         :: dt
  real(amrex_real), intent(inout) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2),3*3*nlayers)
  real(amrex_real), intent(in) :: beta(blo(1):bhi(1), blo(2):bhi(2),3*nlayers)
  real(amrex_real), intent(inout) :: alpha(alo(1):ahi(1), alo(2):ahi(2), nlayers)

  ! local variables
  integer i,j,k
  real(amrex_real) :: dtdx(2)

  real(amrex_real) fluxx(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, Ncomp*nlayers)
  real(amrex_real) fluxy(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, Ncomp*nlayers)

  dtdx = dt/dx

  call compute_swe_flux (dataold, polo, pohi, lo-4, hi+4, Ncomp, nlayers, &
       fluxx, fluxy, lo-4, hi+4, dx, dt, gamma_up, glo, ghi, beta, alpha)

   do   j = lo(2)-4, hi(2)+4
     do i = lo(1)-4, hi(1)+4
        datanew(i,j,:) = dataold(i,j,:) &
             + fluxx(i,j,:) * dtdx(1) + fluxy(i,j,:) * dtdx(2)
        do k = 0, nlayers-1
            alpha(i,j,k+1) = exp(-datanew(i,j,k*Ncomp+1))
            gamma_up(i,j,9*k+1) = alpha(i,j,k+1)**2
        end do
     end do
  end do

  call buoyancy(datanew, pnlo, pnhi, lo-4, hi+4, Ncomp, nlayers, &
       fluxx, lo-4, hi+4, dx)

   do       k = 0, nlayers-1
       do   j = lo(2)-4, hi(2)+4
         do i = lo(1)-4, hi(1)+4
            datanew(i,j,k*Ncomp+1:(k+1)*Ncomp) = &
                datanew(i,j,k*Ncomp+1:(k+1)*Ncomp) + &
                alpha(i,j,k+1) * dt * fluxx(i,j,k*Ncomp+1:(k+1)*Ncomp)
            alpha(i,j,k+1) = exp(-datanew(i,j,k*Ncomp+1))
            gamma_up(i,j,9*k+1) = alpha(i,j,k+1)**2
         end do
      end do
  end do

  call compute_swe_flux (datanew, pnlo, pnhi, lo-2, hi+2, Ncomp, nlayers, &
       fluxx, fluxy, lo-4, hi+4, dx, dt, gamma_up, glo, ghi, beta, alpha)

   do   j = lo(2)-2, hi(2)+2
     do i = lo(1)-2, hi(1)+2
        datanew(i,j,:) = 0.25d0 * (3.d0 * dataold(i,j,:) + &
             datanew(i,j,:) + fluxx(i,j,:) * dtdx(1) + fluxy(i,j,:) * dtdx(2))
         do k = 0, nlayers-1
             alpha(i,j,k+1) = exp(-datanew(i,j,k*Ncomp+1))
             gamma_up(i,j,9*k+1) = alpha(i,j,k+1)**2
         end do
     end do
  end do

  call buoyancy(datanew, pnlo, pnhi, lo-2, hi+2, Ncomp, nlayers, &
       fluxx, lo-4, hi+4, dx)

   do       k = 0, nlayers-1
       do   j = lo(2)-2, hi(2)+2
         do i = lo(1)-2, hi(1)+2
            datanew(i,j,k*Ncomp+1:(k+1)*Ncomp) = &
                datanew(i,j,k*Ncomp+1:(k+1)*Ncomp) + &
                0.25d0 * alpha(i,j,k+1) * dt * fluxx(i,j,k*Ncomp+1:(k+1)*Ncomp)
            alpha(i,j,k+1) = exp(-datanew(i,j,k*Ncomp+1))
            gamma_up(i,j,9*k+1) = alpha(i,j,k+1)**2
         end do
      end do
  end do

  call compute_swe_flux (datanew, pnlo, pnhi, lo, hi, Ncomp, nlayers, &
       fluxx, fluxy, lo-4, hi+4, dx, dt, gamma_up, glo, ghi, beta, alpha)

  do    j = lo(2), hi(2)
     do i = lo(1), hi(1)
        datanew(i,j,:) = 1.0d0 / 3.0d0 * (dataold(i,j,:) + &
             2.0d0 * (datanew(i,j,:) + fluxx(i,j,:) * dtdx(1) + &
             fluxy(i,j,:) * dtdx(2)))
         do k = 0, nlayers-1
             alpha(i,j,k+1) = exp(-datanew(i,j,k*Ncomp+1))
             gamma_up(i,j,9*k+1) = alpha(i,j,k+1)**2
         end do
     end do
  end do

  call buoyancy(datanew, pnlo, pnhi, lo, hi, Ncomp, nlayers, &
       fluxx, lo-4, hi+4, dx)

   do       k = 0, nlayers-1
       do   j = lo(2), hi(2)
         do i = lo(1), hi(1)
            datanew(i,j,k*Ncomp+1:(k+1)*Ncomp) = &
                datanew(i,j,k*Ncomp+1:(k+1)*Ncomp) + &
                alpha(i,j,k+1) * dt * fluxx(i,j,k*Ncomp+1:(k+1)*Ncomp) * &
                2.0d0/3.0d0
            alpha(i,j,k+1) = exp(-datanew(i,j,k*Ncomp+1))
            gamma_up(i,j,9*k+1) = alpha(i,j,k+1)**2
         end do
      end do
  end do

  ! update height for printing
  call height(alpha(pnlo(1):pnhi(1),pnlo(2):pnhi(2),:), datanew, &
              alpha0, M, R, pnlo, pnhi, nlayers, Ncomp)

end subroutine update_swe_data
