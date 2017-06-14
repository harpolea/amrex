subroutine W_swe(q, lo, hi, Ncomp, gamma_up, glo, ghi, W)
    ! Calculate lorentz factor
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: lo(3), hi(3),Ncomp, glo(3), ghi(3)
    real(amrex_real), intent(in) :: q(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    real(amrex_real), intent(in) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9)
    real(amrex_real), intent(out) :: W(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

    integer i,j,l

    do l = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                W(i,j,l) = sqrt((q(i,j,l,2)**2 * gamma_up(i,j,l,1)+&
                      2.0d0 * q(i,j,l,2) * q(i,j,l,3) * &
                      gamma_up(i,j,l,2) + q(i,j,l,3)**2 * &
                      gamma_up(i,j,l,5)) / q(i,j,l,1)**2 + 1.0d0)
                ! nan check
                if (W(i,j,l) /= W(i,j,l)) then
                    W(i,j,l) = 1.0d0
                end if
          end do
      end do
    end do

end subroutine W_swe

subroutine swe_flux(U, f, lo, hi, Ncomp, x_dir,&
                    gamma_up, glo, ghi, beta, alpha)
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: Ncomp
    integer, intent(in) :: lo(3), hi(3), glo(3), ghi(3)
    real(amrex_real), intent(in)  :: U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    real(amrex_real), intent(out) :: f(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    logical, intent(in) :: x_dir
    real(amrex_real), intent(in) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 3*3)
    real(amrex_real), intent(in) :: beta(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 3)
    real(amrex_real), intent(in) :: alpha(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3))

    integer i, j, l
    real(amrex_real) v(2)
    real(amrex_real) W(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real), parameter :: half = 0.5d0

    call W_swe(U, lo, hi, Ncomp, gamma_up, glo, ghi, W)

    do l = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)

                v(1) = U(i,j,l,2) / (U(i,j,l,1) * W(i,j,l))
                v(2) = U(i,j,l,3) / (U(i,j,l,1) * W(i,j,l))

                if (x_dir) then
                  f(i,j,l,1) = U(i,j,l,1) * &
                        (v(1)*gamma_up(i,j,l,1) + v(2)*gamma_up(i,j,l,2) -&
                         beta(i,j,l,1) / alpha(i,j,l))
                  f(i,j,l,2) = U(i,j,l,2) * &
                        (v(1)*gamma_up(i,j,l,1) + v(2)*gamma_up(i,j,l,2) -&
                         beta(i,j,l,1) / alpha(i,j,l)) + half * U(i,j,l,1)**2 / W(i,j,l)**2
                  f(i,j,l,3) = U(i,j,l,3) * &
                        (v(1)*gamma_up(i,j,l,1) + v(2)*gamma_up(i,j,l,2) -&
                         beta(i,j,l,1) / alpha(i,j,l))
                else
                  f(i,j,l,1) = U(i,j,l,1) * &
                        (v(1)*gamma_up(i,j,l,2) + v(2)*gamma_up(i,j,l,5) -&
                         beta(i,j,l,2) / alpha(i,j,l))
                  f(i,j,l,2) = U(i,j,l,2) * &
                        (v(1)*gamma_up(i,j,l,2) + v(2)*gamma_up(i,j,l,5) -&
                         beta(i,j,l,2) / alpha(i,j,l))
                  f(i,j,l,3) = U(i,j,l,3) * &
                        (v(1)*gamma_up(i,j,l,2) + v(2)*gamma_up(i,j,l,5) -&
                         beta(i,j,l,2) / alpha(i,j,l)) +  half * U(i,j,l,1)**2 / W(i,j,l)**2
                end if
            end do
        end do
    end do

end subroutine swe_flux

subroutine compute_swe_flux (U, datalo, datahi, lo, hi, Ncomp, &
     fluxx, fluxy, flo, fhi, dx, dt, gamma_up, glo, ghi, beta, alpha)! bind(C, name="compute_flux")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(3), hi(3), Ncomp, datalo(3), datahi(3), flo(3), fhi(3), glo(3), ghi(3)
  real(amrex_real), intent(in)    :: U(datalo(1):datahi(1), datalo(2):datahi(2), datalo(3):datahi(3), Ncomp)
  real(amrex_real), intent(inout) :: fluxx( flo(1): fhi(1), flo(2): fhi(2), flo(3):fhi(3), Ncomp)
  real(amrex_real), intent(inout) :: fluxy( flo(1): fhi(1), flo(2): fhi(2), flo(3):fhi(3), Ncomp)
  real(amrex_real), intent(in)    :: dt, dx(3)

  ! local variables
  integer i,j,k,l
  real(amrex_real) S_upwind(Ncomp), S_downwind(Ncomp), S(Ncomp), r(Ncomp), ph(Ncomp), f_p(Ncomp), f_m(Ncomp)
  real(amrex_real), parameter :: half = 0.5d0, alph = 0.5d0
  real(amrex_real), intent(in) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 3*3)
  real(amrex_real), intent(in) :: beta(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 3)
  real(amrex_real), intent(in) :: alpha(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3))

  real(amrex_real) Up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3):hi(3), Ncomp)
  real(amrex_real) Um(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3):hi(3), Ncomp)
  real(amrex_real) fp(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3):hi(3), Ncomp)
  real(amrex_real) fm(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3):hi(3), Ncomp)
  real(amrex_real) dxdt(3)

  dxdt = dx/dt

  ! x fluxes
  do        k = lo(3), hi(3)
      do    j = lo(2)-1, hi(2)+1
         do i = lo(1)-1, hi(1)+1
             S_upwind = U(i+1,j,k,:) - U(i,j,k,:)
             S_downwind = U(i,j,k,:) - U(i-1,j,k,:)
             S = half * (S_upwind + S_downwind)

             r = S_upwind / S_downwind

             call phi(r, ph, Ncomp)

             S = S * ph

             Up(i,j,k,:) = U(i,j,k,:) + half * S
             Um(i,j,k,:) = U(i,j,k,:) - half * S
         end do
     end do
 end do

 !write(*,*) g, dt
 !return

 call swe_flux(Up, fp, lo-1, hi+1, Ncomp, .true., gamma_up, glo, ghi, beta, alpha)
 call swe_flux(Um, fm, lo-1, hi+1, Ncomp, .true., gamma_up, glo, ghi, beta, alpha)


 do         k = lo(3), hi(3)
     do     j = lo(2), hi(2)
         do i = lo(1), hi(1)
            f_p = half * (fp(i,j,k,:) + fm(i+1,j,k,:) + &
                alph * dxdt(1) * (Up(i,j,k,:) - Um(i+1,j,k,:)))
            f_m = half * (fp(i-1,j,k,:) + fm(i,j,k,:) + &
                alph * dxdt(1) * (Up(i-1,j,k,:) - Um(i,j,k,:)))

            fluxx(i,j,k,:) = -(f_p - f_m)
         end do
    end do
 end do


  ! y fluxes
  do        l = lo(3), hi(3)
      do    j = lo(2)-1, hi(2)+1
         do i = lo(1)-1, hi(1)+1
             S_upwind = U(i,j+1,l,:) - U(i,j,l,:)
             S_downwind = U(i,j,l,:) - U(i,j-1,l,:)
             S = half * (S_upwind + S_downwind)

             r = S_upwind / S_downwind

             call phi(r, ph, Ncomp)

             S = S * ph

             Up(i,j,l,:) = U(i,j,l,:) + half * S
             Um(i,j,l,:) = U(i,j,l,:) - half * S
         end do
     end do
 end do

 call swe_flux(Up, fp, lo-1, hi+1, Ncomp, .false., gamma_up, glo, ghi, beta, alpha)
 call swe_flux(Um, fm, lo-1, hi+1, Ncomp, .false., gamma_up, glo, ghi, beta, alpha)

 do           k = lo(3), hi(3)
       do     j = lo(2), hi(2)
           do i = lo(1), hi(1)
             f_p = half * (fp(i,j,k,:) + fm(i,j+1,k,:) + &
                alph * dxdt(2) * (Up(i,j,k,:) - Um(i,j+1,k,:)))
             f_m = half * (fp(k,j-1,k,:) + fm(i,j,k,:) + &
                alph * dxdt(2) * (Up(i,j-1,k,:) - Um(i,j,k,:)))

             fluxy(i,j,k,:) = -(f_p - f_m)
           end do
       end do
   end do

  !write(*,*) fluxx(5,5,1)

end subroutine compute_swe_flux

subroutine buoyancy(U, datalo, datahi, lo, hi, Ncomp, buoy, flo, fhi, dx)

     use amrex_fort_module, only : amrex_real
     implicit none

     integer lo(3), hi(3), Ncomp, datalo(3), datahi(3), flo(3), fhi(3)
     real(amrex_real), intent(in)    :: U(datalo(1):datahi(1), datalo(2):datahi(2), datalo(3):datahi(3), Ncomp)
     real(amrex_real), intent(inout) :: buoy( flo(1): fhi(1), flo(2): fhi(2), flo(3):fhi(3), Ncomp)
     real(amrex_real), intent(in)    :: dx(3)

     integer i,j,l,m

     buoy = 0.0d0

     do        l = lo(3), hi(3)
         do    j = lo(2), hi(2)
            do i = lo(1), hi(1)
                do m = l+1, hi(3)-1
                   buoy(i,j,l,2) = buoy(i,j,l,2) + &
                        U(i,j,l,1) * (U(i+1,j,m,1) - &
                        U(i-1,j,m,1)) / (2.0d0 * dx(1))

                   buoy(i,j,l,3) = buoy(i,j,l,3) + &
                        U(i,j,l,1) * (U(i,j+1,m,1) -&
                        U(i,j-1,m,1)) / (2.0d0 * dx(2))
                end do

                ! TODO: need to define densities somewhere
                do m = 1, l-1
                   buoy(i,j,l,2) = buoy(i,j,l,2) + &
                        U(i,j,l,1) * (U(i+1,j,m,1) - &
                        U(i-1,j,m,1)) / (2.0d0 * dx(1))

                   buoy(i,j,l,3) = buoy(i,j,l,3) + &
                        U(i,j,l,1) * (U(i,j+1,m,1) - &
                        U(i,j-1,m,1)) / (2.0d0 * dx(2))
                end do
            end do
         end do
     end do

end subroutine buoyancy

subroutine height(alpha, U, alpha0, M, R, lo, hi, Ncomp)

    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) ::  lo(3), hi(3), Ncomp
    real(amrex_real), intent(in)    :: alpha(lo(1):hi(1),lo(2):hi(2), lo(3):hi(3))
    real(amrex_real), intent(inout)    :: U(lo(1):hi(1),lo(2):hi(2), lo(3):hi(3), Ncomp)
    real(amrex_real), intent(in) :: alpha0, M, R

    U(:,:,:,4) = alpha0 * R**2 / M * (alpha - alpha0)

end subroutine height

subroutine update_swe_data (lo, hi, Ncomp, dataold, polo, pohi, &
     datanew, pnlo, pnhi, dx, dt, gamma_up, glo, ghi, beta, blo, bhi, alpha, alo, ahi, alpha0, M, R) bind(C, name="update_swe_data")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(3), hi(3), Ncomp, polo(3), pohi(3), pnlo(3), pnhi(3), glo(3), ghi(3), blo(3), bhi(3), alo(3), ahi(3)
  real(amrex_real), intent(in)    :: dataold(polo(1):pohi(1),polo(2):pohi(2), polo(3):pohi(3), Ncomp)
  real(amrex_real), intent(inout) :: datanew(pnlo(1):pnhi(1),pnlo(2):pnhi(2), pnlo(3):pnhi(3), Ncomp)
  real(amrex_real), intent(in)    :: dx(3), alpha0, M, R
  real(amrex_real), value         :: dt
  real(amrex_real), intent(inout) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 3*3)
  real(amrex_real), intent(in) :: beta(blo(1):bhi(1), blo(2):bhi(2), blo(3):bhi(3), 3)
  real(amrex_real), intent(inout) :: alpha(alo(1):ahi(1), alo(2):ahi(2), alo(3):ahi(3))

  ! local variables
  integer i,j,k
  real(amrex_real) :: dtdx(3)

  real(amrex_real) fluxx(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, lo(3):hi(3), Ncomp)
  real(amrex_real) fluxy(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, lo(3):hi(3), Ncomp)

  dtdx = dt/dx

  call compute_swe_flux (dataold, polo, pohi, lo-4, hi+4, Ncomp, &
       fluxx, fluxy, lo-4, hi+4, dx, dt, gamma_up, glo, ghi, beta, alpha)

   do   j = lo(2)-4, hi(2)+4
     do i = lo(1)-4, hi(1)+4
        datanew(i,j,:,:) = dataold(i,j,:,:) &
             + fluxx(i,j,:,:) * dtdx(1) + fluxy(i,j,:,:) * dtdx(2)
        do k = lo(3), hi(3)
            alpha(i,j,k) = exp(-datanew(i,j,k,1))
            gamma_up(i,j,k,1) = alpha(i,j,k)**2
        end do
     end do
  end do

  call buoyancy(datanew, pnlo, pnhi, lo-4, hi+4, Ncomp, &
       fluxx, lo-4, hi+4, dx)

   do       k = lo(3), hi(3)
       do   j = lo(2)-4, hi(2)+4
         do i = lo(1)-4, hi(1)+4
            datanew(i,j,k,:) = datanew(i,j,k,:) + &
                               alpha(i,j,k) * dt * fluxx(i,j,k,:)
            alpha(i,j,k) = exp(-datanew(i,j,k,1))
            gamma_up(i,j,k,1) = alpha(i,j,k)**2
         end do
      end do
  end do

  call compute_swe_flux (datanew, pnlo, pnhi, lo-2, hi+2, Ncomp, &
       fluxx, fluxy, lo-4, hi+4, dx, dt, gamma_up, glo, ghi, beta, alpha)

   do   j = lo(2)-2, hi(2)+2
     do i = lo(1)-2, hi(1)+2
        datanew(i,j,:,:) = 0.25d0 * (3.d0 * dataold(i,j,:,:) + &
             datanew(i,j,:,:) + fluxx(i,j,:,:) * dtdx(1) + fluxy(i,j,:,:) * dtdx(2))
         do k = 1,hi(3)-1
             alpha(i,j,k) = exp(-datanew(i,j,k,1))
             gamma_up(i,j,k,1) = alpha(i,j,k)**2
         end do
     end do
  end do

  call buoyancy(datanew, pnlo, pnhi, lo-2, hi+2, Ncomp, &
       fluxx, lo-4, hi+4, dx)

   do       k = lo(3), hi(3)
       do   j = lo(2)-2, hi(2)+2
         do i = lo(1)-2, hi(1)+2
            datanew(i,j,k,:) = &
                datanew(i,j,k,:) + &
                0.25d0 * alpha(i,j,k) * dt * fluxx(i,j,k,:)
            alpha(i,j,k) = exp(-datanew(i,j,k,1))
            gamma_up(i,j,k,1) = alpha(i,j,k)**2
         end do
      end do
  end do

  call compute_swe_flux (datanew, pnlo, pnhi, lo, hi, Ncomp, &
       fluxx, fluxy, lo-4, hi+4, dx, dt, gamma_up, glo, ghi, beta, alpha)

  do    j = lo(2), hi(2)
     do i = lo(1), hi(1)
        datanew(i,j,:,:) = 1.0d0 / 3.0d0 * (dataold(i,j,:,:) + &
             2.0d0 * (datanew(i,j,:,:) + fluxx(i,j,:,:) * dtdx(1) + &
             fluxy(i,j,:,:) * dtdx(2)))
         do k = lo(3), hi(3)
             alpha(i,j,k) = exp(-datanew(i,j,k,1))
             gamma_up(i,j,k,1) = alpha(i,j,k)**2
         end do
     end do
  end do

  call buoyancy(datanew, pnlo, pnhi, lo, hi, Ncomp, &
       fluxx, lo-4, hi+4, dx)

   do       k = lo(3), hi(3)
       do   j = lo(2), hi(2)
         do i = lo(1), hi(1)
            datanew(i,j,k,:) = &
                datanew(i,j,k,:) + &
                alpha(i,j,k) * dt * fluxx(i,j,k,:) * &
                2.0d0/3.0d0
            alpha(i,j,k) = exp(-datanew(i,j,k,1))
            gamma_up(i,j,k,1) = alpha(i,j,k)**2
         end do
      end do
  end do

  ! update height for printing
  call height(alpha(pnlo(1):pnhi(1),pnlo(2):pnhi(2),pnlo(3):pnhi(3)), &
              datanew, alpha0, M, R, pnlo, pnhi, Ncomp)

end subroutine update_swe_data
