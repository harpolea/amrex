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
        ! Van Leer MC
        if (r(i) > 0.d0) then
            ph(i) = min(2.0d0 * r(i) / (1.0d0 + r(i)), 2.0d0 / (1.0d0 + r(i)))
        else
            ph(i) = 0.0d0
        end if

        ! superbee
        !if (r(i) >= 1.0d0) then
        !    ph(i) = min(2.0d0, min(r(i), 2.0d0 / (1.0d0 + r(i))))
        !elseif ( r(i) >= 0.5d0 ) then
        !    ph(i) = 1.0d0
        !elseif ( r(i) > 0.0d0 ) then
        !    ph(i) = 2.0d0 * r(i)
        !end if
    end do

end subroutine phi


subroutine comp_flux(U, f, lo, hi, Ncomp, dir, gamma)
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, intent(in) :: Ncomp, dir
    integer, intent(in) :: lo(3), hi(3)
    real(amrex_real), intent(in)  :: U(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    real(amrex_real), intent(out) :: f(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), Ncomp)
    real(amrex_real), intent(in)  :: gamma

    integer i, j
    real(amrex_real), parameter :: half = 0.5d0
    real(amrex_real) :: p(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

    p = (gamma - 1.0d0) * (U(:,:,:,5) - &
        half * (U(:,:,:,2)**2 + U(:,:,:,3)**2 + U(:,:,:,4)**2) / U(:,:,:,1))

    if (p(lo(1), lo(2), lo(3)) /= p(lo(1), lo(2), lo(3))) then
        write (*,*) "p is nan", p(lo(1), lo(2), lo(3)), U(lo(1), lo(2), lo(3),:)
    end if

    if (dir == 0) then
        f(:,:,:,1) = U(:,:,:,2)
        f(:,:,:,2) = U(:,:,:,2)**2 / U(:,:,:,1) + p
        f(:,:,:,3) = U(:,:,:,3) * U(:,:,:,2) / U(:,:,:,1)
        f(:,:,:,4) = U(:,:,:,4) * U(:,:,:,2) / U(:,:,:,1)
        f(:,:,:,5) = (U(:,:,:,5) + p) * U(:,:,:,2) / U(:,:,:,1)
    else if (dir == 1) then
        f(:,:,:,1) = U(:,:,:,3)
        f(:,:,:,2) = U(:,:,:,2) * U(:,:,:,3) / U(:,:,:,1)
        f(:,:,:,3) = U(:,:,:,3)**2 / U(:,:,:,1) + p
        f(:,:,:,4) = U(:,:,:,4) * U(:,:,:,3) / U(:,:,:,1)
        f(:,:,:,5) = (U(:,:,:,5) + p) * U(:,:,:,3) / U(:,:,:,1)
    else
        f(:,:,:,1) = U(:,:,:,4)
        f(:,:,:,2) = U(:,:,:,2) * U(:,:,:,4) / U(:,:,:,1)
        f(:,:,:,3) = U(:,:,:,3) * U(:,:,:,4) / U(:,:,:,1)
        f(:,:,:,4) = U(:,:,:,4)**2 / U(:,:,:,1) + p
        f(:,:,:,5) = (U(:,:,:,5) + p) * U(:,:,:,4) / U(:,:,:,1)
    end if

end subroutine comp_flux

subroutine compute_flux (U, datalo, datahi, lo, hi, Ncomp,&
     fluxx, fluxy, fluxz, flo, fhi, dx, dt, gamma)

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(3), hi(3), Ncomp, datalo(3), datahi(3), flo(3), fhi(3)
  real(amrex_real), intent(in)    :: U(datalo(1):datahi(1), datalo(2):datahi(2), datalo(3):datahi(3), Ncomp)
  real(amrex_real), intent(inout) :: fluxx(flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3), Ncomp)
  real(amrex_real), intent(inout) :: fluxy(flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3), Ncomp)
  real(amrex_real), intent(inout) :: fluxz(flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3), Ncomp)
  real(amrex_real), intent(in)    :: dt, dx(3), gamma

  ! local variables
  integer i,j,k
  real(amrex_real) S_upwind(Ncomp), S_downwind(Ncomp), S(Ncomp), r(Ncomp), ph(Ncomp), f_p(Ncomp), f_m(Ncomp)
  real(amrex_real), parameter :: half = 0.5d0, alph = 0.1d0

  real(amrex_real) Up(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, Ncomp)
  real(amrex_real) Um(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, Ncomp)
  real(amrex_real) fp(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, Ncomp)
  real(amrex_real) fm(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, Ncomp)
  real(amrex_real) dxdt(3)

  dxdt = dx/dt

  ! x fluxes
  do        k = lo(3)-1, hi(3)+1
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

 !write(*,*) "U", U(lo(1),lo(2),lo(3),:)

 call comp_flux(Up, fp, lo-1, hi+1, Ncomp, 0, gamma)
 call comp_flux(Um, fm, lo-1, hi+1, Ncomp, 0, gamma)

 !write(*,*) "Um", Um(lo(1),lo(2),lo(3),:)
 !write(*,*) "fm", fm(lo(1),lo(2),lo(3),:)

  do        k = lo(3), hi(3)
      do    j = lo(2), hi(2)
         do i = lo(1), hi(1)
            f_p = half * (fp(i,j,k,:) + fm(i+1,j,k,:) + alph * dxdt(1) * (Up(i,j,k,:) - Um(i+1,j,k,:)))
            f_m = half * (fp(i-1,j,k,:) + fm(i,j,k,:) + alph * dxdt(1) * (Up(i-1,j,k,:) - Um(i,j,k,:)))

            fluxx(i,j,k,:) = -(f_p - f_m)
         end do
      end do
  end do

  ! y fluxes
  do        k = lo(3)-1, hi(3)+1
      do    j = lo(2)-1, hi(2)+1
         do i = lo(1)-1, hi(1)+1
             S_upwind = U(i,j+1,k,:) - U(i,j,k,:)
             S_downwind = U(i,j,k,:) - U(i,j-1,k,:)
             S = half * (S_upwind + S_downwind)

             r = S_upwind / S_downwind

             call phi(r, ph, Ncomp)

             S = S * ph

             Up(i,j,k,:) = U(i,j,k,:) + half * S
             Um(i,j,k,:) = U(i,j,k,:) - half * S
         end do
     end do
 end do

 call comp_flux(Up, fp, lo-1, hi+1, Ncomp, 1, gamma)
 call comp_flux(Um, fm, lo-1, hi+1, Ncomp, 1, gamma)

   do        k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             f_p = half * (fp(i,j,k,:) + fm(i,j+1,k,:) + alph * dxdt(2) * (Up(i,j,k,:) - Um(i,j+1,k,:)))
             f_m = half * (fp(i,j-1,k,:) + fm(i,j,k,:) + alph * dxdt(2) * (Up(i,j-1,k,:) - Um(i,j,k,:)))

             fluxy(i,j,k,:) = -(f_p - f_m)
          end do
       end do
   end do

   ! z fluxes
   do        k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
              S_upwind = U(i,j,k+1,:) - U(i,j,k,:)
              S_downwind = U(i,j,k,:) - U(i,j,k-1,:)
              S = half * (S_upwind + S_downwind)

              r = S_upwind / S_downwind

              call phi(r, ph, Ncomp)

              S = S * ph

              Up(i,j,k,:) = U(i,j,k,:) + half * S
              Um(i,j,k,:) = U(i,j,k,:) - half * S
          end do
      end do
  end do

  call comp_flux(Up, fp, lo-1, hi+1, Ncomp, 2, gamma)
  call comp_flux(Um, fm, lo-1, hi+1, Ncomp, 2, gamma)

    do        k = lo(3), hi(3)
        do    j = lo(2), hi(2)
           do i = lo(1), hi(1)
              f_p = half * (fp(i,j,k,:) + fm(i,j,k+1,:) + alph * dxdt(3) * (Up(i,j,k,:) - Um(i,j,k+1,:)))
              f_m = half * (fp(i,j,k-1,:) + fm(i,j,k,:) + alph * dxdt(3) * (Up(i,j,k-1,:) - Um(i,j,k,:)))

              fluxz(i,j,k,:) = -(f_p - f_m)
           end do
        end do
    end do

  !write(*,*) fluxx(5,5,1)

end subroutine compute_flux


subroutine update_data (lo, hi, Ncomp, dataold, polo, pohi, &
     datanew, pnlo, pnhi, dx, dt, gamma) bind(C, name="update_data")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(3), hi(3), Ncomp, polo(3), pohi(3), pnlo(3), pnhi(3)
  real(amrex_real), intent(in)    :: dataold(polo(1):pohi(1),polo(2):pohi(2),polo(3):pohi(3), Ncomp)
  real(amrex_real), intent(inout) :: datanew(pnlo(1):pnhi(1),pnlo(2):pnhi(2),pnlo(3):pnhi(3), Ncomp)
  real(amrex_real), intent(in)    :: dx(3), gamma
  real(amrex_real), value         :: dt

  ! local variables
  integer i,j,k
  real(amrex_real) :: dtdx(3)

  real(amrex_real) fluxx(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, lo(3)-4:hi(3)+4, Ncomp)
  real(amrex_real) fluxy(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, lo(3)-4:hi(3)+4, Ncomp)
  real(amrex_real) fluxz(lo(1)-4:hi(1)+4, lo(2)-4:hi(2)+4, lo(3)-4:hi(3)+4, Ncomp)

  dtdx = dt/dx

  !write(*,*) "old data", dataold(10,10,10,:)

  call compute_flux (dataold, polo, pohi, lo-4, hi+4, Ncomp, &
       fluxx, fluxy, fluxz, lo-4, hi+4, dx, dt, gamma)

   do       k = lo(3)-4, hi(3)+4
       do   j = lo(2)-4, hi(2)+4
         do i = lo(1)-4, hi(1)+4
            datanew(i,j,k,:) = dataold(i,j,k,:) &
                 + fluxx(i,j,k,:) * dtdx(1) + fluxy(i,j,k,:) * dtdx(2) +&
                  fluxz(i,j,k,:) * dtdx(3)
         end do
      end do
  end do

  !write(*,*) "fluxx", fluxx(10,10,10,:)
  !write(*,*) "new data", datanew(10,10,10,:)

  call compute_flux (datanew, pnlo, pnhi, lo-2, hi+2, Ncomp, &
       fluxx, fluxy, fluxz, lo-4, hi+4, dx, dt, gamma)

   do       k = lo(3)-2, hi(3)+2
       do   j = lo(2)-2, hi(2)+2
         do i = lo(1)-2, hi(1)+2
            datanew(i,j,k,:) = 0.25d0 * (3.d0 * dataold(i,j,k,:) + &
                 datanew(i,j,k,:) + fluxx(i,j,k,:) * dtdx(1) + &
                 fluxy(i,j,k,:) * dtdx(2) + fluxz(i,j,k,:) * dtdx(3))
         end do
      end do
  end do

  call compute_flux (datanew, pnlo, pnhi, lo, hi, Ncomp, &
       fluxx, fluxy, fluxz, lo-4, hi+4, dx, dt, gamma)

  do        k = lo(3), hi(3)
      do    j = lo(2), hi(2)
         do i = lo(1), hi(1)
            datanew(i,j,k,:) = 1.0d0 / 3.0d0 * (dataold(i,j,k,:) + &
                 2.0d0 * (datanew(i,j,k,:) + fluxx(i,j,k,:) * dtdx(1) +&
                  fluxy(i,j,k,:) * dtdx(2) + fluxz(i,j,k,:) * dtdx(3)))
         end do
      end do
  end do

end subroutine update_data
