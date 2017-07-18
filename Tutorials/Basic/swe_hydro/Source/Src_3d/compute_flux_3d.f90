module compute_flux_module

  implicit none

  private

  public :: compute_flux_3d

contains

  subroutine compute_flux_3d(lo, hi, dt, dx, &
                             phi, ph_lo, ph_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi, &
                             flxz, fz_lo, fz_hi, &
                             phi_p, phi_m, fp, fm, &
                             slope, glo, ghi, Ncomp, gr, &
                             alpha0, M, R)

    use slope_module, only: slopex, slopey, slopez

    integer, intent(in) :: lo(3), hi(3), glo(3), ghi(3), Ncomp
    double precision, intent(in) :: dt, dx(3)
    integer, intent(in) :: ph_lo(3), ph_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    integer, intent(in) :: fz_lo(3), fz_hi(3)
    double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3), Ncomp)
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),Ncomp)
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),Ncomp)
    double precision, intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),Ncomp)
    double precision, dimension(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3), Ncomp) :: &
         phi_p, phi_m, fp, fm, slope
    logical, intent(in) :: gr
    double precision, intent(in) :: alpha0, M, R

    integer :: i, j, k, l
    double precision :: dxdt(3), f_p(Ncomp), f_m(Ncomp)
    double precision, parameter :: g = 1.0d-3, gamma = 5.0d0/3.0d0, alph = 0.5d0
    double precision :: alpha(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3)), gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9), U_prim(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3), Ncomp), p(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3))

    dxdt = dx/dt
    flxx(:,:,:,:) = 0.0d0
    flxy(:,:,:,:) = 0.0d0
    flxz(:,:,:,:) = 0.0d0

    write(*,*) "phi: ", phi(lo(1)+3, lo(2)+3, lo(3)+3, :)

    call slopex(lo-1, hi+1, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi, Ncomp)

    ! compute phi on x faces using umac to upwind; ignore transverse terms
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
              phi_p(i,j,k,:) = phi(i,j,k,:) + 0.5d0 * slope(i,j,k,:)
              phi_m(i,j,k,:) = phi(i,j,k,:) - 0.5d0 * slope(i,j,k,:)
          end do
       end do
    end do

    if (Ncomp < 5) then
        if (gr) then
            call gr_swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, 0, alpha)
            call gr_swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, 0, alpha)
        else
            call swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, g, 0)
            call swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, g, 0)
        end if
    else
        if (gr) then
            call gr_comp_flux(phi_p, fp, lo, hi, Ncomp, 0, gamma, glo, ghi, alpha, dx, alpha0, M, R)
            call gr_comp_flux(phi_m, fm, lo, hi, Ncomp, 0, gamma, glo, ghi, alpha, dx, alpha0, M, R)
        else
            call comp_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, gamma, 0)
            call comp_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, gamma, 0)
        end if
    end if

    do        k = lo(3), hi(3)
        do    j = lo(2), hi(2)
           do i = lo(1), hi(1)
              f_p = 0.5d0 * (fp(i,j,k,:) + fm(i+1,j,k,:) + alph * dxdt(1) * (phi_p(i,j,k,:) - phi_m(i+1,j,k,:)))
              f_m = 0.5d0 * (fp(i-1,j,k,:) + fm(i,j,k,:) + alph * dxdt(1) * (phi_p(i-1,j,k,:) - phi_m(i,j,k,:)))

              flxx(i,j,k,:) = -(f_p - f_m)

              if (gr) then
                  flxx(i,j,k,:) = flxx(i,j,k,:) * alpha(i,j,k)
              end if
              do l = 1, Ncomp
                  if (flxx(i,j,k,l) /= flxx(i,j,k,l)) then
                      flxx(i,j,k,l) = 0.d0
                  end if
              end do
           end do
        end do
    end do

    call slopey(lo-1, hi+1, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi, Ncomp)

    ! compute phi on y faces using vmac to upwind; ignore transverse terms
    do       k = lo(3), hi(3)
       do    j = lo(2)-1  , hi(2)+1
          do i = lo(1), hi(1)

              phi_p(i,j,k,:) = phi(i,j,k,:) + 0.5d0 * slope(i,j,k,:)
              phi_m(i,j,k,:) = phi(i,j,k,:) - 0.5d0 * slope(i,j,k,:)

          end do
       end do
    end do

    if (Ncomp < 5) then
        if (gr) then
            call gr_swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, 1, alpha)
            call gr_swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, 1, alpha)
        else
            call swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, g, 1)
            call swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, g, 1)
        end if
    else
        if (gr) then
            call gr_comp_flux(phi_p, fp, lo, hi, Ncomp, 1, gamma, glo, ghi, alpha, dx, alpha0, M, R)
            call gr_comp_flux(phi_m, fm, lo, hi, Ncomp, 1, gamma, glo, ghi, alpha, dx, alpha0, M, R)
        else
            call comp_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, gamma, 1)
            call comp_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, gamma, 1)
        end if
    end if

    do        k = lo(3), hi(3)
        do    j = lo(2), hi(2)
           do i = lo(1), hi(1)
              f_p = 0.5d0 * (fp(i,j,k,:) + fm(i,j+1,k,:) + alph * dxdt(2) * (phi_p(i,j,k,:) - phi_m(i,j+1,k,:)))
              f_m = 0.5d0 * (fp(i,j-1,k,:) + fm(i,j,k,:) + alph * dxdt(2) * (phi_p(i,j-1,k,:) - phi_m(i,j,k,:)))

              flxy(i,j,k,:) = -(f_p - f_m)

              if (gr) then
                  flxy(i,j,k,:) = flxy(i,j,k,:) * alpha(i,j,k)
              end if
              do l = 1, Ncomp
                  if (flxy(i,j,k,l) /= flxy(i,j,k,l)) then
                      flxy(i,j,k,l) = 0.d0
                  end if
              end do
           end do
        end do
    end do

    call slopez(lo-1, hi+1, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi, Ncomp)

    ! compute phi on z faces using wmac to upwind; ignore transverse terms
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)

              phi_p(i,j,k,:) = phi(i,j,k,:) + 0.5d0 * slope(i,j,k,:)
              phi_m(i,j,k,:) = phi(i,j,k,:) - 0.5d0 * slope(i,j,k,:)

          end do
       end do
    end do

    if (Ncomp < 5) then
        !if (gr) then
        !    call gr_swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, 2, alpha)
        !    call gr_swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, 2, alpha)
        !else
        !    call swe_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, g, 2)
        !    call swe_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, g, 2)
        !end if
        flxz(:,:,:,:) = 0.0d0
    else
        if (gr) then
            call gr_comp_flux(phi_p, fp, lo, hi, Ncomp, 2, gamma, glo, ghi, alpha, dx, alpha0, M, R)
            call gr_comp_flux(phi_m, fm, lo, hi, Ncomp, 2, gamma, glo, ghi, alpha, dx, alpha0, M, R)
        else
            call comp_flux(phi_p, fp, glo, ghi, lo, hi, Ncomp, gamma, 2)
            call comp_flux(phi_m, fm, glo, ghi, lo, hi, Ncomp, gamma, 2)
        end if

        do        k = lo(3), hi(3)
            do    j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  f_p = 0.5d0 * (fp(i,j,k,:) + fm(i,j,k+1,:) + alph * dxdt(3) * (phi_p(i,j,k,:) - phi_m(i,j,k+1,:)))
                  f_m = 0.5d0 * (fp(i,j,k-1,:) + fm(i,j,k,:) + alph * dxdt(3) * (phi_p(i,j,k-1,:) - phi_m(i,j,k,:)))

                  flxz(i,j,k,:) = -(f_p - f_m)

                  if (gr) then
                      flxz(i,j,k,:) = flxz(i,j,k,:) * alpha(i,j,k)
                  end if
                  do l = 1, Ncomp
                      if (flxz(i,j,k,l) /= flxz(i,j,k,l)) then
                          flxz(i,j,k,l) = 0.d0
                      end if
                  end do
               end do
            end do
        end do

        if (gr) then
            write(*,*) "fz before: ", flxz(lo(1), lo(2), lo(3), 4)
            call calc_gamma_up(gamma_up, glo, ghi, lo, hi, alpha0, M, R, dx)

            call cons_to_prim(phi, ph_lo, ph_hi, U_prim, ph_lo, ph_hi, p, ph_lo, ph_hi, lo, hi, Ncomp, gamma, alpha0, M, R, dx)

            call gr_sources(flxz(:,:,:,4), fz_lo, fz_hi, phi, ph_lo, ph_hi, &
                p, ph_lo, ph_hi, alpha, glo, ghi, gamma_up, glo, ghi, &
                M, R, gamma, Ncomp, lo, hi, dx)

            write(*,*) "fz after: ", flxz(lo(1), lo(2), lo(3), 4)
        end if

        !write(*,*) "slope", slope(lo(1), lo(2), lo(3), :)
        !write(*,*) "phi", phi(lo(1), lo(2), lo(3), :)
        !write(*,*) "phi_p", phi_p(lo(1), lo(2), lo(3), :)
        !write(*,*) "phi_m", phi_m(lo(1), lo(2), lo(3), :)
        !write(*,*) "fp", fp(lo(1), lo(2), lo(3), :)
        !write(*,*) "fm", fm(lo(1), lo(2), lo(3), :)
    end if

  end subroutine compute_flux_3d

  subroutine swe_flux(U, f, glo, ghi, lo, hi, Ncomp, g, x_dir)
      implicit none

      integer, intent(in) :: Ncomp
      double precision, intent(in) :: g
      integer, intent(in) :: glo(3), ghi(3), lo(3), hi(3)
      double precision, intent(in)  :: U(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), Ncomp)
      double precision, intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), Ncomp)
      integer, intent(in) :: x_dir

      integer i, j, k

      f(:,:,:,:) = 0.0d0

      if (x_dir == 0) then
          do k = lo(3), hi(3)
              do    j = lo(2), hi(2)
                 do i = lo(1)-1, hi(1)+1
                    f(i,j,k,1) =  U(i,j,k,2)
                    f(i,j,k,2) = U(i,j,k,2)**2/U(i,j,k,1) + 0.5d0 * g * U(i,j,k,1)**2
                    f(i,j,k,3) =  U(i,j,k,2) * U(i,j,k,3) / U(i,j,k,1)
                 end do
              end do
          end do
      else if (x_dir == 1) then
          do k = lo(3), hi(3)
              do    j = lo(2)-1, hi(2)+1
                 do i = lo(1), hi(1)
                    f(i,j,k,1) =  U(i,j,k,3)
                    f(i,j,k,2) = U(i,j,k,2) * U(i,j,k,3) / U(i,j,k,1)
                    f(i,j,k,3) =  U(i,j,k,3)**2/U(i,j,k,1) + 0.5d0 * g * U(i,j,k,1)**2
                 end do
              end do
          end do
      end if

  end subroutine swe_flux

  subroutine W_swe(q, lo, hi, Ncomp, gamma_up, glo, ghi, W)
      ! Calculate lorentz factor
      implicit none

      integer, intent(in) :: lo(3), hi(3), Ncomp, glo(3), ghi(3)
      double precision, intent(in) :: q(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), Ncomp)
      double precision, intent(in) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9)
      double precision, intent(out) :: W(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3))

      integer i, j, k

      do k = lo(3), hi(3)
          do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                  W(i,j,k) = sqrt((q(i,j,k,2)**2 * gamma_up(i,j,k,1)+ &
                        2.0d0 * q(i,j,k,2) * q(i,j,k,3) * gamma_up(i,j,k,2) + &
                        q(i,j,k,3)**2 * gamma_up(i,j,k,5)) / q(i,j,k,1)**2 + 1.0d0)
                  ! nan check
                  if (W(i,j,k) /= W(i,j,k)) then
                      W(i,j,k) = 1.0d0
                  end if
            end do
        end do
    end do

  end subroutine W_swe

  subroutine gr_swe_flux(U, f, glo, ghi, lo, hi, Ncomp, x_dir, alpha)
      implicit none

      integer, intent(in) :: Ncomp
      integer, intent(in) :: glo(3), ghi(3), lo(3), hi(3)
      double precision, intent(in)  :: U(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), Ncomp)
      double precision, intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), Ncomp)
      integer, intent(in) :: x_dir
      double precision, intent(out) :: alpha(glo(1):ghi(1),glo(2):ghi(2), glo(3):ghi(3))

      integer :: i, j, k
      double precision :: v(2), W(glo(1):ghi(1),glo(2):ghi(2), glo(3):ghi(3))
      double precision :: beta(glo(1):ghi(1),glo(2):ghi(2), glo(3):ghi(3),3)
      double precision :: gamma_up(glo(1):ghi(1),glo(2):ghi(2), glo(3):ghi(3),9)

      alpha = exp(-U(:,:,:,1))
      beta = 0.0d0

      gamma_up = 0.0d0
      gamma_up(:,:,:,1) = 1.0d0
      gamma_up(:,:,:,5) = 1.0d0

      call W_swe(U, lo-1, hi+1, Ncomp, gamma_up, glo, ghi, W)

      f(:,:,:,:) = 0.0d0

      if (x_dir == 0) then
          do k = lo(3), hi(3)
              do    j = lo(2), hi(2)
                 do i = lo(1)-1, hi(1)+1
                    if (U(i,j,k,1) < 1.d-20) then
                        v(1) = U(i,j,k,2)
                        v(2) = U(i,j,k,3)
                    else
                        v(1) = U(i,j,k,2) / (U(i,j,k,1) * W(i,j,k))
                        v(2) = U(i,j,k,3) / (U(i,j,k,1) * W(i,j,k))
                    end if

                    f(i,j,k,1) =  U(i,j,k,1) * &
                          (v(1)*gamma_up(i,j,k,1) + v(2)*gamma_up(i,j,k,2) -&
                           beta(i,j,k,1) / alpha(i,j,k))
                    f(i,j,k,2) = U(i,j,k,2) * &
                          (v(1)*gamma_up(i,j,k,1) + v(2)*gamma_up(i,j,k,2) -&
                           beta(i,j,k,1) / alpha(i,j,k)) + 0.5d0 * U(i,j,k,1)**2 / W(i,j,k)**2
                    f(i,j,k,3) =  U(i,j,k,3) * &
                          (v(1)*gamma_up(i,j,k,1) + v(2)*gamma_up(i,j,k,2) -&
                           beta(i,j,k,1) / alpha(i,j,k))

                    !f(i,j,k,:) = f(i,j,k,:)
                 end do
              end do
          end do
      else if (x_dir == 1) then
          do k = lo(3), hi(3)
              do    j = lo(2)-1, hi(2)+1
                 do i = lo(1), hi(1)
                     if (U(i,j,k,1) < 1.d-20) then
                         v(1) = U(i,j,k,2)
                         v(2) = U(i,j,k,3)
                     else
                         v(1) = U(i,j,k,2) / (U(i,j,k,1) * W(i,j,k))
                         v(2) = U(i,j,k,3) / (U(i,j,k,1) * W(i,j,k))
                     end if

                    f(i,j,k,1) = U(i,j,k,1) * &
                          (v(1)*gamma_up(i,j,k,2) + v(2)*gamma_up(i,j,k,5) -&
                           beta(i,j,k,2) / alpha(i,j,k))
                    f(i,j,k,2) = U(i,j,k,2) * &
                          (v(1)*gamma_up(i,j,k,2) + v(2)*gamma_up(i,j,k,5) -&
                           beta(i,j,k,2) / alpha(i,j,k))
                    f(i,j,k,3) =  U(i,j,k,3) * &
                          (v(1)*gamma_up(i,j,k,2) + v(2)*gamma_up(i,j,k,5) -&
                           beta(i,j,k,2) / alpha(i,j,k)) +  0.5d0 * U(i,j,k,1)**2 / W(i,j,k)**2

                    !f(i,j,k,:) = f(i,j,k,:)
                 end do
              end do
          end do
      end if

  end subroutine gr_swe_flux

  subroutine comp_flux(U, f, glo, ghi, lo, hi, Ncomp, gamma, x_dir)
      implicit none

      integer, intent(in) :: Ncomp
      double precision, intent(in) :: gamma
      integer, intent(in) :: glo(3), ghi(3), lo(3), hi(3)
      double precision, intent(in)  :: U(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), Ncomp)
      double precision, intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), Ncomp)
      integer, intent(in) :: x_dir

      integer i, j, k
      double precision :: p(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3))

      p = (gamma - 1.d0) * (U(:,:,:,5) - 0.5d0 * (U(:,:,:,2)**2 + U(:,:,:,3)**2 + U(:,:,:,4)**2) / U(:,:,:,1))

      if (x_dir == 0) then
          do k = lo(3), hi(3)
              do    j = lo(2), hi(2)
                 do i = lo(1)-1, hi(1)+1
                    f(i,j,k,1) =  U(i,j,k,2)
                    f(i,j,k,2) = U(i,j,k,2)**2 / U(i,j,k,1) + p(i,j,k)
                    f(i,j,k,3) =  U(i,j,k,2) * U(i,j,k,3) / U(i,j,k,1)
                    f(i,j,k,4) =  U(i,j,k,2) * U(i,j,k,4) / U(i,j,k,1)
                    f(i,j,k,5) = (U(i,j,k,5) + p(i,j,k)) * U(i,j,k,2) / U(i,j,k,1)
                 end do
              end do
          end do
      else if (x_dir == 1) then
          do k = lo(3), hi(3)
              do    j = lo(2)-1, hi(2)+1
                 do i = lo(1), hi(1)
                    f(i,j,k,1) =  U(i,j,k,3)
                    f(i,j,k,2) = U(i,j,k,3) * U(i,j,k,2) / U(i,j,k,1)
                    f(i,j,k,3) =  U(i,j,k,3)**2 / U(i,j,k,1) + p(i,j,k)
                    f(i,j,k,4) = U(i,j,k,3) * U(i,j,k,2) / U(i,j,k,1)
                    f(i,j,k,5) = (U(i,j,k,5) + p(i,j,k)) * U(i,j,k,3) / U(i,j,k,1)
                 end do
              end do
          end do
      else
          do k = lo(3)-1, hi(3)+1
              do    j = lo(2), hi(2)
                 do i = lo(1), hi(1)
                    f(i,j,k,1) =  U(i,j,k,4)
                    f(i,j,k,2) = U(i,j,k,4) * U(i,j,k,2) / U(i,j,k,1)
                    f(i,j,k,3) = U(i,j,k,4) * U(i,j,k,3) / U(i,j,k,1)
                    f(i,j,k,4) = U(i,j,k,4)**2 / U(i,j,k,1) + p(i,j,k)
                    f(i,j,k,5) = (U(i,j,k,5) + p(i,j,k)) * U(i,j,k,4) / U(i,j,k,1)
                 end do
              end do
          end do
      end if

  end subroutine comp_flux

  subroutine f_of_p(f, p, U, Ncomp, gamma, gamma_up)
      implicit none

      integer, intent(in) :: Ncomp
      double precision, intent(in)  :: U(Ncomp), p, gamma, gamma_up(9)
      double precision, intent(out) :: f

      double precision :: sq

      sq = sqrt((U(5) + p + U(1))**2 - U(2)**2*gamma_up(1)- &
          2.0d0 * U(2) * U(3) * gamma_up(2) - &
          2.0d0 * U(2) * U(4) * gamma_up(3) - &
          U(3)**2 * gamma_up(5) - &
          2.0d0 * U(3) * U(4) * gamma_up(6) -&
          U(4)**2 * gamma_up(9))

      f = (gamma - 1.0d0) * sq / (U(5) + p + U(1)) * &
          (sq - p * (U(5) + p + U(1)) / sq - U(1)) - p

  end subroutine f_of_p

  subroutine zbrent(p, x1, b, U, Ncomp, gamma, gamma_up)
      ! route finder using brent's method
      implicit none

      integer, intent(in) :: Ncomp
      double precision, intent(out) :: p
      double precision, intent(in)  :: U(Ncomp), gamma, gamma_up(9), x1
      double precision, intent(inout) :: b

      double precision, parameter :: TOL = 1.0d-12
      integer, parameter :: ITMAX = 100

      double precision a, c, d, fa, fb, fc, fs, s
      logical mflag, con1, con2, con3, con4, con5
      integer i

      a = x1
      c = 0.0d0
      d = 0.0d0
      call f_of_p(fa, a, U, Ncomp, gamma, gamma_up)
      call f_of_p(fb, b, U, Ncomp, gamma, gamma_up)
      fc = 0.0d0

      if (fa * fb >= 0.0d0) then
          p = b
          return
      end if

      if (abs(fa) < abs(fb)) then
          d = a
          a = b
          b = d

          d = fa
          fa = fb
          fb = d
      end if

      c = a
      fc = fa

      mflag = .true.

      do i = 1, ITMAX
          if (fa /= fc .and. fb /= fc) then
              s = a*fb*fc / ((fa-fb) * (fa-fc)) + b*fa*fc / ((fb-fa)*(fb-fc)) +&
                  c*fa*fb / ((fc-fa)*(fc-fb))
          else
              s = b - fb * (b-a) / (fb-fa)
          end if

          con1 = .false.

          if (0.25d0 * (3.0d0 * a + b) < b) then
              if ( s < 0.25d0 * (3.0d0 * a + b) .or. s > b) then
                  con1 = .true.
              end if
          else if (s < b .or. s > 0.25d0  * (3.0d0 * a + b)) then
              con1 = .true.
          end if

          con2 = mflag .and. abs(s - b) >= 0.5d0 * abs(b-c)

          con3 = (.not. mflag) .and. abs(s-b) >= 0.5d0 * abs(c-d)

          con4 = mflag .and. abs(b-c) < TOL

          con5 = (.not. mflag) .and. abs(c-d) < TOL

          if (con1 .or. con2 .or. con3 .or. con4 .or. con5) then
              s = 0.5d0 * (a + b)
              mflag = .true.
          else
              mflag = .false.
          end if

          call f_of_p(fs, s, U, Ncomp, gamma, gamma_up)

          if (abs(fa) < abs(fb)) then
              d = a
              a = b
              b = d

              d = fa
              fa = fb
              fb = d
          end if

          d = c
          c = b
          fc = fb

          if (fa * fs < 0.0d0) then
              b = s
              fb = fs
          else
              a = s
              fa = fs
          end if

          if (fb == 0.0d0 .or. fs == 0.0d0 .or. abs(b-a) < TOL) then
              p = b
              return
          end if

      end do

      p = x1

  end subroutine zbrent

  subroutine cons_to_prim(U, clo, chi, U_prim, prlo, prhi, p, plo, phi, lo, hi, Ncomp, gamma, alpha0, M, R, dx) bind(C, name="cons_to_prim")
      ! convert from conserved variables (D, Sx, Sy, tau) to primitive variables (rho, v^x, v^y, eps). Also outputs the pressure
      implicit none

      integer, intent(in) :: Ncomp
      integer, intent(in) :: clo(3), chi(3), prlo(3), prhi(3), plo(3), phi(3), lo(3), hi(3)
      double precision, intent(in)  :: U(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3), Ncomp)
      double precision, intent(out) :: U_prim(prlo(1):prhi(1), prlo(2):prhi(2), prlo(3):prhi(3), Ncomp)
      double precision, intent(out) :: p(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
      double precision, intent(in)  :: gamma, alpha0, M, R, dx(3)

      double precision gamma_up(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 9)
      double precision :: pmin, pmax, ssq, q(Ncomp), fmin, fmax, sq, h, W2
      integer :: i, j, k, l

      call calc_gamma_up(gamma_up, lo, hi, lo, hi, alpha0, M, R, dx)

      write(*,*) "cons_to_prim"

      do k = lo(3), hi(3)
          do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                  q = U(i,j,k,:)

                  !HACK
                  if (any(q /= q)) then
                      write(*,*) "q is nan: ", q, i, j, k
                      stop
                  end if
                  if (q(1) /= q(1)) then
                      write(*,*) "q is nan: ", q, i, j, k
                      stop
                      q(1) = 1.0d0
                  end if
                  do l = 2, 4
                      if (q(l) /= q(l)) then
                          q(l) = 0.0d0
                      end if
                  end do
                  if (q(5) == 1.0d0 .or. q(5) /= q(5)) then
                      q(5) = 1.5d0
                  end if
                  ssq = q(2)**2 * gamma_up(i,j,k,1) + &
                      2.0d0 * q(2) * q(3) * gamma_up(i,j,k,2) + &
                      2.0d0 * q(2) * q(4) * gamma_up(i,j,k,3) + &
                      q(3)**2 * gamma_up(i,j,k,5) + &
                      2.0d0 * q(3) * q(4) * gamma_up(i,j,k,6) + &
                      q(4)**2 * gamma_up(i,j,k,9)

                  if (ssq /= ssq) then
                      ssq = 0.d0
                  end if

                  pmin = (1.0d0 - ssq)**2 * q(5) * (gamma - 1.0d0)
                  pmax = (gamma - 1.0d0) * (q(5) + q(1)) / (2.0d0 - gamma)

                  if (pmin < 0.0d0) then
                      pmin = 0.d0
                  end if

                  if (pmax < 0.d0 .or. pmax < pmin) then
                      pmax = 1.0d0
                  end if

                  call f_of_p(fmin, pmin, q, Ncomp, gamma, gamma_up(i,j,k,:))
                  call f_of_p(fmax, pmax, q, Ncomp, gamma, gamma_up(i,j,k,:))

                  if (fmin * fmax > 0.0d0) then
                      pmin = 0.d0
                  end if

                  call f_of_p(fmin, pmin, q, Ncomp, gamma, gamma_up(i,j,k,:))

                  if (fmin * fmax > 0.0d0) then
                      pmax = pmax * 10.d0
                  end if

                  call zbrent(p(i,j,k), pmin, pmax, q, Ncomp, gamma, gamma_up(i,j,k,:))

                  if (p(i,j,k) /= p(i,j,k) .or. p(i,j,k) < 0.0d0 .or. p(i,j,k) > 1.0d0) then
                      p(i,j,k) = abs((gamma - 1.0d0) * (q(5) + q(1)) / (2.0d0 - gamma))

                      if (p(i,j,k) > 1.0d0) then
                          p(i,j,k) = 1.0d0
                      end if
                  end if

                  sq = sqrt((q(5) + p(i,j,k) + q(1))**2 - ssq)

                  if (sq /= sq) then
                      sq = q(5) + p(i,j,k) + q(1)
                  end if

                  h = 1.0d0 + gamma * (sq - p(i,j,k) * (q(5) + p(i,j,k) + q(1)) / sq - q(1)) / q(1)
                  W2 = 1.0d0 + ssq / (q(1) * h)**2

                  !write(*,*) "p, sq", p(i,j,k), sq

                  U_prim(i,j,k,1) = q(1) * sq / (q(5) + p(i,j,k) + q(1))
                  U_prim(i,j,k,2) = (gamma_up(i,j,k,1) * q(2) + &
                      gamma_up(i,j,k,2) * q(3) + gamma_up(i,j,k,3) * q(4)) /&
                      (W2 * h * U_prim(i,j,k,1))
                  U_prim(i,j,k,3) = (gamma_up(i,j,k,4) * q(2) + &
                      gamma_up(i,j,k,5) * q(3) + gamma_up(i,j,k,6) * q(4)) /&
                      (W2 * h * U_prim(i,j,k,1))
                  U_prim(i,j,k,4) = (gamma_up(i,j,k,3) * q(2) + &
                      gamma_up(i,j,k,6) * q(3) + gamma_up(i,j,k,9) * q(4)) /&
                      (W2 * h * U_prim(i,j,k,1))
                  U_prim(i,j,k,5) = (h - 1.0d0) / gamma

                  if (U_prim(i,j,k,2) /= U_prim(i,j,k,2)) then
                      write(*,*) "p, q(1), sq, W, h", p(i,j,k), q(2), sq, W2, h
                  end if

              end do
          end do
      end do

      write(*,*) "U_comp", U(lo(1)+4, lo(2)+4, lo(3)+4, :)
      write(*,*) "U_prim", U_prim(lo(1), lo(2), lo(3), :)
      !write(*,*) "p", p(lo(1)+3, lo(2)+3, lo(3)+3)
      !write(*,*) "gamma_up", gamma_up(lo(1), lo(2), lo(3), 7:9)
      !write(*,*) "alpha0, R, M", alpha0, M, R

  end subroutine cons_to_prim

  subroutine gr_comp_flux(U, f, lo, hi, Ncomp, dir, gamma, glo, ghi, alpha, dx, alpha0, M, R)
      implicit none

      integer, intent(in) :: Ncomp, dir
      integer, intent(in) :: lo(3), hi(3), glo(3), ghi(3)
      double precision, intent(inout)  :: U(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), Ncomp)
      double precision, intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), Ncomp)
      double precision, intent(in)  :: gamma, dx(3), alpha0, M, R
      double precision, intent(out) :: alpha(glo(1):ghi(1),glo(2):ghi(2), glo(3):ghi(3))

      integer :: i,j,k
      double precision :: p(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1)
      double precision :: U_prim(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, Ncomp)
      double precision :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9)
      double precision :: beta(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 3)
      !double precision, parameter :: M = 1.0d0, R = 100.0d0
      !double precision alpha0

      alpha = sqrt(1.0d0 - 2.0d0 * M / (R+1.0))
      beta = 0.0d0
      gamma_up = 0.0d0
      gamma_up(:,:,:,1) = 1.0d0
      gamma_up(:,:,:,5) = 1.0d0
      gamma_up(:,:,:,9) = alpha**2

      !alpha0 = sqrt(1.0d0 - 2.0d0 * M / R)

      write(*,*) "gr_comp_flux"
      write(*,*) "U: ", U(lo(1)+3, lo(2)+3, lo(3)+3, :)

      ! HACK: while the boundaries don't work
      !do k = lo(3)-1, hi(3)+1
        !  do j = lo(2)-1, hi(2)+1
        !      do i = lo(1)-1, hi(1)+1
        !          if (U(i,j,k,1) == 0.0d0) then
        !              U(i,j,k,1) = 1.0d0
        !          end if
        !          if (U(i,j,k,5) == 0.0d0 .or. U(i,j,k,5) == 1.0d0) then
        !              U(i,j,k,5) = 1.5d0
        !          end if
        !      end do
         ! end do
      !end do


      call cons_to_prim(U, glo, ghi, U_prim, lo-1, hi+1, p, lo-1, hi+1, lo-1, hi+1, Ncomp, gamma, alpha0, M, R, dx)

      f(:,:,:,:) = 0.0d0

      if (dir == 0) then
          do k = lo(3), hi(3)
              do j = lo(2), hi(2)
                  do i = lo(1)-1, hi(1)+1
                      f(i,j,k,1) = U(i,j,k,1) * (U_prim(i,j,k,2) - &
                           beta(i,j,k,1) / alpha(i,j,k))
                      f(i,j,k,2) = U(i,j,k,2) * (U_prim(i,j,k,2) - &
                           beta(i,j,k,1) / alpha(i,j,k)) + p(i,j,k)
                      f(i,j,k,3) = U(i,j,k,3) * (U_prim(i,j,k,2) - &
                           beta(i,j,k,1) / alpha(i,j,k))
                      f(i,j,k,4) = U(i,j,k,4) * (U_prim(i,j,k,2) - &
                           beta(i,j,k,1) / alpha(i,j,k))
                      f(i,j,k,5) = U(i,j,k,5) * (U_prim(i,j,k,2) - &
                           beta(i,j,k,1) / alpha(i,j,k)) + p(i,j,k) * U_prim(i,j,k,2)

                      !f(i,j,k,:) = f(i,j,k,:)
                   end do
               end do
           end do
      else if (dir == 1) then
          do k = lo(3), hi(3)
              do j = lo(2)-1, hi(2)+1
                  do i = lo(1), hi(1)
                      !write(*,*) "y, alpha", alpha(lo(1),lo(2),lo(3))
                      f(i,j,k,1) = U(i,j,k,1) * (U_prim(i,j,k,3) - &
                           beta(i,j,k,2) / alpha(i,j,k))
                      f(i,j,k,2) = U(i,j,k,2) * (U_prim(i,j,k,3) - &
                          beta(i,j,k,2) / alpha(i,j,k))
                      f(i,j,k,3) = U(i,j,k,3) * (U_prim(i,j,k,3) - &
                           beta(i,j,k,2) / alpha(i,j,k)) + p(i,j,k)
                      f(i,j,k,4) = U(i,j,k,4) * (U_prim(i,j,k,3) - &
                           beta(i,j,k,2) / alpha(i,j,k))
                      f(i,j,k,5) = U(i,j,k,5) * (U_prim(i,j,k,3) - &
                           beta(i,j,k,2) / alpha(i,j,k)) + p(i,j,k) * U_prim(i,j,k,3)

                      !f(i,j,k,:) = f(i,j,k,:)
                   end do
               end do
           end do
       else
           do k = lo(3)-1, hi(3)+1
               do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                       !write(*,*) "y, alpha", alpha(lo(1),lo(2),lo(3))
                       f(i,j,k,1) = U(i,j,k,1) * (U_prim(i,j,k,4) - &
                            beta(i,j,k,3) / alpha(i,j,k))
                       f(i,j,k,2) = U(i,j,k,2) * (U_prim(i,j,k,4) - &
                           beta(i,j,k,3) / alpha(i,j,k))
                       f(i,j,k,3) = U(i,j,k,3) * (U_prim(i,j,k,4) - &
                            beta(i,j,k,3) / alpha(i,j,k))
                       f(i,j,k,4) = U(i,j,k,4) * (U_prim(i,j,k,4) - &
                            beta(i,j,k,3) / alpha(i,j,k)) + p(i,j,k)
                       f(i,j,k,5) = U(i,j,k,5) * (U_prim(i,j,k,4) - &
                            beta(i,j,k,2) / alpha(i,j,k)) + p(i,j,k) * U_prim(i,j,k,4)

                       !f(i,j,k,:) = f(i,j,k,:)
                    end do
                end do
            end do

      end if
  end subroutine gr_comp_flux

end module compute_flux_module
