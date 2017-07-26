
subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  double precision, intent(in) :: problo(*), probhi(*)

  ! nothing needs to be done here

end subroutine amrex_probinit


subroutine initdata(level, time, lo, hi, &
     phi, phi_lo, phi_hi, &
     dx, prob_lo, Ncomp, alpha0, M, R, p, nlayers, gamma, rho_swe) bind(C, name="initdata")

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), phi_lo(3), phi_hi(3), Ncomp, nlayers
  double precision, intent(in) :: time
  double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1), &
                                        phi_lo(2):phi_hi(2), &
                                        phi_lo(3):phi_hi(3), Ncomp)
  double precision, intent(in) :: dx(3), prob_lo(3), alpha0, M, R, p(nlayers), gamma
  double precision, intent(inout) :: rho_swe(nlayers)

  integer          :: dm
  integer          :: i,j,k
  double precision :: x, y, z, r2, alpha, z_surface, gamma_surf
  double precision :: Kk, rho_ref, rho, rhoh, cs, h, gamma_z
  double precision :: rho_z, rhoh_temp, cs_temp, W, u

  if (phi_lo(3) .eq. 0 .and. phi_hi(3) .eq. 0) then
     dm = 2
  else
     dm = 3
  end if

  phi(:,:,:,:) = 0.0d0

  z_surface = 1.0d0
  Kk = 1.0d0
  rho_ref = 1.0d0

  gamma_surf = (1.0d0 - M * z_surface / (R*alpha0)**2) / alpha0

  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k = phi_lo(3),phi_hi(3)
     do j = phi_lo(2),phi_hi(2)
        z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
        do i = phi_lo(1),phi_hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

           if ( dm.eq. 2) then
              r2 = ((x-0.5d0)**2 + (y-0.75d0)**2) / 0.01d0
           else
              r2 = ((x-0.5d0)**2 + (y-0.75d0)**2 + (z-0.5d0)**2) / 0.01d0
           end if
           alpha = sqrt(1.0d0 - 2.0d0 * M / (1.0d0 + R))!exp(-r2) + R)) !alpha0 + M * (1.0d0 + exp(-r2)) / (R**2 * alpha0)
           phi(i,j,k,1) = 1.d0! + exp(-r2)
           if (Ncomp == 4 .and. dm == 2) then ! compressible
               phi(i,j,k,1) = 1.d0! + exp(-r2)
               phi(i,j,k,4) = phi(i,j,k,1)
           else if (Ncomp == 5) then ! compressible
               !phi(i,j,k,1) = 1.d0 + exp(-r2)
               !phi(i,j,k,5) = phi(i,j,k,1)

               gamma_z = (1.0d0 - M * z / (R*alpha0)**2) / alpha0

               rho = rho_ref * (gamma_z - gamma_surf)**(1.0d0/gamma)
               rhoh = rho + gamma * Kk * rho**gamma / (gamma - 1.0d0)
               W = 1.0d0
               phi(i,j,k,1) = rho * W
               phi(i,j,k,2:4) = 0.0d0
               phi(i,j,k,5) = rhoh * W**2 - Kk * rho**gamma - phi(i,j,k,1)
           else !swe
               ! height
               if (k < 1) then
                   h = R*alpha0**2 / M * (1.0d0 - alpha0/R * (p(1) / (Kk * rho_ref**gamma) + gamma_surf))
               else if (k > nlayers) then
                   h = R*alpha0**2 / M * (1.0d0 - alpha0/R * gamma_surf)
               else
                   h = R*alpha0**2 / M * (1.0d0 - alpha0/R * (p(k) / (Kk * rho_ref**gamma) + gamma_surf))
               end if

               !write(*,*) "alpha0 = ", alpha0, "gamma = ", gamma, "gamma_surf = ", gamma_surf, "h = ", h

               if (k > 0 .and. k <= nlayers) then
                   rho_swe(k) = p(k)**(1.0d0 / gamma) / Kk

                   gamma_z = p(k) + gamma_surf

                   z = (1.0d0 - gamma_z * alpha0) * (R*alpha0)**2 / M
               else if (k < 1) then
                   gamma_z = p(1) + gamma_surf

                   z = (1.0d0 - gamma_z * alpha0) * (R*alpha0)**2 / M
               else
                   gamma_z = gamma_surf

                   z = (1.0d0 - gamma_z * alpha0) * (R*alpha0)**2 / M
               end if

               phi(i,j,k,1) = -log(z * M / (alpha0 * R) + alpha0)
               !-0.5d0 * log(1.0d0 - 2.0d0 / h)
               !-log(alpha)

               !rho = (p(k) / K)**(1.0d0/gamma)
               !rhoh = rho + gamma * p(k) / (gamma - 1.0d0)
               !cs = sqrt(gamma * p(k) / rhoh)
               !gamma_z = (1.0d0 - M * h / (R*alpha0)**2) / alpha0
               !rho_z = rho_ref * (gamma_z - gamma_surf)**(1.0d0/gamma)
               !rhoh_temp = rho_z + gamma * Kk * rho_z**gamma / (gamma - 1.0d0)
               !cs_temp = sqrt(gamma * Kk *rho_z**gamma / rhoh_temp)
               !J = -log((sqrt(gamma - 1.0d0) + cs_temp) / (sqrt(gamma - 1.0d0) - cs_temp)) / sqrt(gamma - 1.0d0)
               !a = log((sqrt(gamma - 1.0d0) + cs) / (sqrt(gamma - 1.0d0) - cs)) / sqrt(gamma - 1.0d0)
               u = 0.0d0 !(exp(2.0d0 * (J + a)) - 1.0d0) / (1.0d0 + exp(2.0d0 * (J + a)))

               W = 1.0d0 !/ sqrt(1.0d0 - u**2)

               phi(i,j,k,1) = phi(i,j,k,1) * W
               phi(i,j,k,2) = phi(i,j,k,1) * u * W
               phi(i,j,k,3) = phi(i,j,k,1) * u * W

           end if
        end do
     end do
  end do
  !$omp end parallel do

  write(*,*) phi(lo(1), lo(2), lo(3), :)

end subroutine initdata
