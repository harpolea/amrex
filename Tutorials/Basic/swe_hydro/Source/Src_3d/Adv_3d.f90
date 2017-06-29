
subroutine advect(time, lo, hi, &
     &            uin , ui_lo, ui_hi, & ! state in
     &            uout, uo_lo, uo_hi, & ! state out
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            flxz, fz_lo, fz_hi, &
     &            dx,dt,Ncomp,gr) bind(C, name="advect")

  use mempool_module, only : bl_allocate, bl_deallocate
  use compute_flux_module, only : compute_flux_3d

  implicit none

  integer, intent(in) :: lo(3), hi(3), Ncomp
  double precision, intent(in) :: dx(3), dt, time
  integer, intent(in) :: ui_lo(3), ui_hi(3)
  integer, intent(in) :: uo_lo(3), uo_hi(3)
  integer, intent(in) :: fx_lo(3), fx_hi(3)
  integer, intent(in) :: fy_lo(3), fy_hi(3)
  integer, intent(in) :: fz_lo(3), fz_hi(3)
  double precision, intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3),Ncomp)
  double precision, intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),Ncomp)
  double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3), Ncomp)
  double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3), Ncomp)
  double precision, intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3), Ncomp)
  logical, intent(in) :: gr

  integer :: i, j, k
  integer :: glo(3), ghi(3)
  double precision :: dtdx(3)

  ! Some compiler may not support 'contiguous'.  Remove it in that case.
  double precision, dimension(:,:,:,:), pointer, contiguous :: &
       phi_p, phi_m, fp, fm, slope, u1

  dtdx = dt/dx

  glo = lo - 4
  ghi = hi + 4

  ! edge states
  call bl_allocate(phi_p, glo, ghi, Ncomp)
  call bl_allocate(phi_m, glo, ghi, Ncomp)
  call bl_allocate(fp   , glo, ghi, Ncomp)
  call bl_allocate(fm   , glo, ghi, Ncomp)
  ! slope
  call bl_allocate(slope  , glo, ghi, Ncomp)

  call bl_allocate(u1  , glo, ghi, Ncomp)

  ! We like to allocate these **pointers** here and then pass them to a function
  ! to remove their pointerness for performance, because normally pointers could
  ! be aliasing.  We need to use pointers instead of allocatable arrays because
  ! we like to use BoxLib's bl_allocate to allocate memeory instead of the intrinsic
  ! allocate.  Bl_allocate is much faster than allocate inside OMP.
  ! Note that one MUST CALL BL_DEALLOCATE.

  ! call a function to compute flux
  call compute_flux_3d(lo-3, hi+3, dt, dx, &
                       uin, ui_lo, ui_hi, &
                       flxx, fx_lo, fx_hi, &
                       flxy, fy_lo, fy_hi, &
                       flxz, fz_lo, fz_hi, &
                       phi_p, phi_m, fp, fm, slope, glo, ghi, Ncomp, gr)

  ! Do a conservative update
  do       k = lo(3)-3, hi(3)+3
     do    j = lo(2)-3, hi(2)+3
        do i = lo(1)-3, hi(1)+3
           u1(i,j,k,:) = uin(i,j,k,:) + &
                  flxx(i,j,k,:) * dtdx(1) + flxy(i,j,k,:) * dtdx(2) &
                + flxz(i,j,k,:) * dtdx(3)
        enddo
     enddo
  enddo

  !write(*,*) "uin:  ", uin(lo(1)+1, lo(2)+1,lo(3)+1,:)
  !write(*,*) "flxx: ", flxx(lo(1)+1, lo(2)+1,lo(3)+1,:)
  !write(*,*) "flxy: ", flxy(lo(1)+1, lo(2)+1,lo(3)+1,:)
  !write(*,*) "flxz: ", flxz(lo(1)+1, lo(2)+1,lo(3)+1,:)
  !write(*,*) "u1:   ", u1(lo(1)+1, lo(2)+1,lo(3)+1,:)

  call compute_flux_3d(lo, hi, dt, dx, &
                       u1, glo, ghi, &
                       flxx, fx_lo, fx_hi, &
                       flxy, fy_lo, fy_hi, &
                       flxz, fz_lo, fz_hi, &
                       phi_p, phi_m, fp, fm, slope, glo, ghi, Ncomp, gr)

   ! Do a conservative update
   do       k = lo(3), hi(3)
      do    j = lo(2), hi(2)
         do i = lo(1), hi(1)
            uout(i,j,k,:) = 0.5d0 * (uin(i,j,k,:) + u1(i,j,k,:) + &
                   flxx(i,j,k,:) * dtdx(1) + flxy(i,j,k,:) * dtdx(2) &
                 + flxz(i,j,k,:) * dtdx(3) )
         enddo
      enddo
   enddo

   !write(*,*) "flxx: ", flxx(lo(1)+1, lo(2)+1,lo(3)+1,:)
   !write(*,*) "flxy: ", flxy(lo(1)+1, lo(2)+1,lo(3)+1,:)
   !write(*,*) "flxz: ", flxz(lo(1)+1, lo(2)+1,lo(3)+1,:)
   !write(*,*) "uout: ", uout(lo(1)+1, lo(2)+1,lo(3)+1,:)

  ! Scale by face area in order to correctly reflx
  !do       k = lo(3), hi(3)
    ! do    j = lo(2), hi(2)
    !    do i = lo(1), hi(1)+1
    !       flxx(i,j,k,:) = flxx(i,j,k,:) * (dt * dx(2)*dx(3))
    !    enddo
     !enddo
  !enddo
  !do       k = lo(3), hi(3)
    ! do    j = lo(2), hi(2)+1
    !    do i = lo(1), hi(1)
    !       flxy(i,j,k,:) = flxy(i,j,k,:) * (dt * dx(1)*dx(3))
    !    enddo
     !enddo
  !enddo
  !do       k = lo(3), hi(3)+1
    ! do    j = lo(2), hi(2)
    !    do i = lo(1), hi(1)
    !       flxz(i,j,k,:) = flxz(i,j,k,:) * (dt * dx(1)*dx(2))
    !    enddo
     !enddo
  !enddo

  call bl_deallocate(phi_p)
  call bl_deallocate(phi_m)
  call bl_deallocate(fp)
  call bl_deallocate(fm)
  call bl_deallocate(slope)
  call bl_deallocate(u1)

end subroutine advect
