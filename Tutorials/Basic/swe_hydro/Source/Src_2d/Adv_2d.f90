
subroutine advect(time, lo, hi, &
     &            uin , ui_lo, ui_hi, &
     &            uout, uo_lo, uo_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            dx,dt,Ncomp,gr) bind(C, name="advect")

  use mempool_module, only : bl_allocate, bl_deallocate
  use compute_flux_module, only : compute_flux_2d

  implicit none

  integer, intent(in) :: lo(2), hi(2), Ncomp
  double precision, intent(in) :: dx(2), dt, time
  integer, intent(in) :: ui_lo(2), ui_hi(2)
  integer, intent(in) :: uo_lo(2), uo_hi(2)
  integer, intent(in) :: fx_lo(2), fx_hi(2)
  integer, intent(in) :: fy_lo(2), fy_hi(2)
  double precision, intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),Ncomp)
  double precision, intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),Ncomp)
  double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),Ncomp)
  double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),Ncomp)
  logical, intent(in) :: gr

  integer :: i, j
  integer :: glo(2), ghi(2)
  double precision :: dtdx(2)

  ! Some compiler may not support 'contiguous'.  Remove it in that case.
  double precision, dimension(:,:,:), pointer, contiguous :: phi_p, phi_m, fp, fm, slope, u1

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

  !write(*,*) "uin", uin(lo(1),lo(2),:)
  !write(*,*) ui_lo, ui_hi, lo, hi, fx_lo, fx_hi

  ! call a function to compute flux
  call compute_flux_2d(lo-3, hi+3, dt, dx, &
                       uin, ui_lo, ui_hi, &
                       flxx, fx_lo, fx_hi, &
                       flxy, fy_lo, fy_hi, &
                       phi_p, phi_m, fp, fm, slope, glo, ghi, Ncomp, gr)


  ! Do a conservative update
  do    j = lo(2)-3,hi(2)+3
     do i = lo(1)-3,hi(1)+3
        u1(i,j,:) = uin(i,j,:) + &
             flxx(i,j,:)  * dtdx(1) + flxy(i,j,:)  * dtdx(2)
     enddo
  enddo

  !write(*,*)  "flxx", flxx(lo(1),lo(2),:), "flxy", flxx(lo(1),lo(2),:)
  !write(*,*)  "u1", u1(lo(1),lo(2),1)

  call compute_flux_2d(lo, hi, dt, dx, &
                       u1, glo, ghi, &
                       flxx, fx_lo, fx_hi, &
                       flxy, fy_lo, fy_hi, &
                       phi_p, phi_m, fp, fm, slope, glo, ghi, Ncomp, gr)

   do    j = lo(2),hi(2)
      do i = lo(1),hi(1)
         uout(i,j,:) = 0.5d0 * (uin(i,j,:) + u1(i,j,:) + &
              (flxx(i,j,:)  * dtdx(1) + flxy(i,j,:)  * dtdx(2) ))
      enddo
   enddo

  ! Scale by face area in order to correctly reflx
  !do    j = lo(2), hi(2)
    ! do i = lo(1), hi(1)+1
    !    flxx(i,j,:) = flxx(i,j,:) * (dt * dx(2))
     !enddo
  !enddo

  ! Scale by face area in order to correctly reflx
  !do    j = lo(2), hi(2)+1
    ! do i = lo(1), hi(1)
    !    flxy(i,j,:) = flxy(i,j,:) * (dt * dx(1))
     !enddo
  !enddo

  !write(*,*) "uout", uout(lo(1),lo(2),:)

  call bl_deallocate(phi_p)
  call bl_deallocate(phi_m)
  call bl_deallocate(fp)
  call bl_deallocate(fm)
  call bl_deallocate(slope)
  call bl_deallocate(u1)

end subroutine advect
