
subroutine get_face_velocity(level, time, &
     q, q_lo, q_hi, &
     vx, vx_l1, vx_l2, vx_h1, vx_h2, &
     vy, vy_l1, vy_l2, vy_h1, vy_h2, &
     Ncomp) bind(C, name="get_face_velocity")

  use mempool_module, only : bl_allocate, bl_deallocate

  implicit none

  integer, intent(in) :: level, Ncomp
  double precision, intent(in) :: time
  integer, intent(in) :: q_lo(2), q_hi(2)
  integer, intent(in) :: vx_l1, vx_l2, vx_h1, vx_h2
  integer, intent(in) :: vy_l1, vy_l2, vy_h1, vy_h2
  double precision, intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2), Ncomp)
  double precision, intent(out) :: vx(vx_l1:vx_h1,vx_l2:vx_h2)
  double precision, intent(out) :: vy(vy_l1:vy_h1,vy_l2:vy_h2)

  integer :: i, j
  !double precision, parameter :: M_PI = 3.141592653589793238462643383279502884197d0

  ! x velocity
  do j = vx_l2, vx_h2
     do i = vx_l1, vx_h1
        vx(i,j) =  0.5d0 * (q(i-1,j,2)/q(i-1,j,1) + q(i,j,2)/q(i,j,1))
     end do
  end do

  ! y velocity
  do j = vy_l2, vy_h2
     do i = vy_l1, vy_h1
        vy(i,j) = 0.5d0 * (q(i,j-1,3)/q(i,j-1,1) + q(i,j,3)/q(i,j,1))
     end do
  end do

end subroutine get_face_velocity
