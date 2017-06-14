
subroutine get_face_velocity(level, time, &
     q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
     vx, vx_l1, vx_l2, vx_l3, vx_h1, vx_h2, vx_h3, &
     vy, vy_l1, vy_l2, vy_l3, vy_h1, vy_h2, vy_h3, &
     vz, vz_l1, vz_l2, vz_l3, vz_h1, vz_h2, vz_h3, &
     Ncomp) bind(C, name="get_face_velocity")

  use mempool_module, only : bl_allocate, bl_deallocate

  implicit none

  integer, intent(in) :: level, Ncomp
  double precision, intent(in) :: time
  integer, intent(in) :: q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  integer, intent(in) :: vx_l1, vx_l2, vx_l3, vx_h1, vx_h2, vx_h3
  integer, intent(in) :: vy_l1, vy_l2, vy_l3, vy_h1, vy_h2, vy_h3
  integer, intent(in) :: vz_l1, vz_l2, vz_l3, vz_h1, vz_h2, vz_h3
  double precision, intent(out) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,Ncomp)
  double precision, intent(out) :: vx(vx_l1:vx_h1,vx_l2:vx_h2,vx_l3:vx_h3)
  double precision, intent(out) :: vy(vy_l1:vy_h1,vy_l2:vy_h2,vy_l3:vy_h3)
  double precision, intent(out) :: vz(vz_l1:vz_h1,vz_l2:vz_h2,vz_l3:vz_h3)

  integer :: i, j, k

  ! x velocity
  do k = vx_l3, vx_h3
      do j = vx_l2, vx_h2
         do i = vx_l1, vx_h1
            vx(i,j,k) =  0.5d0 * (q(i-1,j,k,2) + q(i,j,k,2))
         end do
      end do
  end do

  ! y velocity
  do k = vy_l3, vy_h3
      do j = vy_l2, vy_h2
         do i = vy_l1, vy_h1
            vy(i,j,k) = 0.5d0 * (q(i,j-1,k,3) + q(i,j,k,3))
         end do
      end do
  end do

  vz = 1.d0

end subroutine get_face_velocity
