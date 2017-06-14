
subroutine get_face_velocity(level, time, &
     q, q_l1, q_l2, q_h1, q_h2, &
     vx, vx_l1, vx_l2, vx_h1, vx_h2, &
     vy, vy_l1, vy_l2, vy_h1, vy_h2, &
     Ncomp) bind(C, name="get_face_velocity")

  use probdata_module, only : adv_vel

  implicit none

  integer, intent(in) :: level, Ncomp
  double precision, intent(in) :: time
  integer, intent(in) :: q_l1, q_l2, q_h1, q_h2
  integer, intent(in) :: vx_l1, vx_l2, vx_h1, vx_h2
  integer, intent(in) :: vy_l1, vy_l2, vy_h1, vy_h2
  double precision, intent(in) :: q(q_l1:q_h1,q_l2:q_h2, Ncomp)
  double precision, intent(out) :: vx(vx_l1:vx_h1,vx_l2:vx_h2)
  double precision, intent(out) :: vy(vy_l1:vy_h1,vy_l2:vy_h2)

  vx = adv_vel(1)
  vy = adv_vel(2)

end subroutine get_face_velocity
