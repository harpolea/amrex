
! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the state
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag        <=  integer tag array
! ::: tag_lo,hi   => index extent of tag array
! ::: state       => state array
! ::: state_lo,hi => index extent of state array
! ::: set         => integer value to tag cell for refinement
! ::: clear       => integer value to untag cell
! ::: lo,hi       => work region we are allowed to change
! ::: dx          => cell size
! ::: problo      => phys loc of lower left corner of prob domain
! ::: time        => problem evolution time
! ::: level       => refinement level of this array
! ::: -----------------------------------------------------------

subroutine state_error(tag,tag_lo,tag_hi, &
                       state,state_lo,state_hi, &
                       set,clear,&
                       lo,hi,&
                       dx,problo,time,phierr,Ncomp) bind(C, name="state_error")

  implicit none

  integer          :: lo(3),hi(3),Ncomp
  integer          :: state_lo(3),state_hi(3)
  integer          :: tag_lo(3),tag_hi(3)
  double precision :: state(state_lo(1):state_hi(1), &
                            state_lo(2):state_hi(2), &
                            state_lo(3):state_hi(3),Ncomp)
  integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
  double precision :: problo(3),dx(3),time, phierr
  integer          :: set,clear

  double precision :: ax, ay, az, r2
  integer          :: i, j, k, n, dim

  if (state_lo(3) .eq. state_hi(3)) then
     dim = 2
  else
     dim = 3
  end if

  ! Tag on regions of high phi
     do       k = lo(3), hi(3)
        do    j = lo(2), hi(2)
           do i = lo(1), hi(1)
               ! NOTE: changed to be less than for gr swe
               if (dim == 2) then
                  if ((Ncomp < 4 .and. state(i,j,k,1) .le. phierr) .or. (Ncomp >= 4 .and. state(i,j,k,1) .ge. phierr)) then
                     tag(i,j,k) = set
                  endif
              else
                  if ((Ncomp < 5 .and. state(i,j,k,1) .le. phierr) .or. (Ncomp >= 5 .and. state(i,j,k,1) .ge. phierr)) then
                     tag(i,j,k) = set
                  endif
              end if
           enddo
        enddo
     enddo
     tag = clear

     do       k = lo(3), hi(3)
        do    j = lo(2), hi(2)
           az = problo(3) + (dble(k)+0.5d0) * dx(3)
           ay = problo(2) + (dble(j)+0.5d0) * dx(2)
           do i = lo(1), hi(1)
               ax = problo(1) + (dble(i)+0.5d0) * dx(1)
               r2 = ((ax-0.5d0)**2 + (ay-0.75d0)**2 + (az-0.5d0)**2) / 0.005d0
               if ( 1.0d0 + exp(-r2) .ge. phierr) then
                   tag(i,j,k) = set
               end if
           enddo
        enddo
    enddo



end subroutine state_error
