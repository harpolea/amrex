module amrex_basefab_nd_module

  use amrex_fort_module, only : amrex_real

  implicit none

contains

  ! dst = src
  subroutine amrex_fort_fab_copy(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
       bind(c,name='amrex_fort_fab_copy')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = sblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_fab_copy
    
  
  ! copy from multi-d array to 1d array
  function amrex_fort_fab_copytomem (lo, hi, dst, src, slo, shi, ncomp) result(nelems) &
       bind(c,name='amrex_fort_fab_copytomem')
    use iso_c_binding, only : c_long
    integer(c_long) :: nelems
    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp
    real(amrex_real)             :: dst(*)
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)

    integer :: i, j, k, n, nx
    integer(c_long) :: offset

    nx = hi(1)-lo(1)+1
    offset = 1-lo(1)
    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(offset+i) = src(i,j,k,n) 
             end do
             offset = offset + nx
          end do
       end do
    end do    

    nelems = offset - (1-lo(1))
  end function amrex_fort_fab_copytomem


  ! copy from 1d array to multi-d array
  function amrex_fort_fab_copyfrommem (lo, hi, dst, dlo, dhi, ncomp, src) result(nelems) &
       bind(c,name='amrex_fort_fab_copyfrommem')
    use iso_c_binding, only : c_long
    integer(c_long) :: nelems
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), ncomp
    real(amrex_real), intent(in   ) :: src(*)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i, j, k, n, nx
    integer(c_long) :: offset

    nx = hi(1)-lo(1)+1
    offset = 1-lo(1)
    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n)  = src(offset+i)
             end do
             offset = offset + nx
          end do
       end do
    end do    

    nelems = offset - (1-lo(1))
  end function amrex_fort_fab_copyfrommem
  

  subroutine amrex_fort_fab_setval(lo, hi, dst, dlo, dhi, ncomp, val) &
       bind(c,name='amrex_fort_fab_setval')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), ncomp
    real(amrex_real), intent(in) :: val
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i, j, k, n

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = val
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_fab_setval


  subroutine amrex_fort_fab_setval_ifnot (lo, hi, dst, dlo, dhi, ncomp, msk, mlo, mhi, val) &
       bind(c,name='amrex_fort_fab_setval_ifnot')
    real(amrex_real), intent(in) :: val
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), mlo(3), mhi(3), ncomp
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    
    integer :: i,j,k,n

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (msk(i,j,k) .eq. 0) then
                   dst(i,j,k,n) = val
                end if
             end do
          end do
       end do
    end do    
  end subroutine amrex_fort_fab_setval_ifnot
    

  function amrex_fort_fab_norminfmask (lo, hi, msk, mlo, mhi, src, slo, shi, ncomp) result(nrm) &
       bind(c,name='amrex_fort_fab_norminfmask')
    integer, intent(in) :: lo(3), hi(3), mlo(3), mhi(3), slo(3), shi(3), ncomp
    integer         , intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real) :: nrm

    integer :: i,j,k,n

    nrm = 0.0_amrex_real
    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (msk(i,j,k).eq.1) then
                   nrm = max(nrm, abs(src(i,j,k,n)))
                end if
             end do
          end do
       end do
    end do
  end function amrex_fort_fab_norminfmask


  function amrex_fort_fab_norm (lo, hi, src, slo, shi, ncomp, p) result(nrm) &
       bind(c,name='amrex_fort_fab_norm')
    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp, p
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real) :: nrm

    integer :: i,j,k,n

    nrm = 0.0_amrex_real
    if (p .eq. 0) then ! max norm
       do n = 1, ncomp
          do       k = lo(3), hi(3)
             do    j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   nrm = max(nrm, abs(src(i,j,k,n)))
                end do
             end do
          end do
       end do
    else if (p .eq. 1) then
       do n = 1, ncomp
          do       k = lo(3), hi(3)
             do    j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   nrm = nrm + abs(src(i,j,k,n))
                end do
             end do
          end do
       end do
    end if
  end function amrex_fort_fab_norm


  function amrex_fort_fab_sum (lo, hi, src, slo, shi, ncomp) result(sm) &
       bind(c,name='amrex_fort_fab_sum')
    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real) :: sm

    integer :: i,j,k,n

    sm = 0.0_amrex_real
    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                sm = sm + src(i,j,k,n)
             end do
          end do
       end do
    end do
  end function amrex_fort_fab_sum


  subroutine amrex_fort_fab_plus(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
       bind(c,name='amrex_fort_fab_plus')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = sblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = dst(i,j,k,n) + src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_fab_plus


  subroutine amrex_fort_fab_minus(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
       bind(c,name='amrex_fort_fab_minus')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = sblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = dst(i,j,k,n) - src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_fab_minus


  subroutine amrex_fort_fab_mult(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
       bind(c,name='amrex_fort_fab_mult')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = sblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = dst(i,j,k,n) * src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_fab_mult


  subroutine amrex_fort_fab_divide(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
       bind(c,name='amrex_fort_fab_divide')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = sblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = dst(i,j,k,n) / src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_fab_divide


  subroutine amrex_fort_fab_protdivide(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
       bind(c,name='amrex_fort_fab_protdivide')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = sblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (src(i+off(1),j+off(2),k+off(3),n) .ne. 0._amrex_real) then
                   dst(i,j,k,n) = dst(i,j,k,n) / src(i+off(1),j+off(2),k+off(3),n)
                end if
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_fab_protdivide


  ! dst = a/src
  subroutine amrex_fort_fab_invert(lo, hi, dst, dlo, dhi, ncomp, a) &
       bind(c,name='amrex_fort_fab_invert')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), ncomp
    real(amrex_real), intent(in   ) :: a
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = a / dst(i,j,k,n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_fab_invert


  ! dst += a*src
  subroutine amrex_fort_fab_saxpy(lo, hi, dst, dlo, dhi, a, src, slo, shi, sblo, ncomp) &
       bind(c,name='amrex_fort_fab_saxpy')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(amrex_real), intent(in   ) :: a
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = sblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = dst(i,j,k,n) + a * src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_fab_saxpy


  ! dst = src + a*dst
  subroutine amrex_fort_fab_xpay(lo, hi, dst, dlo, dhi, a, src, slo, shi, sblo, ncomp) &
       bind(c,name='amrex_fort_fab_xpay')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(amrex_real), intent(in   ) :: a
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = sblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = src(i+off(1),j+off(2),k+off(3),n) + a * dst(i,j,k,n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_fab_xpay


  ! dst = a*x + b*y
  subroutine amrex_fort_fab_lincomb(lo, hi, dst, dlo, dhi, a, x, xlo, xhi, xblo, &
       b, y, ylo, yhi, yblo, ncomp) bind(c,name='amrex_fort_fab_lincomb')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), xlo(3), xhi(3), xblo(3), &
         ylo(3), yhi(3), yblo(3), ncomp
    real(amrex_real), intent(in   ) :: a, b
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    real(amrex_real), intent(in   ) ::   x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3),ncomp)
    real(amrex_real), intent(in   ) ::   y(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3),ncomp)
    
    integer :: i,j,k,n,xoff(3),yoff(3)

    xoff = xblo - lo
    yoff = yblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = a * x(i+xoff(1),j+xoff(2),k+xoff(3),n) &
                     +         b * y(i+yoff(1),j+yoff(2),k+yoff(3),n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_fab_lincomb

  ! dst = dst + src1*src2
  subroutine amrex_fort_fab_addproduct(lo, hi, dst, dlo, dhi, src1, s1lo, s1hi, src2, s2lo, s2hi,ncomp) &
       bind(c,name='amrex_fort_fab_addproduct')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), s1lo(3), s1hi(3), s2lo(3), s2hi(3), ncomp
    real(amrex_real), intent(in   ) :: src1(s1lo(1):s1hi(1),s1lo(2):s1hi(2),s1lo(3):s1hi(3),ncomp)
    real(amrex_real), intent(in   ) :: src2(s2lo(1):s2hi(1),s2lo(2):s2hi(2),s2lo(3):s2hi(3),ncomp)
    real(amrex_real), intent(inout) ::  dst( dlo(1): dhi(1), dlo(2): dhi(2), dlo(3): dhi(3),ncomp)
    
    integer :: i,j,k,n

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = src1(i,j,k,n) * src2(i,j,k,n) + dst(i,j,k,n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_fab_addproduct
  
  ! dot_product
  function amrex_fort_fab_dot(lo, hi, x, xlo, xhi, y, ylo, yhi, yblo, ncomp) result(dp) &
       bind(c,name='amrex_fort_fab_dot')
    integer, intent(in) :: lo(3), hi(3), xlo(3), xhi(3), ylo(3), yhi(3), yblo(3), ncomp
    real(amrex_real), intent(in) :: x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3),ncomp)
    real(amrex_real), intent(in) :: y(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3),ncomp)
    real(amrex_real) :: dp

    integer :: i,j,k,n, off(3)

    dp = 0.0_amrex_real

    off = yblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dp = dp + x(i,j,k,n)*y(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end function amrex_fort_fab_dot

  ! dot_product
  function amrex_fort_fab_dot_mask(lo, hi, x, xlo, xhi, y, ylo, yhi, yblo, m, mlo, mhi, ncomp) result(dp) &
       bind(c,name='amrex_fort_fab_dot_mask')
    integer, intent(in) :: lo(3), hi(3), xlo(3), xhi(3), ylo(3), yhi(3), yblo(3), mlo(3), mhi(3), ncomp
    real(amrex_real), intent(in) :: x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3),ncomp)
    real(amrex_real), intent(in) :: y(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3),ncomp)
    integer         , intent(in) :: m(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    real(amrex_real) :: dp

    integer :: i,j,k,n, off(3)

    dp = 0.0_amrex_real

    off = yblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dp = dp + x(i,j,k,n)*y(i+off(1),j+off(2),k+off(3),n)*m(i,j,k)
             end do
          end do
       end do
    end do
  end function amrex_fort_fab_dot_mask

  ! dst = src
  subroutine amrex_fort_ifab_copy(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
       bind(c,name='amrex_fort_ifab_copy')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    integer, intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    integer, intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = sblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_ifab_copy

  ! copy from multi-d array to 1d array
  function amrex_fort_ifab_copytomem (lo, hi, dst, src, slo, shi, ncomp) result(nelems) &
       bind(c,name='amrex_fort_ifab_copytomem')
    use iso_c_binding, only : c_long
    integer(c_long) :: nelems
    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp
    integer             :: dst(*)
    integer, intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)

    integer :: i, j, k, n, nx
    integer(c_long) :: offset

    nx = hi(1)-lo(1)+1
    offset = 1-lo(1)
    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(offset+i) = src(i,j,k,n) 
             end do
             offset = offset + nx
          end do
       end do
    end do    

    nelems = offset - (1-lo(1))
  end function amrex_fort_ifab_copytomem


  ! copy from 1d array to multi-d array
  function amrex_fort_ifab_copyfrommem (lo, hi, dst, dlo, dhi, ncomp, src) result(nelems) &
       bind(c,name='amrex_fort_ifab_copyfrommem')
    use iso_c_binding, only : c_long
    integer(c_long) :: nelems
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), ncomp
    integer, intent(in   ) :: src(*)
    integer, intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i, j, k, n, nx
    integer(c_long) :: offset

    nx = hi(1)-lo(1)+1
    offset = 1-lo(1)
    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n)  = src(offset+i)
             end do
             offset = offset + nx
          end do
       end do
    end do    

    nelems = offset - (1-lo(1))
  end function amrex_fort_ifab_copyfrommem
  

  subroutine amrex_fort_ifab_setval(lo, hi, dst, dlo, dhi, ncomp, val) &
       bind(c,name='amrex_fort_ifab_setval')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), ncomp
    integer, intent(in) :: val
    integer, intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i, j, k, n

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = val
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_ifab_setval


  subroutine amrex_fort_ifab_plus(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
       bind(c,name='amrex_fort_ifab_plus')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    integer, intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    integer, intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = sblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = dst(i,j,k,n) + src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_ifab_plus


  subroutine amrex_fort_ifab_minus(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
       bind(c,name='amrex_fort_ifab_minus')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    integer, intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    integer, intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = sblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = dst(i,j,k,n) - src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_ifab_minus

end module amrex_basefab_nd_module
