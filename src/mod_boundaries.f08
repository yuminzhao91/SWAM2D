!------------------------------------------------------------------------------
! MODULE: boundaries
!------------------------------------------------------------------------------
!> \brief subroutines to read, to extend and calculate physical parameter
!> models.
!------------------------------------------------------------------------------
!> \author Damien Pageot
!> \date 09 Jan 2017
!------------------------------------------------------------------------------
! Revision history
! 09 Jan 2017: Add doxygen documentation
!------------------------------------------------------------------------------

module boundaries

  implicit none

contains

  subroutine sponges(nsp, fac, h, spg)
    
    integer :: i, nsp
    real :: fac, h, spg(3, nsp+1), q

    fac = acos(fac)/(float(nsp)*h)
    
    spg(:, : ) = 0.

    do i=1,nsp+1
       spg(1, i) = cos(fac*(float(nsp-i+1)*h))
       spg(2, i) = cos(fac*(float(nsp-i+1)*h-h/2.))
       spg(3, i) = cos(fac*(float(nsp-i+1)*h+h/2.))
       !q = float(i-1)*h
       !spg(1, i) = (0.5*5*(300./(float(nsp+1)*h))*(q/(float(nsp+1)*h))*(q/(float(nsp+1)*h)))
       !spg(2, i) = (0.5*5*(300./(float(nsp+1)*h))*(q/(float(nsp+1)*h))*(q/(float(nsp+1)*h)))
       !spg(3, i) = (0.5*5*(300./(float(nsp+1)*h))*(q/(float(nsp+1)*h))*(q/(float(nsp+1)*h)))
    end do
    
    spg(2, nsp) = 1.

  end subroutine sponges

  subroutine addpmlx(n1e, n2e, nsp, spga, spgb, u)

    integer :: i, n1e, n2e, nsp
    real :: spga(nsp+1), spgb(nsp+1), u(n1e, n2e)

    do i=1,nsp+1
       u(:, i) = u(:, i) * spga(i)
       u(:, n2e-i) = u(:, n2e-i) * spgb(i)
    end do

  end subroutine addpmlx

  subroutine addpmlz(n1e, n2e, nsp, spga, spgb, u, isurf)

    integer :: i, n1e, n2e, nsp, isurf
    real :: spga(nsp+1), spgb(nsp+1), u(n1e, n2e)

    do i=1,nsp+1
       if(isurf == 0) u(i, :) = u(i, :) *spga(i)
       u(n1e-i, :) = u(n1e-i, :) * spgb(i)
    end do

  end subroutine addpmlz

end module boundaries
