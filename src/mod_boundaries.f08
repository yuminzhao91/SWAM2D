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

  subroutine pmlmod(vpe, n1e, n2e, h, npml, apml, ppml, &
       pmlx0, pmlx1, pmlz0, pmlz1)

    integer :: n1e, n2e, npml, ppml, i
    real :: apml, val0, val1, val2
    real, dimension(n1e, n2e) :: vpe, pmlx0, pmlx1, pmlz0, pmlz1

    real :: r, vpmax, d0, l, h
    
    ! >> Initialize pml arrays
    pmlx0(:, :) = 0.
    pmlx1(:, :) = 0.
    pmlz0(:, :) = 0.
    pmlz1(:, :) = 0.

    r = 0.001
    vpmax = maxval(vpe)
    l = float(npml-1)*h
    d0 = apml*vpmax*log(1./r)/(2*l)
    
    do i=1,npml+1
       val0 = float(npml-i+1)*h
       val1 = float(npml-i+1)*h-(h/2.)
       val2 = float(npml-i+1)*h+(h/2.)
       pmlz0(i, :) = d0*(val0/l)**2
       pmlz0(n1e+1-i,:) = d0*(val0/l)**ppml
       pmlz1(i, :) = d0*(val1/l)**2
       pmlz1(n1e+1-i, :) = d0*(val2/l)**ppml
       pmlx0(:, i) = d0*(val0/l)**2
       pmlx0(:, n2e+1-i) = d0*(val0/l)**ppml
       pmlx1(:, i) = d0*(val1/l)**2
       pmlx1(:, n2e+1-i) = d0*(val2/l)**ppml
    enddo

    pmlx0(:, :) = pmlx0(:, :)/2.
    pmlx1(:, :) = pmlx1(:, :)/2.
    pmlz0(:, :) = pmlz0(:, :)/2.
    pmlz1(:, :) = pmlz1(:, :)/2.
    
  end subroutine pmlmod

  
  subroutine sponges(nsp, fac, h, spg)
    
    integer :: i, nsp
    real :: fac, h, spg(3, nsp+1), q

    fac = acos(fac)/(float(nsp)*h)

    spg(:, : ) = 0.

    do i=1,nsp+1
       spg(1, i) = cos(fac*(float(nsp-i+1)*h))
       spg(2, i) = cos(fac*(float(nsp-i+1)*h-h/2.))
       spg(3, i) = cos(fac*(float(nsp-i+1)*h+h/2.))
       !q = float(nsp-i+1)*h
       !spg(1, i) = 1.-(5*(600./(float(nsp+1)*h))*(q/(float(nsp+1)*h))*(q/(float(nsp+1)*h)))
       !spg(2, i) = 1.-(5*(600./(float(nsp+1)*h))*(q/(float(nsp+1)*h))*(q/(float(nsp+1)*h)))
       !spg(3, i) = 1.-(5*(600./(float(nsp+1)*h))*(q/(float(nsp+1)*h))*(q/(float(nsp+1)*h)))
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
