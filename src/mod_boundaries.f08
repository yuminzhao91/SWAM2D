!------------------------------------------------------------------------------
! LICENSE
!------------------------------------------------------------------------------
! This file is part of SWAM2D.
!
! SWAM2D is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! SWAM2D is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with SWAM2D. If not, see <http://www.gnu.org/licenses/>.
!------------------------------------------------------------------------------
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
       pmlx0, pmlx1, pmlz0, pmlz1, isurf)

    ! ADD IF ISURF HERE
    ! DO NOT CONSTRUCT TOP PML IN CASE OF FREE SURFACE.
    integer :: n1e, n2e, npml, ppml, i, isurf
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
       if(isurf == 1)then
          pmlz0(i, :) = d0*(val0/l)**2
          pmlz1(i, :) = d0*(val1/l)**2
       endif
       pmlz0(n1e+1-i,:) = d0*(val0/l)**ppml
       !pmlz1(i, :) = d0*(val1/l)**2
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

  subroutine dirichlet(n1e, n2e, uxx, uxz, uzx, uzz)
    ! implement Dirichlet boundary conditions on the four edges of the grid
    integer :: n1e, n2e
    real, dimension(n1e,n2e) :: uxx, uxz, uzx, uzz

    uxx(1, :) = 0.
    uxx(n1e, :)= 0.
    uxx(:, 1) = 0.
    uxx(:, n2e)= 0.

    uxz(1, :) = 0.
    uxz(n1e, :)= 0.
    uxz(:, 1) = 0.
    uxz(:, n2e)= 0.
    
    uzx(1, :) = 0.
    uzx(n1e, :)= 0.
    uzx(:, 1) = 0.
    uzx(:, n2e)= 0.
    
    uzz(1, :) = 0.
    uzz(n1e, :)= 0.
    uzz(:, 1) = 0.
    uzz(:, n2e)= 0.
    
  end subroutine dirichlet
  
end module boundaries
