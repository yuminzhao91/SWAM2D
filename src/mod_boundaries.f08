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
  use types
  
  implicit none

contains

  subroutine pmlmod (tmod, tbnd)

    type(typemod) :: tmod
    type(typebnd) :: tbnd
    
    ! ADD IF ISURF HERE
    ! DO NOT CONSTRUCT TOP PML IN CASE OF FREE SURFACE.
    integer :: i 
    real :: val0, val1, val2

    real :: r, vpmax, d0, l !, h
    
    ! >> Initialize pml arrays
    tbnd%pmlx0(:, :) = 0.
    tbnd%pmlx1(:, :) = 0.
    tbnd%pmlz0(:, :) = 0.
    tbnd%pmlz1(:, :) = 0.

    r = 0.0001 !0.9
    vpmax = maxval(tmod%vpe)
    l = float(tbnd%npml-1)*tmod%h
    !d0 = tbnd%apml*vpmax*log(1./r)/(2*l)
    d0 = float(tbnd%ppml+1)*tbnd%apml*log(1./r)/(2.*l)
    
    do i=1,tbnd%npml+1
       val0 = float(tbnd%npml-i+1)*tmod%h
       val1 = float(tbnd%npml-i+1)*tmod%h-(tmod%h/2.)
       val2 = float(tbnd%npml-i+1)*tmod%h+(tmod%h/2.)
       if (tbnd%isurf == 1) then
          tbnd%pmlz0(i, :) = 0. !d0*(val0/l)**2
          tbnd%pmlz1(i, :) = 0. !d0*(val1/l)**2
       else
          tbnd%pmlz0(i, :) = d0*(val0/l)**tbnd%ppml !2
          tbnd%pmlz1(i, :) = d0*(val1/l)**tbnd%ppml !2
       endif
       !tbnd%pmlz0(i, :) = d0*(val0/l)**ppml
       !tbnd%pmlz1(i, :) = d0*(val1/l)**ppml
       tbnd%pmlz0(tmod%n1e+1-i,:) = d0*(val0/l)**tbnd%ppml
       tbnd%pmlz1(tmod%n1e+1-i, :) = d0*(val2/l)**tbnd%ppml
       tbnd%pmlx0(:, i) = d0*(val0/l)**tbnd%ppml !2 !ppml
       tbnd%pmlx0(:, tmod%n2e+1-i) = d0*(val0/l)**tbnd%ppml
       tbnd%pmlx1(:, i) = d0*(val1/l)**tbnd%ppml !2 !ppml
       tbnd%pmlx1(:, tmod%n2e+1-i) = d0*(val2/l)**tbnd%ppml
    enddo

    tbnd%pmlx0(:, :) = tbnd%pmlx0(:, :)/2.
    tbnd%pmlx1(:, :) = tbnd%pmlx1(:, :)/2.
    tbnd%pmlz0(:, :) = tbnd%pmlz0(:, :)/2.
    tbnd%pmlz1(:, :) = tbnd%pmlz1(:, :)/2.
    
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
