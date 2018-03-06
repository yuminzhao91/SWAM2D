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
! MODULE: types
!------------------------------------------------------------------------------
!> \brief types
!------------------------------------------------------------------------------
!> \author Damien Pageot
!> \date 09 Jan 2017
!------------------------------------------------------------------------------
! Revision history
! 09 Jan 2017: Add doxygen documentation
!------------------------------------------------------------------------------

module types
  
  implicit none

  !* Run parameters
  type typerun
     integer              :: isnap
     real                 :: tmax, dt, dtsnap
     character(len=80)    :: frun
  end type typerun
  
  !* Parameters related to physical parameter and grid size
  type typemod
     integer              :: n1, n2
     integer              :: n1e, n2e
     real                 :: h
     real,    allocatable :: vp(:, :), vs(:, :), ro(:, :)
     real,    allocatable :: vpe(:, :), vse(:, :), roe(:, :)
     real,    allocatable :: bux(:, :), buz(:, :)
     real,    allocatable :: mu0(:, :), mue(:, :)
     real,    allocatable :: lb0(:, :), lbmu(:, :)
     character(len=80)    :: fvp, fvs, fro
  end type typemod

  !* Parameters related to boundary conditions 
  type typebnd
     integer              :: isurf
     integer              :: npml, ppml
     real                 :: apml
     real, allocatable    :: pmlx0(:, :), pmlx1(:, :)
     real, allocatable    :: pmlz0(:, :), pmlz1(:, :)
  end type typebnd

  !* Parameters related to acquisition
  type typeacq
     integer              :: srctype, srcfunc
     integer              :: ixs, izs, nrec
     real                 :: dts, sigma, f0, t0, xs, zs
     real, allocatable    :: tsrc(:), gsrc(:, :)
     character(len=80)    :: facqui
  end type typeacq
  
end module types
