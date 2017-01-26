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

  type typemod
     real,    allocatable :: vp(:, :), vs(:, :), ro(:, :)
     real,    allocatable :: vpe(:, :), vse(:, :), roe(:, :)
     real,    allocatable :: bux(:, :), buz(:, :)
     real,    allocatable :: mu0(:, :), mue(:, :)
     real,    allocatable :: lb0(:, :), lbmu(:, :)
  end type typemod

end module types
