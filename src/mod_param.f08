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
! MODULE: param
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

module param
  
  use types

  implicit none
  
contains

  subroutine parread (trun, tmod, tbnd, tacq)

    type(typerun) :: trun
    type(typemod) :: tmod
    type(typebnd) :: tbnd
    type(typeacq) :: tacq

    open (101, file='param.input', status='old')
    read (101, *) !#[run]
    read (101, *) trun%frun
    read (101, *) trun%tmax, trun%dt
    read (101, *) !#[materials]
    read (101, *) tmod%fvp, tmod%fvs, tmod%fro
    read (101, *) tmod%n1, tmod%n2, tmod%h
    read (101, *) !#[boundaries]
    read (101, *) tbnd%isurf, tbnd%npml, tbnd%apml, tbnd%ppml
    read (101, *) ! #[source]
    read (101, *) tacq%srctype, tacq%srcfunc, tacq%sigma
    read (101, *) tacq%f0, tacq%t0
    read (101, *) tacq%xs, tacq%zs 
    read (101, *) !#[receiver]
    read (101, *) tacq%facqui
    read (101, *) tacq%dts
    read (101, *) !#[snapshot]
    read (101, *) trun%isnap
    read (101, *) trun%dtsnap
    close (101)

  end subroutine parread

end module param
