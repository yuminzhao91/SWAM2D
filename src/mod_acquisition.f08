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
! MODULE: acquisition
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

module acquisition

  implicit none

contains

  subroutine acqread(fname, nsp, h, nrec, recp)
    
    integer              :: irec, nrec, nsp, stat
    character(len=*)     :: fname
    real :: xr, zr, h
    integer,    allocatable :: recp(:, :)

    nrec = 0
    open(11, file=fname, status='old')
    do while(.true.)
       read(11, *, iostat=stat)
       !# if EOF: exit while loop
       if(stat /= 0) exit
       !# else: increment receiver number
       nrec = nrec+1
    end do

    rewind(11)

    allocate(recp(nrec, 2))
    
    do irec=1,nrec
       read(11, *) xr, zr
       recp(irec, 1) = nint(xr/h)+nsp+1
       recp(irec, 2) = nint(zr/h)+nsp+1
    end do

    close(11)

  end subroutine acqread
    
  subroutine acqfree(recp)

    integer, allocatable :: recp(:, :)

    deallocate(recp)

  end subroutine acqfree

end module acquisition
