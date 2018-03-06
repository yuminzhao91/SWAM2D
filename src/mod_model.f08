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
! MODULE: model
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

module model

  use types
  
  implicit none
  
contains

  subroutine modread (iunit, fname, n1, n2, v)
    !! subroutine: modread
    !> \brief read an input binary file (single precision).
    !> \param[in] fname name of the input physical parameter file
    !> \param[in] n1 The number of grid points in the first direction (z)
    !> \param[in] n2 The number of grid points in the second direction (x)
    !> \param[out] v single precision array of size (n1, n2) containing
    !> parameter values.
    integer :: n1, n2, iunit
    real :: v(n1, n2)
    character(len=*) :: fname

    open(iunit, file=fname, access='direct', recl=n1*n2*4)
    read(iunit, rec=1) v
    close(iunit)

  end subroutine modread

  subroutine modext(v, n1, n2, nsp, ve)
    !! subroutine: modext
    !> \brief extend parameter model with absording boundary condition
    !> (ABC) layers.
    !> \param[in] v physical parameter array of size (n1, n2)
    !> \param[in] n1 The number of grid points in the first direction (z)
    !> \param[in] n2 The number of grid points in the second direction (x)
    !> \param[in] nsp The number of grid points added for ABC layers 
    !> \param[out] ve physical parameter array of size (n1+2*nsp, n2+2*nsp)
    integer :: i, n1, n2, nsp
    real :: ve(n1+2*nsp, n2+2*nsp)
    real :: v(n1,n2)
    
    ve(nsp+1:n1+nsp,nsp+1:n2+nsp) = v(1:n1,1:n2)
  
    do i=1,nsp
       ve(:,i) = ve(:,nsp+1)
       ve(:,n2+nsp+i) = ve(:,n2+nsp)
    end do
  
    do i=1,nsp
       ve(i, :) = ve(nsp+1, :)
       ve(n1+nsp+i, :) = ve(n1+nsp, :)
    end do

    
  end subroutine modext

  
  subroutine modbuo (tmod)
    !! subroutine: modbuo
    !> \brief calculate buoyancies bux and buz from extend density
    !> model to fulfill staggered-grid 4th order finite-differences
    !> requierments.
    !> \param[in] roe extend density model (array of size (n1e, n2e))
    !> \param[in] n1e The number of grid points in the first direction (z)
    !> of the extended grid
    !> \param[in] n2 The number of grid points in the second direction (x)
    !> of the extended grid
    !> \param[out] bux buoyancy model for the x-direction
    !> \param[out] buz buoyancy model for the z-direction

    type(typemod) :: tmod
    
    integer :: i1, i2

    do i2=1,tmod%n2e-1
       tmod%bux(:, i2) = (1./2.)*(1./tmod%roe(:, i2)+1./tmod%roe(:, i2+1))
    end do
    tmod%bux(:,tmod%n2e) = 1./tmod%roe(:,tmod%n2e)

    do i1=1,tmod%n1e-1
       tmod%buz(i1, :) = (1./2.)*(1./tmod%roe(i1, :)+1./tmod%roe(i1+1, :))
    end do
    tmod%buz(tmod%n1e,:) = 1./tmod%roe(tmod%n1e,:)

  end subroutine modbuo

  subroutine modlame (tmod) !(vpe, vse, roe, n1e, n2e, mu0, mue, lb0, lbmu)
    !! subroutine: modlame
    !> \brief calculate Lame's parameters \f$ \lambda \f$, \$f \mu \f$
    !> and \f$ mu_{e} \f$ from extend parameter models to fulfill
    !> staggered-grid 4th order finite-differences requierments.
    !> \param[in] vpe extend P-wave velocity model
    !> \param[in] vse extend S-wave velocity model
    !> \param[in] roe extend density model
    !> \param[in] n1e The number of grid points in the first direction (z)
    !> of the extended grid
    !> \param[in] n2 The number of grid points in the second direction (x)
    !> of the extended grid
    !> \param[out] mu0 \f$ \mu \f$ on normal grid (grid points)
    !> \param[out] mue \f$ \mu \f$ on staggered grid (intermediate gird points)
    !> \param[out] lb0 \f$ \lambda \f$ on normal grid (grid points)
    !> \param[out] lbmu \f$ \lambda+2\mu \f$ on normal grid (grid points)
    type(typemod) :: tmod
    integer :: i1, i2 !, n1e, n2e
    !real, dimension(n1e, n2e) :: vpe, vse, roe, mu0, mue, lb0, lbmu

    tmod%mu0(:, :) = tmod%vse(:, :)*tmod%vse(:, :)*tmod%roe(:, :)
  
    do i2=1,tmod%n2e-1
       do i1=1, tmod%n1e-1
          tmod%mue(i1, i2) = 1./((1./4.)* &
               (1./tmod%mu0(i1,i2)+ &
               1./tmod%mu0(i1+1,i2)+ &
               1./tmod%mu0(i1,i2+1)+1./tmod%mu0(i1+1,i2+1)))
       end do
    end do
    tmod%mue(:,tmod%n2e) = tmod%mue(:,tmod%n2e-1)
    tmod%mue(tmod%n1e,:) = tmod%mue(tmod%n1e-1,:)

    tmod%lb0(:, :) = tmod%vpe(:, :)*tmod%vpe(:, :)*tmod%roe(:, :) &
         -2.*tmod%mu0(:, :)
    
    tmod%lbmu(:,:) = tmod%lb0(:, :)+2.*tmod%mu0(:,:)

  end subroutine modlame

end module model
