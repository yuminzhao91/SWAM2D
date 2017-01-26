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
program main

  use types

  use param
  use model
  use source
  use boundaries
  use marching
  use acquisition

  implicit none
  
  integer :: n1e, n2e, is1, is2
  integer :: nt, nts, ntsnap
  integer :: nrec
  real    :: tmax, t0, f0, sigma, xs, zs
  real, allocatable :: tsrc(:), gsrc(:, :)
  real, allocatable :: pmlx0(:, :), pmlx1(:, :)
  real, allocatable :: pmlz0(:, :), pmlz1(:, :)
  real, allocatable :: recx(:, :), recz(:, :), recp(:, :)
  integer, allocatable :: recpos(:, :)

  integer :: n1, n2, isurf, npml, srctype, srcfunc, isnap
  real    :: dt, h, dts, apml, dtsnap
  integer :: ppml
  character(len=80) :: frun
  character(len=80) :: fvp, fvs, fro, facqui

  type(typemod) :: tmod

  !# >> Read input parameter file
  call parread(frun, tmax, dt, fvp, fvs, fro, n1, n2, &
       h, isurf, npml, apml, ppml, srctype, srcfunc, sigma, &
       f0, t0, xs, zs, facqui, dts, isnap, dtsnap)
  
  !# >> Get source position index on extend grid
  call srcindex(xs, zs, npml, h, is1, is2)

  nt = nint(tmax/dt)+1
  nts = nint(tmax/dts)+1
  ntsnap = nint(tmax/dtsnap)+1

  write(*, *) nt, nts, nt/nts+1, ntsnap

  !# >> Acquisition
  call acqread(facqui, npml, h, nrec, recpos)

  !# >> Allocate recorded seismogram arrays
  allocate(recx(nts, nrec), recz(nts, nrec), recp(nts, nrec))

  !# >> Extend model dimensions with PML
  n1e = n1+2*npml
  n2e = n2+2*npml
  
  !# >> Allocate model matrices
  allocate(tmod%vp(n1, n2), tmod%vs(n1, n2), tmod%ro(n1, n2))
  
  !# >> Read model files
  call modread(fvp, n1, n2, tmod%vp)
  call modread(fvs, n1, n2, tmod%vs)
  call modread(fro, n1, n2, tmod%ro)

  !# >> Extend models with pmls
  allocate(tmod%vpe(n1e, n2e), tmod%vse(n1e, n2e), tmod%roe(n1e, n2e))

  call modext(tmod%vp, n1, n2, npml, tmod%vpe)
  call modext(tmod%vs, n1, n2, npml, tmod%vse)
  call modext(tmod%ro, n1, n2, npml, tmod%roe)

  !# >> Buoyancy and Lame parameters
  allocate(tmod%bux(n1e, n2e), tmod%buz(n1e, n2e))
  allocate(tmod%mu0(n1e, n2e), tmod%mue(n1e, n2e))
  allocate(tmod%lb0(n1e, n2e), tmod%lbmu(n1e,n2e))
  
  call modbuo(tmod%roe, n1e, n2e, tmod%bux, tmod%buz)
  call modlame(tmod%vpe, tmod%vse, tmod%roe, n1e, n2e, tmod%mu0, tmod%mue, tmod%lb0, tmod%lbmu)

  !# >> Stability condition
  write(*, *) 'Courant::', dt*maxval(tmod%vp)/h

  allocate(pmlx0(n1e, n2e), pmlx1(n1e, n2e))
  allocate(pmlz0(n1e, n2e), pmlz1(n1e, n2e))

  !# >> Calculate PMLs
  call pmlmod(tmod%vpe, n1e, n2e, h, npml, apml, ppml, &
       pmlx0, pmlx1, pmlz0, pmlz1)
  
  !# >> Deallocate model matrices
  deallocate(tmod%vp, tmod%vs, tmod%ro)
  deallocate(tmod%vpe, tmod%vse, tmod%roe)

  !# >> Source
  allocate(tsrc(nt), gsrc(n1e, n2e)) 
  call ricker(nt, dt, f0, t0, tsrc)

  !# srcspread strongly affect results remove and replace with single point gsrc
  !call srcspread(n1e, n2e, npml, is1, is2, h, gsrc, sigma)
  gsrc(:, :) = 0.
  gsrc(is1, is2) = 1.

  write(*, * ) 'EVOLUTION'
  call evolution(n1e, n2e, npml, h, dt, nt, nts, ntsnap, nrec, srctype, tsrc, &
       gsrc, recx, recz, recp, recpos, tmod, isurf, &
       pmlx0, pmlx1, pmlz0, pmlz1, isnap)

  !# write seismos
  write(*, * ) nts, nrec
  open(41, file='recx.bin', access='direct', recl=nts*nrec*4)
  write(41, rec=1) recx
  close(41)
  open(42, file='recz.bin', access='direct', recl=nts*nrec*4)
  write(42, rec=1) recz
  close(42)
  open(43, file='recp.bin', access='direct', recl=nts*nrec*4)
  write(43, rec=1) recp
  close(43)

  call acqfree(recpos)
  deallocate(recx, recz, recp)

  deallocate(tsrc, gsrc)
  deallocate(tmod%bux, tmod%buz)
  deallocate(tmod%mu0, tmod%mue)
  deallocate(tmod%lb0, tmod%lbmu)

  deallocate(pmlx0, pmlx1, pmlz0, pmlz1)
  
end program main
