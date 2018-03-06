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

  integer :: nt, nts, ntsnap
  real, allocatable :: recx(:, :), recz(:, :), recp(:, :)
  integer, allocatable :: recpos(:, :)

  type(typerun) :: trun
  type(typemod) :: tmod
  type(typebnd) :: tbnd
  type(typeacq) :: tacq
  
  !# >> Read input parameter file
  call parread (trun, tmod, tbnd, tacq)
  
  !# >> Get source position index on extend grid
  call srcindex (tmod, tbnd, tacq)

  nt = nint(trun%tmax/trun%dt)+1
  nts = nint(trun%tmax/tacq%dts)+1
  ntsnap = nint(trun%tmax/trun%dtsnap)+1

  write(*, *) nt, nts, nt/nts+1, ntsnap

  !# >> Acquisition
  write(*,*) 'acqusition'
  call acqread (tmod, tbnd, tacq, recpos)
  !(tacq%facqui, tbnd%npml, tmod%h, tacq%nrec, recpos)

  !# >> Allocate recorded seismogram arrays
  write(*,*) 'allocate 1'
  allocate (recx(nts, tacq%nrec), recz(nts, tacq%nrec), recp(nts, tacq%nrec))

  !# >> Extend model dimensions with PML
  tmod%n1e = tmod%n1+2*tbnd%npml
  tmod%n2e = tmod%n2+2*tbnd%npml
  
  !# >> Allocate model matrices
  write(*,*) 'allocate 2'
  allocate (tmod%vp(tmod%n1, tmod%n2))
  allocate (tmod%vs(tmod%n1, tmod%n2))
  allocate (tmod%ro(tmod%n1, tmod%n2))
  
  !# >> Read model files
  write(*,*) 'read model files'
  call modread (21, tmod%fvp, tmod%n1, tmod%n2, tmod%vp)
  call modread (22, tmod%fvs, tmod%n1, tmod%n2, tmod%vs)
  call modread (23, tmod%fro, tmod%n1, tmod%n2, tmod%ro)

  !# >> Extend models with pmls
  write(*,*) 'allocate 3'
  allocate (tmod%vpe(tmod%n1e, tmod%n2e))
  allocate (tmod%vse(tmod%n1e, tmod%n2e))
  allocate (tmod%roe(tmod%n1e, tmod%n2e))

  call modext (tmod%vp, tmod%n1, tmod%n2, tbnd%npml, tmod%vpe)
  call modext (tmod%vs, tmod%n1, tmod%n2, tbnd%npml, tmod%vse)
  call modext (tmod%ro, tmod%n1, tmod%n2, tbnd%npml, tmod%roe)
  
  !# >> Buoyancy and Lame parameters
  write(*,*) 'allocate 4'
  allocate (tmod%bux(tmod%n1e, tmod%n2e), tmod%buz(tmod%n1e, tmod%n2e))
  allocate (tmod%mu0(tmod%n1e, tmod%n2e), tmod%mue(tmod%n1e, tmod%n2e))
  allocate (tmod%lb0(tmod%n1e, tmod%n2e), tmod%lbmu(tmod%n1e,tmod%n2e))

  write(*,*) 'lame'
  call modbuo (tmod)
  call modlame (tmod)

  !# >> Stability condition
  write(*, *) trun%dt, maxval(tmod%vpe), tmod%h
  write(*, *) 'Courant::', trun%dt*maxval(tmod%vpe)/tmod%h
  
  !# >> Deallocate model matrices
  write(*,*) 'deallocate'
  deallocate (tmod%vp, tmod%vs, tmod%ro)
  deallocate (tmod%vpe, tmod%vse, tmod%roe)
  
  !# >> Allocate PMLs
  allocate (tbnd%pmlx0(tmod%n1e, tmod%n2e))
  allocate (tbnd%pmlx1(tmod%n1e, tmod%n2e))
  allocate (tbnd%pmlz0(tmod%n1e, tmod%n2e))
  allocate (tbnd%pmlz1(tmod%n1e, tmod%n2e))

  !# >> Calculate PMLs
  call pmlmod (tmod, tbnd)

  !# >> Source
  allocate (tacq%tsrc(nt))
  allocate (tacq%gsrc(tmod%n1e, tmod%n2e)) 
  call ricker (nt, trun%dt, tacq%f0, tacq%t0, tacq%tsrc)

  !# >> Spread source over several grid points
  call srcspread (tmod%n1e, tmod%n2e, tbnd%npml, tacq%izs, tacq%ixs, tmod%h, &
       tacq%gsrc, tacq%sigma)

  write( *, * ) 'ENTERING EVOLUTION'
  call evolution (nt, nts, ntsnap, &
       tacq%nrec, tacq%srctype, tacq%tsrc, &
       tacq%gsrc, recx, recz, recp, recpos, trun, tmod, tbnd, tacq, tbnd%isurf, trun%isnap)

  !# write seismos
  write(*, * ) nts, tacq%nrec
  open(41, file='recx.bin', access='direct', recl=nts*tacq%nrec*4)
  write(41, rec=1) recx
  close(41)
  open(42, file='recz.bin', access='direct', recl=nts*tacq%nrec*4)
  write(42, rec=1) recz
  close(42)
  open(43, file='recp.bin', access='direct', recl=nts*tacq%nrec*4)
  write(43, rec=1) recp
  close(43)

  call acqfree (recpos)
  deallocate (recx, recz, recp)

  deallocate (tacq%tsrc, tacq%gsrc)
  deallocate (tmod%bux, tmod%buz)
  deallocate (tmod%mu0, tmod%mue)
  deallocate (tmod%lb0, tmod%lbmu)

  deallocate (tbnd%pmlx0, tbnd%pmlx1)
  deallocate (tbnd%pmlz0, tbnd%pmlz1)
  
end program main
