program main

  use types

  use param
  use model
  use source
  use boundaries
  use marching
  use acquisition

  implicit none
  
  integer :: i1, i2, n1e, n2e, is1, is2
  integer :: nt, nts, ntsnap
  integer :: nrec
  real    :: tmax, t0, f0, sigma, xs, zs
  real, allocatable :: tsrc(:), gsrc(:, :)
  real, allocatable :: pmlx0(:, :), pmlx1(:, :)
  real, allocatable :: pmlz0(:, :), pmlz1(:, :)
  real, allocatable :: spg(:,:)
  real, allocatable :: recx(:, :), recz(:, :)
  integer, allocatable :: recp(:, :)

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
  call acqread(facqui, npml, h, nrec, recp)
  allocate(recx(nts, nrec), recz(nts, nrec))

  !# >> Extend dimensions
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
  call pmlmod(tmod%vpe, n1e, n2e, h, npml, apml, ppml, &
       pmlx0, pmlx1, pmlz0, pmlz1)
  
  !# >> Deallocate model matrices
  deallocate(tmod%vp, tmod%vs, tmod%ro)
  deallocate(tmod%vpe, tmod%vse, tmod%roe)

  !# >> Source
  allocate(tsrc(nt), gsrc(n1e, n2e)) 
  call ricker(nt, dt, f0, t0, tsrc)
  call srcspread(n1e, n2e, npml, is1, is2, h, gsrc, sigma)

  allocate(spg(3, npml+1))

  spg(:, : ) = 1.

  write(*, * ) 'SPONGES'
  call sponges(npml, apml, h, spg)
  
  write(*, * ) 'EVOLUTION'
  call evolution(n1e, n2e, npml, h, dt, dts, dtsnap, nt, nts, ntsnap, nrec, srctype, tsrc, &
       gsrc, spg, recx, recz, recp, tmod, isurf, &
       pmlx0, pmlx1, pmlz0, pmlz1, isnap)

  !# write seismos
  write(*, * ) nts, nrec
  open(41, file='recx.bin', access='direct', recl=nt*nrec*4)
  write(41, rec=1) recx
  close(41)
  open(42, file='recz.bin', access='direct', recl=nt*nrec*4)
  write(42, rec=1) recz
  close(42)

  call acqfree(recp)
  deallocate(recx, recz)

  deallocate(tsrc, gsrc)
  deallocate(tmod%bux, tmod%buz)
  deallocate(tmod%mu0, tmod%mue)
  deallocate(tmod%lb0, tmod%lbmu)

  deallocate(pmlx0, pmlx1, pmlz0, pmlz1)
  
end program main
