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
  integer :: nt, nts
  integer :: nrec
  real    :: tmax, t0, f0, sigma, xs, zs
  real, allocatable :: tsrc(:), gsrc(:, :), spg(:,:)
  real, allocatable :: recx(:, :), recz(:, :)
  integer, allocatable :: recp(:, :)

  integer :: n1, n2, isurf, npml, srctype, srcfunc
  real    :: dt, h, dts, apml
  character(len=80) :: frun
  character(len=80) :: fvp, fvs, fro, facqui

  type(typemod) :: tmod

  !# >> Read input parameter file
  call parread(frun, tmax, dt, fvp, fvs, fro, n1, n2, &
       h, isurf, npml, apml, srctype, srcfunc, sigma, &
       f0, t0, xs, zs, facqui, dts)
  
  !# Read acquisition parameters
  is1 = nint(xs/h)+npml+1
  is2 = nint(zs/h)+npml+1

  nt = nint(tmax/dt)+1
  nts = nint(tmax/dts)+1

  write(*, *) nt, nts, nt/nts+1

  !# Acquisition
  call acqread(facqui, npml, h, nrec, recp)
  allocate(recx(nts, nrec), recz(nts, nrec))

  !# Extend dimensions
  n1e = n1+2*npml
  n2e = n2+2*npml
  
  !# Allocate model matrices
  allocate(tmod%vp(n1, n2), tmod%vs(n1, n2), tmod%ro(n1, n2))
  
  !# Read model files
  call modread(fvp, n1, n2, tmod%vp)
  call modread(fvs, n1, n2, tmod%vs)
  call modread(fro, n1, n2, tmod%ro)

  !# Extend models with pmls
  allocate(tmod%vpe(n1e, n2e), tmod%vse(n1e, n2e), tmod%roe(n1e, n2e))

  call modext(tmod%vp, n1, n2, npml, tmod%vpe)
  call modext(tmod%vs, n1, n2, npml, tmod%vse)
  call modext(tmod%ro, n1, n2, npml, tmod%roe)

  !# Buoyancy and Lame parameters
  allocate(tmod%bux(n1e, n2e), tmod%buz(n1e, n2e))
  allocate(tmod%mu0(n1e, n2e), tmod%mue(n1e, n2e))
  allocate(tmod%lb0(n1e, n2e), tmod%lbmu(n1e,n2e))
  
  do i2=1,n2e-1
     tmod%bux(:, i2) = (1./2.)*(1./tmod%roe(:, i2)+1./tmod%roe(:, i2+1))
  end do
  tmod%bux(:,n2e) = 1./tmod%roe(:,n2e)

  do i1=1,n1e-1
     tmod%buz(i1, :) = (1./2.)*(1./tmod%roe(i1, :)+1./tmod%roe(i1+1, :))
  end do
  tmod%buz(n1e,:) = 1./tmod%roe(n1e,:)

  tmod%mu0(:, :) = tmod%vse(:, :)*tmod%vse(:, :)*tmod%roe(:, :)
  
  do i2=1,n2e-1
     do i1=1, n2e-1
        tmod%mue(i1, i2) = (1./4.)*(tmod%mu0(i1,i2)+tmod%mu0(i1+1,i2)+tmod%mu0(i1,i2+1)+tmod%mu0(i1+1,i2+1))
     end do
  end do
  tmod%mue(:,n2e) = tmod%mue(:,n2e-1)
  tmod%mue(n1e,:) = tmod%mue(n1e-1,:)

  tmod%lb0(:, :) = tmod%vpe(:, :)*tmod%vpe(:, :)*tmod%roe(:, :)-2.*tmod%mu0(:, :)

  tmod%lbmu(:,:) = tmod%lb0(:, :)+2.*tmod%mu0(:,:)


  !# Stability condition
  write(*, *) 'Courant::', dt*maxval(tmod%vp)/h

  !# Deallocate model matrices
  deallocate(tmod%vp, tmod%vs, tmod%ro)
  deallocate(tmod%vpe, tmod%vse, tmod%roe)

  !# Source
  allocate(tsrc(nt), gsrc(n1e, n2e)) 
  call ricker(nt, dt, f0, t0, tsrc)
  call srcspread(n1e, n2e, npml, is1, is2, h, gsrc, sigma)

  allocate(spg(3, npml+1))

  spg(:, : ) = 1.

  write(*, * ) 'SPONGES'
  call sponges(npml, apml, h, spg)

  write(*, * ) 'EVOLUTION'
  call evolution(n1e, n2e, npml, h, dt, dts, nt, nts, nrec, srctype, tsrc, gsrc, spg, recx, recz, recp, tmod, isurf)

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

end program main
