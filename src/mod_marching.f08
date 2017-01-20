!------------------------------------------------------------------------------
! MODULE: marching
!------------------------------------------------------------------------------
!> @brief
!> ...
!!
!> @details
!> \latexonly
!> \begin{small}
!> \begin{eqnarray}
!> D^{+}_{t}u_{t}(m,n,l-1/2) =\frac{1}{\rho(m,n)}[
!> D^{-}_{x}\tau_{xx}(m+1/2,n,l) +
!> D^{-}_{z}\tau_{xz}(m,n+1/2,l)]
!> \\
!> D^{+}_{t}w_{t}(m+1/2,n+1/2,l-1/2)=\frac{1}{\rho(m+1/2,n+1/2)}[
!> D^{+}_{x}\tau_{xz}(m,n+1/2,l) +
!> D^{+}_{z}\tau_{zz}(m+1/2,n,l)]
!> \\
!> D^{+}_{t}\tau_{xx}(m+1/2,n,l) = 
!> [\lambda(m+1/2,n)+2\mu(m+1/2,n)]D^{+}_{x}u_{t}(m,n,l+1/2) \nonumber\\
!> + \lambda(m+1/2,n)D^{-}_{z}w_{t}(m+1/2,n+1/2,l+1/2)]&
!> \\
!> D^{+}_{t}\tau_{xz}(m,n+1/2,l) = 
!> \mu(m,n+1/2)[D^{+}_{z}u_{t}(m,n,l+1/2)+
!> D^{-}_{x}w_{t}(m+1/2,n+1/2,l+1/2)]
!> \\
!> D^{+}_{t}\tau_{zz} = 
!> [\lambda(m+1/2,n)+2\mu(m+1/2,n)]D^{-}_{z}w_{t}(m+1/2,n+1/2,l+1/2)\nonumber \\
!> + \lambda(m+1/2,n)D^{+}_{x}u_{t}(m,n,l+1/2)&
!> \end{eqnarray}
!> \end{small}
!> \endlatexonly
!!
!> \cite levander1988fourth
!!
!> @author Damien Pageot
!> @date 09 Jan 2017
!------------------------------------------------------------------------------
! Revision history

module marching
  
  use types

  use deriv
  use boundaries

  implicit none

contains
  
  subroutine evolution(n1e, n2e, nsp, h,  dt, dts, dtsnap, nt, nts, ntsnap, nrec, srctype, &
       tsrc, gsrc, spg, recx, recz, recp, tmod, isurf, &
       pmlx0, pmlx1, pmlz0, pmlz1, isnap)

    type(typemod) :: tmod

    integer :: it, n1e, n2e, nsp, nt, i1, i2, nts, nrec, its, ets, itsnap, ntsnap, isnap
    integer :: etsnap, ix, iz, irec, j, itt, srctype, isurf
    real :: start, finish, dtsnap
    real :: dt, dts, tsrc(nt), gsrc(n1e, n2e), spg(3, nsp+1), h, full
    real :: dth
    
    real :: recx(nts, nrec), recz(nts, nrec)
    integer :: recp(nrec,2)

    real, allocatable :: ux(:, :), uz(:, :), txx(:, :), tzz(:, :), txz(:, :)
    real, allocatable :: uxx(:, :), uxz(:, :)
    real, allocatable :: uzx(:, :), uzz(:, :)
    real, allocatable :: txxx(:, :), txxz(:, :)
    real, allocatable :: tzzx(:, :), tzzz(:, :)
    real, allocatable :: txzx(:, :), txzz(:, :)

    real, allocatable :: d1(:, :), d2(:, :)

    character(len=80) :: snapfile

    ! >> ADDING PML TEMP. IN MARCHING MOD
    real, dimension(n1e, n2e) :: pmlx0, pmlx1, pmlz0, pmlz1


    open(901, file='test_pmlz0.bin', access='direct', recl=n1e*n2e*4)
    write(901, rec=1) pmlz0
    close(901)

    open(902, file='test_pmlz1.bin', access='direct', recl=n1e*n2e*4)
    write(902, rec=1) pmlz1
    close(902)

    ! >> END PML


    dth = dt/h
    
    write(*, *) 'EVOLUTION ALLOCATION'
    allocate(d1(n1e, n2e), d2(n1e, n2e))
    allocate(ux(n1e, n2e), uz(n1e, n2e), txx(n1e, n2e), tzz(n1e, n2e), txz(n1e, n2e))
    allocate(uxx(n1e, n2e), uxz(n1e, n2e))
    allocate(uzx(n1e, n2e), uzz(n1e, n2e))
    allocate(txxx(n1e, n2e), txxz(n1e, n2e))
    allocate(tzzx(n1e, n2e), tzzz(n1e, n2e))
    allocate(txzx(n1e, n2e), txzz(n1e, n2e))

    uxx(:, :) = 0.
    uxz(:, :) = 0.
    txxx(:, :) = 0.
    txxz(:, :) = 0.
    tzzx(:, :) = 0.
    tzzz(:, :) = 0.
    txzx(:, :) = 0.
    txzz(:, :) = 0.

    full = 0.
    its = 1
    itt = 1
    ets = (nt-1)/(nts-1)
    etsnap = (nt-1)/(ntsnap-1)

    do it=1,nt
       call cpu_time(start)

       !# UX
       
       call dxforward(txx, n1e, n2e, d2)
       call dzbackward(txz, n1e, n2e, d1)

       !call addpmlx(n1e, n2e, nsp, spg(3,:), spg(3,:), uxx)
       !call addpmlz(n1e, n2e, nsp, spg(1,:), spg(1,:), uxz, isurf)

       uxx(:, :) = (((1./dt-pmlx1(:,:))*uxx(:, :)+(1./h)*tmod%bux(:, :)*d2(:, :))/(1./dt+pmlx1(:, :)))
       uxz(:, :) = (((1./dt-pmlz0(:,:))*uxz(:, :)+(1./h)*tmod%bux(:, :)*d1(:, :))/(1./dt+pmlz0(:, :))) !uxz(:, :)+dth*tmod%bux(:, :)*d1(:, :)
       ux(:, :) = uxx(:, :)+uxz(:, :)
       if(isurf == 1)then
          ux(1:nsp,:) = 0.
       end if
       
       !# UZ
       call dxbackward(txz, n1e, n2e, d2)
       call dzforward(tzz, n1e, n2e, d1)

       !call addpmlx(n1e, n2e, nsp, spg(1,:), spg(1,:), uzx)
       !call addpmlz(n1e, n2e, nsp, spg(3,:), spg(3,:), uzz, isurf)

       if(srctype== 2)then
          uzx(:, :) = (((1./dt-pmlx0(:,:))*uzx(:, :)+(1./h)*tmod%buz(:, :)*d2(:, :) &
               +(tsrc(it)*gsrc(:,:))*dt/(h*h))/(1./dt+pmlx0(:, :)))
       else
          uzx(:, :) = (((1./dt-pmlx0(:,:))*uzx(:, :)+(1./h)*tmod%buz(:, :)*d2(:, :))/(1./dt+pmlx0(:, :)))
       end if
       uzz(:, :) = (((1./dt-pmlz1(:,:))*uzz(:, :)+(1./h)*tmod%buz(:, :)*d1(:, :))/(1./dt+pmlz1(:, :)))

       uz(:, :) = uzx(:, :)+uzz(:, :)
       if(isurf == 1)then
          uz(1:nsp,:) = 0.
       end if
       
       ! implement Dirichlet boundary conditions on the four edges of the grid
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

       !# TXX -- TZZ
       call dxbackward(ux, n1e, n2e, d2)
       call dzbackward(uz, n1e, n2e, d1)

       !call addpmlx(n1e, n2e, nsp, spg(1,:), spg(1,:), txxx)
       !call addpmlz(n1e, n2e, nsp, spg(1,:), spg(1,:), txxz, isurf)

       if(srctype == 1)then
          txxx(:, :) = (((1./dt-pmlx0(:,:))*txxx(:, :)+(1./h)*tmod%lbmu(:, :)*d2(:, :)&
               +(tsrc(it)*gsrc(:,:))/(h*h))/(1./dt+pmlx0(:, :)))
       else
          txxx(:, :) = (((1./dt-pmlx0(:,:))*txxx(:, :)+(1./h)*tmod%lbmu(:, :)*d2(:, :))/(1./dt+pmlx0(:, :)))
       end if
       txxz(:, :) = (((1./dt-pmlz0(:,:))*txxz(:, :)+(1./h)*tmod%lb0(:, :)*d1(:, :))/(1./dt+pmlz0(:, :)))

       txx(:, :) = txxx(:, :) + txxz(:, :)
       if(isurf == 1)then
          txx(nsp+1, :) = txxx(nsp+1, :)
       end if
       
       !call addpmlx(n1e, n2e, nsp, spg(1,:), spg(1,:), tzzx)
       !call addpmlz(n1e, n2e, nsp, spg(1,:), spg(1,:), tzzz, isurf)

       if(srctype == 0 .or. srctype == 1)then
          tzzx(:, :) = (((1./dt-pmlx0(:,:))*tzzx(:, :)+(1./h)*tmod%lb0(:, :)*d2(:, :)&
               +(tsrc(it)*gsrc(:,:))/(h*h))/(1./dt+pmlx0(:, :)))
       else
          tzzx(:, :) = (((1./dt-pmlx0(:,:))*tzzx(:, :)+(1./h)*tmod%lb0(:, :)*d2(:, :))/(1./dt+pmlx0(:, :)))
       end if
       tzzz(:, :) = (((1./dt-pmlz0(:,:))*tzzz(:, :)+(1./h)*tmod%lbmu(:, :)*d1(:, :))/(1./dt+pmlz0(:, :)))
       if(isurf == 1)then
          tzzx(nsp+1,:) = 0.
          tzzz(nsp+1,:) = 0.
          tzzx(nsp,:) = -tzzx(nsp+2,:)
          tzzz(nsp,:) = -tzzz(nsp+2,:)
          tzzx(nsp-1,:) = -tzzx(nsp+3,:)
          tzzz(nsp-1,:) = -tzzz(nsp+3,:)
       end if
       tzz(:, :) = tzzx(:, :) + tzzz(:, :)

       !# TXZ
       call dxforward(uz, n1e, n2e, d2)
       call dzforward(ux, n1e, n2e, d1)

       !call addpmlx(n1e, n2e, nsp, spg(1,:), spg(1,:), txzx)
       !call addpmlz(n1e, n2e, nsp, spg(3,:), spg(3,:), txzz, isurf)

       txzx(:, :) = (((1./dt-pmlx1(:,:))*txzx(:, :)+(1./h)*tmod%mue(:, :)*d2(:, :))/(1./dt+pmlx1(:, :)))
       txzz(:, :) = (((1./dt-pmlz1(:,:))*txzz(:, :)+(1./h)*tmod%mue(:, :)*d1(:, :))/(1./dt+pmlz1(:, :)))
       if(isurf == 1)then
          txzx(nsp+1,:) = 0.
          txzz(nsp+1,:) = 0.
          txzx(nsp,:) = -txzx(nsp+2,:)
          txzz(nsp,:) = -txzz(nsp+2,:)
          txzx(nsp-1,:) = -txzx(nsp+3,:)
          txzz(nsp-1,:) = -txzz(nsp+3,:)
       end if
       txz(:, :) = txzx(:, :)+txzz(:, :)

       call cpu_time(finish)
       full = full+(finish-start)
       write(*, * ) it, nt, finish-start, full, sqrt(maxval(ux)**2+maxval(uz)**2)

       if((itsnap == etsnap .or. it == 1) .and. isnap == 1)then
          itsnap = 1
          if(it < 10)then
             write (snapfile, "(A8,I1)") "snap0000", it
             open(31, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(31, rec=1) uz
             close(31)
          else if(it >= 10 .and. it < 100)then
             write (snapfile, "(A7,I2)") "snap000", it
             open(31, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(31, rec=1) uz
             close(31)
          else if(it >= 100 .and. it < 1000)then
             write (snapfile, "(A6,I3)") "snap00", it
             open(31, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(31, rec=1) uz
             close(31)
          else if(it >= 1000 .and. it < 10000)then
             write (snapfile, "(A5,I4)") "snap0", it
             open(31, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(31, rec=1) uz
             close(31)
          end if
       else
          itsnap = itsnap+1
       endif

       if(its == ets .or. it == 1)then
          its = 1
          do irec=1,nrec
             ix = recp(irec, 1)
             iz = recp(irec, 2) 
             recx(itt, irec) = ux(iz, ix)
             recz(itt, irec) = uz(iz, ix)
          end do
          itt = itt + 1
       else
          its = its + 1
       end if

    end do

    deallocate(d1, d2)
    deallocate(ux, uz, txx, tzz, txz)
    deallocate(uxx, uzz)
    deallocate(txxx, txxz)
    deallocate(tzzx, tzzz)
    deallocate(txzx, txzz)

  end subroutine evolution

end module marching
