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
  
  subroutine evolution(n1e, n2e, nsp, h,  dt, dts, nt, nts, nrec, srctype, tsrc, gsrc, spg, recx, recz, recp, tmod, isurf)

    type(typemod) :: tmod

    integer :: it, n1e, n2e, nsp, nt, i1, i2, nts, nrec, its, ets, ix, iz, irec, j, itt, srctype, isurf
    real :: start, finish
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

    do it=1,nt
       call cpu_time(start)

       !# UX
       
       call dxbackward(txx, n1e, n2e, d2)
       call dzbackward(txz, n1e, n2e, d1)

       call addpmlz(n1e, n2e, nsp, spg(1,:), spg(1,:), uxz, isurf)
       call addpmlx(n1e, n2e, nsp, spg(1,:), spg(1,:), uxx)

       uxx(:, :) = uxx(:, :)+dth*tmod%bux(:, :)*d2(:, :)
       uxz(:, :) = uxz(:, :)+dth*tmod%bux(:, :)*d1(:, :)
       ux(:, :) = uxx(:, :)+uxz(:, :)
       if(isurf == 1)then
          ux(1:nsp,:) = 0.
       end if
       
       !# UZ
       call dxforward(txz, n1e, n2e, d2)
       call dzforward(tzz, n1e, n2e, d1)

       call addpmlz(n1e, n2e, nsp, spg(2,:), spg(3,:), uzz, isurf)
       call addpmlx(n1e, n2e, nsp, spg(2,:), spg(3,:), uzx)

       if(srctype== 2)then
          uzx(:, :) = uzx(:, :)+dth*tmod%buz(:, :)*d2(:, :)+(tsrc(it)*gsrc(:,:))*dt*dt/(h*h)
       else
          uzx(:, :) = uzx(:, :)+dth*tmod%buz(:, :)*d2(:, :)
       end if
       uzz(:, :) = uzz(:, :)+dth*tmod%buz(:, :)*d1(:, :)

       uz(:, :) = uzx(:, :)+uzz(:, :)
       if(isurf == 1)then
          uz(1:nsp,:) = 0.
       end if
       
       !# TXX -- TZZ
       call dxforward(ux, n1e, n2e, d2)
       call dzbackward(uz, n1e, n2e, d1)

       call addpmlx(n1e, n2e, nsp, spg(2,:), spg(3,:), txxx)
       call addpmlz(n1e, n2e, nsp, spg(1,:), spg(1,:), txxz, isurf)

       if(srctype == 1)then
          txxx(:, :) = txxx(:, :)+dth*tmod%lbmu(:, :)*d2(:, :)+(tsrc(it)*gsrc(:,:))*dt/(h*h)
       else
          txxx(:, :) = txxx(:, :)+dth*tmod%lbmu(:, :)*d2(:, :)
       end if
       txxz(:, :) = txxz(:, :)+dth*tmod%lb0(:, :)*d1(:, :)

       txx(:, :) = txxx(:, :) + txxz(:, :)
       if(isurf == 1)then
          txx(nsp+1, :) = txxx(nsp+1, :)
       end if
       
       call addpmlx(n1e, n2e, nsp, spg(2,:), spg(3,:), tzzx)
       call addpmlz(n1e, n2e, nsp, spg(1,:), spg(1,:), tzzz, isurf)

       if(srctype == 0 .or. srctype == 1)then
          tzzx(:, :) = tzzx(:, :)+dth*tmod%lb0(:, :)*d2(:, :)+(tsrc(it)*gsrc(:,:))*dt/(h*h)
       else
          tzzx(:, :) = tzzx(:, :)+dth*tmod%lb0(:, :)*d2(:, :)
       end if
       tzzz(:, :) = tzzz(:, :)+dth*tmod%lbmu(:, :)*d1(:, :)
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
       call dxbackward(uz, n1e, n2e, d2)
       call dzforward(ux, n1e, n2e, d1)

       call addpmlx(n1e, n2e, nsp, spg(1,:), spg(1,:), txzx)
       call addpmlz(n1e, n2e, nsp, spg(2,:), spg(3,:), txzz, isurf)

       txzx(:, :) = txzx(:, :)+dth*tmod%mue(:, :)*d2(:, :)
       txzz(:, :) = txzz(:, :)+dth*tmod%mue(:, :)*d1(:, :)
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
