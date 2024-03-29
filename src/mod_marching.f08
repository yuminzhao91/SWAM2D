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
  
  subroutine evolution(n1e, n2e, nsp, h,  dt, nt, nts, ntsnap, nrec, srctype, &
       tsrc, gsrc, recx, recz, recp, recpos, tmod, isurf, &
       pmlx0, pmlx1, pmlz0, pmlz1, isnap)

    type(typemod) :: tmod

    integer :: i1, i2, it, n1e, n2e, nsp, nt, nts, nrec, its, ets, itsnap, ntsnap, isnap
    integer :: etsnap, ix, iz, irec, itt, srctype, isurf
    real :: start, finish
    real :: dt, tsrc(nt), gsrc(n1e, n2e), h, full
    real :: dth

    real :: recx(nts, nrec), recz(nts, nrec), recp(nts, nrec)
    integer :: recpos(nrec, 2)

    real, allocatable :: ux(:, :), uz(:, :), txx(:, :), tzz(:, :), txz(:, :)
    real, allocatable :: uxx(:, :), uxz(:, :)
    real, allocatable :: uzx(:, :), uzz(:, :)
    real, allocatable :: txxx(:, :), txxz(:, :)
    real, allocatable :: tzzx(:, :), tzzz(:, :)
    real, allocatable :: txzx(:, :), txzz(:, :)

    real, allocatable :: press(:, :)

    real, allocatable :: d1(:, :), d2(:, :)

    character(len=80) :: snapfile


    ! >> ADDING PML TEMP. IN MARCHING MOD
    real, dimension(n1e, n2e) :: pmlx0, pmlx1, pmlz0, pmlz1


    !open(901, file='test_pmlz0.bin', access='direct', recl=n1e*n2e*4)
    !write(901, rec=1) pmlz0
    !close(901)

    !open(902, file='test_pmlz1.bin', access='direct', recl=n1e*n2e*4)
    !write(902, rec=1) pmlz1
    !close(902)

    ! >> END PML
   
    dth = dt/h
    
    write(*, *) 'EVOLUTION ALLOCATION'
    allocate(d1(n1e, n2e), d2(n1e, n2e))

    ! >> Allocate velocity fields
    allocate(ux(n1e, n2e), uz(n1e, n2e))

    ! >> Allocate stress fields
    allocate(txx(n1e, n2e), tzz(n1e, n2e), txz(n1e, n2e))

    ! >> Allocate splited velocity fields
    allocate(uxx(n1e, n2e), uxz(n1e, n2e))
    allocate(uzx(n1e, n2e), uzz(n1e, n2e))

    ! >> Allocate splited stress fields 
    allocate(txxx(n1e, n2e), txxz(n1e, n2e))
    allocate(tzzx(n1e, n2e), tzzz(n1e, n2e))
    allocate(txzx(n1e, n2e), txzz(n1e, n2e))

    ! >> Allocate pressure field
    allocate(press(n1e, n2e))

    ! >> Initialize velocity fields (including splitted)
    ux(:, :) = 0.
    uz(:, :) = 0.
    uxx(:, :) = 0.
    uxz(:, :) = 0.
    uzx(:, :) = 0.
    uzz(:, : ) = 0.

    ! >> Initialize stress fields (including splitted)
    txx(:, :) = 0.
    tzz(:, :) = 0.
    txz(:, :) = 0.
    txxx(:, :) = 0.
    txxz(:, :) = 0.
    tzzx(:, :) = 0.
    tzzz(:, :) = 0.
    txzx(:, :) = 0.
    txzz(:, :) = 0.

    ! >> Initialize marching and sampling parameters
    full = 0.
    its = 1
    itsnap = 1
    itt = 1
    ets = (nt-1)/(nts-1)
    etsnap = (nt-1)/(ntsnap-1)

    ! >> Start marching
    do it=1,nt
       call cpu_time(start)

       !# UX
       !# dux/dut = (1/rho)[dtxx/dx+dtxz/dz]
       
       call dxforward(txx, n1e, n2e, d2)
       call dzbackward(txz, n1e, n2e, d1)
       
       uxx(:, :) = (((1./dt-pmlx1(:,:))*uxx(:, :)+(1./h)*tmod%bux(:, :)*d2(:, :))/(1./dt+pmlx1(:, :)))
       uxz(:, :) = (((1./dt-pmlz0(:,:))*uxz(:, :)+(1./h)*tmod%bux(:, :)*d1(:, :))/(1./dt+pmlz0(:, :)))
       ux(:, :) = uxx(:, :)+uxz(:, :)
       if(isurf == 1)then
          ux(1:nsp,:) = 0.
       end if
       
       !# UZ
       call dxbackward(txz, n1e, n2e, d2)
       call dzforward(tzz, n1e, n2e, d1)

       if(srctype== 2)then
          uzx(:, :) = (((1./dt-pmlx0(:,:))*uzx(:, :)+(1./h)*tmod%buz(:, :)*d2(:, :))/(1./dt+pmlx0(:, :))) &
               +tmod%buz*(tsrc(it)*gsrc(:, :)*dt/(h*h))
       else
          uzx(:, :) = (((1./dt-pmlx0(:,:))*uzx(:, :)+(1./h)*tmod%buz(:, :)*d2(:, :))/(1./dt+pmlx0(:, :)))
       end if
       uzz(:, :) = (((1./dt-pmlz1(:,:))*uzz(:, :)+(1./h)*tmod%buz(:, :)*d1(:, :))/(1./dt+pmlz1(:, :)))

       uz(:, :) = uzx(:, :)+uzz(:, :)
       if(isurf == 1)then
          uz(1:nsp,:) = 0.
       end if
       
       call dirichlet( n1e, n2e, uxx, uxz, uzx, uzz)
       ! implement Dirichlet boundary conditions on the four edges of the grid
       
       !# PRESSURE
       do i2=1,n2e-1
          do i1=2,n1e-1
             press(i1, i2) = (-tmod%lbmu(i1, i2)/h)*(ux(i1,i2)-ux(i1,i2-1)+uz(i1,i2)-uz(i1-1,i2))
          enddo
       enddo

       !# TXX -- TZZ
       call dxbackward(ux, n1e, n2e, d2)
       call dzbackward(uz, n1e, n2e, d1)

       if(srctype == 1)then
          txxx(:, :) = (((1./dt-pmlx0(:,:))*txxx(:, :)+(1./h)*tmod%lbmu(:, :)*d2(:, :)&
               )/(1./dt+pmlx0(:, :)))+(tsrc(it)*gsrc(:,:))/(h*h)*dt
       else
          txxx(:, :) = (((1./dt-pmlx0(:,:))*txxx(:, :)+(1./h)*tmod%lbmu(:, :)*d2(:, :))/(1./dt+pmlx0(:, :)))
       end if
       txxz(:, :) = (((1./dt-pmlz0(:,:))*txxz(:, :)+(1./h)*tmod%lb0(:, :)*d1(:, :))/(1./dt+pmlz0(:, :)))

       txx(:, :) = txxx(:, :) + txxz(:, :)
       if(isurf == 1)then
          txx(nsp+1, :) = txxx(nsp+1, :)
       end if
       
       if(srctype == 0 .or. srctype == 1)then
          tzzx(:, :) = (((1./dt-pmlx0(:,:))*tzzx(:, :)+(1./h)*tmod%lb0(:, :)*d2(:, :)&
               )/(1./dt+pmlx0(:, :)))+(tsrc(it)*gsrc(:,:))/(h*h)*dt
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
             write (snapfile, "(A9,I1)") "snapz0000", it
             open(31, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(31, rec=1) uz
             close(31)
             write (snapfile, "(A9,I1)") "snapx0000", it
             open(32, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(32, rec=1) ux
             close(32)
             write (snapfile, "(A9,I1)") "snapp0000", it
             open(33, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(33, rec=1) press
             close(33)
          else if(it >= 10 .and. it < 100)then
             write (snapfile, "(A8,I2)") "snapz000", it
             open(31, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(31, rec=1) uz
             close(31)
             write (snapfile, "(A8,I2)") "snapx000", it
             open(32, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(32, rec=1) ux
             close(32)
             write (snapfile, "(A8,I2)") "snapp000", it
             open(33, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(33, rec=1) press
             close(33)
          else if(it >= 100 .and. it < 1000)then
             write (snapfile, "(A7,I3)") "snapz00", it
             open(31, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(31, rec=1) uz
             close(31)
             write (snapfile, "(A7,I3)") "snapx00", it
             open(32, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(32, rec=1) ux
             close(32)
             write (snapfile, "(A7,I3)") "snapp00", it
             open(33, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(33, rec=1) press
             close(33)
          else if(it >= 1000 .and. it < 10000)then
             write (snapfile, "(A6,I4)") "snapz0", it
             open(31, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(31, rec=1) uz
             close(31)
             write (snapfile, "(A6,I4)") "snapx0", it
             open(32, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(32, rec=1) ux
             close(32)
             write (snapfile, "(A6,I4)") "snapp0", it
             open(33, file=snapfile, access='direct', recl=n1e*n2e*4)
             write(33, rec=1) press
             close(33)
          end if
       else
          itsnap = itsnap+1
       endif

       if(its == ets .or. it == 1)then
          its = 1
          do irec=1,nrec
             ix = recpos(irec, 1)
             iz = recpos(irec, 2) 
             recx(itt, irec) = (ux(iz, ix)+ux(iz,ix-1))/2.
             recz(itt, irec) = (uz(iz, ix)+uz(iz-1,ix))/2.
             recp(itt, irec) = press(iz, ix)
          end do
          itt = itt + 1
       else
          its = its + 1
       end if

    end do

    deallocate(d1, d2)

    ! >> Free velocity fields 
    deallocate(ux, uz)

    ! >> Free splitted velocity fields
    deallocate(uxx, uxz, uzx, uzz)

    ! >> Free stress fields
    deallocate(txx, tzz, txz)

    ! >> Free splitted stress fields
    deallocate(txxx, txxz, tzzx, tzzz, txzx, txzz)

    ! >> Free pressure field
    deallocate(press)

  end subroutine evolution

end module marching
