!------------------------------------------------------------------------------
! MODULE: source
!------------------------------------------------------------------------------
!> \brief subroutines related to source time function generation, source
!> positioning  and source spreading on grid
!------------------------------------------------------------------------------
!> \author Damien Pageot
!> \date 09 Jan 2017
!------------------------------------------------------------------------------
! Revision history
! 09 Jan 2017: Add doxygen documentation
!------------------------------------------------------------------------------

module source
  
  use types

  implicit none

contains

  subroutine srcaddexp(ve, n1e, n2e, dt, dx, src, g)
    !> \brief Add explosive source. 
    integer:: n1e, n2e
    real :: dt, dx, src
    real :: ve(n1e, n2e), g(n1e, n2e)

    ve(:, :) = ve(:, :)+(dt/dx)*g(:, :)*src
    
  end subroutine srcaddexp

  subroutine srcaddforce(ve, n1e, n2e, dt, dx, src, g)
    !> \brief Add vertical or horizontal force source.
    integer:: i1, i2, n1e, n2e
    real :: dt, dx, src
    real :: ve(n1e, n2e), g(n1e, n2e)

    ve(:, :) = ve(:, :)+(dt*dt/dx)*g(:, :)*src
    
  end subroutine srcaddforce
  
  subroutine srcindex(xs, zs, npml, h, is1, is2)
    
    integer :: npml, is1, is2
    real :: xs, zs, h

    is1 = nint(zs/h)+npml+1
    is2 = nint(xs/h)+npml+1

  end subroutine srcindex

  subroutine ricker(nt, dt, f0, t0, tsrc)
    
    integer :: it, nt
    real :: dt, f0, t0, tsrc(nt)
    real :: sigma, pi, t

    pi = 4.*atan(1.)

    do it=1,nt
       t = float(it-1)*dt-t0
       sigma = (pi*f0*(t))*(pi*f0*(t)) 
       tsrc(it) = -1.*(1.-2.*sigma)*exp(-1.*sigma)
    end do

    open(11, file='fricker.bin', access='direct', recl=nt*4)
    write(11, rec=1) tsrc
    close(11)

  end subroutine ricker

  subroutine srcspread(n1e, n2e, nsp, is1, is2, h, gsrc, sigma)
    
    integer :: i1, i2, is1, is2, n1e, n2e, nsp
    real :: gsrc(n1e, n2e), betasum
    real :: x, z, xs, zs, p1, p2, sigma, h

    gsrc(:, :) = 0.

    betasum = 0.

    if( sigma .lt. 0.)then
       gsrc(is1, is2) = 1.
    else
    xs = float(is2-1)*h
    zs = float(is1-1)*h
    do i2=nsp+1,n2e-nsp
       x = float(i2-1)*h 
       do i1=nsp+1,n1e-nsp
          z = float(i1-1)*h
          p2 = (x-xs)*(x-xs)/(sigma*sigma)
          p1 = (z-zs)*(z-zs)/(sigma*sigma)
          gsrc(i1, i2) = exp(-1.*p1-p2)
          betasum = betasum+exp(-1.*p1-p2)
       end do
    end do

    do i2=1,n2e
       do i1=1,n1e
          gsrc(i1, i2) = gsrc(i1, i2)/betasum
       end do
    end do
    endif

  end subroutine srcspread

end module source
