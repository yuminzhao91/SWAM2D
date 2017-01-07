module param
  implicit none

contains

  subroutine parread(frun, tmax, dt, fvp, fvs, fro, n1, n2, h, &
       isurf, npml, apml, srctype, srcfunc, sigma, f0, t0, &
       xs, zs, facqui, dts)

    integer :: n1, n2, isurf, npml, srctype, srcfunc
    real :: tmax, dt, h, apml, sigma, f0, t0, xs, zs, dts
    character(len=*) :: frun, fvp, fvs, fro, facqui

    open(101, file='param.input', status='old')
    read(101, *) !#[run]
    read(101, *) frun
    read(101, *) tmax, dt
    read(101, *) !#[materials]
    read(101, *) fvp, fvs, fro
    read(101, *) n1, n2, h
    read(101, *) !#[boundaries]
    read(101, *) isurf, npml, apml
    read(101, *) ! #[source]
    read(101, *) srctype, srcfunc, sigma
    read(101, *) f0, t0
    read(101, *) xs, zs 
    read(101, *) !#[receiver]
    read(101, *) facqui
    read(101, *) dts
    close(101)

  end subroutine parread

end module param
