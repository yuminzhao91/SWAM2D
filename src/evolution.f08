module marching

  use deriv

  implicit none

contains
  
  subroutine evolution(n1e, n2e, nsp, ux, uz, dt, dts, nt, tsrc, gsrc, spg)

    integer :: it, n1e, n2e, nsp, nt
    real :: dt, dts, tsrc(nt), gsrc(n1e, n2e), spg(3, nsp+1)
    real :: ux(n1e, n2e), uz(n1e, n2e)
    
    real, allocatable :: uxx(:, :), uxz(:, :), uzx(:, :), uzz(:, :)
    real, allocatable :: txx(:, :), txxx(:, :), txxz(:, :)
    real, allocatable :: tzz(:, :), tzzx(:, :), tzzz(:, :)
    real, allocatable :: txz(:, :), txzx(:, :), txzz(:, :)

    real, allocatable :: d1(:, :), d2(:, :)

    allocate(uxx(n1e, n2e), uxz(n1e, n2e))
    allocate(uzx(n1e, n2e), uzz(n1e, n2e))
    allocate(txx(n1e, n2e), txxx(n1e, n2e), txxz(n1e, n2e))
    allocate(tzz(n1e, n2e), tzzx(n1e, n2e), tzzz(n1e, n2e))
    allocate(txz(n1e, n2e), txzx(n1e, n2e), txzz(n1e, n2e))

    

    deallocate(uxx, uzz)
    deallocate(txx, txxx, txxz)
    deallocate(tzz, tzzx, tzzz)
    deallocate(txz, txzx, txzz)

end module marching
