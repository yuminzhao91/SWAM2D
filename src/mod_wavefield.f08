module wavefield

  use types
  use deriv
  
  implicit none
  
contains

  subroutine field_ux (trun, tmod, tbnd, ux, txx, txz)

    implicit none

    type (typerun) :: trun
    type (typemod) :: tmod
    type (typebnd) :: tbnd
    
    real :: ux(tmod%n1e, tmod%n2e)
    real :: txx(tmod%n1e, tmod%n2e)
    real :: txz(tmod%n1e, tmod%n2e)
    
    real, allocatable    :: d1(:, :), d2(:, :)
    real, allocatable    :: uxx(:, :), uxz(:, :)
    
    !# Allocate local arrays
    allocate (d1(tmod%n1e, tmod%n2e), d2(tmod%n1e, tmod%n2e))
    allocate (uxx(tmod%n1e, tmod%n2e), uxz(tmod%n1e, tmod%n2e))

    !# Calculate derivatives
    call dxforward (txx, tmod%n1e, tmod%n2e, d2)
    call dzbackward (txz, tmod%n1e, tmod%n2e, tbnd%npml, d1, tbnd%isurf)

    !# Calculate wavefields
    uxx(:, :) = (((1./trun%dt-tbnd%pmlx1(:,:))*uxx(:, :) &
         +(1./tmod%h)*tmod%bux(:, :)*d2(:, :))/(1./trun%dt+tbnd%pmlx1(:, :)))
    uxz(:, :) = (((1./trun%dt-tbnd%pmlz0(:,:))*uxz(:, :) &
         +(1./tmod%h)*tmod%bux(:, :)*d1(:, :))/(1./trun%dt+tbnd%pmlz0(:, :)))

    ux(:, :) = uxx(:, :) + uxz(:, :)
    
    !# Deallocate local arrays
    deallocate (uxx, uxz)
    deallocate (d1, d2)
    
  end subroutine field_ux
  
end module wavefield
