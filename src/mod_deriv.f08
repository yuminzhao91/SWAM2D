!------------------------------------------------------------------------------
! MODULE: deriv
!------------------------------------------------------------------------------
!> \brief subroutines to calculate 4th order forward and backward derivative
!> operators for the x-direction and z-direction:
!> \f$ D^{+}_{x}\f$, \f$ D^{-}_{x}\f$, \f$ D^{+}_{z}\f$ and
!> \f$ D^{-}_{z}\f$ (see \cite levander1988fourth).
!------------------------------------------------------------------------------
!> \author Damien Pageot
!> \date 09 Jan 2017
!------------------------------------------------------------------------------
! Revision history
! 09 Jan 2017: Add doxygen documentation
!------------------------------------------------------------------------------

module deriv

  implicit none

  real, parameter :: c1=(9./8.), c2=(-1./24.)
  
contains

  subroutine dxforward(f, n1, n2, d)
    !> \brief 4th order forward derivative in the x-direction.\n
    !> \f$ D^{+}_{x} = c1*[f(i,j+1)-f(i,j)]+c2*[f(i,j+2)-f(i,j-1)]\f$.
    !> \param[out] d  derivative
    !> \param[in]  f  array of size n1*n2 to derive
    !> \param[in]  n1 The number of grid points in the first direction (z)
    !> \param[in]  n2 The number of grid points in the second direction (x)
    integer :: i2
    integer :: n1, n2
    real :: f(n1, n2)

    real :: d(n1, n2)

    d(:, : ) = 0.

    !! >> 4th order derivative
    do i2=2,n2-2
       d(:,i2) = c1*(f(:,i2+1)-f(:,i2))+c2*(f(:,i2+2)-f(:,i2-1))
    end do

    !! >> 2nd order derivative
    d(:,1) = f(:,2)-f(:,1)
    d(:,n2-1) = f(:,n2)-f(:,n2-1)

  end subroutine dxforward

  subroutine dxbackward(f, n1, n2, d)
    !> @brief 4th order backward derivative in the x-direction.\n
    !> \f$ D^{-}_{x} = c1*[f(i,j)-f(i,j-1)]+c2*[f(i,j+1)-f(i,j-2)]\f$.
    !> @param[out] d  derivative
    !> @param[in]  f  array of size n1*n2 to derive
    !> @param[in]  n1 The number of grid points in the first direction (z)
    !> @param[in]  n2 The number of grid points in the second direction (x)
    integer :: i2
    integer :: n2, n1
    real :: f(n1, n2)

    real :: d(n1, n2)

    d(:, : ) = 0.

    !! >> 4th order derivative
    do i2=3,n2-1
       d(:,i2) = c1*(f(:,i2)-f(:,i2-1))+c2*(f(:,i2+1)-f(:,i2-2))
    end do
    
    !! >> 2nd order derivative
    d(:,2) = f(:,2)-f(:,1)
    d(:,n2) = f(:,n2)-f(:,n2-1)

  end subroutine dxbackward

  subroutine dzforward(f, n1, n2, d)
    !> @brief 4th order forward derivative in the z-direction.\n
    !> \f$ D^{+}_{z} = c1*[f(i+1,j)-f(i,j)]+c2*[f(i+2,j)-f(i-1,j)]\f$.
    !> @param[out] d  derivative
    !> @param[in]  f  array of size n1*n2 to derive
    !> @param[in]  n1 The number of grid points in the first direction (z)
    !> @param[in]  n2 The number of grid points in the second direction (x)
    integer :: i1
    integer :: n2, n1
    real :: f(n1, n2)

    real :: d(n1, n2)

    d(:, : ) = 0.

    !! >> 4th order derivative
    do i1=2,n1-2
       d(i1,:) = c1*(f(i1+1,:)-f(i1,:))+c2*(f(i1+2,:)-f(i1-1,:))
    end do
    
    !! >> 2nd order derivative
    d(1,:) = f(2,:)-f(1,:)
    d(n1-1,:) = f(n1,:)-f(n1-1,:)

  end subroutine dzforward

  subroutine dzbackward(f, n1, n2, d)
    !> @brief 4th order backward derivative in the z-direction.\n
    !> \f$ D^{-}_{z} = c1*[f(i,j)-f(i-1,j)]+c2*[f(i+1,j)-f(i-2,j)]\f$.
    !> @param[out] d  derivative
    !> @param[in]  f  array of size n1*n2 to derive
    !> @param[in]  n1 The number of grid points in the first direction (z)
    !> @param[in]  n2 The number of grid points in the second direction (x)
    integer :: i1
    integer :: n2, n1
    real :: f(n1, n2)

    real :: d(n1, n2)

    d(:, : ) = 0.

    !! >> 4th order derivative
    do i1=3,n1-1
       d(i1,:) = c1*(f(i1,:)-f(i1-1,:))+c2*(f(i1+1,:)-f(i1-2,:))
    end do

    !! >> 2nd order derivative
    d(2,:) = f(2,:)-f(1,:)
    d(n1,:) = f(n1,:)-f(n1-1,:)

  end subroutine dzbackward
  
end module deriv
