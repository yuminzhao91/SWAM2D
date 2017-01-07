module deriv

  implicit none

  real, parameter :: c1=(9./8.), c2=(-1./24.)
  
contains

  subroutine dxforward(f, n1, n2, d)
    ! >> dxforward
    !    D(i) = f(i+1)-f(i)
    integer :: i2
    integer :: n1, n2
    real :: f(n1, n2)

    real :: d(n1, n2)

    d(:, : ) = 0.

    ! >> 4th order derivative
    do i2=2,n2-2
       d(:,i2) = c1*(f(:,i2+1)-f(:,i2))+c2*(f(:,i2+2)-f(:,i2-1))
    end do

    ! >> 2nd order derivative
    d(:,1) = f(:,2)-f(:,1)
    d(:,n2-1) = f(:,n2)-f(:,n2-1)

  end subroutine dxforward

  subroutine dxbackward(f, n1, n2, d)
    ! >> dxbackward
    !    D(i) = f(i)-f(i-1)
    integer :: i2
    integer :: n2, n1
    real :: f(n1, n2)

    real :: d(n1, n2)

    d(:, : ) = 0.

    ! >> 4th order derivative
    do i2=3,n2-1
       d(:,i2) = c1*(f(:,i2)-f(:,i2-1))+c2*(f(:,i2+1)-f(:,i2-2))
    end do
    
    ! >> 2nd order derivative
    d(:,2) = f(:,2)-f(:,1)
    d(:,n2) = f(:,n2)-f(:,n2-1)

  end subroutine dxbackward

  subroutine dzforward(f, n1, n2, d)
    ! >> dzforward
    !    D(i) = f(i+1)-f(i)
    integer :: i1
    integer :: n2, n1
    real :: f(n1, n2)

    real :: d(n1, n2)

    d(:, : ) = 0.

    ! >> 4th order derivative
    do i1=2,n1-2
       d(i1,:) = c1*(f(i1+1,:)-f(i1,:))+c2*(f(i1+2,:)-f(i1-1,:))
    end do
    
    ! >> 2nd order derivative
    d(1,:) = f(2,:)-f(1,:)
    d(n1-1,:) = f(n1,:)-f(n1-1,:)

  end subroutine dzforward

  subroutine dzbackward(f, n1, n2, d)
    ! >> dzbackward
    !    D(i) = f(i)-f(i-1)
    integer :: i1
    integer :: n2, n1
    real :: f(n1, n2)

    real :: d(n1, n2)

    d(:, : ) = 0.

    ! >> 4th order derivative
    do i1=3,n1-1
       d(i1,:) = c1*(f(i1,:)-f(i1-1,:))+c2*(f(i1+1,:)-f(i1-2,:))
    end do

    ! >> 2nd order derivative
    d(2,:) = f(2,:)-f(1,:)
    d(n1,:) = f(n1,:)-f(n1-1,:)

  end subroutine dzbackward
  
end module deriv
