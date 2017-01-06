module model

  implicit none
  
contains

  subroutine modread(fname, n1, n2, v)

    integer :: n1, n2
    real :: v(n1, n2)
    character(len=*) :: fname

    open(11, file=fname, access='direct', recl=n1*n2*4)
    read(11, rec=1) v
    close(11)

  end subroutine modread

  subroutine modext(v, n1, n2, nsp, ve)

    integer :: i, n1, n2, n1e, n2e, nsp
    real :: v(n1, n2), ve(n1+2*nsp, n2+2*nsp)

    ve(nsp+1:n1+nsp,nsp+1:n2+nsp) = v(1:n1,1:n2)

    do i=1,nsp
       ve(:,i) = ve(:,nsp+1)
       ve(:,n2+nsp+i) = ve(:,n2+nsp)
    end do

    do i=1,nsp
       ve(i, :) = ve(nsp+1, :)
       ve(n1+nsp+i, :) = ve(n1+nsp, :)
    end do

  end subroutine modext
  
end module model
