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
  
  subroutine modbuo(roe, n1e, n2e, bux, buz)

    integer :: i1, i2, n1e, n2e
    real, dimension(n1e, n2e) :: roe, bux, buz

    do i2=1,n2e-1
       bux(:, i2) = (1./2.)*(1./roe(:, i2)+1./roe(:, i2+1))
    end do
    bux(:,n2e) = 1./roe(:,n2e)

    do i1=1,n1e-1
       buz(i1, :) = (1./2.)*(1./roe(i1, :)+1./roe(i1+1, :))
    end do
    buz(n1e,:) = 1./roe(n1e,:)

  end subroutine modbuo

  subroutine modlame(vpe, vse, roe, n1e, n2e, mu0, mue, lb0, lbmu)

    integer :: i1, i2, n1e, n2e
    real, dimension(n1e, n2e) :: vpe, vse, roe, mu0, mue, lb0, lbmu

    mu0(:, :) = vse(:, :)*vse(:, :)*roe(:, :)
  
    do i2=1,n2e-1
       do i1=1, n2e-1
          mue(i1, i2) = (1./4.)*(mu0(i1,i2)+mu0(i1+1,i2)+mu0(i1,i2+1)+mu0(i1+1,i2+1))
       end do
    end do
    mue(:,n2e) = mue(:,n2e-1)
    mue(n1e,:) = mue(n1e-1,:)

    lb0(:, :) = vpe(:, :)*vpe(:, :)*roe(:, :)-2.*mu0(:, :)
    
    lbmu(:,:) = lb0(:, :)+2.*mu0(:,:)

  end subroutine modlame

end module model
