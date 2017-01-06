module acquisition

  implicit none

contains

  subroutine acqread(fname, nsp, h, nrec, recp)
    
    integer              :: irec, nrec, nsp, stat
    character(len=*)     :: fname
    real :: xr, zr, h
    integer,    allocatable :: recp(:, :)

    nrec = 0
    open(11, file=fname, status='old')
    do while(.true.)
       read(11, *, iostat=stat)
       !# if EOF: exit while loop
       if(stat /= 0) exit
       !# else: increment receiver number
       nrec = nrec+1
    end do

    rewind(11)

    allocate(recp(nrec, 2))
    
    do irec=1,nrec
       read(11, *) xr, zr
       recp(irec, 1) = nint(xr/h)+nsp+1
       recp(irec, 2) = nint(zr/h)+nsp+1
    end do

    close(11)

  end subroutine acqread
    
  subroutine acqfree(recp)

    integer, allocatable :: recp(:, :)

    deallocate(recp)

  end subroutine acqfree

end module acquisition
