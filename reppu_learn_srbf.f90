program reppuana
  use reservoir
  use reppu_par
  use polemap
  use mtmod
!$  use omp_lib

  implicit none

  integer, parameter :: nlearn=nt-576
  integer, parameter :: nspinup=10

  integer, parameter :: nhist=12
  integer, parameter :: nhistref=1

  real(8) :: r, rexp
  integer, dimension(nhist) :: iearr
  integer :: iesum
  integer :: nls

  real(8) :: al,au
  real(8), dimension(nnodes) :: xnodes,xdnodes
  real(8), dimension(nparmax,nhist) :: zarrhist
  real(8), dimension(nparmax) :: zarr
  real(8), dimension(nnodes,nnodes) :: xcov, xcovinv
  real(8), dimension(nnodes) :: xmean
  real(8), dimension(nnodes,npole) :: xyvec
  real(8) :: yy,xmon,day,xhr
  real(8) :: xdoy, xdoypr

  real(8), dimension(npole) :: yobs
  real(8), dimension(npole,nt) :: Ymat
  real(8), dimension(npole) :: ymean

  integer :: i,j,k,kk
  integer :: jopt,kcount,imon

  real(8), dimension(nnodes) :: sarr
  real(8), dimension(nnodes,npole) :: warr
  real(8), parameter :: elim=1.0e-6

  integer, parameter :: lwork=300330001 ! 3*nnodes**2+(5+2k)*nnodes+1 (k=log(n)/log(2))
  integer, parameter :: liwork=50002 ! 5n+2

  real(8), dimension(lwork) :: work
  integer, dimension(liwork) :: iwork
  integer :: info
  integer :: iopt
  integer, parameter :: iecri=6

  real(8), parameter :: rsig=1.0e-1 !! For nparms=4

  character(3) :: optch
  character(4) :: yeararr
  character(6) :: skipdir
  integer :: iskipyear

  integer, dimension(12) :: nmday = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
  integer, dimension(12) :: nmdoy

  real, parameter :: f0=-1.0
  real, parameter :: pfact=1.0e-3
  real :: dum

  nmdoy(1)=0
  do i=2,12
    nmdoy(i)=nmdoy(i-1)+nmday(i-1)
  end do

  iopt=0 ! AU:0, AE:1, AUlat:2, ALlat:3, AUMLT:4, ALMLT:5
  jopt=0 ! Linear obs:0, Nonlinear obs:1

  xmean(:)=0.0
  xcov(:,:)=0.0
  xcovinv(:,:)=0.0
  zarr(:)=0.0
  zarrhist(:,:)=0.0
  xyvec(:,:)=0.0
  ymean(:)=0.0
  iearr(:)=0

  kk=0
  kcount=0

  call init_reservoir(nparams)

  call sgrnd(111)
  do i=1,nnodes
    xnodes(i)=tanh(normalrnd())
    dum=normalrnd()

    do j=1,nparams
      r=grnd()
      if(r < 0.2) then
        r=grnd()
        dum=normalrnd()
      end if
    end do

    do j=1,nnodes
      r=grnd()
      if(r < 0.2) then
        r=grnd()
        dum=normalrnd()
      end if
    end do
  end do

  open(25,file='weightSECS_reppu.dat',form='unformatted',status='old')
  read(25) Ymat
  close(25)

  open(15,file='swall3_5min.txt',status='old')

  xdoy=0.0

  do k=1,nlearn
    kk=kk+1

    read(15,*) yy,xmon,day,xhr,(zarr(i),i=1,4)

    do j=1,nlon
      yobs(:)=Ymat(:,k)
    end do

    zarr(1)=0.2*zarr(1)
    zarr(2)=0.2*zarr(2)
    zarr(3)=2.0*zarr(3)

    xdoypr=xdoy

    xdoy=1.0*((yy-2001)*365+int(nint(yy-2001)/4))
    imon=nint(xmon)
    if(mod(nint(yy),4)==0.and.imon > 2) then
      xdoy=xdoy+nmdoy(imon)+day+xhr/24.0
    else
      xdoy=xdoy+nmdoy(imon)+day+xhr/24.0-1
    end if

    zarr(5)=cos(2*pi*xdoy/365.24)
    zarr(6)=sin(2*pi*xdoy/365.24)
    zarr(7)=cos(2*pi*xhr/24.0)
    zarr(8)=sin(2*pi*xhr/24.0)


    call forward(zarr,xnodes)

    if(xdoy-xdoypr > 0.1) then
      write(6,'(a,i8,3i6)') 'Reset', kcount,nint(yy),nint(xmon),nint(day)
      kk=0
    end if

    if(kk > nspinup) then
      kcount=kcount+1
!$omp parallel do private(i)
      do j=1,nnodes
        xmean(j)=xmean(j)+xnodes(j)
        do i=1,nnodes
          xcov(i,j)=xcov(i,j)+xnodes(i)*xnodes(j)
        end do
        do i=1,npole
          xyvec(j,i)=xyvec(j,i)+xnodes(j)*yobs(i)
        end do
      end do
!$omp end parallel do

      ymean(:)=ymean(:)+yobs(:)
    end if
  end do
  close(15)

  write(6,*) 'Data read.'

!  nls = nlearn - nspinup
  nls = kcount
  do j=1,nnodes
    xmean(j)=xmean(j)/nls
  end do

  write(6,*) ymean(1),nls

  ymean(:)=ymean(:)/nls
  write(6,*) ymean(1),nls

  do j=1,nnodes
    do i=1,nnodes
      xcov(i,j)=xcov(i,j)/nls-xmean(i)*xmean(j)
    end do
    do i=1,npole
      xyvec(j,i)=xyvec(j,i)/nls-xmean(j)*ymean(i)
    end do
  end do

  write(6,*) 'Comincio EVD.'
  call dsyevd('v','u',nnodes,xcov,nnodes,sarr,work,lwork,iwork,liwork,info)
  write(6,*) 'EVD Done.'

  if(info /= 0) write(6,*) info

  write(6,*) sarr(1),sarr(nnodes)

  do i=1,nnodes
    sarr(i)=1.0/sqrt(sarr(i)+rsig*rsig)
  end do

  do j=1,nnodes
    do i=1,nnodes
      xcov(i,j)=sarr(j)*xcov(i,j)
    end do
  end do

  call dsyrk('u','n',nnodes,nnodes,1.0d0,xcov,nnodes,0.0d0,xcovinv,nnodes)
  call dsymm('l','u',nnodes,npole,1.0d0,xcovinv,nnodes,xyvec,nnodes,0.0d0,warr,nnodes)


  open(56,file='wsv_reppu.dat',status='unknown',form='unformatted')

  write(56) ymean
  write(56) xmean
  write(56) warr

  close(56)

  stop

end program reppuana
