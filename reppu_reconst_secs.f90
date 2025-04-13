program reppureconst
  use reservoir
  use reppu_par
  use SECS
  use polemap
  use mtmod

  implicit none

!  integer, parameter :: nlearn=664
  integer, parameter :: nlearn=nt-576
  ! integer, parameter :: nt=756

  ! integer, parameter :: nparams=4
  ! integer, parameter :: ncomp=10

!  integer, parameter :: nparams=10

  integer, parameter :: nspinup=10

!  real(8), parameter :: xi=0.75
  real(8), parameter :: xi=0.5
!    real(8), parameter :: xi=1.0
  integer, parameter :: nhist=12
  integer, parameter :: nhistref=1

  real(8) :: r, rexp
  integer :: iesum
  integer :: nls

!  integer, parameter :: nx=nlat*nlon

  real(8) :: al,au
  real(8), dimension(nnodes) :: xnodes,xdnodes
  real(8), dimension(nparmax,nhist) :: zarrhist
  real(8), dimension(nparmax) :: zarr
  ! real(8), dimension(nnodes,nnodes) :: xcov, xcovinv
  real(8), dimension(nnodes) :: xmean
  ! real(8), dimension(nnodes,nx) :: xyvec
  real(8) :: yy,xmon,day,xhr
  real(8) :: xdoy

  real(8), dimension(npole) :: wpred
  real(8), dimension(nlat*nlon) :: ppred
  real(8), dimension(npole) :: ymean

!  real(8), parameter :: xlambda=10000.0
!  real(8), parameter :: xlambda=2500.0
!  real(8), parameter :: xlambda=400.0
!  real(8), parameter :: xlambda=100.0
!  real(8), parameter :: xlambda=9.0
  real(8), parameter :: xlambda=1.0
  integer :: i,j,k,kk,l
!  real(8) :: check
!  real(8) :: alpred, aupred
!  real(8) :: pd,xnsw,vsw,tsw,bx,by,bz,alobs,auobs,sym
!  integer :: iday
!  real(8) :: hour
  integer :: jopt,kcount,imon

  real(8), dimension(nnodes) :: sarr
  real(8), dimension(nnodes,npole) :: warr
  real(8), parameter :: elim=1.0e-6

  integer, parameter :: lwork=300330001 ! 3*nnodes**2+(5+2k)*nnodes+1 (k=log(n)/log(2))
  integer, parameter :: liwork=50002 ! 5n+2

!   integer, parameter :: lwork=50450433 ! 3*nnodes**2+(5+2k)*nnodes+1 (k=log(n)/log(2))
!   integer, parameter :: liwork=20482 ! 5n+2
! !  integer, parameter :: lwork=3171400 ! 3*nnodes**2+(5+2k)*nnodes+1 (k=log(n)/log(2))
  real(8), dimension(lwork) :: work
  integer, dimension(liwork) :: iwork
  integer :: info
  integer :: iopt
  integer, parameter :: iecri=6
  character(3) :: optch
  character(4) :: yeararr
  character(6) :: skipdir
  integer :: iskipyear

  integer, dimension(12) :: nmday = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
  integer, dimension(12) :: nmdoy

  real(8), dimension(nlon*nlat,npole) :: Hmat
  real(8) :: ophi,olambda
  real(8), dimension(nlon) :: pharr = (/1.688, 6.188, 10.688, 15.188, 19.688, 24.188, 28.688, &
       33.188, 37.688, 42.188, 46.688, 51.188, 55.688, 60.188, 64.688, 69.188, 73.688, 78.188, &
       82.688, 87.188, 91.688, 96.188, 100.688, 105.188, 109.688, 114.188, 118.688, 123.188, &
       127.688, 132.188, 136.688, 141.188, 145.687, 150.188, 154.688, 159.188, 163.688, 168.187, &
       172.687, 177.187, 181.687, 186.187, 190.687, 195.187, 199.687, 204.187, 208.687, 213.187, &
       217.687, 222.187, 226.687, 231.187, 235.687, 240.187, 244.687, 249.187, 253.687, 258.187, &
       262.687, 267.187, 271.687, 276.187, 280.687, 285.187, 289.687, 294.187, 298.687, 303.187, &
       307.687, 312.187, 316.687, 321.187, 325.687, 330.187, 334.687, 339.187, 343.687, 348.187, &
       352.687, 357.187/)
!       352.687, 357.187, 367.688/)

  real(8), dimension(nlat) :: xlarr = (/53.109, 55.172, 57.234, 59.297, 61.359, 63.422, &
       65.484, 66.783, 67.826, 68.870, 69.913, 70.957, 72.000, 73.044, 74.087, 75.130, &
       76.174, 77.217, 78.261, 79.304, 80.348, 81.391, 82.435, 83.478, 84.522, 85.565, &
       86.609, 87.652, 88.696, 89.739/)


  nmdoy(1)=0
  do i=2,12
    nmdoy(i)=nmdoy(i-1)+nmday(i-1)
  end do


  iopt=0 ! AU:0, AE:1, AUlat:2, ALlat:3, AUMLT:4, ALMLT:5
  jopt=0 ! Linear obs:0, Nonlinear obs:1

!  call get_command_argument(1,optch)
!  read(optch,*,iostat=ios) iopt
!  if(ios/=0) iopt=0

!   call get_command_argument(1,yeararr)
!   read(yeararr,*,iostat=ios) iskipyear

! !  write(skipdir,'("l",i4.4,"_")') iskipyear
!   write(skipdir,'("l",i4.4,"/")') iskipyear

  call genpoles

  xmean(:)=0.0
  zarr(:)=0.0
  zarrhist(:,:)=0.0

  kk=0
  kcount=0

  call init_reservoir(nparams)
  do i=1,nnodes
    xnodes(i)=tanh(normalrnd())
  end do

  open(55,file='wsv_reppu.dat',form='unformatted')

  read(55) ymean
  read(55) xmean
  read(55) warr

  close(55)


!  open(15,file='sw756.txt')
  open(15,file='swall3_5min.txt')
  do k=1,nt
    kk=kk+1

!    read(15,*) (zarr(i),i=1,nparams)
!    read(15,*) yy,xmon,day,xhr,(zarr(i),i=1,nparams)
    read(15,*) yy,xmon,day,xhr,(zarr(i),i=1,4)
    zarr(1)=0.2*zarr(1)
    zarr(2)=0.2*zarr(2)
    ! zarr(1)=0.1*zarr(1)
    ! zarr(2)=0.1*zarr(2)
    zarr(3)=2.0*zarr(3)
!    zarr(3)=5.0*zarr(3)
    ! zarr(4)=20.0*zarr(4)

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
    ! ! zarr(3)=10.0*zarr(3)
    ! ! zarr(4)=20.0*zarr(4)


    call forward(zarr,xnodes)

    wpred(:)=ymean(:)
    do i=1,nnodes
      xdnodes(i)=xnodes(i)-xmean(i)
    end do
    call dgemv('t',nnodes,npole,1.0d0,warr,nnodes,xdnodes,1,1.0d0,wpred,1)

    if(k==45673) then
!$omp parallel do private(ophi,olambda,i,k)
      do j=1,nlon
        ophi=pharr(j)*D2R
        do i=1,nlat
          olambda=xlarr(i)*D2R

          do l=1,npole
            Hmat(i+nlat*(j-1),l)=garea(l)*SECSstream(olambda,ophi,sgridrlat(l),sgridphi(l))
            Hmat(i+nlat*(j-1),l)=Hmat(i+nlat*(j-1),l)+garea(l)*SECSstream(olambda,ophi,-sgridrlat(l),sgridphi(l))
          end do
        end do
      end do
!$omp end parallel do

      call dgemv('n',nlat*nlon,npole,1.0d0,Hmat,nlat*nlon,wpred,1,1.0d0,ppred,1)

      open(36,file='predmap_secs.dat')
      do j=1,nlon
        write(36,'(30f16.6)') (ppred(i+nlat*(j-1)),i=1,nlat)
      end do
      close(36)
    end if

  end do
  close(15)

  stop
end program reppureconst
  
