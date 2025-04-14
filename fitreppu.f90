program fitPCA
  use reppu_par
  use SRBF
  use polemap

  implicit none

  integer, parameter :: nlearn=nt
  integer, parameter :: nspinup=10

  integer, parameter :: ihstart=0
  integer, parameter :: ihend=21

  real(8), dimension(nlon*nlat,nlearn) :: Ymat
  real(8), dimension(nlon*nlat,1) :: Pmean
  real(8), dimension(npole) :: alpha
  real(8), dimension(npole,nlearn) :: Wmat
  real(8), dimension(npole,1) :: Wmean
  real(8), dimension(npole,npole) :: covrt
  real(8), dimension(nlon*nlat,npole) :: Hmat
  real(8), dimension(nlon*nlat,npole) :: HVmat
  real(8), dimension(nlon*nlat,nlon*nlat) :: HVHmat
  real(8), parameter :: rsig = 1.0d-8, xi=1.0d2
  real(8), parameter :: eps = 1.0d-8
  real(8) :: ophi,olambda

  integer, parameter :: mlat=20
  integer, parameter :: mlon=72

  real(8), dimension(nlon,nlat,0:ncomp) :: potmap

  real(8), dimension(nlon) :: pharr = (/1.688, 6.188, 10.688, 15.188, 19.688, 24.188, 28.688, &
       33.188, 37.688, 42.188, 46.688, 51.188, 55.688, 60.188, 64.688, 69.188, 73.688, 78.188, &
       82.688, 87.188, 91.688, 96.188, 100.688, 105.188, 109.688, 114.188, 118.688, 123.188, &
       127.688, 132.188, 136.688, 141.188, 145.687, 150.188, 154.688, 159.188, 163.688, 168.187, &
       172.687, 177.187, 181.687, 186.187, 190.687, 195.187, 199.687, 204.187, 208.687, 213.187, &
       217.687, 222.187, 226.687, 231.187, 235.687, 240.187, 244.687, 249.187, 253.687, 258.187, &
       262.687, 267.187, 271.687, 276.187, 280.687, 285.187, 289.687, 294.187, 298.687, 303.187, &
       307.687, 312.187, 316.687, 321.187, 325.687, 330.187, 334.687, 339.187, 343.687, 348.187, &
       352.687, 357.187/)

  real(8), dimension(nlat) :: xlarr = (/53.109, 55.172, 57.234, 59.297, 61.359, 63.422, &
       65.484, 66.783, 67.826, 68.870, 69.913, 70.957, 72.000, 73.044, 74.087, 75.130, &
       76.174, 77.217, 78.261, 79.304, 80.348, 81.391, 82.435, 83.478, 84.522, 85.565, &
       86.609, 87.652, 88.696, 89.739/)

  integer, parameter :: lwork=75155001, liwork=25003

  real(8), dimension(lwork) :: work
  integer, dimension(liwork) :: iwork
  character(9) :: potfile

  integer :: i,j,k,m
  integer :: iyear,idoy
  integer :: info,info1,info2,info3
  real(8) :: sfunc

  real :: dlon, dlat

  real, parameter :: f0=-1.0
  real, parameter :: pfact=1.0e-3

! theta=theta-90.
! theta=theta/360*np.pi*2

! #north

  call genpoles
  call calccov(xi)
  covrt(:,:)=covinv(:,:)
  
  call dsyevd('v','u',npole,covrt,npole,alpha,work,lwork,iwork,liwork,info)
  write(6,*) alpha(1),alpha(npole)
  do i=1,m
    if(alpha(i) > eps*alpha(m)) then
      alpha(i)=1.0/sqrt(alpha(i))
    end if
  end do

  do j=1,m
    do i=1,m
!      covrt(i,j)=zeta(j)*covrt(i,j)
      covrt(i,j)=alpha(j)*covrt(i,j)
    end do
  end do

  call dsyrk('u','n',npole,npole,1.0d0,covrt,npole,0.0d0,cov,npole)

  ! date='20150317'
  ! iyear=2015
  ! idoy=76
!  date='20170327'
  iyear=2017
  idoy=86


  open(15,file='fpot.dat',status='old')
  do k=1,nlearn
!    kk=kk+1

    do j=1,nlon
      read(15,'(30f16.6)') (Ymat(i+nlat*(j-1),k),i=1,nlat)

      do i=1,nlat
        Ymat(i+nlat*(j-1),k)=f0*pfact*Ymat(i+nlat*(j-1),k)
      end do

    end do
    read(15,*)
  end do
  close(15)

!$omp parallel do private(ophi,olambda,i,k)
  do j=1,nlon
    ophi=pharr(j)*D2R
    do i=1,nlat
      olambda=xlarr(i)*D2R

      do k=1,npole
        Hmat(i+nlat*(j-1),k)=garea(k)*SRBFstream(olambda,ophi,sgridrlat(k),sgridphi(k))
        Hmat(i+nlat*(j-1),k)=Hmat(i+nlat*(j-1),k)+garea(k)*SRBFstream(olambda,ophi,-sgridrlat(k),sgridphi(k))
      end do
    end do
  end do
!$omp end parallel do


!! HV
  call dsymm('r','u',nlon*nlat,npole,1.0d0,cov,npole,Hmat,nlon*nlat,&
       0.0d0,HVmat,nlon*nlat)
!! HVH^T
  call dgemm('n','t',nlon*nlat,nlon*nlat,npole,1.0d0,HVmat,nlon*nlat,&
       Hmat,nlon*nlat,0.0d0,HVHmat,nlon*nlat)

  write(6,*) Hmat(1,1),HVHmat(1,1)

!! HVH^T+R
  do i=1,nlon*nlat
    HVHmat(i,i)=HVHmat(i,i)+rsig
  end do

  call dpotrf('u',nlon*nlat,HVHmat,nlon*nlat,info1)
  call dpotrs('u',nlon*nlat,nlearn,HVHmat,nlon*nlat,Ymat,nlon*nlat,info2)
  call dgemm('t','n',npole,nlearn,nlon*nlat,1.0d0,HVmat,nlon*nlat,&
       Ymat,nlon*nlat,0.0d0,Wmat,npole)

  open(16,file='weightSRBF_reppu.dat',form='unformatted')
  write(16) Wmat
  close(16)
 
900 continue
  
  call deletepoles

end program fitPCA
