program fitPCA
  use reppu_par
  use SECS
  use polemap

  implicit none

!  integer, parameter :: nlearn=nt-576
  integer, parameter :: nlearn=nt
  integer, parameter :: nspinup=10

!  real(8), parameter :: PI=3.14159265358979323846, D2R=0.017453292519943295769

  integer, parameter :: ihstart=0
  integer, parameter :: ihend=21
!  integer, parameter :: nlon=32, nlat=15

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
  ! ! integer, parameter :: mlat=15
  ! ! integer, parameter :: mlon=32
!  real(8), dimension(mlon,mlat,ncomp+1) :: potmap
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
!       352.687, 357.187, 367.688/)

  real(8), dimension(nlat) :: xlarr = (/53.109, 55.172, 57.234, 59.297, 61.359, 63.422, &
       65.484, 66.783, 67.826, 68.870, 69.913, 70.957, 72.000, 73.044, 74.087, 75.130, &
       76.174, 77.217, 78.261, 79.304, 80.348, 81.391, 82.435, 83.478, 84.522, 85.565, &
       86.609, 87.652, 88.696, 89.739/)

!  integer, parameter :: lwork=27087010, liwork=15002
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

  write(6,*) cov(1,1),alpha(1),covrt(1,1)

!   do j=1,mlat
!     olambda=(49.0+2*j)*D2R
!     do i=1,mlon
!       ophi=(5*i-2.5)*D2R

! !      yvec(i+nlon*(j-1))=Ymat(i+nlon*(j-1),1)
!       call set_obs(olambda,ophi)
        
!       do k=1,5
!         sfunc=SECSstream(sgridrlat(k),sgridphi(k))
!         write(6,*) sfunc
!       end do
!     end do
!   end do

!   goto 900

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

  ! open(15,file='PCApattern.dat')
  ! do j=1,nlon
  !   read(15,'(30f16.6)') (Pmean(i+nlat*(j-1),1),i=1,nlat)
  ! end do

  ! do m=1,ncomp
  !   do j=1,nlon
  !     read(15,'(30f16.6)') (Ymat(i+nlat*(j-1),m),i=1,nlat)
  !   end do
  ! end do
  ! close(15)

!$omp parallel do private(ophi,olambda,i,k)
  do j=1,nlon
    ophi=pharr(j)*D2R
    do i=1,nlat
      olambda=xlarr(i)*D2R

! !      yvec(i+nlon*(j-1))=Ymat(i+nlon*(j-1),1)
!       call set_obs(olambda,ophi)

      do k=1,npole
        Hmat(i+nlat*(j-1),k)=garea(k)*SECSstream(olambda,ophi,sgridrlat(k),sgridphi(k))
        Hmat(i+nlat*(j-1),k)=Hmat(i+nlat*(j-1),k)+garea(k)*SECSstream(olambda,ophi,-sgridrlat(k),sgridphi(k))
        ! Hmat(i+nlat*(j-1),k)=garea(k)*SECSstream(sgridrlat(k),sgridphi(k))
        ! Hmat(i+nlat*(j-1),k)=Hmat(i+nlat*(j-1),k)+garea(k)*SECSstream(-sgridrlat(k),sgridphi(k))
      end do
    end do
  end do
!$omp end parallel do

  ! call dgemm('n','n',nlon*nlat,npole,npole,1.0,Hmat,nlon*nlat,covrt,npole,&
  !      0.0,HXmat,nlon*nlat)
  ! call dsyrk('u','n',nlon*nlat,npole,1.0,HXmat,nlon*nlat,0.0,HVHmat,nlon*nlat)

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
  ! call dpotrs('u',nlon*nlat,1,HVHmat,nlon*nlat,Pmean,nlon*nlat,info3)
  ! call dgemm('t','n',npole,1,nlon*nlat,1.0d0,HVmat,nlon*nlat,&
  !      Pmean,nlon*nlat,0.0d0,Wmean,npole)

  call dpotrs('u',nlon*nlat,nlearn,HVHmat,nlon*nlat,Ymat,nlon*nlat,info2)
  call dgemm('t','n',npole,nlearn,nlon*nlat,1.0d0,HVmat,nlon*nlat,&
       Ymat,nlon*nlat,0.0d0,Wmat,npole)

!  write(6,*) HVmat(1,1),Wmat(1,1)

  open(16,file='weightSECS_reppu.dat',form='unformatted')
  write(16) Wmat
  close(16)
 
!   potmap(:,:,:)=0.0
!   dlat = 40.0/nlat
!   dlon = 360.0/nlon
!   do j=1,nlat
!     olambda=(50.0+dlat*(j-0.5))*D2R
! !$omp parallel do private(ophi,sfunc,k,m)
!     do i=1,nlon
! !      ophi=(5*i-2.5)*D2R
!       ophi=(dlon*(i-0.5))*D2R
! !      ophi=pharr(i)*D2R

! ! !      yvec(i+nlon*(j-1))=Ymat(i+nlon*(j-1),1)
! !       call set_obs(olambda,ophi)
        
!       ! do k=1,5
!       !   sfunc=SECSstream(sgridrlat(k),sgridphi(k))
!       !   write(6,*) Wmat(k,1),sfunc
!       ! end do

!       do k=1,npole
!         sfunc=garea(k)*(SECSstream(olambda,ophi,sgridrlat(k),sgridphi(k))+SECSstream(olambda,ophi,-sgridrlat(k),sgridphi(k)))
! !        sfunc=garea(k)*(SECSstream(sgridrlat(k),sgridphi(k))+SECSstream(-sgridrlat(k),sgridphi(k)))

!         potmap(i,j,0)=potmap(i,j,0)+wmean(k,1)*sfunc
!         do m=1,ncomp
!           potmap(i,j,m)=potmap(i,j,m)+Wmat(k,m)*sfunc
!         end do
!       end do
!     end do
! !$omp end parallel do
!   end do

!   do m=0,ncomp
!     write(potfile,'("pot",i2.2,".dat")') m
!     open(36,file=potfile)
!     do i=1,nlon
!       do j=1,nlat
! !      olambda=(49.0+2*j)*D2R
! !        ophi=(5*i-2.5)*D2R
!         write(36,*) potmap(i,j,m)
!       end do
!     end do
!     close(36)
!   end do

900 continue
  
  call deletepoles

end program fitPCA
