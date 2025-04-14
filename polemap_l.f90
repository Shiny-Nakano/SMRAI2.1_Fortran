module polemap
! !$  use omp_lib
  implicit none

  real(8), private, parameter :: pi=3.14159265358979323846, d2r=0.017453292519943295769

  integer, parameter :: npole=2000

  real(8), allocatable, dimension(:) :: sgriddlat,sgridmlt
  real(8), allocatable, dimension(:) :: sgridrlat,sgridphi
  real(8), allocatable, dimension(:) :: sgridglat,sgridglon
  real(8), allocatable, dimension(:) :: smirglat,smirglon
  real(8), allocatable, dimension(:) :: garea

  integer, private :: istat

contains
  subroutine genpoles
    integer :: k

    allocate(sgriddlat(npole),stat=istat)
    allocate(sgridmlt(npole),stat=istat)
    allocate(sgridrlat(npole),stat=istat)
    allocate(sgridphi(npole),stat=istat)
    allocate(sgridglat(npole),stat=istat)
    allocate(sgridglon(npole),stat=istat)
    allocate(smirglat(npole),stat=istat)
    allocate(smirglon(npole),stat=istat)
    allocate(garea(npole),stat=istat)

    open(98,file='polepoints_sp.dat',status='old')
    do k=1,npole
      read(98,*) sgridrlat(k),sgridphi(k)
      sgriddlat(k)=sgridrlat(k)/d2r

      sgridmlt(k)=12.0/pi*sgridphi(k)
      if(sgridmlt(k) < 0.0) then
        sgridmlt(k)=sgridmlt(k)+24.0
      end if

      garea(k)=2.0*pi/npole
    end do
    close(98)

    return
  end subroutine genpoles


!   subroutine recalc_geogra
!     real(8), dimension(3,3) :: Cmat
!     real(8), dimension(3) :: rsm,rgeo
!     real(8) :: xmag,ymag,zmag
!     real(8) :: riono,rho
!     integer :: i,k
!     real(8) :: theta, thetam, phi2

!     do i=1,3
!       rsm(:)=0.0
!       rsm(i)=1.0
!       call magsm(xmag,ymag,zmag,rsm(1),rsm(2),rsm(3),-1)
!       call geomag(Cmat(1,i),Cmat(2,i),Cmat(3,i),xmag,ymag,zmag,-1)
!     end do

!     riono=1.0+100.0/6.4e3

! !$omp parallel do private(theta,thetam,rho,phi2,rsm,rgeo)
!     do k=1,npole
!       theta=0.5*pi-sgridrlat(k)
!       thetam=0.5*pi+sgridrlat(k)
!       phi2=sgridphi(k)-pi

!       call sphcar(riono,theta,phi2,rsm(1),rsm(2),rsm(3),1)
!       call dgemv('n',3,3,1.0d0,Cmat,3,rsm,1,0.0d0,rgeo,1)
!       call sphcar(rho,theta,sgridglon(k),rgeo(1),rgeo(2),rgeo(3),-1)
!       sgridglat(k)=0.5*pi-theta

!       call sphcar(riono,thetam,phi2,rsm(1),rsm(2),rsm(3),1)
!       call dgemv('n',3,3,1.0d0,Cmat,3,rsm,1,0.0d0,rgeo,1)
!       call sphcar(rho,thetam,smirglon(k),rgeo(1),rgeo(2),rgeo(3),-1)
!       smirglat(k)=0.5*pi-thetam
!     end do
! !$omp end parallel do

!     return
!   end subroutine recalc_geogra


!   ! subroutine recalc_geogra
!   !   real(8) :: xgeo,ygeo,zgeo
!   !   real(8) :: xmag,ymag,zmag
!   !   real(8) :: xsm,ysm,zsm
!   !   real(8) :: riono,rho
!   !   integer :: k
!   !   real(8) :: theta, thetam, phi2

!   !   riono=1.0+100.0/6.4e3

!   !   k=1
!   !   do k=1,npole
!   !     theta=0.5*pi-sgridrlat(k)
!   !     thetam=0.5*pi+sgridrlat(k)
!   !     phi2=sgridphi(k)-pi

!   !     call sphcar(riono,theta,phi2,xsm,ysm,zsm,1)
!   !     call magsm(xmag,ymag,zmag,xsm,ysm,zsm,-1)
!   !     call geomag(xgeo,ygeo,zgeo,xmag,ymag,zmag,-1)
!   !     call sphcar(rho,theta,sgridglon(k),xgeo,ygeo,zgeo,-1)
!   !     sgridglat(k)=0.5*pi-theta

!   !     call sphcar(riono,thetam,phi2,xsm,ysm,zsm,1)
!   !     call magsm(xmag,ymag,zmag,xsm,ysm,zsm,-1)
!   !     call geomag(xgeo,ygeo,zgeo,xmag,ymag,zmag,-1)
!   !     call sphcar(rho,thetam,smirglon(k),xgeo,ygeo,zgeo,-1)
!   !     smirglat(k)=0.5*pi-thetam
!   !   end do

!   !   return
!   ! end subroutine recalc_geogra


!   subroutine geo2smmat(Cmat)
!     real(8), intent(inout), dimension(3,3) :: Cmat
!     real(8), dimension(3) :: rgeo
!     real(8) :: xmag,ymag,zmag
!     integer :: i

!     do i=1,3
!       rgeo(:)=0.0
!       rgeo(i)=1.0
!       call geomag(rgeo(1),rgeo(2),rgeo(3),xmag,ymag,zmag,1)
!       call magsm(xmag,ymag,zmag,Cmat(1,i),Cmat(2,i),Cmat(3,i),1)
!     end do

!     return
!   end subroutine geo2smmat


!   subroutine calcvar(fs)
!     real(8) :: fs
!     integer :: i

!     cov(:,:)=0.0
!     covinv(:,:)=0.0

!     do i=1,npole
!       cov(i,i)=1.0
! !      cov(i,i)=sin(sgridrlat(i))**(2*fs)
! !      cov(i,i) = exp(-fs*cos(sgridrlat(i)**2)

! !      cov(i,i)=exp(8.0*sin(sgridrlat(i))**2)
! !      cov(i,i)=exp(-0.5*((sgriddlat(i)-90.0)/20.0)**2)

! !      cov(i,i)=exp(-0.5*((sgriddlat(i)-90.0)/200.0)**2)
! !      cov(i,i)=exp(-0.5*((sgriddlat(i)-90.0)/20.0)**2)
!       covinv(i,i)=1.0/cov(i,i)
!     end do

! !     cov(npole,npole)=1.0
! !     covinv(npole,npole)=1.0

!     return
!   end subroutine calcvar

  
!   subroutine calccov(fs)
!     real(8) :: fs
!     real(8), parameter :: xkappac = 32.92 !! 1.0/(2.0*(1-cos(10.0*pi/180.0)))
! !    real(8), parameter :: xkappac = 14.674 !! 1.0/(2.0*(1-cos(15.0*pi/180.0)))
! !    real(8), parameter :: xkappac = 3.73 !! 1.0/(2.0*(1-cos(30.0*pi/180.0)))
! !    real(8), parameter :: xkappac = 131.4 !! 1.0/(2.0*(1-cos(5.0*pi/180.0)))
!     real(8) :: dist
!     real(8), dimension(3) :: r1, r2
! !    real(8), parameter :: xi2=100.0
!     real(8), parameter :: xi2=0.01
! !    integer, parameter :: lwork=50030001
!     integer, parameter :: lwork=75155001
!     real(8), dimension(lwork) :: work
!     integer, parameter :: liwork=25003
!     integer, dimension(liwork) :: iwork
!     real(8), parameter :: elim1=1.0e-8
! !    real(8) :: a1,a2,a3
!     real(8), dimension(3) :: alpha, beta
!     real(8) :: a1, a2, a3
! !    real(8), parameter :: delta=0.3
! !    real(8), parameter :: delta=0.5
!     real(8), parameter :: delta=1.0
!     real(8) :: gamma1, gamma2, q
!     real(8) :: sblat, cfact
!     integer :: i,j
!     integer :: info

!     sblat = sin(boundlat*d2r)
!     cov(:,:)=0.0

!     do j=1,npole
!       r2(1)=delta*cos(sgridphi(j))*cos(sgridrlat(j))
!       r2(2)=delta*sin(sgridphi(j))*cos(sgridrlat(j))
!       r2(3)=sin(sgridrlat(j))
!       gamma2=sqrt( (delta*cos(sgridrlat(j)))**2 + sin(sgridrlat(j))**2 )

!       do i=1,npole
!         r1(1)=delta*cos(sgridphi(i))*cos(sgridrlat(i))
!         r1(2)=delta*sin(sgridphi(i))*cos(sgridrlat(i))
!         r1(3)=sin(sgridrlat(i))
!         gamma1=sqrt( (delta*cos(sgridrlat(i)))**2 + sin(sgridrlat(i))**2 )

!         q = (r1(1)*r2(1)+r1(2)*r2(2)+r1(3)*r2(3)) / (gamma1*gamma2)

! !        cov(i,j)=exp(xkappac*(r1(1)*r2(1)+r1(2)*r2(2)+r1(3)*r2(3)-1.0))
!         cfact=(sin(sgridrlat(i))-sblat)*(sin(sgridrlat(j))-sblat)/(1-sblat)**2
!         cov(i,j)=cfact*exp( xkappac*(q-1.0) )
!       end do

!       cov(j,j)=cov(j,j)+xi2
!     end do

!     covinv(:,:)=cov(:,:)

!     call dsyevd('v','u',npole,covinv,npole,zeta,work,lwork,iwork,liwork,info)
! !    write(6,*) work(1),iwork(1),info

!     do i=1,npole
!       if(zeta(i)/zeta(npole) < elim1) then
!         if(zeta(i) < 0.0) then
!           write(6,*) 'neg eig', zeta(i)
!         end if

!         zeta(i)=0.0
!       else
!         zeta(i)=1.0/sqrt(zeta(i))
!       end if
!     end do

!     do j=1,npole
!       do i=1,npole
!         hhrt(i,j)=zeta(j)*covinv(i,j)
!       end do
!     end do

!     call dsyrk('u','n',npole,npole,1.0d0,hhrt,npole,0.0d0,covinv,npole)

!     do j=1,npole-1
!       do i=j+1,npole
!         covinv(i,j)=covinv(j,i)
!       end do
!     end do

!     return
!   end subroutine calccov

  
  subroutine deletepoles
    deallocate(sgriddlat,stat=istat)
    deallocate(sgridmlt,stat=istat)
    deallocate(sgridrlat,stat=istat)
    deallocate(sgridphi,stat=istat)
    deallocate(sgridglat,stat=istat)
    deallocate(sgridglon,stat=istat)
    deallocate(smirglat,stat=istat)
    deallocate(smirglon,stat=istat)

    deallocate(garea,stat=istat)

    return
  end subroutine deletepoles
end module polemap
