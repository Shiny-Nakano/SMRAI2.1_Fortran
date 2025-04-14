!!
!! SMRAI 2.1 in Fortran
!! 
!!  An emulator of REPPU MHD model based on SMRAI2 (Kataoka et al., 2024)
!!   
!!  coded by S. Nakano (Apr. 2025)
!! 
program reppuassim
  use reppu_par
  use reservoir
  use SRBF
  use polemap
  use mtmod
!$  use omp_lib

  implicit none

  real(8), parameter :: Rearth=6.4e3
  real(8), parameter :: Riono=1.0+100.0/Rearth

  real(8) :: r

  real(8), dimension(nnodes) :: xnodes,dxnodes
  real(8), dimension(nparmax) :: zarr
  real(8), dimension(nnodes) :: xmean

  real(8), dimension(npole) :: pmean, pmap
  integer :: i,j,k,kk

  real(8), dimension(nnodes,npole) :: warr

  character(18) :: swdfile

  character(8) :: date
  character(14) :: assimout

  integer, parameter :: mlat=20
  integer, parameter :: mlon=72
  real(8), dimension(npole,mlat*mlon) :: Atmat
  real(8), dimension(mlat*mlon) :: potmap
  character(14) :: potfile
  real(8), dimension(3,3) :: Convmat
!  real(8), dimension(3) :: rgeo,rsm

  integer :: ihour,imin,iimin
  integer :: iyy, idy, ihh, imm

  real(8) :: xnsw,vsw,tsw,pd
  real(8) :: bx,by,bz
  real(8) :: xhr,xdoy


  call get_command_argument(1,swdfile)

  call genpoles
  call reconst_mat(Atmat)

  zarr(:)=0.0

  call init_reservoir(nparams)

  do i=1,nnodes
    xnodes(i)=tanh(normalrnd())  !! Initializing the state vector
  end do

  open(55,file='wsv_reppu.dat',form='unformatted')
  read(55) pmean
  read(55) xmean
  read(55) warr
  close(55)

!200 format(2i4,2i3,3f9.2,f12.1,3f8.2)

  open(25,file=swdfile,status='old')
  
 do ihour=0,23
    write(assimout,'("pcemltr_",i2.2,".dat")') ihour
    open(16,file=assimout,form='unformatted')

    hourloop: do imin=0,55,5
      read(25,*) iyy,idy,ihh,imm,pd,xnsw,vsw,tsw,bx,by,bz

      zarr(1)=0.2*by
      zarr(2)=0.2*bz
      zarr(3)=2.0*(log10(vsw)-2.5)
      zarr(4)=log10(xnsw)-1.0

      xhr=ihh+imm/60.0
      xdoy=(iyy-2001)+idy+xhr/24.0
      zarr(5)=cos(2*pi*xdoy/365.24)
      zarr(6)=sin(2*pi*xdoy/365.24)
      zarr(7)=cos(2*pi*xhr/24.0)
      zarr(8)=sin(2*pi*xhr/24.0)

      call forward(zarr,xnodes)

      pmap(:)=pmean(:)
      do i=1,nnodes
        dxnodes(i)=xnodes(i)-xmean(i)
      end do

      call dgemv('t',nnodes,npole,1.0d0,warr,nnodes,dxnodes,1,1.0d0,pmap,1)
      call dgemv('t',npole,mlat*mlon,1.0d0,Atmat,npole,pmap,1,0.0d0,potmap,1)

      write(16) pmap

      write(potfile,'("potest",2i2.2,".dat")') ihour, imin
      open(36,file=potfile)
      do j=1,mlat
        do i=1,mlon
          write(36,*) potmap(i+(j-1)*mlon)
        end do
      end do
      close(36)

    end do hourloop
    close(16)
  end do
  close(25)


  call deletepoles

  stop
contains
  subroutine reconst_mat(At)
    real(8), intent(inout), dimension(npole,mlat*mlon) :: At
    real(8) :: dellon, dellat
    real(8) :: olambda, ophi

    dellon = 360.0/mlon
    dellat = 40.0/mlat

    do j=1,mlat
      olambda=(49.0+dellat*j)*D2R
!$omp parallel do private(ophi,k)
      do i=1,mlon
        ophi=(dellon*i-2.5)*D2R

        do k=1,npole
          At(k,i+(j-1)*mlon)=garea(k)*(SRBFstream(olambda,ophi,sgridrlat(k),sgridphi(k)) &
               +SRBFstream(olambda,ophi,-sgridrlat(k),sgridphi(k)))
        end do
      end do
!$omp end parallel do
    end do

    return
  end subroutine reconst_mat

end program reppuassim

