!!
!! Echo state network model
!!
!!  coded by S. Nakano (May 2021)
!!
module reservoir
  use mtmod
!$  use omp_lib

  implicit none
  private

  real(8), parameter :: pi=3.141592653589793238463

  integer, public, parameter :: nparmax=10
  integer :: npars

  integer, public, parameter :: nnodes=1500

  integer, dimension(nnodes) :: nedge
  integer, dimension(nnodes,nnodes) :: iedge

  real(8), parameter :: delta=0.999
!  real(8), parameter :: delta=0.98

  real(8), parameter :: winrate=0.2

  real(8), dimension(nnodes,nnodes) :: wmat, wmatsv
  real(8), dimension(nnodes,nparmax) :: winput
  real(8), dimension(nnodes) :: gmm
  real(8), dimension(nnodes) :: bias


  public init_reservoir
  public forward
  public normalrnd
contains
  subroutine init_reservoir(np)
    integer, intent(in) :: np 
    integer :: i,j,k
    real(8) :: r

    integer, parameter :: lwork=5*nnodes
    real(8), dimension(lwork) :: work,rwork
    real(8), dimension(nnodes) :: svec
    real(8), dimension(nnodes,nnodes) :: umat,vtmat
    integer :: info
    integer :: nconnect
    real(8) :: dum

    npars=np

    call sgrnd(111)

    winput(:,:)=0.0
    wmat(:,:)=0.0

    do i=1,nnodes
      dum=tanh(normalrnd())
      bias(i)=0.3*normalrnd()

      do j=1,npars
        r=grnd()
        if(r < winrate) winput(i,j)=weight_val()
      end do

      k=0

      do j=1,nnodes
        r=grnd()
        if(r < winrate) wmat(i,j)=weight_val()
      end do
    end do

    wmatsv(:,:)=wmat(:,:)

    call dgesvd('A','A',nnodes,nnodes,wmatsv,nnodes,svec,umat,nnodes,&
         vtmat,nnodes,work,lwork,info)

    do j=1,nnodes
      do i=1,nnodes
        wmat(i,j)=delta*wmat(i,j)/svec(1)
      end do
    end do

    return
  end subroutine init_reservoir


  real(8) function weight_val()
    real(8) :: r
    real(8), parameter :: wscale=1.0

    r=grnd()
    weight_val=wscale*normalrnd()

    return
  end function weight_val


  subroutine forward(zinput,xnds)
    real(8), dimension(nparmax) :: zinput
    real(8), dimension(nnodes) :: xnds
    real(8), dimension(nnodes) :: xnew
    real(8) :: xs
    integer :: i,j
!    real(8), parameter :: alpha=0.5
    real(8), parameter :: alpha=0.0

    call dgemv('n',nnodes,npars,1.0d0,winput,nnodes,zinput,1,0.0d0,xnew,1)
    call dgemv('n',nnodes,nnodes,1.0d0,wmat,nnodes,xnds,1,1.0d0,xnew,1)

!$omp parallel do
    do i=1,nnodes
      xnds(i)=alpha*xnds(i)+(1.0-alpha)*tanh(xnew(i)+bias(i))
    end do
!$omp end parallel do

  end subroutine forward

  real(8) function normalrnd()
    real(8) :: a, b
!    real(8), parameter :: pi=3.141592653589793238463

    a=(1.0d0-grnd())
    b=(1.0d0-grnd())

    normalrnd=sqrt(-2.0d0*log(a))*sin(2.0d0*pi*b)

    return
  end function normalrnd

end module reservoir
