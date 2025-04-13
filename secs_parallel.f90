module SECS
  implicit none

!  real(8), dimension(3) :: ob

  private
  real(8), parameter :: pi=3.14159265358979323846, d2r=0.017453292519943295769
! !  real(8) :: fazm,flam
!   private : product
  real(8), dimension(3) :: obs_azim, obs_lam

!  real(8), parameter :: xkappa = 32.92 !! 1.0/(2.0*(1-cos(10.0*pi/180.0)))
!  real(8), parameter :: xkappa = 51.4 !! 1.0/(2.0*(1-cos(8.0*pi/180.0)))
!  real(8), parameter :: xkappa = 131.4 !! 1.0/(2.0*(1-cos(5.0*pi/180.0)))
!   real(8), parameter :: xkappa = 205.26 !! 1.0/(2.0*(1-cos(4.0*pi/180.0)))
  real(8), parameter :: xkappa = 525.33 !! 1.0/(2.0*(1-cos(2.5*pi/180.0)))
!  real(8), parameter :: xkappa = 820.785 !! 1.0/(2.0*(1-cos(2.0*pi/180.0)))

!  public set_obs
  public SECSfact
  public SECSfact2
  public SECSstream
!  public set_secs_param
contains
!   subroutine set_secs_param(xk)
!     xkappa=xk
!     return
!   end subroutine set_secs_param

  subroutine product(a,b,c)
    real(8), dimension(3) :: a,b
    real(8), dimension(3) :: c

    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)

    return
  end subroutine product


  ! subroutine set_obs(olat,ophi,ob)
  !   real(8), intent(in) :: olat,ophi
  !   real(8), dimension(3), intent(inout) :: ob
  !   real(8) :: cop,sop,col,sol

  !   cop = cos(ophi)
  !   sop = sin(ophi)
  !   col = cos(olat)
  !   sol = sin(olat)

  !   ob(1) = cop*col
  !   ob(2) = sop*col
  !   ob(3) = sol
  !   ! oabs = sqrt(ob(1)*ob(1)+ob(2)*ob(2)+ob(3)*ob(3))
  !   ! write(6,*) oabs

  !   obs_azim(1) = - sop
  !   obs_azim(2) = cop
  !   obs_azim(3) = 0.0

  !   obs_lam(1) = - cop * sol
  !   obs_lam(2) = - sop * sol
  !   obs_lam(3) = col

  !   return
  ! end subroutine set_obs


  real(8) function SECSstream(olat,ophi,slat,sphi)
    real(8), intent(in) :: olat,ophi
    real(8), intent(in) :: slat,sphi
    real(8), dimension(3) :: sc
    real(8) :: soabs
    real(8), dimension(3) :: ob
    real(8) :: col

    col = cos(olat)
    ob(1) = cos(ophi)*col
    ob(2) = sin(ophi)*col
    ob(3) = sin(olat)

    sc(1)=cos(sphi)*cos(slat)
    sc(2)=sin(sphi)*cos(slat)
    sc(3)=sin(slat)

    soabs = sqrt(sc(1)*sc(1)+sc(2)*sc(2)+sc(3)*sc(3))
!     soabs = sqrt(sc(1)*sc(1)+sc(2)*sc(2)+sc(3)*sc(3))*oabs
! !    theta = acos(ddot(3,sc,1,ob,1)/soabs)

    SECSstream = exp(xkappa*((sc(1)*ob(1)+sc(2)*ob(2)+sc(3)*ob(3))/soabs-1.0))
!    SECSstream = exp(xkappa*(cos(theta)-1.0))

    return
  end function SECSstream



  subroutine SECSfact(olat,ophi,slat,sphi,fazm,flam)
    real(8) :: slat,sphi
!    real(8) :: olat,ophi
    real(8), dimension(3) :: sc,secs_azim
!    real(8), dimension(3) :: obs_azim, obs_lam
    real(8), intent(inout) :: fazm,flam
!    real(8), parameter :: xkappa = 500.0
    real(8) :: sfact
    real(8) :: soabs
    real(8), dimension(3) :: ob
    real(8), intent(in) :: olat,ophi
    real(8) :: cop,sop,col,sol

    sc(1)=cos(sphi)*cos(slat)
    sc(2)=sin(sphi)*cos(slat)
    sc(3)=sin(slat)

    cop = cos(ophi)
    sop = sin(ophi)
    col = cos(olat)
    sol = sin(olat)

    ob(1) = cop*col
    ob(2) = sop*col
    ob(3) = sol

    ! obs_azim(1) = - sop
    ! obs_azim(2) = cop
    ! obs_azim(3) = 0.0

    ! obs_lam(1) = - cop * sol
    ! obs_lam(2) = - sop * sol
    ! obs_lam(3) = col

    soabs = sqrt(sc(1)*sc(1)+sc(2)*sc(2)+sc(3)*sc(3))
    sfact=exp(xkappa*((sc(1)*ob(1)+sc(2)*ob(2)+sc(3)*ob(3))/soabs-1.0))

    call product(sc,ob,secs_azim)
    secs_azim(:) = secs_azim(:)*xkappa*sfact

    ! fazm = ddot(3,obs_azim,1,secs_azim,1)
    ! flam = ddot(3,obs_lam,1,secs_azim,1)

    ! fazm=obs_azim(1)*secs_azim(1)+obs_azim(2)*secs_azim(2)+obs_azim(3)*secs_azim(3)
    ! flam=obs_lam(1)*secs_azim(1)+obs_lam(2)*secs_azim(2)+obs_lam(3)*secs_azim(3)

    fazm=-sop*secs_azim(1)+cop*secs_azim(2)
    flam=-sol*(cop*secs_azim(1)+sop*secs_azim(2))+col*secs_azim(3)


!!  kappa**2 * sin(theta) * exp(kappa*(cos(theta)-1.0))
!!   * sin(theta) is included in the length of the cross product
!!   *   ( |a x b| = |a||b| sin(theta) )

!     if(abs(theta) > 1.0e-6) then
!       call product(sc,ob,secs_azim)
!       secs_azim(:) = secs_azim(:) * cos(theta/2) / (sin(theta) * sin(theta/2))
! !!  cot(theta/2) / sin(theta)
! !!   * sin(theta) is for normalizing the length of the cross product
! !!   *   ( |a x b| = |a||b| sin(theta) )

!       fazm = ddot(3,obs_azim,1,secs_azim,1)
!       flam = ddot(3,obs_lam,1,secs_azim,1)
!     else
!       fazm = 0.0
!       flam = 0.0
!     end if

    return
  end subroutine SECSfact


  real(8) function SECSfact2(olat,ophi,slat,sphi,angle)
    real(8), intent(in) :: slat,sphi,angle
    real(8), intent(in) :: olat,ophi
    real(8) :: faz, fla

    call SECSfact(olat,ophi,slat,sphi,faz,fla)
!    return fp*cos(angle)+fl*sin(angle)


!    SECSfact2=-faz*cos(angle)-fla*sin(angle)
    SECSfact2=faz*cos(angle)+fla*sin(angle)

    return 
  end function SECSfact2
end module SECS
