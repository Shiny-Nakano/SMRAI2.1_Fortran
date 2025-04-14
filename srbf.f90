module SRBF
  implicit none

  private
  real(8), parameter :: pi=3.14159265358979323846, d2r=0.017453292519943295769
  real(8), dimension(3) :: obs_azim, obs_lam

!  real(8), parameter :: xkappa = 131.4 !! 1.0/(2.0*(1-cos(5.0*pi/180.0)))
  real(8), parameter :: xkappa = 525.33 !! 1.0/(2.0*(1-cos(2.5*pi/180.0)))
!  real(8), parameter :: xkappa = 820.785 !! 1.0/(2.0*(1-cos(2.0*pi/180.0)))

  public SRBFfact
  public SRBFfact2
  public SRBFstream
!  public set_secs_param
contains
  subroutine product(a,b,c)
    real(8), dimension(3) :: a,b
    real(8), dimension(3) :: c

    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)

    return
  end subroutine product


  real(8) function SRBFstream(olat,ophi,slat,sphi)
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
    SRBFstream = exp(xkappa*((sc(1)*ob(1)+sc(2)*ob(2)+sc(3)*ob(3))/soabs-1.0))

    return
  end function SRBFstream



  subroutine SRBFfact(olat,ophi,slat,sphi,fazm,flam)
    real(8) :: slat,sphi
    real(8), dimension(3) :: sc,secs_azim
    real(8), intent(inout) :: fazm,flam
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

    soabs = sqrt(sc(1)*sc(1)+sc(2)*sc(2)+sc(3)*sc(3))
    sfact=exp(xkappa*((sc(1)*ob(1)+sc(2)*ob(2)+sc(3)*ob(3))/soabs-1.0))

    call product(sc,ob,secs_azim)
    secs_azim(:) = secs_azim(:)*xkappa*sfact

    fazm=-sop*secs_azim(1)+cop*secs_azim(2)
    flam=-sol*(cop*secs_azim(1)+sop*secs_azim(2))+col*secs_azim(3)

    return
  end subroutine SRBFfact


  real(8) function SRBFfact2(olat,ophi,slat,sphi,angle)
    real(8), intent(in) :: slat,sphi,angle
    real(8), intent(in) :: olat,ophi
    real(8) :: faz, fla

    call SRBFfact(olat,ophi,slat,sphi,faz,fla)

    SRBFfact2=faz*cos(angle)+fla*sin(angle)

    return 
  end function SRBFfact2
end module SRBF
