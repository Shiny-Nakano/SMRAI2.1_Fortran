module reppu_par
  integer, parameter :: nt=45792
!  integer, parameter :: nt=756
!  integer, parameter :: nlon=15, nlat=32
!  integer, parameter :: mlat=15, mlon=32
  integer, parameter :: nlat=30, nlon=80
  integer, parameter :: npca=50

  integer, parameter :: nparams=8
!  integer, parameter :: nparams=9

!  integer, parameter :: nparams=4
!  integer, parameter :: ncomp=10
  integer, parameter :: ncomp=15

  real(8), parameter :: pi=3.14159265358979323846, d2r=0.017453292519943295769
!  real(8), parameter :: pi=3.141592653589793238463
end module reppu_par
