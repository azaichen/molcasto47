
module definitions

  implicit none

  ! define the byte length of integer and real variables:
  integer, parameter       :: KINT = kind(1)
  integer(KINT), parameter :: KREAL = kind(1.0d0)
  integer(KINT), parameter :: LCHARS = 160

  ! input/output units
  integer(KINT), parameter :: inp=5, out=6, err=0, iuo=10, iuh = 11, i47=47

end module definitions

