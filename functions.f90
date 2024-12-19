!This modules contains functions and normalization factors

module functions

contains

real(KIND=8) function norm(alpha,l)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: alpha
REAL(KIND=8), parameter :: pi = 3.14159265358979323
INTEGER, INTENT(IN) :: l

norm = sqrt((2.0D0*alpha/pi)**(3.0D0/2.0D0)*(4.0D0*alpha)**l/dfact(2*l-1))

!write(*,*) alpha, l, norm

end function norm

real(KIND=8) function dfact(n)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
INTEGER :: i

if (n<-1) then
stop 'Error: Negative factorial!'
endif

dfact = 1
if (n == -1 .or. n == 0 .or. n == 1) then
dfact = 1
else
do i = n, 1, -2
dfact = dfact * i
end do
end if

end function dfact

end module functions
