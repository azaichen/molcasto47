module orbitals
contains

subroutine getorbitals(nbas, coef, occ)

  ! read orbital data from a Molcas-formatted orbital data file (ascii)

  use definitions
  implicit none

  integer(KINT), intent(in) :: nbas

  !real(KREAL), dimension(nbas) :: nrow
  real(KREAL), dimension(nbas,nbas), intent(out) :: coef, occ

  integer(KINT) :: ios,i,j

  character*(LCHARS) :: line

  real(KREAL), parameter :: zero=0.0d0

  ! ============================================================================

  write(out,'(/1x,a)') 'Processing orbital data file ...'

  occ = zero
  coef = zero

  ! read data from file unit iuo (opened in main program)

  do
    read(iuo,*,iostat=ios) line

    ! we'll terminate the loop when the end of the file is reached
    if (ios .ne. 0) exit

    !	  write(*,*) line
    if (trim(line) == '#ORB') then
      write(out,*) ' orbitals found in orbitals file'
      do i = 1, Nbas
        read(iuo,*)
        read(iuo,*) (coef(j,i), j = 1, nbas)
      enddo
    endif

    if (trim(line) == '#OCC') then
      read(iuo,*)
      write(out,*) ' occupation numbers found in orbital file'
      read(iuo,*) (occ(i,i), i = 1, nbas)
    endif

  enddo

  write(out,'(1x,a/)') '... finished processing orbital data file'

  return

end subroutine getorbitals

end module orbitals
