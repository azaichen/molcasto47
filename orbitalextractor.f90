module orbitals
contains

subroutine getorbitals(nbas, coef, occ, sym, basis)

  ! read orbital data from a Molcas-formatted orbital data file (ascii)

  use definitions
  implicit none

!  integer, intent(in) :: irreps 
  integer(KINT), intent(in) :: nbas 
  logical, intent(in) :: sym
  integer(8), dimension(:), intent(in) :: basis

  !real(KREAL), dimension(nbas) :: nrow
  real(KREAL), dimension(nbas,nbas), intent(out) :: coef, occ

  integer(KINT) :: ios,i,j,m,irrep

  character*(LCHARS) :: line

  real(KREAL), parameter :: zero=0.0d0

  ! ============================================================================

  write(out,'(/1x,a)') 'Processing orbital data file ...'

  occ = zero
  coef = zero
  m = 0
  !write(*,*) basis
  ! read data from file unit iuo (opened in main program)
  if (SYM) then 
    do
    read(iuo,*,iostat=ios) line

    ! we'll terminate the loop when the end of the file is reached
    if (ios .ne. 0) exit

    !	  write(*,*) line
    if (trim(line) == '#ORB') then
      write(out,*) ' orbitals found in orbitals file'
      do irrep = 1,size(basis)
        do i = 1,basis(irrep)
          read(iuo,*)
          read(iuo,*) (coef(j+m,i+m), j = 1, basis(irrep))
        enddo
        m = m + basis(irrep)
      enddo
    endif

    if (trim(line) == '#OCC') then
      read(iuo,*)
      write(out,*) ' occupation numbers found in orbital file'
      read(iuo,*) (occ(i,i), i = 1, nbas)
    endif

  enddo

  else
  
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
  endif

  write(out,'(1x,a/)') '... finished processing orbital data file'

  return

end subroutine getorbitals

end module orbitals
