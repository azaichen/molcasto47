!Read data and Mayer and Wiberg Bond Indexes from NBO 47 files
!Not completed

program test_read
implicit none 
integer, parameter :: nb=11
character(len=100) :: F47
real(8), dimension(:,:), allocatable :: Pf, Sf
integer, dimension(:), allocatable :: cent
real(8), dimension(:,:), allocatable :: Ma, Wib
real(8) :: pop
integer :: k, p

call getarg(1,F47)

call READPS(F47, Pf, Sf, nb)

!call CENTERS(F47, cent, nb)

!write(*,*) cent

!call MAYERBI(Pf, Sf, cent, Ma, Wib, nb)

write(*,'(4E20.12)') Sf 

!write(*,*) 'Mayer BI:'

!do k = 1,4
!   write(*,'(4F8.4)') (Ma(k, p), p = 1,4)
!end do

!write(*,*) 

!write(*,*) 'Wiberg BI:'

!do k = 1,4
!   write(*,'(4F8.4)') (Wib(k, p), p = 1,4)
!end do

contains 

subroutine MAYERBI(P, S, C, M, W, nbas)

implicit none

integer, intent(in) ::  nbas
integer :: i, j
integer :: d
integer, dimension(nbas), intent(in) :: C
real(8), dimension(nbas,nbas), intent(in) :: P,S
real(8), dimension(nbas,nbas) :: PS
real(8), dimension(:,:), allocatable, intent(out) :: M, W
real(8), parameter :: ALPHA=1.0D0, BETA=0.0D0 
!real(8) :: pop

!d = Nbas

!call DGEMM('N','N',nb,nb,nb,ALPHA,P,nb,S,nb,BETA,PS,nb)

d = MAXVAL(C)

!pop = 0
!do i = 1,nb 
!   pop = pop+PS(i,i)
!enddo

!write(*,*) pop

!do i = 1,nb
!   write(*,'(38F6.2)') (PS(i,j), j = 1,nb)
!end do

ALLOCATE(M(d,d))
ALLOCATE(W(d,d))

M = 0.0D0
W = 0.0D0


do i = 1, nbas 
   do j = 1, nbas
    if (C(i) .ne. C(j)) then
        W(C(i),C(j)) = W(C(i),C(j)) + P(i,j)*P(j,i)
 !       M(C(i),C(j)) = M(C(i),C(j)) + PS(i,j)*PS(j,i)
    endif
   enddo
enddo

end subroutine

subroutine CENTERS(file47, C, nbas)

implicit none

integer, intent(in) ::  nbas
integer :: i, j, ios, l
character(len=100), intent(in) :: file47
character(len=300) :: input_line
integer, allocatable, dimension(:), intent(out) :: C
integer, dimension(11) :: cent10
character(10), parameter :: basind='$BASIS' 
character(10) :: junk

ALLOCATE(C(nbas))

open(unit=10, file=file47, status="old", action='read') ! open file
    do
     read(10, '(a)',iostat=ios) input_line   ! read the lines
     if (ios /=0) exit
     if (trim(input_line) == basind) then  
         i = 0
         l = 0
         do 
         i = i+1
         if (i == 1) then
         read(10, *, iostat=ios) junk, junk, (cent10(j), j = 1,11)
		 else
         read(10, *, iostat=ios) (cent10(j), j = 1,11)
		 endif
		    do j = 1,10
              l = l+1
                if (l<=nbas) then
                   C(l)=cent10(j)
                endif
            enddo
         if (ios /=0) exit
         enddo
     end if
     enddo

close(10)

end subroutine

subroutine READPS(file47, P, S, nbas)

implicit none

integer, intent(in) ::  nbas
integer :: i, j, ios, dtriang, l
character(len=100), intent(in) :: file47
character(len=300) :: input_line
real(8), allocatable, dimension(:,:), intent(out) :: P, S
real(8), allocatable, dimension(:) :: triangp, triangs
real(8), dimension(4) :: triang4
character(10), parameter :: pind=' $DENSITY', sind=' $OVERLAP' 

dtriang=(nbas*nbas-nbas)/2+nbas 

ALLOCATE(S(nbas,nbas))
ALLOCATE(P(nbas,nbas))
ALLOCATE(triangs(dtriang))
ALLOCATE(triangp(dtriang))

S = 0.0D0 
P = 0.0D0
pop = 0.0D0

open(unit=10, file=file47, status="old", action='read') ! open file
    do
     read(10, '(a)',iostat=ios) input_line   ! read the lines
     if (ios /=0) exit
     if (trim(input_line) == sind) then  ! check if it is a comment line or data line
	     write(*,*) 'Overlap is found!'
          i = 0
           l = 0
          do 
           i = i+1
         read(10, *, iostat=ios) (triang4(j), j = 1,4)
           do j = 1,4
              l = l+1 
                if (l<=dtriang) then 
                  triangs(l)=triang4(j)
     !             write(*,*) l, triangs(l)
              endif
           enddo
           if (ios /=0) exit
           enddo
      end if
     if (trim(input_line) == pind) then  ! check if it is a comment line or data line
          i = 0
           l = 0
          do 
           i = i+1
         read(10, *, iostat=ios) (triang4(j), j = 1,4)
           do j = 1,4
              l = l+1 
                if (l<=dtriang) then 
                  triangp(l)=triang4(j)
     !             write(*,*) l, triangp(l)
              endif
           enddo
           if (ios /=0) exit
           enddo 
     end if
     enddo

close(10)

l = 0
      do i = 1, nbas
            do j = 1, i
             l = l + 1
             if (l <= dtriang) then
              S(i,j) = triangs(l) 
              S(j,i) = S(i,j) 
             endif
!             write(*,*) i, j, S(i,j)
             enddo
       enddo


l  = 0      
      do i = 1, nbas
           do j = 1, i
           l = l+1
           if (l <= dtriang) then
            P(i,j) = triangp(l)
            P(j,i) = P(i,j)  
           endif
           if (i==j) then
           pop = pop + P(i,i)
           endif
!           write(*,*) i, j, P(i,j)
           enddo
      enddo

DEALLOCATE(triangs)
DEALLOCATE(triangp)

end subroutine 

! subroutine keyword(line, word, val1, val2, val3)
! implicit none
! integer :: ierr
! integer, parameter :: iu = 301
! character (len=1000), intent(in) :: line
! character (len=15) :: word 
! character (len=15) :: first 
! character  (len=15) :: word1, word2, word3, word4, word5, word5
! integer, intent(out) :: val1, val2, val3 
! 
! 
! open (unit=iu,file=,action="read")
! do
!    read (iu,"(a)",iostat=ierr) line ! read line into character variable
!    if (ierr /= 0) exit
!    read (line,*) first                   ! read first word of line
!    if (trim(word) == '$GENNBO') then     ! found search string at beginning of line
!       read (line,*) word1, val1, word2, val2, word3, val3  
! 	  write (*,*) word1, val1, word2, val2, word3, val3
! 	  exit
!    end if
! end do
! end subroutine


end
