!> \brief
!! This subroutine fills the entries of the psp1 array.
!!
!! Comments:
!!
!! Versions:
!!   - 1.0 Katz 08/16/2010
!!
!! Additional notes:\par
!!  none
!!
!! Source code:
!!   \include fillPsp.f90


SUBROUTINE fillpsp(iHi,iLo,iLoN,iHiN,ig,igN,ncsp1,npsp1,face, &
                   csp1,csp2,psp1,psp2)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: iHi,iLo,iLoN,iHiN,ig,igN,ncsp1,npsp1
INTEGER,INTENT(IN   ),DIMENSION(2,iLo:iHi+ig) :: face
INTEGER,INTENT(IN   ),DIMENSION(ncsp1) :: csp1
INTEGER,INTENT(INOUT),DIMENSION(npsp1) :: psp1
INTEGER,INTENT(IN   ),DIMENSION(iLoN:iHiN+igN+1) :: csp2
INTEGER,INTENT(INOUT),DIMENSION(iLoN:iHiN+igN+1) :: psp2
INTEGER :: n,m,j,k,l,n1,n2
INTEGER,ALLOCATABLE,DIMENSION(:) :: flag


ALLOCATE(flag(iLoN:iHiN+igN))
flag = 0
DO n=iLoN,iHiN
   flag(n) = n
DO m=csp2(n)+1,csp2(n+1)
   j  = csp1(m)
   DO k=1,2
      l  = face(k,j)
      n1 = flag(l)
      IF (n1 /= n) THEN
         n2       = psp2(n)+1
         psp2(n)  = n2
         psp1(n2) = l
         flag(l)  = n
      END IF
   END DO
END DO
END DO
DO n=iHiN+igN+1,iLoN+1,-1
   psp2(n) = psp2(n-1)
END DO
psp2(iLoN) = 0

DEALLOCATE(flag)

!do n=iLoN,iHiN
!   write(*,'(i12,10i8)')n,(psp1(j),j=psp2(n)+1,psp2(n+1))
!end do


END SUBROUTINE fillpsp
