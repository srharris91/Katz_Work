!> \brief
!! This subroutine finds the dimension of the psp1 array.
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
!!   \include findCspDim.f90


SUBROUTINE findpspdim(iHi,iLo,iLoN,iHiN,ig,igN,ncsp1,psp2,face, &
                      csp1,csp2,npsp1)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: iHi,iLo,iLoN,iHiN,ig,igN,ncsp1
INTEGER,INTENT(  OUT) :: npsp1
INTEGER,INTENT(IN   ),DIMENSION(2,iLo:iHi+ig) :: face
INTEGER,INTENT(IN   ),DIMENSION(ncsp1) :: csp1
INTEGER,INTENT(IN   ),DIMENSION(iLoN:iHiN+igN+1) :: csp2
INTEGER,INTENT(INOUT),DIMENSION(iLoN:iHiN+igN+1) :: psp2
INTEGER :: n,m,j,k,l,n1
INTEGER,ALLOCATABLE,DIMENSION(:) :: flag


ALLOCATE(flag(iLoN:iHiN+igN))
psp2 = 0
flag = 0
DO n=iLoN,iHiN
   flag(n) = n
DO m=csp2(n)+1,csp2(n+1)
   j  = csp1(m)
   DO k=1,2
      l  = face(k,j)
      n1 = flag(l)
      IF (n1 /= n) THEN
         psp2(n+1) = psp2(n+1)+1
         flag(l)   = n
      END IF
   END DO
END DO
END DO

DO n=iLoN+1,iHiN+igN+1
   psp2(n) = psp2(n)+psp2(n-1)
END DO

npsp1 = psp2(iHiN+igN+1)

DEALLOCATE(flag)


END SUBROUTINE findpspdim
