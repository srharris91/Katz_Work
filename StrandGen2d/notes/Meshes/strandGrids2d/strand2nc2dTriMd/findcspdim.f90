!> \brief
!! This subroutine finds the dimension of the csp1 array.
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


SUBROUTINE findcspdim(iHi,iLo,iLoN,iHiN,ig1,ig2,ig,igN,csp2,face,ncsp1)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: iHi,iLo,iLoN,iHiN,ig1,ig2,ig,igN
INTEGER,INTENT(  OUT) :: ncsp1
INTEGER,INTENT(IN   ),DIMENSION(2,iLo:iHi+ig) :: face
INTEGER,INTENT(INOUT),DIMENSION(iLoN:iHiN+igN+1) :: csp2
INTEGER :: n,j,k


csp2 = 0
DO n=iLo,iHi+ig1+ig2
DO k=1,2
   j       = face(k,n)+1
   csp2(j) = csp2(j)+1
END DO
END DO
IF (MINVAL(csp2(iLoN+1:iHiN+1)) == 0) THEN
   WRITE(*,*)'*** hanging node detected in grid in edgeExtract ***'
   STOP
END IF

DO n=iLoN+1,iHiN+igN+1
   csp2(n) = csp2(n)+csp2(n-1)
END DO

ncsp1 = csp2(iHiN+igN+1)


END SUBROUTINE findcspdim
