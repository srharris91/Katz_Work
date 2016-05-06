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


SUBROUTINE findcspdim(nFaces, &
                      nNodes, &
                      nGfaces, &
                      nGnodes, &
                      ncsp1, &
                      face, &
                      csp2)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: nFaces, &
                         nNodes, &
                         nGfaces, &
                         nGnodes
INTEGER,INTENT(  OUT) :: ncsp1
INTEGER,INTENT(IN   ),DIMENSION(2,nFaces) :: face
INTEGER,INTENT(INOUT),DIMENSION(nNodes+1) :: csp2
INTEGER :: n,j,k


csp2 = 0

DO n=1,nFaces
DO k=1,2
   j       = face(k,n)+1
   csp2(j) = csp2(j)+1
END DO
END DO

!do n=2,nNodes-nGnodes+1
!   write(*,*)n,csp2(n)
!end do

IF (MINVAL(csp2(2:nNodes-nGnodes+1)) == 0) THEN
   WRITE(*,*)'*** hanging node detected in grid in edgeExtract ***'
   STOP
END IF

DO n=2,nNodes+1
   csp2(n) = csp2(n)+csp2(n-1)
END DO

ncsp1 = csp2(nNodes+1)


END SUBROUTINE findcspdim
