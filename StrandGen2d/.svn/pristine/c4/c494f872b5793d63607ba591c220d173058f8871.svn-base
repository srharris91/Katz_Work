!> \brief
!! This subroutine determines the strand clipping index.
!!
!! Comments:
!!
!! Versions:
!!   - 1.0 Katz 07/23/2010
!!
!! Additional notes:\par
!!   none
!!
!! Source code:
!!   \include strandClipSimple.F90


SUBROUTINE strandclipsimple(nFaces, &
		            nNodes, &
                            nGfaces, &
                            nGnodes, &
                            nPstr, &
                            face, &
                            fClip, &
                            nClip)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: nFaces, &
		         nNodes, &
                         nGfaces, &
                         nGnodes, &
                         nPstr
INTEGER,INTENT(IN   ),DIMENSION(2,nFaces) :: face
INTEGER,INTENT(  OUT),DIMENSION(nNodes) :: nClip
INTEGER,INTENT(  OUT),DIMENSION(nFaces) :: fClip
INTEGER :: n,i,j,k


WRITE(*,*)
WRITE(*,*)'*** only simple clipping in strandClip is implemented ***'
WRITE(*,*)

fClip(1:nFaces-nGfaces) = nPStr

nClip = 0
DO n=1,nFaces
DO k=1,2
   i = face(k,n)
   IF (i > 0) THEN
      j = fClip(n)
      nClip(i) = MAX(nClip(i),j)
   END IF
END DO
END DO


END SUBROUTINE strandclipsimple
