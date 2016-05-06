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


SUBROUTINE findpspdim(nFaces, &
                      nNodes, &
                      nGfaces, &
                      nGnodes, &
                      ncsp1, &
                      npsp1, &
                      face, &
                      csp1, &
                      csp2, &
                      psp2)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: nFaces, &
                         nNodes, &
                         nGfaces, &
                         nGnodes, &
                         ncsp1
INTEGER,INTENT(  OUT) :: npsp1
INTEGER,INTENT(IN   ),DIMENSION(2,nFaces) :: face
INTEGER,INTENT(IN   ),DIMENSION(ncsp1) :: csp1
INTEGER,INTENT(IN   ),DIMENSION(nNodes+1) :: csp2
INTEGER,INTENT(  OUT),DIMENSION(nNodes+1) :: psp2
INTEGER :: n,m,j,k,l,n1
INTEGER,DIMENSION(nNodes) :: flag


psp2 = 0
flag = 0
DO n=1,nNodes-nGnodes
   flag(n) = n
DO m=csp2(n)+1,csp2(n+1)
   j  = csp1(m)
   DO k=1,2
      l  = face(k,j)
      IF (l > 0) THEN
         n1 = flag(l)
         IF (n1 /= n) THEN
            psp2(n+1) = psp2(n+1)+1
            flag(l)   = n
         END IF
      END IF
   END DO
END DO
END DO

DO n=2,nNodes+1
   psp2(n) = psp2(n)+psp2(n-1)
END DO

npsp1 = psp2(nNodes+1)


END SUBROUTINE findpspdim
