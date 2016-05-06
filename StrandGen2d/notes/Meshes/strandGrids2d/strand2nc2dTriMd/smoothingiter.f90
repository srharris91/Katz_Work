!> \brief
!! This subroutine generates the strand pointing vectors.
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
!!   \include smoothingIter.F90


SUBROUTINE smoothingiter(iLo,iHi,iLoN,iHiN,ig,igN,nNodes,nBNodes,npsp1, &
                         bNode,bNorm,pv,pv0,ang,psp1,psp2,rL2)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: iLo,iHi,iLoN,iHiN,ig,igN,nNodes,nBNodes,npsp1
INTEGER,INTENT(IN   ),DIMENSION(  nBNodes) :: bNode
REAL,   INTENT(IN   ),DIMENSION(2,nBNodes) :: bNorm
REAL,   INTENT(INOUT),DIMENSION(2,iLoN:iHiN+igN) :: pv
REAL,   INTENT(IN   ),DIMENSION(2,iLoN:iHiN+igN) :: pv0
REAL,   INTENT(IN   ),DIMENSION(  iLoN:iHiN+igN) :: ang
INTEGER,INTENT(IN   ),DIMENSION(npsp1) :: psp1
INTEGER,INTENT(IN   ),DIMENSION(  iLoN:iHiN+igN+1) :: psp2
REAL,   INTENT(  OUT) :: rL2
INTEGER :: n,i,m,j
INTEGER,ALLOCATABLE,DIMENSION(:) :: flag
REAL,   ALLOCATABLE,DIMENSION(:,:) :: pv2
REAL :: ax,ay,a,ds,nx,ny


! flag boundary nodes with the edgeTags of the 2 edges which meet at it

ALLOCATE(flag(iLoN:iHiN+igN))
flag = 0
DO n=1,nBNodes
   i       = bNode(n)
   flag(i) = n
END DO


! iteratively smooth pointing vectors

ALLOCATE(pv2(2,iLoN:iHiN+igN))
pv2 = pv
DO n=iLoN,iHiN
   ax = 0.
   ay = 0.
   j  = flag(n)
   DO m=psp2(n)+1,psp2(n+1)
      i  = psp1(m)
      ax = ax+pv(1,i)
      ay = ay+pv(2,i)
   END DO
   IF (j > 0) THEN
      nx = bNorm(1,j)
      ny = bNorm(2,j)
      a  = ax*nx+ay*ny
      ax = ax-a*nx
      ay = ay-a*ny
   END IF
   ds = 1./SQRT(ax*ax+ay*ay)
   ax = ax*ds
   ay = ay*ds
   ds = ax*pv0(1,n)+ay*pv0(2,n)
   IF (ds >= ang(n)) THEN
      pv2(1,n) = ax
      pv2(2,n) = ay
   END IF
END DO


! check for convergence

rL2  = 0.
DO n=iLoN,iHiN
   a   = 1.-pv(1,n)*pv2(1,n)-pv(2,n)*pv2(2,n)
   rL2 = rL2+a*a
END DO
rL2 = SQRT(rL2/REAL(iHiN-iLoN+1))


! set vectors to new smoothed values

pv = pv2

DEALLOCATE(flag,pv2)


END SUBROUTINE smoothingiter
