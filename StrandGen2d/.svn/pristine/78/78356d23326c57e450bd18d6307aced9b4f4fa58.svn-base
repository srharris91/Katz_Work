!> \brief
!! This subroutine finds sharp corner nodes.
!!
!! Comments:
!!
!! Versions:
!!   - 1.0 Workman 09/21/2011
!!   - 1.1 Katz    09/18/2012
!!
!! Additional notes:
!!   
!! Source code:
!!   \include findsharpcorners.f90


SUBROUTINE findsharpcorners(&
     nSurfNodes, &
     nSurfFaces, &
     nSharp, &
     nBndEdges, &
     surfFaces, &
     bndEdges, &
     xSurf, &
     angle, &
     smax, &
     nNewNodes, &
     addNode)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: &
     nSurfNodes, &
     nSurfFaces, &
     nBndEdges, &
     smax
INTEGER,INTENT(  OUT) :: &
     nSharp, &
     nNewNodes
REAL,   INTENT(IN   ) :: angle
INTEGER,INTENT(IN   ),DIMENSION(2,nSurfFaces) :: surfFaces
INTEGER,INTENT(IN   ),DIMENSION(  nBndEdges ) :: bndEdges
REAL,   INTENT(IN   ),DIMENSION(2,nSurfNodes) :: xSurf
INTEGER,INTENT(  OUT),DIMENSION(  nSurfNodes) :: addNode
INTEGER :: n,n1,n2,i
REAL    :: dx,dy,ds,nx,ny,pi,a,b,cp
REAL,DIMENSION(  nSurfNodes) :: flag
REAL,DIMENSION(2,nSurfNodes) :: Rnorm,Lnorm


! find left and right normal vectors at each node
DO n=1,nSurfFaces
   n1          = surfFaces(1,n)
   n2          = surfFaces(2,n)
   dx          = xSurf(1,n2)-xSurf(1,n1)
   dy          = xSurf(2,n2)-xSurf(2,n1)
   ds          = 1./SQRT(dx*dx+dy*dy)
   nx          =-dy*ds
   ny          = dx*ds
   Rnorm(1,n1) = nx
   Rnorm(2,n1) = ny
   Lnorm(1,n2) = nx
   Lnorm(2,n2) = ny
END DO


! mark boundary nodes
flag = 0
DO n=1,nBndEdges
   flag(bndEdges(n)) = 1
END DO


! for sharp convex nodes, determine the number of new nodes to add
IF (angle > 180.) THEN
   WRITE(*,*)'angle must be less than 180. in inputs'
   STOP
END IF
nSharp    = 0
nNewNodes = 0
addNode   = 1
pi        = 4.*ATAN(1.)
a         = angle*pi/180.
DO n=1,nSurfNodes
IF (flag(n) == 0) THEN
   nx         = Rnorm(1,n)
   ny         = Rnorm(2,n)
   dx         = Lnorm(1,n)
   dy         = Lnorm(2,n)
   ds         = nx*dx+ny*dy
   ds         = ACOS(ds)
   cp         = nx*dy-ny*dx
IF (ds > a .AND. cp > 0.) THEN
   IF (smax < 0) THEN
      i       = abs(smax)
   ELSE
      b       =(ds-a)/(pi-a)*REAL(smax)
      i       = MAX(CEILING(b),2) ! add at least 3 nodes at a sharp corner
   END IF
   nSharp     = nSharp+1
   nNewNodes  = nNewNodes+i
   addNode(n) = 1+i
END IF
END IF
END DO


END SUBROUTINE findsharpcorners
