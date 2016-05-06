!> \brief
!! This subroutine forms the cells surrounding points structures.
!! \param pid Process ID.
!! \param nFaces Number of faces.
!! \param nGfaces Number of ghost faces on the partition.
!! \param nNodes Number of Nodes.
!! \param nGnodes Number of ghost nodes on the partition.
!! \param nEdges Number of edges.
!! \param nBedges Number of boundary edges.
!! \param nPstr Number of cells along a strand.
!! \param face List of faces.
!! \param edge List of edges.
!! \param edgn List of nodes corresponding to each edge.
!! \param x Node coordinates
!! \param facs Face areas in the structured direction.
!! \param facu Face areas in the unstructured direction.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-9-20
!! \par Further Documentation:
!! \par Source code:
!! \include src/Numerics/faceArea.F90


SUBROUTINE facearea( &
     pid, &
     ndim, &
     nFaces, &
     nGfaces, &
     nNodes, &
     nGnodes, &
     nEdges, &
     nBedges, &
     nPstr, &
     face, &
     edge, &
     edgn, &
     x, &
     facs, &
     facu)


IMPLICIT NONE

INTEGER,INTENT(IN   ) :: &
     pid, &
     ndim, &
     nFaces, &
     nGfaces, &
     nNodes, &
     nGnodes, &
     nEdges, &
     nBedges, &
     nPstr
INTEGER,INTENT(INOUT),DIMENSION(2,           nFaces) :: face
INTEGER,INTENT(INOUT),DIMENSION(2,           nEdges) :: edge
INTEGER,INTENT(INOUT),DIMENSION(             nEdges) :: edgn
REAL,   INTENT(INOUT),DIMENSION(ndim,0:nPstr,nNodes) :: x
REAL,   INTENT(INOUT),DIMENSION(ndim,0:nPstr,nEdges) :: facs
REAL,   INTENT(INOUT),DIMENSION(ndim,0:nPstr,nFaces) :: facu
INTEGER :: n,i,j,n1,n2,c1,c2
REAL :: xm,ym,ax,ay,xp1,yp1,xp2,yp2
REAL,                 DIMENSION(ndim,0:nPStr+1,nFaces+nBedges) :: sumA


! make arrays 1-based
face = face+1
edge = edge+1
edgn = edgn+1


! Compute structured face areas
facs = 0.
DO n=1,nEdges
   n1          = edgn(n)
DO j=1,nPstr
   xp1         = x(1,j-1,n1)
   yp1         = x(2,j-1,n1)
   xp2         = x(1,j  ,n1)
   yp2         = x(2,j  ,n1)
   facs(1,j,n) = yp2-yp1
   facs(2,j,n) = xp1-xp2
END DO
END DO


! potentially flip the sign on the boundary edge face areas
j = nEdges-nBedges
DO n=1,nBedges
   j  = j+1
   c1 = edge(1,j)
   c2 = edge(2,j)
   n1 = face(1,c1)
   IF (edgn(j) == n1) facs(:,:,j) =-facs(:,:,j)
END DO

! Compute unstructured face areas
facu = 0.
DO n=1,nFaces
   n1          = face(1,n)
   n2          = face(2,n)
DO j=0,nPStr
   xp1         = x(1,j,n2)
   yp1         = x(2,j,n2)
   xp2         = x(1,j,n1)
   yp2         = x(2,j,n1)
   facu(1,j,n) = yp2-yp1
   facu(2,j,n) = xp1-xp2
END DO
END DO


! Check that face areas sum to 0.
sumA = 0.
DO n=1,nEdges
   c1        = edge(1,n)
   c2        = edge(2,n)
DO j=1,nPStr
   sumA(:,j,c1) = sumA(:,j,c1)+facs(:,j,n)
   sumA(:,j,c2) = sumA(:,j,c2)-facs(:,j,n)
END DO
END DO
DO n=1,nFaces
DO j=0,nPstr
   sumA(:,j  ,n) = sumA(:,j  ,n)+facu(:,j,n)
   sumA(:,j+1,n) = sumA(:,j+1,n)-facu(:,j,n)
END DO
END DO
xm = 0.
ym = 0.
DO n=1,nFaces-nGfaces
DO j=1,nPstr
   xm = MAX(xm,ABS(sumA(1,j,n)))
   ym = MAX(ym,ABS(sumA(2,j,n)))
END DO
END DO

WRITE(*,*)
WRITE(*,'(A19,ES24.14,A1,ES24.14,A1)') &
     'max face area sum (',xm,',',ym,')'
WRITE(*,*)


! make arrays 0-based
face = face-1
edge = edge-1
edgn = edgn-1


END SUBROUTINE facearea
