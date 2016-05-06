!> \brief
!! This subroutine computes cell-center coordinates.
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
!! \param fClip Facewise clipping index.
!! \param x Node coordinates
!! \param facs Face areas in the structured direction.
!! \param facu Face areas in the unstructured direction.
!! \param xc Cell-center coordinates.
!! \param v Cell volumes.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-9-20
!! \par Further Documentation:
!! \par Source code:
!! \include src/Numerics/volume.F90


SUBROUTINE volume( &
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
     fClip, &
     x, &
     facs, &
     facu, &
     xc, &
     v)


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
INTEGER,INTENT(INOUT),DIMENSION(2,             nFaces        ) :: face
INTEGER,INTENT(INOUT),DIMENSION(2,             nEdges        ) :: edge
INTEGER,INTENT(INOUT),DIMENSION(               nEdges        ) :: edgn
INTEGER,INTENT(INOUT),DIMENSION(               nFaces+nBedges) :: fClip
REAL,   INTENT(INOUT),DIMENSION(ndim,0:nPstr  ,nNodes        ) :: x
REAL,   INTENT(INOUT),DIMENSION(ndim,0:nPstr  ,nEdges        ) :: facs
REAL,   INTENT(INOUT),DIMENSION(ndim,0:nPstr  ,nFaces        ) :: facu
REAL,   INTENT(INOUT),DIMENSION(ndim,0:nPstr+1,nFaces+nBedges) :: xc
REAL,   INTENT(INOUT),DIMENSION(     0:nPstr+1,nFaces+nBedges) :: v
INTEGER :: n,i,j,n1,n2,c1,c2,m
REAL :: lx,xa,sum,xp1,yp1,xp2,yp2,xct,yct,vt


! make arrays 1-based
face = face+1
edge = edge+1
edgn = edgn+1


! compute cell volumes
v = 0.
DO n=1,nEdges
   n1      = edgn(n)
   c1      = edge(1,n)
   c2      = edge(2,n)
DO j=1,nPStr
   xp1     = x(1,j-1,n1)
   xp2     = x(1,j  ,n1)
   lx      = facs(1,j,n)
   xa      = .5*(xp1+xp2)
   v(j,c1) = v(j,c1)+xa*lx
   v(j,c2) = v(j,c2)-xa*lx
END DO
END DO

DO n=1,nFaces
   n1       = face(1,n)
   n2       = face(2,n)
DO j=0,nPStr
   lx       = facu(1,j,n)
   xa       = .5*(x(1,j,n1)+x(1,j,n2))
   v(j  ,n) = v(j  ,n)+xa*lx
   v(j+1,n) = v(j+1,n)-xa*lx
END DO
END DO


! Compute total volume for testing purposes
sum = 0.
DO n=1,nFaces-nGfaces
DO j=1,fclip(n)
   sum = sum+v(j,n)
END DO
END DO
WRITE(*,*)
WRITE(*,*)'Fine level total volume: ',sum
WRITE(*,*)


! Compute boundary value "volumes"
DO n=1,nFaces
   v(0,n)       = v(1,n)
   v(nPstr+1,n) = v(nPstr,n)
END DO

m = nEdges-nBedges
DO n=1,nBedges
   m       = m+1
   c1      = edge(1,m)
   c2      = edge(2,m)
DO j=1,nPStr
   v(j,c2) = v(j,c1)
END DO
   v(0,c2) = v(1,c2)
   v(nPstr+1,c2) = v(nPstr,c2)
END DO


! compute centroid location for the entire grid for testing purposes
xct = 0.
yct = 0.
vt  = 0.
DO n=1,nFaces-nGfaces
DO j=1,nPstr
   xct = xct+xc(1,j,n)*v(j,n)
   yct = yct+xc(2,j,n)*v(j,n)
   vt  = vt+v(j,n)
END DO
END DO
xct = xct/vt
yct = yct/vt

WRITE(*,*)
WRITE(*,'(A15,ES24.14,A1,ES24.14,A1)') &
     'grid centroid (',xct,',',yct,')'
WRITE(*,*)


! make arrays 0-based
face = face-1
edge = edge-1
edgn = edgn-1


END SUBROUTINE volume
