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
!! \param x Node coordinates
!! \param facs Face areas in the structured direction.
!! \param facu Face areas in the unstructured direction.
!! \param xc Cell-center coordinates
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-9-20
!! \par Further Documentation:
!! \par Source code:
!! \include src/Numerics/cellCoord.F90


SUBROUTINE cellcoord( &
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
     facu, &
     xc)


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
REAL,   INTENT(INOUT),DIMENSION(ndim,0:nPstr  ,nNodes        ) :: x
REAL,   INTENT(INOUT),DIMENSION(ndim,0:nPstr  ,nEdges        ) :: facs
REAL,   INTENT(INOUT),DIMENSION(ndim,0:nPstr  ,nFaces        ) :: facu
REAL,   INTENT(INOUT),DIMENSION(ndim,0:nPstr+1,nFaces+nBedges) :: xc
INTEGER :: n,n1,n2,n3,n4,j,c1,c2,m,i
REAL :: third,sixth,lx,ly,lz,xp1,yp1,xp2,yp2,xp3,yp3,dxa,dya,dxb,dyb, &
     v1,v2,pi,a,ds,dx,dy


! make arrays 1-based
face = face+1
edge = edge+1
edgn = edgn+1

third = 1./3.


DO n=1,nFaces
   n1 = face(1,n)
   n2 = face(2,n)
DO j=1,nPStr
   xp1 = x(1,j-1,n1)
   yp1 = x(2,j-1,n1)
   xp2 = x(1,j-1,n2)
   yp2 = x(2,j-1,n2)
   xp3 = x(1,j  ,n2)
   yp3 = x(2,j  ,n2)
   dxa = xp2-xp1
   dya = yp2-yp1
   dxb = xp3-xp1
   dyb = yp3-yp1
   v1  = .5*(dxa*dyb-dya*dxb)
   xc(1,j,n) = third*(xp1+xp2+xp3)*v1
   xc(2,j,n) = third*(yp1+yp2+yp3)*v1
   xp1 = x(1,j-1,n1)
   yp1 = x(2,j-1,n1)
   xp2 = x(1,j  ,n2)
   yp2 = x(2,j  ,n2)
   xp3 = x(1,j  ,n1)
   yp3 = x(2,j  ,n1)
   dxa = xp2-xp1
   dya = yp2-yp1
   dxb = xp3-xp1
   dyb = yp3-yp1
   v2  = .5*(dxa*dyb-dya*dxb)
   xc(1,j,n) = xc(1,j,n)+third*(xp1+xp2+xp3)*v2
   xc(2,j,n) = xc(2,j,n)+third*(yp1+yp2+yp3)*v2
   xc(1,j,n) = xc(1,j,n)/(v1+v2)
   xc(2,j,n) = xc(2,j,n)/(v1+v2)
END DO
END DO


! set surface and boundary locations to be the quadrature points
DO n=1,nFaces
   n1              = face(1,n)
   n2              = face(2,n)
   xc(:,0,n)       = .5*(x(:,0,    n1)+x(:,0,    n2))
   xc(:,nPstr+1,n) = .5*(x(:,nPstr,n1)+x(:,nPstr,n2))
END DO
DO n=nEdges-nBedges+1,nEdges
   n1         = edgn(n)
   c2         = edge(2,n)
DO j=1,nPStr
   xc(:,j,c2) = .5*(x(:,j-1,n1)+x(:,j,n1))
END DO
   xc(:,0      ,c2) = x(:,0    ,n1)
   xc(:,nPstr+1,c2) = x(:,nPstr,n1)
END DO


! make arrays 0-based
face = face-1
edge = edge-1
edgn = edgn-1


END SUBROUTINE cellcoord
