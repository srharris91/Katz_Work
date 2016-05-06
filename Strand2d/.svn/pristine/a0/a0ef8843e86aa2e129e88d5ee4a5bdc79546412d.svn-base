!> \brief
!! This subroutine sets the edge data structure for fine MG levels.
!! \param pid Process ID.
!! \param nFaces Number of faces.
!! \param nGfaces Number of ghost faces on the partition.
!! \param nNodes Number of nodes.
!! \param nGnodes Number of ghost nodes on the partition.
!! \param nEdges Number of edges, including boundary and partition edges.
!! \param nBedges Number of boundary edges.
!! \param nPedges Number of partition edges.
!! \param face List of faces.
!! \param bEdge List of boundary edges.
!! \param pEdge List of partition edges.
!! \param bTag Tag for each boundary edge.
!! \param edge List of edges (interior first, then partition, then boundary).
!! \param edgp Edge extension cells for limiting purposes.
!! \param edgn The node to which each edge corresponds.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-9-20
!! \par Further Documentation:
!! \par Source code:
!! \include src/Numerics/edgeExtract.F90


SUBROUTINE edgeextract( &
     pid, &
     nFaces, &
     nGfaces, &
     nNodes, &
     nGnodes, &
     nEdges, &
     nBedges, &
     nPedges, &
     face, &
     bEdge, &
     pEdge, &
     bTag, &
     edge, &
     edgp, &
     edgn);


IMPLICIT NONE

INTEGER,INTENT(IN   ) :: &
     pid, &
     nFaces, &
     nGfaces, &
     nNodes, &
     nGnodes, &
     nEdges, &
     nBedges, &
     nPedges
INTEGER,INTENT(INOUT),DIMENSION(2,nFaces ) :: face
INTEGER,INTENT(IN   ),DIMENSION(  nBedges) :: bEdge
INTEGER,INTENT(IN   ),DIMENSION(  nPedges) :: pEdge
INTEGER,INTENT(INOUT),DIMENSION(  nBedges) :: bTag
INTEGER,INTENT(  OUT),DIMENSION(2,nEdges ) :: edge
INTEGER,INTENT(  OUT),DIMENSION(2,nEdges ) :: edgp
INTEGER,INTENT(  OUT),DIMENSION(  nEdges ) :: edgn

INTEGER :: i,j,k,l,m,n,mm,kk,ne,n1,k1,ncsp1,c1,c2,nesf1
REAL    :: xa,ya,xb,yb,xe,ye,xd,yd,dxa,dya,dxb,dyb,cpr,third
INTEGER,ALLOCATABLE,DIMENSION(:) :: bflag,nflag,csp1,csp2,esf1,esf2


! make connectivity 1-based
face = face+1
bTag = bTag+1

ALLOCATE(&
     bflag(nBedges  ), &
     csp2 (nNodes +1))


! Form cells surrounding points (csp) arrays
csp2 = 0
DO n=1,nFaces
DO k=1,2
   j       = face(k,n)+1
   csp2(j) = csp2(j)+1
END DO
END DO

IF (MINVAL(csp2(2:nNodes+1)) == 0) THEN !test for hanging nodes
   WRITE(*,*)'*** hanging node detected in grid in edgeExtract ***'
   STOP
END IF
DO n=2,nNodes+1 !add csp2 to previous values
   csp2(n) = csp2(n)+csp2(n-1)
END DO
ncsp1 = csp2(nNodes+1)
ALLOCATE(csp1(ncsp1))

DO n=1,nFaces
DO k=1,2
   j        = face(k,n)
   i        = csp2(j)+1
   csp2(j) = i
   csp1(i) = n
END DO
END DO
DO n=nNodes+1,2,-1
   csp2(n) = csp2(n-1)
END DO
csp2(1) = 0

!DO n=1,nNodes
!   write(*,*)n,(csp1(m),m=csp2(n)+1,csp2(n+1))
!END DO


! flag boundary nodes with bTag. If a partition boundary, set to -1.
! IF ghost node, set to -2.
ALLOCATE(nflag(nNodes))
nflag = 0
DO n=1,nBedges
   i        = bEdge(n)
   nflag(i) = bTag(n)
END DO
DO n=1,nPedges
   i        = bEdge(n)
   nflag(i) =-1
END DO
DO n=nNodes-nGnodes+1,nNodes
   nflag(n) =-2
END DO

!do n=1,nnodes
!   write(*,*)n,nflag(n)
!end do


! extract interior edges
edge = 0
edgn = 0
ne   = 0
DO n=1,nNodes-nGnodes
IF (nflag(n) == 0) THEN !interior edge
   ne       = ne+1
   edgn(ne) = n
   DO m=csp2(n)+1,csp2(n+1)
      j = csp1(m)
      DO k=1,2
      IF (n == face(k,j)) THEN
         k1 = k
         EXIT
      END IF
      END DO
   n1 = 2
   IF (k1 == 2) n1 = 1
   edge(n1,ne) = j
   END DO
END IF
END DO


! extract partition edges
DO n=1,nNodes-nGnodes
IF (nflag(n) == -1) THEN !partition edge
   ne       = ne+1
   edgn(ne) = n
   DO m=csp2(n)+1,csp2(n+1)
      j = csp1(m)
      IF (j <= nFaces) THEN ! put interior cell first, and boundary dof second
         edge(1,ne) = j
      ELSE
         edge(2,ne) = j
      END IF
   END DO
END IF
END DO


! extract boundary edges, adding boundary dofs
k = nFaces
DO n=1,nNodes-nGnodes
IF (nflag(n) > 0) THEN !boundary edge
   ne       = ne+1
   edgn(ne) = n
   DO m=csp2(n)+1,csp2(n+1)
      j = csp1(m)
      IF (j <= nFaces) THEN ! put interior cell first, and boundary dof second
         edge(1,ne) = j
      ELSE
         edge(2,ne) = j
      END IF
   END DO
   k          = k+1
   edge(2,ne) = k
END IF
END DO

IF (ne /= nEdges) THEN
   WRITE(*,*)'*** edge dimension miscounted in edgeExtract ***'
   STOP
END IF
IF (k  /= nFaces+nBedges) THEN
   WRITE(*,*)'*** DOFs miscounted in edgeExtract ***'
   STOP
END IF

!do n=1,nEdges
!   write(*,*)nedges,n,edge(:,n),edgn(n)
!end do


! form edges surrounding faces
ALLOCATE(esf2(nFaces+nBedges+1))
esf2 = 0
DO n=1,nEdges
DO k=1,2
   j       = edge(k,n)+1
   esf2(j) = esf2(j)+1
END DO
END DO
DO n=2,nFaces+nBedges+1
   esf2(n) = esf2(n)+esf2(n-1)
END DO
nesf1 = esf2(nFaces+nBedges+1)
ALLOCATE(esf1(nesf1))
DO n=1,nEdges
DO k=1,2
   j       = edge(k,n)
   i       = esf2(j)+1
   esf2(j) = i
   esf1(i) = n
END DO
END DO
DO n=nFaces+nBedges+1,2,-1
   esf2(n) = esf2(n-1)
END DO
esf2(1) = 0

!do n=1,nfaces+nbedges
!   write(*,*)n,(esf1(j),j=esf2(n)+1,esf2(n+1))
!end do


! find extension cells (note: I do not find extension cells for boundary
! edges; I'll set lim = 1 here)
edgp = 1
DO n=1,nEdges-nBedges
DO k=1,2
   c1 = edge(k,n)
   oppEdge:DO j=esf2(c1)+1,esf2(c1+1)
      i = esf1(j)
      IF (i /= n) EXIT oppEdge
   END DO oppEdge
   c2 = edge(1,i)
   IF (c2 == c1) c2 = edge(2,i)
   edgp(k,n) = c2
END DO
END DO

!do n=1,nedges
!   write(*,*)n,edge(:,n),edgp(:,n)
!end do


DEALLOCATE( &
     nflag, &
     bflag, &
     csp1, &
     csp2, &
     esf1, &
     esf2)


! make connectivity 0-based
face = face-1
bTag = bTag-1
edge = edge-1
edgp = edgp-1
edgn = edgn-1


END SUBROUTINE edgeextract
