!> \brief
!! This subroutine performs stage 1 of the surface mesh partitioning.
!!
!! Comments:
!!
!! Versions:
!!   - 1.0 Katz 08/03/2010
!!
!! Additional notes:\par
!!  none
!!
!! Source code:
!!   \include partitionstage1.f90


SUBROUTINE partitionstage1(nStrandBlocks, &
                           nSurfFacesD, &
                           nSurfNodesD, &
                           nBndEdgesD, &
                           nPrtEdgesD, &
                           globalFaceMap, &
                           globalNodeMap, &
                           globalEdgeMap, &
                           prtEdges, &
                           nPrtEdges, &
                           ncsp1, &
                           csp1, &
                           csp2, &
                           nPtsPerStrand, &
                           nFringeD, &
                           surfaceOnly, &
                           nSurfNodes, &
                           nSurfFaces, &
                           nBndEdges, &
                           nSurfPatches, &
                           nEdgePatches, &
                           nSharpG,&
                           surfFaces, &
                           bndEdges, &
                           nodeClip, &
                           faceClip, &
                           faceTag, &
                           bndTag, &
                           sFlag, &
                           xSurf, &
                           pointingVec, &
                           xStrand, &
                           bndNormal, &
                           id, &
                           pid, &
                           nFaces, &
                           nNodes, &
                           nGfaces, &
                           nGnodes, &
                           nBedges, &
                           nPedges, &
                           nPstr, &
                           nFringe, &
                           nScomp, &
                           nBcomp, &
                           nSharp)
              
IMPLICIT NONE


INTEGER,INTENT(IN   ) :: nStrandBlocks, &
                         nSurfFacesD, &
                         nSurfNodesD, &
                         nBndEdgesD, &
                         nPrtEdgesD, &
                         ncsp1, &
                         nPtsPerStrand, &
                         nFringeD, &
                         surfaceOnly, &
                         nSurfNodes, &
                         nSurfFaces, &
                         nBndEdges, &
                         nSurfPatches, &
                         nEdgePatches, &
                         nSharpG
INTEGER,INTENT(INOUT) :: nPrtEdges
INTEGER,INTENT(INOUT),DIMENSION(  nSurfFacesD) :: globalFaceMap
INTEGER,INTENT(INOUT),DIMENSION(  nSurfNodesD) :: globalNodeMap
INTEGER,INTENT(INOUT),DIMENSION(2,nPrtEdgesD ) :: globalEdgeMap
INTEGER,INTENT(INOUT),DIMENSION(2,nPrtEdgesD ) :: prtEdges
INTEGER,INTENT(INOUT),DIMENSION(ncsp1        ) :: csp1
INTEGER,INTENT(INOUT),DIMENSION(nSurfNodes+1 ) :: csp2
INTEGER,INTENT(IN   ),DIMENSION(2,nSurfFaces) :: surfFaces
INTEGER,INTENT(IN   ),DIMENSION(  nBndEdges ) :: bndEdges
INTEGER,INTENT(IN   ),DIMENSION(  nSurfNodes) :: nodeClip
INTEGER,INTENT(IN   ),DIMENSION(  nSurfFaces) :: faceClip
INTEGER,INTENT(IN   ),DIMENSION(  nSurfFaces) :: faceTag
INTEGER,INTENT(IN   ),DIMENSION(  nBndEdges ) :: bndTag
INTEGER,INTENT(IN   ),DIMENSION(  nSurfNodes) :: sFlag
REAL   ,INTENT(IN   ),DIMENSION(2,nSurfNodes) :: xSurf
REAL   ,INTENT(IN   ),DIMENSION(2,nSurfNodes) :: pointingVec
REAL   ,INTENT(IN   ),DIMENSION(0:nPtsPerStrand) :: xStrand
REAL   ,INTENT(IN   ),DIMENSION(2,nEdgePatches) :: bndNormal
INTEGER,INTENT(INOUT),DIMENSION(nStrandBlocks) :: id
INTEGER,INTENT(INOUT),DIMENSION(nStrandBlocks) :: pid
INTEGER,INTENT(INOUT),DIMENSION(nStrandBlocks) :: nFaces
INTEGER,INTENT(INOUT),DIMENSION(nStrandBlocks) :: nNodes
INTEGER,INTENT(INOUT),DIMENSION(nStrandBlocks) :: nGfaces
INTEGER,INTENT(INOUT),DIMENSION(nStrandBlocks) :: nGnodes
INTEGER,INTENT(INOUT),DIMENSION(nStrandBlocks) :: nBedges
INTEGER,INTENT(INOUT),DIMENSION(nStrandBlocks) :: nPedges
INTEGER,INTENT(INOUT),DIMENSION(nStrandBlocks) :: nPstr
INTEGER,INTENT(INOUT),DIMENSION(nStrandBlocks) :: nFringe
INTEGER,INTENT(INOUT),DIMENSION(nStrandBlocks) :: nScomp
INTEGER,INTENT(INOUT),DIMENSION(nStrandBlocks) :: nBcomp
INTEGER,INTENT(INOUT),DIMENSION(nStrandBlocks) :: nSharp

INTEGER :: ihim,dFaces,dFacesE,dFacesm,i,j,k,l,m,n,n1,n2
INTEGER,ALLOCATABLE,DIMENSION(:) :: nflag,nflag2,cflag,bflag,ilo,ihi
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: edge


k = (2*nSurfFaces+nBndEdges)/2
ALLOCATE(nflag(nSurfNodes), &
     nflag2(nSurfNodes), &
     cflag(nSurfFaces), &
     bflag(nSurfNodes), &
     ilo(nStrandBlocks), &
     ihi(nStrandBlocks), &
     edge(2,k))

! flag the boundary nodes with the surface patch number
bflag = 0
DO i=1,nBndEdges
   k        = bndEdges(i)
   bflag(k) = bndTag(i)
END DO

! fill in easy dimensions
ihim    = 0
dFaces  =(nSurfFaces)/nStrandBlocks
dFacesE = nSurfFaces-nStrandBlocks*dFaces+dFaces
DO n=1,nStrandBlocks
   dFacesm    = dFaces
   IF (n==nStrandBlocks) dFacesm = dFacesE
   id(n)      = n-1
   pid(n)     = n-1
   ilo(n)     = ihim+1
   ihi(n)     = ilo(n)+dFacesm-1
   ihim       = ihi(n)
   nPstr(n)   = nPtsPerStrand
   nFringe(n) = nFringeD
   nBcomp(n)  = nEdgePatches
   nScomp(n)  = nSurfPatches
END DO


! this assumes 1 strand block!
nSharp(1) = nSharpG


! form cells surrounding points structures
csp2 = 0
DO n=1,nSurfFaces
DO k=1,2
   j       = surfFaces(k,n)+1
   csp2(j) = csp2(j)+1
END DO
END DO
DO n=2,nSurfNodes+1
   csp2(n) = csp2(n)+csp2(n-1)
END DO
IF (csp2(nSurfNodes+1) > ncsp1) THEN
   WRITE(*,*)'*** insufficient dimension ncsp1 in partitionstage1.f90 ***'
   STOP
END IF
DO n=1,nSurfFaces
DO k=1,2
   j       = surfFaces(k,n)
   i       = csp2(j)+1
   csp2(j) = i
   csp1(i) = n
END DO
END DO
DO n=nSurfNodes+1,2,-1
   csp2(n) = csp2(n-1)
END DO
csp2(1) = 0


! fill in global face map array
DO n=1,nStrandBlocks
DO i=ilo(n),ihi(n)
   globalFaceMap(i) = id(n)
END DO
END DO

! find the number of partitions edges on each partition
edge = 0
DO n=1,nSurfNodes
DO m=csp2(n)+1,csp2(n+1)
   j = csp1(m)
   IF (surfFaces(1,j) == n) edge(2,n) = j
   IF (surfFaces(2,j) == n) edge(1,n) = j
END DO
END DO

nPedges = 0
k       = 0
DO n=1,nSurfNodes
   i  = edge(1,n)
   j  = edge(2,n)
IF (MIN(i,j) > 0) THEN
   n1 = globalFaceMap(i)
   n2 = globalFaceMap(j)
IF (n1 /= n2) THEN
   nPedges(n1)   = nPedges(n1)+1
   nPedges(n2)   = nPedges(n2)+1
   k             = k+1
   prtEdges(:,k) = edge(:,n)
   globalEdgeMap(1,k) = n1
   globalEdgeMap(2,k) = n2
END IF
END IF
END DO
nPrtEdges = k
IF (nPrtEdges > nPrtEdgesD) THEN
   WRITE(*,*)'*** nPrtEdges counted incorrectly in partitionstage1.f90 ***'
   STOP
END IF

! for each partition...
ihim = 0
DO n=1,nStrandBlocks

   ! find the number of nodes
   nflag = 0
   DO i=ilo(n),ihi(n)
   DO k=1,2
      nflag(surfFaces(k,i)) = 1
   END DO
   END DO
   k = 0
   DO i=1,nSurfNodes
   IF (nflag(i) == 1) THEN
      k                = k+1
      globalNodeMap(i) = id(n) !there will be nodes on more than one id
   END IF
   END DO
   nNodes(n) = k

   ! find the number of ghost cells
   cflag = 0
   DO i=1,nSurfNodes
   IF (nflag(i) == 1) THEN
      DO m=csp2(i)+1,csp2(i+1)
         j        = csp1(m)
         cflag(j) = 1
      END DO
   END IF
   END DO
   DO i=ilo(n),ihi(n)
      cflag(i) = 0
   END DO
   k = 0
   DO i=1,nSurfFaces
      k = k+cflag(i)
   END DO
   nGfaces(n) = k
   nFaces(n)  = ihi(n)-ilo(n)+1+nGfaces(n)

   ! find the number of ghost nodes
   nflag2 = 0
   DO i=1,nSurfFaces
   IF (cflag(i) == 1) THEN
      DO k=1,2
         nflag2(surfFaces(k,i)) = 1
      END DO
   END IF
   END DO
   DO i=1,nSurfNodes
      nflag(i) = MAX(0,nflag2(i)-nflag(i))
   END DO
   k = 0
   DO i=1,nSurfNodes
      k = k+nflag(i)
   END DO
   nGnodes(n) = k
   nNodes(n)  = nNodes(n)+nGnodes(n)

   ! find the number of boundary edges
   k = 0
   DO i=ilo(n),ihi(n)
      n1 = surfFaces(1,i)
      n2 = surfFaces(2,i)
      IF (bflag(n1) > 0) k = k+1
      IF (bflag(n2) > 0) k = k+1
   END DO
   nBedges(n) = k

END DO

DEALLOCATE(cflag, &
     nflag, &
     nflag2, &
     bflag, &
     ilo, &
     ihi, &
     edge)

! check to make sure partition numbers sum correctly
k = 0
DO n=1,nStrandBlocks
   k = k+nFaces(n)-nGfaces(n)
END DO
IF (k /= nSurfFaces) THEN
   WRITE(*,*)'*** faces partitioned incorrectly in partitionstage1.f90 ***'
   STOP
END IF
k = 0
DO n=1,nStrandBlocks
   k = k+nBedges(n)
END DO
IF (k /= nBndEdges) THEN
   WRITE(*,*)'*** bnd. edges partitioned incorrectly in partitionstage1.f90 ***'
   STOP
END IF

!do n=1,nstrandblocks
!   write(*,*)n
!   write(*,*)id(n)
!   write(*,*)pid(n)
!   write(*,*)nFaces(n)
!   write(*,*)nNodes(n)
!   write(*,*)nGfaces(n)
!   write(*,*)nGnodes(n)
!   write(*,*)nBedges(n)
!   write(*,*)nPedges(n)
!   write(*,*)nPstr(n)
!   write(*,*)nFringe(n)
!   write(*,*)nScomp(n)
!   write(*,*)nBcomp(n)
!   write(*,*)nSharp(n)
!   write(*,*)
!end do


END SUBROUTINE partitionstage1
