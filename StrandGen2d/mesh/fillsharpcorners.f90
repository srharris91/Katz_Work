!> \brief
!! Refills global mesh arrays to account for new sharp corners
!!
!! Comments:
!!
!! Versions:
!!   - 1.0 Workman 10/02/2011
!!   - 1.1 Katz    09/18/2012
!!
!! Additional notes:
!!   
!! Source code:
!!   \include fillsharpcorners.f90


SUBROUTINE fillsharpcorners(&
     nSurfNodes, &
     nSurfFaces, &
     nEdgePatches, &
     nSharp, &
     nBndEdges, &
     surfFaces, &
     bndEdges, &
     faceTag, &
     bndTag, &
     nodeClip, &
     xSurf, &
     bndNormal, &
     nNewNodes, &
     addNode, &
     nSurfNodesT, &
     nSurfFacesT, &
     surfFacesT, &
     nodeClipT, &
     faceClipT, &
     faceTagT, &
     sFlagT, &
     xSurfT, &
     pointingVecT)
		
IMPLICIT NONE

INTEGER,INTENT(IN   ) :: &
     nSurfNodes, &
     nSurfFaces, &
     nEdgePatches, &
     nSharp, &
     nBndEdges, &
     nNewNodes, &
     nSurfNodesT, &
     nSurfFacesT
INTEGER,INTENT(IN   ),DIMENSION(2,nSurfFaces  ) :: surfFaces
INTEGER,INTENT(INOUT),DIMENSION(  nBndEdges   ) :: bndEdges
INTEGER,INTENT(IN   ),DIMENSION(  nSurfFaces  ) :: faceTag
INTEGER,INTENT(IN   ),DIMENSION(  nBndEdges   ) :: bndTag
INTEGER,INTENT(IN   ),DIMENSION(  nSurfNodes  ) :: nodeClip
REAL,   INTENT(IN   ),DIMENSION(2,nSurfNodes  ) :: xSurf
REAL,   INTENT(IN   ),DIMENSION(2,nEdgePatches) :: bndNormal
INTEGER,INTENT(IN   ),DIMENSION(  nSurfNodes  ) :: addNode
INTEGER,INTENT(  OUT),DIMENSION(2,nSurfFacesT ) :: surfFacesT
INTEGER,INTENT(  OUT),DIMENSION(  nSurfNodesT ) :: nodeClipT
INTEGER,INTENT(  OUT),DIMENSION(  nSurfFacesT ) :: faceClipT
INTEGER,INTENT(  OUT),DIMENSION(  nSurfFacesT ) :: faceTagT
INTEGER,INTENT(  OUT),DIMENSION(  nSurfNodesT ) :: sFlagT
REAL,   INTENT(  OUT),DIMENSION(2,nSurfNodesT ) :: xSurfT
REAL,   INTENT(  OUT),DIMENSION(2,nSurfNodesT ) :: pointingVecT
INTEGER :: j,k,m,n,i,n1,n2,nScomp,cL,cR,nL,nR
REAL :: dx,dy,ds,nx,ny,ax,ay,a1,nxL,nyL,nxR,nyR,dt,tL,t
REAL,PARAMETER :: eps = 1.D-14
INTEGER,DIMENSION(nSurfNodesT) :: nmap,Rface,Lface
INTEGER,DIMENSION(nSurfNodes ) :: nTag,Rnode,Lnode


! establish new node numbering and fill in sFlag
j      = 0
k      = 0
sFlagT = 0
DO n=1,nSurfNodes
   i = addNode(n)
   IF (i > 1) j = j+1
DO m=1,i
   k       = k+1
   nmap(k) = n
   IF (i > 1) sFlagT(k) = j
END DO
END DO
IF (j /= nSharp .OR. k /= nSurfNodesT) THEN
   WRITE(*,*)'*** Problem adding sharp nodes in fillsharpcorners.f90 ***'
   STOP
END IF


! fill in an inverse map which gives rightmost and leftmost nodes
Rnode = 0
Lnode = nSurfNodesT+1
DO n=1,nSurfNodesT
   i        = nmap(n)
   Rnode(i) = MAX(Rnode(i),n)
   Lnode(i) = MIN(Lnode(i),n)
END DO


! update boundary edges (assume no sharp corners on boundaries)
DO n=1,nBndEdges
   n1          = bndEdges(n)
   bndEdges(n) = Rnode(n1)
END DO


! set node clipping and surface coordinates
DO n=1,nSurfNodesT
   i            = nmap(n)
   nodeClipT(n) = nodeClip(i)
   xSurfT(:,n)  = xSurf(:,i)
END DO


! insert new faces into the connectivity
DO n=1,nSurfFaces
   n1              = surfFaces(1,n)
   n2              = surfFaces(2,n)
   surfFacesT(1,n) = Rnode(n1)
   surfFacesT(2,n) = Lnode(n2)
   faceTagT(n)     = faceTag(n)
END DO
k = nSurfFaces
DO n=1,nSurfNodesT-1
   i               = sFlagT(n  )
   j               = sFlagT(n+1)
IF (i == j .AND. i /= 0) THEN
   k               = k+1
   surfFacesT(1,k) = n
   surfFacesT(2,k) = n+1
END IF
END DO
IF (k /= nSurfFacesT) THEN
   WRITE(*,*)'*** Problem adding faces in fillsharpcorners.f90 ***'
   STOP
END IF


! set face clipping (this needs to be done here?)
DO n=1,nSurfFacesT
   n1           = surfFacesT(1,n)
   n2           = surfFacesT(2,n)
   faceClipT(n) = MIN(nodeClipT(n1),nodeClipT(n2))
END DO


! set face tags (may have problems where two boundary types meet)
DO n=1,nSurfFaces
   n1          = surfFaces(1,n)
   n2          = surfFaces(2,n)
   nTag(n1)    = faceTag(n)
   nTag(n2)    = faceTag(n)
END DO
DO n=nSurfFaces+1,nSurfFacesT
   n1          = surfFacesT(1,n)
   n2          = nmap(n1)
   faceTagT(n) = nTag(n2)
END DO


! initialize pointing vectors
! first, find right and left cells for each node
Rface = 0
Lface = 0
DO n=1,nSurfFacesT
   n1        = surfFacesT(1,n)
   n2        = surfFacesT(2,n)
   Rface(n1) = n
   Lface(n2) = n
END DO


! next find the node normal as the average of the cell normals
pointingVecT = 0.
DO n=1,nSurfFacesT
   n1                 = surfFacesT(1,n)
   n2                 = surfFacesT(2,n)
   dx                 = xSurfT(1,n2)-xSurfT(1,n1)
   dy                 = xSurfT(2,n2)-xSurfT(2,n1)
   pointingVecT(1,n1) = pointingVecT(1,n1)-dy
   pointingVecT(2,n1) = pointingVecT(2,n1)+dx
   pointingVecT(1,n2) = pointingVecT(1,n2)-dy
   pointingVecT(2,n2) = pointingVecT(2,n2)+dx
END DO
DO n=1,nSurfNodesT
!~    dx                 = pointingVecT(1,n) !Shaun change, uncomment these lines and comment the next 2 lines to get the actual pointing vector instead of just straight vertical
!~    dy                 = pointingVecT(2,n)
   dx                 = 0.
   dy                 = 1.
   ds                 = 1./MAX(eps,SQRT(dx*dx+dy*dy))
   pointingVecT(1,n)  = dx*ds
   pointingVecT(2,n)  = dy*ds
END DO

! constrain boundary vectors to lie in their planes
nScomp = MAXVAL(faceTag)
DO n=1,nBndEdges
   i                 = bndEdges(n)
   j                 = bndTag(n)-nScomp
   nx                = bndNormal(1,j)
   ny                = bndNormal(2,j)
   k                 = MAX(Rface(i),Lface(i))
   n1                = surfFacesT(1,k)
   n2                = surfFacesT(2,k)
   dx                = xSurfT(1,n2)-xSurfT(1,n1)
   dy                = xSurfT(2,n2)-xSurfT(2,n1)
   ds                = 1./SQRT(dx*dx+dy*dy)
   ax                =-dy*ds
   ay                = dx*ds
   a1                = ax*nx+ay*ny
   ax                = ax-a1*nx
   ay                = ay-a1*ny
   ds                = 1./SQRT(ax*ax+ay*ay)
   pointingVecT(1,i) = ax*ds
   pointingVecT(2,i) = ay*ds
END DO

! fix sharp nodes
DO n=1,nSurfNodes
   nL                = Lnode(n)
   nR                = Rnode(n)
IF (nL /= nR) THEN
   cL                = Lface(nL)
   cR                = Rface(nR)
   n1                = surfFacesT(1,cL)
   n2                = surfFacesT(2,cL)
   dx                = xSurfT(1,n2)-xSurfT(1,n1)
   dy                = xSurfT(2,n2)-xSurfT(2,n1)
   ds                = 1./SQRT(dx*dx+dy*dy)
   nxL               =-dy*ds
   nyL               = dx*ds
   n1                = surfFacesT(1,cR)
   n2                = surfFacesT(2,cR)
   dx                = xSurfT(1,n2)-xSurfT(1,n1)
   dy                = xSurfT(2,n2)-xSurfT(2,n1)
   ds                = 1./SQRT(dx*dx+dy*dy)
   nxR               =-dy*ds
   nyR               = dx*ds
   dt                = nxR*nxL+nyR*nyL
   dt                = ACOS(dt)/REAL(nR-nL)
   tL                = ATAN2(nyL,nxL)
DO m=nL,nR
   t                 = tL-REAL(m-nL)*dt
   pointingVecT(1,m) = COS(t)
   pointingVecT(2,m) = SIN(t)
END DO
END IF
END DO

!write(*,*)'bndEdges:'
!do n=1,nBndEdges
!   write(*,*)n,bndEdges(n)
!end do
!write(*,*)
!write(*,*)'surfFaces:'
!do n=1,nSurfFacesT
!   write(*,*)n,surfFacesT(:,n)
!end do
!write(*,*)
!write(*,*)'faceClip, faceTag:'
!do n=1,nSurfFacesT
!   write(*,*)n,faceClipT(n),faceTagT(n)
!end do
!write(*,*)
!write(*,*)'nodeClip, sFlag:'
!do n=1,nSurfNodesT
!   write(*,*)n,nodeClipT(n),sFlagT(n)
!end do
!write(*,*)
!write(*,*)'xSurf, pointingVec:'
!do n=1,nSurfNodesT
!   write(*,*)n,xSurfT(:,n),pointingVecT(:,n)
!end do
!stop


END SUBROUTINE fillsharpcorners

