!> \brief
!! This subroutine reads the mesh data.
!!
!! Comments:
!!  surfFaces   - dimensioned (2,nSurfFaces), the 2 nodes forming each surface
!!                face, in order such that the domain to be gridded is on
!!                the left.
!!  bndNodes    - dimensioned (nBndNodes), a list of all boundary nodes
!!  bndNormal   - dimensioned (2,nBndNodes), a list of (nx,ny) unit normals
!!                for each plane at a boundary node
!!  faceTag     - dimensioned (nSurfFaces), the surface patch number for each
!!                surface face
!!  nodeTag     - dimensioned (nBndNodes), the boundary patch number for each
!!                boundary node
!!  xSurf       - dimensioned (2,nSurfNodes), the coordinates for each surface
!!                node
!!  nodeClip    - dimensioned (nSurfNodes), the clipping index for each surface
!!                node. Note that this can be derived from faceClip, but is
!!                alread included here.
!!  faceClip    - dimensioned (nSurfFaces), the clipping index for each surface
!!                face.
!!  pointingVec - dimensioned (2,nSurfNodes), the strand pointing vector at each
!!                surface node
!!  xStrand     - dimensioned (0:nPtsPerStrand), the 1d distribution of
!!                strand coordinates.
!!
!! Versions:
!!   - 1.0 Katz 08/03/2010
!!
!! Additional notes:\par
!!  none
!!
!! Source code:
!!   \include readmeshdata.f90


SUBROUTINE readmeshdata(iunit, &
		        surfaceOnly0, &
                        nMesh0, &
                        nSurfNodes0, &
                        nSurfFaces0, &
                        nBndEdges0, &
                        nSurfPatches0, &
                        nEdgePatches0, &
                        nSharp0, &
                        nPtsPerStrand0, &
                        surfFaces0, &
                        bndEdges0, &
                        sFlag0, &
                        nodeClip0, &
                        faceClip0, &
                        faceTag0, &
                        bndTag0, &
                        xSurf0, &
                        pointingVec0, &
                        xStrand0, &
                        bndNormal0)

USE avdefs

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: iunit, &
		         surfaceOnly0, &
                         nMesh0, &
                         nSurfNodes0, &
                         nSurfFaces0, &
                         nBndEdges0, &
                         nSurfPatches0, &
                         nEdgePatches0, &
                         nPtsPerStrand0
INTEGER,INTENT(  OUT) :: nSharp0
INTEGER,INTENT(  OUT),DIMENSION(2,nSurfFaces0   ) :: surfFaces0
INTEGER,INTENT(  OUT),DIMENSION(  nBndEdges0    ) :: bndEdges0
INTEGER,INTENT(  OUT),DIMENSION(  nSurfNodes0   ) :: sFlag0
INTEGER,INTENT(  OUT),DIMENSION(  nSurfNodes0   ) :: nodeClip0
INTEGER,INTENT(  OUT),DIMENSION(  nSurfFaces0   ) :: faceClip0
INTEGER,INTENT(  OUT),DIMENSION(  nSurfFaces0   ) :: faceTag0
INTEGER,INTENT(  OUT),DIMENSION(  nBndEdges0    ) :: bndTag0
REAL   ,INTENT(  OUT),DIMENSION(2,nSurfNodes0   ) :: xSurf0
REAL   ,INTENT(  OUT),DIMENSION(2,nSurfNodes0   ) :: pointingVec0
REAL   ,INTENT(  OUT),DIMENSION(0:nPtsPerStrand0) :: xStrand0
REAL   ,INTENT(  OUT),DIMENSION(2,nEdgePatches0 ) :: bndNormal0
INTEGER :: k,n
TYPE(strand2d_surf_patch), &
     ALLOCATABLE,DIMENSION(:,:) :: strandSurfPatch
TYPE(strand2d_edge_patch), &
     ALLOCATABLE,DIMENSION(:,:) :: strandEdgePatch


! assume no sharp corners initially
nSharp0 = 0
sFlag0  = 0


ALLOCATE(strandSurfPatch(nSurfPatches0,nMesh0), &
         strandEdgePatch(nEdgePatches0,nMesh0))

DO n=1,nSurfPatches0
   READ(iunit)strandSurfPatch(n,1)%surfPatchID  ,&
           strandSurfPatch(n,1)%surfPatchBody   ,&
           strandSurfPatch(n,1)%surfPatchComp   ,&
           strandSurfPatch(n,1)%surfPatchBCType
END DO

DO n=1,nEdgePatches0
   READ(iunit)strandEdgePatch(n,1)%edgePatchID  ,&
           strandEdgePatch(n,1)%edgePatchBody   ,&
           strandEdgePatch(n,1)%edgePatchComp   ,&
           strandEdgePatch(n,1)%edgePatchBCType ,&
           strandEdgePatch(n,1)%nx              ,&
           strandEdgePatch(n,1)%ny
END DO

READ(iunit)((surfFaces0  (k,n),k=1,2),n=1,nSurfFaces0   )
READ(iunit)( bndEdges0   (  n),       n=1,nBndEdges0    )
READ(iunit)( nodeClip0   (  n),       n=1,nSurfNodes0   )
READ(iunit)( faceClip0   (  n),       n=1,nSurfFaces0   )
READ(iunit)( faceTag0    (  n),       n=1,nSurfFaces0   )
READ(iunit)( bndTag0     (  n),       n=1,nBndEdges0    )
READ(iunit)((xSurf0      (k,n),k=1,2),n=1,nSurfNodes0   )
READ(iunit)((pointingVec0(k,n),k=1,2),n=1,nSurfNodes0   )
IF (surfaceOnly0 == 0) &
READ(iunit)( xStrand0    (  n),       n=0,nPtsPerStrand0)
DO n=1,nEdgePatches0
   bndNormal0(1,n) = strandEdgePatch(n,1)%nx
   bndNormal0(2,n) = strandEdgePatch(n,1)%ny
END DO


CLOSE(iunit)

DEALLOCATE(strandSurfPatch)
DEALLOCATE(strandEdgePatch)

END SUBROUTINE readmeshdata
