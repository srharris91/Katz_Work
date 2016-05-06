PROGRAM airfoilgrid

IMPLICIT NONE

LOGICAL :: surfaceOnly
CHARACTER(40) :: gfile,ptsFile
INTEGER :: iu,n,k
INTEGER :: nSurfNodes,nSurfFaces,nBndNodes,nPtsPerStrand,nSurfPatches, &
           nNodePatches
INTEGER,ALLOCATABLE,DIMENSION(:  ) :: bndNodes,faceTag,nodeTag,nodeClip,faceClip
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: surfFaces
REAL   ,ALLOCATABLE,DIMENSION(:  ) :: xStrand
REAL   ,ALLOCATABLE,DIMENSION(:,:) :: bndNormal,xSurf,pointingVec


WRITE(*,*)'surface points file:'
READ(*,*)ptsFile
WRITE(*,*)'output file name:'
READ(*,*)gfile
iu = 7
OPEN(iu,FILE=ptsFile,STATUS='old')
READ(iu,*)
READ(iu,*)
READ(iu,*)
READ(iu,*)
READ(iu,*)nSurfNodes
READ(iu,*)
READ(iu,*)
READ(iu,*)
READ(iu,*)
READ(iu,*)
READ(iu,*)
READ(iu,*)

! preliminaries

surfaceOnly   = .true.
nSurfNodes    = nSurfNodes-1 !end points are duplicated
nSurfFaces    = nSurfNodes !closed shapes only
nBndNodes     = 0
nPtsPerStrand = 0
nSurfPatches  = 1
nNodePatches  = 0

ALLOCATE(surfFaces  (2,nSurfFaces), &
         bndNodes   (  nBndNodes ), &
         bndNormal  (2,nBndNodes ), &
         faceTag    (  nSurfFaces), &
         nodeTag    (  nBndNodes ), &
         xSurf      (2,nSurfNodes), &
         nodeClip   (  nSurfNodes), &
         faceClip   (  nSurfFaces), &
         pointingVec(2,nSurfNodes), &
         xStrand    (0:nPtsPerStrand))


! read coordinates

DO n=1,nSurfNodes
   READ(iu,*)(xSurf(k,n),k=1,2)
END DO


! form connectivity

DO n=1,nSurfFaces
   surfFaces(1,n) = n
   surfFaces(2,n) = n+1
   faceTag(n)     = 1
END DO
surfFaces(2,nSurfFaces) = 1
nodeClip = 0
faceClip = 0
pointingVec = 0.
xStrand = 0.


! write to file

iu = 12
OPEN(iu,FILE=TRIM(gfile),STATUS='replace',FORM='unformatted', &
     CONVERT='big_endian')

WRITE(iu)surfaceOnly
WRITE(iu)nSurfNodes
WRITE(iu)nSurfFaces
WRITE(iu)nBndNodes
WRITE(iu)nPtsPerStrand
WRITE(iu)nSurfPatches
WRITE(iu)nNodePatches

WRITE(iu)((surfFaces  (k,n),k=1,2),n=1,nSurfFaces   )
WRITE(iu) (bndNodes   (  n)       ,n=1,nBndNodes    )
WRITE(iu)((bndNormal  (k,n),k=1,2),n=1,nBndNodes    )
WRITE(iu) (faceTag    (  n)       ,n=1,nSurfFaces   )
WRITE(iu) (nodeTag    (  n)       ,n=1,nBndNodes    )
WRITE(iu)((xSurf      (k,n),k=1,2),n=1,nSurfNodes   )
WRITE(iu) (nodeClip   (  n)       ,n=1,nSurfNodes   )
WRITE(iu) (faceClip   (  n)       ,n=1,nSurfFaces   )
WRITE(iu)((pointingVec(k,n),k=1,2),n=1,nSurfNodes   )
WRITE(iu) (xStrand(n),n=0,nPtsPerStrand)

CLOSE(iu)

DEALLOCATE(surfFaces,   &
           bndNodes,    &
           bndNormal,   &
           faceTag,     &
           nodeTag,     &
           xSurf,       &
           nodeClip,    &
           faceClip,    &
           pointingVec, &
           xStrand)


END PROGRAM airfoilgrid
