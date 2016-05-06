PROGRAM plateGrid

IMPLICIT NONE

LOGICAL :: surfaceOnly
CHARACTER(40) :: gfile
INTEGER :: iu,i,k,n,npl,nen,isf,nmax
INTEGER :: nSurfNodes,nSurfFaces,nBndNodes,nPtsPerStrand,nSurfPatches, &
           nNodePatches
REAL    :: x0,x1,dx0,dx1,rmax,dx0a,dx1a
INTEGER,ALLOCATABLE,DIMENSION(:  ) :: bndNodes,faceTag,nodeTag,nodeClip,faceClip
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: surfFaces
REAL   ,ALLOCATABLE,DIMENSION(:  ) :: xStrand,xsb,xsp
REAL   ,ALLOCATABLE,DIMENSION(:,:) :: bndNormal,xSurf,pointingVec


WRITE(*,*)'number of points on the plate (even):'
READ(*,*)npl
WRITE(*,*)'output file name:'
READ(*,*)gfile


! Compute 1-d stretching function before plate

nen  = npl/2
isf  = 2
x0   =-2.
x1   = 0.
dx0  = .48
dx1  = .005
n    = nen+1
rmax = 1.1
nmax = 200
dx0a = 1.
ALLOCATE(xsb(nmax))
CALL sfuns(isf,x0,x1,dx0,dx1,n,rmax,nmax,dx0a,dx1a,xsb)


! Compute 1-d stretching function on plate

isf  = 2
x0   = 0.
x1   = 2.
dx0  = .005
dx1  = .08
n    = npl+1
rmax = 1.1
nmax = 200
dx0a = 1.
ALLOCATE(xsp(nmax))
CALL sfuns(isf,x0,x1,dx0,dx1,n,rmax,nmax,dx0a,dx1a,xsp)


! preliminaries

surfaceOnly   = .true.
nSurfNodes    = nen+npl+1
nSurfFaces    = nen+npl
nBndNodes     = 2
nPtsPerStrand = 0
nSurfPatches  = 2
nNodePatches  = 2

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


! Generate coordinates

DO i=1,nSurfNodes
   IF (i <= nen) THEN
      xSurf(1,i) = xsb(i)
      xSurf(2,i) = 0.
   ELSE
      xSurf(1,i) = xsp(i-nen)
      xSurf(2,i) = 0.
   END IF
END DO


! form connectivity

DO i=1,nSurfFaces
   surfFaces(1,i) = i
   surfFaces(2,i) = i+1
   IF (i <= nen) THEN
      faceTag(i) = 1
   ELSE
      faceTag(i) = 2
   END IF
END DO
bndNodes(1) = 1
bndNodes(2) = nSurfNodes
bndNormal(1,1) =-1.
bndNormal(2,1) = 0.
bndNormal(1,2) = 1.
bndNormal(2,2) = 0.
nodeTag(1) = 3
nodeTag(2) = 4
nodeClip = 0
faceClip = 0
pointingVec = 0.
xStrand = 0.


! write to file

iu = 12
OPEN(iu,FILE=TRIM(gfile),STATUS='replace',FORM='unformatted')

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
DEALLOCATE(xsp,xsb)


END PROGRAM plateGrid
