PROGRAM curve

IMPLICIT NONE

LOGICAL :: surfaceOnly
CHARACTER(40) :: gfile
INTEGER :: iu,i,k,n,npl,nen,isf,nmax,perturb
INTEGER :: nSurfNodes,nSurfFaces,nBndNodes,nPtsPerStrand,nSurfPatches, &
           nNodePatches
REAL    :: dx,r0,pi,al,al0,alN,ds
REAL,PARAMETER :: fact = 0.5
INTEGER,ALLOCATABLE,DIMENSION(:  ) :: bndNodes,faceTag,nodeTag,nodeClip, &
                                      faceClip
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: surfFaces
REAL   ,ALLOCATABLE,DIMENSION(:  ) :: xStrand,xsb,xsp
REAL   ,ALLOCATABLE,DIMENSION(:,:) :: bndNormal,xSurf,pointingVec


WRITE(*,*)'number of cells on the curve:'
READ(*,*)npl
WRITE(*,*)'perturb? (0=no, 1=yes):'
READ(*,*)perturb
WRITE(*,*)'output file name:'
READ(*,*)gfile


! preliminaries

surfaceOnly   = .true.
nSurfNodes    = npl+1
nSurfFaces    = npl
nBndNodes     = 2
nPtsPerStrand = 0
nSurfPatches  = 1
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

pi = 4.*ATAN(1.)
r0 = .5
al0 = .01*pi
alN = .49*pi
dx = (alN-al0)/REAL(npl)
DO i=1,nSurfNodes
   al         = alN-REAL(i-1)*dx
IF (perturb == 1 .AND. i > 1 .AND. i < nSurfNodes) THEN
   ds         = fact*dx*(RAND(0)-.5)
   al         = al+ds
END IF
   xSurf(1,i) = r0*COS(al)
   xSurf(2,i) = r0*SIN(al)
END DO


! form connectivity

DO i=1,nSurfFaces
   surfFaces(1,i) = i
   surfFaces(2,i) = i+1
   faceTag(i) = 1
END DO
bndNodes(1) = 1
bndNodes(2) = nSurfNodes
bndNormal(1,1) =-SIN(alN)
bndNormal(2,1) = COS(alN)
bndNormal(1,2) = SIN(al0)
bndNormal(2,2) =-COS(al0)
nodeTag(1) = 2
nodeTag(2) = 3
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


END PROGRAM curve
