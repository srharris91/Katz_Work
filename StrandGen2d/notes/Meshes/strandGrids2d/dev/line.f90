PROGRAM line

IMPLICIT NONE

LOGICAL :: surfaceOnly
CHARACTER(40) :: gfile,arg
INTEGER :: iu,i,k,n,npl,nen,isf,nmax,perturb,m,narg
INTEGER :: nSurfNodes,nSurfFaces,nBndNodes,nPtsPerStrand,nSurfPatches, &
           nNodePatches
REAL    :: dx,angle,pi
REAL,PARAMETER :: fact = 0.5
INTEGER,ALLOCATABLE,DIMENSION(:  ) :: bndNodes,faceTag,nodeTag,nodeClip,faceClip
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: surfFaces
REAL   ,ALLOCATABLE,DIMENSION(:  ) :: xStrand,xsb,xsp
REAL   ,ALLOCATABLE,DIMENSION(:,:) :: bndNormal,xSurf,pointingVec


!WRITE(*,*)'number of cells on the line:'
!READ(*,*)npl
!WRITE(*,*)'perturb? (0=no, 1=yes):'
!READ(*,*)perturb
!WRITE(*,*)'angle of boundary planes (deg):'
!READ(*,*)angle
!WRITE(*,*)'output file name:'
!READ(*,*)gfile

! read inputs

narg = iargc()
m    = 0
DO WHILE (m < narg)
   m    = m+1
   CALL GETARG(m,arg)
   IF (arg(1:2) == '-n') THEN
      m = m+1
      CALL GETARG(m,arg)
      READ(arg,*)npl
   ELSE IF (arg(1:2) == '-p') THEN
      m = m+1
      CALL GETARG(m,arg)
      READ(arg,*)perturb
   ELSE IF (arg(1:2) == '-a') THEN
      m = m+1
      CALL GETARG(m,arg)
      READ(arg,*)angle
   ELSE IF (arg(1:2) == '-m') THEN
      m = m+1
      CALL GETARG(m,arg)
      gfile = arg
   END IF
END DO
pi = 4.*ATAN(1.)
angle = angle*pi/180.


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

dx = 1./REAL(nSurfFaces)
DO i=1,nSurfNodes
   xSurf(1,i) = dx*REAL(i-1)
   xSurf(2,i) = 0.
IF (perturb == 1 .AND. i > 1 .AND. i < nSurfNodes) THEN
   xSurf(1,i) = xSurf(1,i)+fact*dx*(RAND(0)-.5)
END IF
END DO


! form connectivity

DO i=1,nSurfFaces
   surfFaces(1,i) = i
   surfFaces(2,i) = i+1
   faceTag(i) = 1
END DO
bndNodes(1) = 1
bndNodes(2) = nSurfNodes
bndNormal(1,1) =-SIN(angle)
bndNormal(2,1) = COS(angle)
bndNormal(1,2) = SIN(angle)
bndNormal(2,2) =-COS(angle)
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


END PROGRAM line
