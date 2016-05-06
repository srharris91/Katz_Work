PROGRAM sine

IMPLICIT NONE

INTEGER :: nSurfElem,nSurfNode,nBndNode,nCompBd
INTEGER :: n,m,narg,iu,perturb
INTEGER,ALLOCATABLE,DIMENSION(:) :: order
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: surfElem
INTEGER,ALLOCATABLE,DIMENSION(:) :: surfElemTag
INTEGER,ALLOCATABLE,DIMENSION(:) :: bndNode
INTEGER,ALLOCATABLE,DIMENSION(:) :: bndNodeTag
REAL :: freq,dx,pi
REAL,ALLOCATABLE,DIMENSION(:,:) :: surfX,bndNodeNormal
CHARACTER(80) :: fileName,arg


! read inputs
narg = iargc()
m    = 0
DO WHILE (m < narg)
   m = m+1
   CALL GETARG(m,arg)
   IF (arg(1:2) == '-n') THEN
      m = m+1
      CALL GETARG(m,arg)
      READ(arg,*)nSurfElem
   ELSE IF (arg(1:2) == '-p') THEN
      m = m+1
      CALL GETARG(m,arg)
      READ(arg,*)perturb
   ELSE IF (arg(1:2) == '-a') THEN
      m = m+1
      CALL GETARG(m,arg)
      READ(arg,*)freq
   ELSE IF (arg(1:2) == '-m') THEN
      m = m+1
      CALL GETARG(m,arg)
      fileName = arg
   END IF
END DO


! compute dimensions of mesh and allocate data
nSurfNode = nSurfElem+1
nBndNode  = 2
nCompBd   = 3

ALLOCATE( &
     order(nSurfElem), &
     surfElem(2,nSurfElem), &
     surfElemTag(nSurfElem), &
     surfX(2,nSurfNode), &
     bndNode(nBndNode), &
     bndNodeTag(nBndNode), &
     bndNodeNormal(2,nBndNode))


! set surface tags and boundary information
order              = 1 !linear elements
surfElemTag        = 0
bndNode(1)         = 0
bndNode(2)         = nSurfNode-1
bndNodeTag(1)      = 1
bndNodeTag(2)      = 2
bndNodeNormal(1,1) =-1.
bndNodeNormal(2,1) = 0.
bndNodeNormal(1,2) = 1.
bndNodeNormal(2,2) = 0.


! generate grid data
CALL SRAND(0)
DO n=1,nSurfElem
   surfElem(1,n) = n-1
   surfElem(2,n) = n
END DO
pi = 4.*ATAN(1.)
dx = 2.*pi/REAL(nSurfElem)
DO n=1,nSurfNode
   surfX(1,n) = dx*REAL(n-1)
   surfX(2,n) = SIN(freq*surfX(1,n))
END DO


! write the mesh to file
iu = 7
OPEN(UNIT=iu,FILE=filename)
WRITE(iu,*)nSurfElem,nSurfNode,nCompBd,nBndNode
DO n=1,nSurfElem
   WRITE(iu,*)order(n),surfElem(:,n),surfElemTag(n)
END DO
DO n=1,nSurfNode
   WRITE(iu,*)surfX(:,n)
END DO
DO n=1,nBndNode
   WRITE(iu,*)bndNode(n),bndNodeTag(n),bndNodeNormal(:,n)
END DO
CLOSE(iu)

WRITE(*,*)
WRITE(*,'(A,A)')  'Successfully wrote: ', filename
WRITE(*,*)


END PROGRAM sine
