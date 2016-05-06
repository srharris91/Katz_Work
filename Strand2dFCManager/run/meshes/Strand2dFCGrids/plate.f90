PROGRAM plate

IMPLICIT NONE

INTEGER :: nSurfElem,nSurfNode,nBndNode,nCompBd
INTEGER :: n,m,narg,iu,isf,nen,npl,npl0,nmax,fact
INTEGER,ALLOCATABLE,DIMENSION(:) :: order
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: surfElem
INTEGER,ALLOCATABLE,DIMENSION(:) :: surfElemTag
INTEGER,ALLOCATABLE,DIMENSION(:) :: bndNode
INTEGER,ALLOCATABLE,DIMENSION(:) :: bndNodeTag
REAL :: pi,angle,dx,x0,x1,dx0,dx1,rmax,dx0a,dx1a
REAL,ALLOCATABLE,DIMENSION(:  ) :: xsb,xsp
REAL,ALLOCATABLE,DIMENSION(:,:) :: surfX,bndNodeNormal
CHARACTER(80) :: fileName,arg


! read inputs
narg = iargc()
m    = 0
DO WHILE (m < narg)
   m    = m+1
   CALL GETARG(m,arg)
   IF (arg(1:2) == '-n') THEN ! -n number of points on plate in total.
      m = m+1
      CALL GETARG(m,arg)
      READ(arg,*)npl
   ELSE IF (arg(1:2) == '-o') THEN ! -o FileName this is the specified output filename.mesh
      m = m+1
      CALL GETARG(m,arg)
      fileName = arg
   END IF
END DO

npl0 = npl
!~ npl  = 128 !commented out by Shaun Feb. 2015
fact = npl/npl0


! Compute 1-d stretching function before plate
nen  = npl/2
isf  = 2
x0   =-1.
x1   = 0.
dx0  = .96
dx1  = .01
n    = nen+1
rmax = 1.1
nmax = 200
dx0a = 1.
ALLOCATE(xsb(nmax))
CALL sfuns(isf,x0,x1,dx0,dx1,n,rmax,nmax,dx0a,dx1a,xsb)


! Compute 1-d stretching function on plate
isf  = 2
x0   = 0.
x1   = 1.
!~ x1   = 4.
dx0  = .01
dx1  = .48
n    = npl+1
rmax = 1.1
nmax = 200
dx0a = 1.
ALLOCATE(xsp(nmax))
CALL sfuns(isf,x0,x1,dx0,dx1,n,rmax,nmax,dx0a,dx1a,xsp)


! compute dimensions of mesh and allocate data
nSurfElem = nen+npl
nSurfNode = nSurfElem+1
nBndNode  = 2
nCompBd   = 4

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
bndNode(1)         = 0
bndNode(2)         =(nSurfNode-1)/fact
bndNodeTag(1)      = 2
bndNodeTag(2)      = 3
bndNodeNormal(1,1) =-1.
bndNodeNormal(2,1) = 0.
bndNodeNormal(1,2) = 1.
bndNodeNormal(2,2) = 0.


! generate grid data
DO n=1,nSurfElem
   surfElem(1,n) = n-1
   surfElem(2,n) = n
   IF (n <= nen) THEN
      surfElemTag(n) = 0
   ELSE
      surfElemTag(n) = 1
   END IF
END DO

DO n=1,nSurfNode
   IF (n <= nen) THEN
      surfX(1,n) = xsb(n)
   ELSE
      surfX(1,n) = xsp(n-nen)
   END IF
   surfX(2,n) = 0.
END DO


! write the mesh to file
iu = 7
OPEN(UNIT=iu,FILE=filename)
WRITE(iu,*)nSurfElem/fact,(nSurfNode-1)/fact+1,nCompBd,nBndNode
DO n=1,nSurfElem,fact
   WRITE(iu,*)order(n),surfElem(1,n)/fact,surfElem(1,n)/fact+1,surfElemTag(n)
END DO
DO n=1,nSurfNode,fact
   WRITE(iu,*)surfX(:,n)
END DO
DO n=1,nBndNode
   WRITE(iu,*)bndNode(n),bndNodeTag(n),bndNodeNormal(:,n)
END DO
CLOSE(iu)

WRITE(*,*)
WRITE(*,'(A,A)')  'Successfully wrote: ', filename
WRITE(*,*)


END PROGRAM plate
