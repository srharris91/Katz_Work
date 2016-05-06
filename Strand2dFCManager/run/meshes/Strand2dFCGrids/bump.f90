PROGRAM bump

IMPLICIT NONE

INTEGER :: n,m,k,narg,i,j,iu,level,skip,nSurfNodeF,nSurfNode,nSurfElem
INTEGER :: nBndNode,nCompBd,isf,nmax
INTEGER,ALLOCATABLE,DIMENSION(:  ) :: surfOrder,surfElemTag,bndNode,bndNodeTag
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: surfElem
REAL :: pi,x0,x1,dx0,dx1,rmax,dx0a,dx1a
REAL,ALLOCATABLE,DIMENSION(:  ) :: xT
REAL,ALLOCATABLE,DIMENSION(:,:) :: surfX,bndNodeNormal
CHARACTER(80) :: meshfile,arg


! read inputs
! level = 1 is 512 x 128
! level = 2 is 264 x 64
! level = 3 is 128 x 32
! level = 4 is 64  x 16
! level = 5 is 32  x 8
! level = 6 is 16  x 4
! level = 7 is 8   x 2
narg = iargc()
m    = 0
DO WHILE (m < narg)
   m    = m+1
   CALL GETARG(m,arg)
   IF (arg(1:2) == '-n') THEN
      m = m+1
      CALL GETARG(m,arg)
      READ(arg,*)level
   ELSE IF (arg(1:2) == '-o') THEN
      m = m+1
      CALL GETARG(m,arg)
      meshfile = arg
   END IF
END DO

IF (level >= 8) THEN
   WRITE(*,*)'*** please choose level 7 or lower ***'
   STOP
END IF


! determine grid dimensions
nSurfNodeF = 513
skip       = 2**(level-1)
nSurfNode  =(nSurfNodeF-1)/skip+1
nSurfElem  = nSurfNode -1
nCompBd    = 3
nBndNode   = 2
pi         = 4.*ATAN(1.)


! allocate data
ALLOCATE(surfOrder(nSurfElem))
ALLOCATE(surfElem(2,nSurfElem))
ALLOCATE(surfElemTag(nSurfElem))
ALLOCATE(surfX(2,nSurfNode))
ALLOCATE(bndNode(nBndNode))
ALLOCATE(bndNodeTag(nBndNode))
ALLOCATE(bndNodeNormal(2,nBndNode))


! determine stretching function along the bump
isf  = 2
x0   = 0.
x1   = 1.
dx0  = .064
dx1  = .0002
rmax = 1.1
nmax = nSurfNodeF*2
dx0a = 1.
n    =(nSurfNodeF-1)/2+1
ALLOCATE(xT(nmax))
CALL sfuns(isf,x0,x1,dx0,dx1,n,rmax,nmax,dx0a,dx1a,xT)
k = 0
DO j=1,(nSurfNodeF-1)/2+1,skip
   k          = k+1
!~    surfX(1,k) = 25.*(xT(j)-1.) !Shaun change, switch with line below for normal
   surfX(1,k) = 5.*(xT(j)-1.)
END DO
k = k-1
DO j=(nSurfNodeF-1)/2+1,1,-skip
   k          = k+1
!~    surfX(1,k) =-25.*(xT(j)-1.) !Shaun change, switch with line below for normal
   surfX(1,k) =-5.*(xT(j)-1.)
END DO
DEALLOCATE(xT)

DO n=1,nSurfNode
   IF (surfX(1,n) < -.5 .OR. surfX(1,n) > .5) THEN
      surfX(2,n) = 0.
   ELSE
      surfX(2,n) = 0.05*(SIN(pi*(.5+surfX(1,n))))**4
   END IF
END DO


! fill in connectivity data
surfOrder   = 1
surfElemTag = 0
!Shauns change to surfElemTag to make 0 before bump and 1 on and after bump
DO n=1,nSurfElem
   IF (surfX(1,n) < -.5) THEN !Shauns change, set inviscid before bump
      surfElemTag(n)=0
   ELSE
      surfElemTag(n)=1 !viscous on and after bump
   END IF
END DO
   
   
DO n=1,nSurfElem
   surfElem(1,n) = n-1
   surfElem(2,n) = n
END DO


! fill in boundary node arrays
bndNode(1)         = 0
bndNode(2)         = nSurfNode-1
bndNodeTag(1)      = 1
bndNodeTag(2)      = 2
bndNodeNormal(1,1) =-1.
bndNodeNormal(2,1) = 0.
bndNodeNormal(1,2) = 1.
bndNodeNormal(2,2) = 0.


! write to mesh file
OPEN(iu,FILE=meshfile,STATUS='replace',FORM='formatted')
WRITE(iu,*)nSurfElem,nSurfNode,nCompBd,nBndNode
DO n=1,nSurfElem
   WRITE(iu,*)surfOrder(n),surfElem(:,n),surfElemTag(n)
END DO
DO n=1,nSurfNode
   WRITE(iu,*)surfX(:,n)
END DO
DO n=1,nBndNode
   WRITE(iu,*)bndNode(n),bndNodeTag(n),bndNodeNormal(:,n)
END DO
CLOSE(iu)


! deallocate data
DEALLOCATE(surfOrder)
DEALLOCATE(surfElem)
DEALLOCATE(surfElemTag)
DEALLOCATE(surfX)
DEALLOCATE(bndNode)
DEALLOCATE(bndNodeTag)
DEALLOCATE(bndNodeNormal)


WRITE(*,*)
WRITE(*,*)'*** successfully wrote ',TRIM(meshfile),' ***'
WRITE(*,*)


END PROGRAM bump
