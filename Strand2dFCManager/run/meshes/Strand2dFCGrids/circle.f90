PROGRAM circle

IMPLICIT NONE

INTEGER :: nSurfElem,nSurfNode,nBndNode,nCompBd
INTEGER :: n,m,narg,iu,i,k,km,orderM,perturb,n0,n1
INTEGER,ALLOCATABLE,DIMENSION(:) :: order
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: surfElem
INTEGER,ALLOCATABLE,DIMENSION(:) :: surfElemTag
INTEGER,ALLOCATABLE,DIMENSION(:) :: bndNode
INTEGER,ALLOCATABLE,DIMENSION(:) :: bndNodeTag
REAL :: pi,t,dt,t0,t1
REAL,ALLOCATABLE,DIMENSION(:  ) :: xsb,xsp
REAL,ALLOCATABLE,DIMENSION(:,:) :: surfX,bndNodeNormal
CHARACTER(80) :: fileName,arg


! read inputs
narg = iargc()
m    = 0
DO WHILE (m < narg)
   m    = m+1
   CALL GETARG(m,arg)
   IF (arg(1:2) == '-n') THEN
      m = m+1
      CALL GETARG(m,arg)
      READ(arg,*)n
   ELSE IF (arg(1:2) == '-p') THEN
      m = m+1
      CALL GETARG(m,arg)
      READ(arg,*)perturb
   ELSE IF (arg(1:2) == '-o') THEN
      m = m+1
      CALL GETARG(m,arg)
      fileName = arg
   END IF
END DO


! compute dimensions of mesh and allocate data
orderM    = 3 !cubic elements
nSurfElem = n
nSurfNode = n*orderM
nBndNode  = 0
nCompBd   = 1

ALLOCATE( &
     order(nSurfElem), &
     surfElem(orderM+1,nSurfElem), &
     surfElemTag(nSurfElem), &
     surfX(2,nSurfNode), &
     bndNode(nBndNode), &
     bndNodeTag(nBndNode), &
     bndNodeNormal(2,nBndNode))


! set surface tags and boundary information
order              = orderM
surfElemTag        = 0


! generate grid data
km = 0
k  = km+orderM
DO n=1,nSurfElem
   surfElem(1,n) = km
   surfElem(2,n) = k
   DO i=1,orderM-1
      surfElem(2+i,n) = km+i
   END DO
   km            = k
   k             = k+orderM
   IF (k == nSurfNode) k = 0
END DO

pi = 4.*ATAN(1.)
dt = 2.*pi/REAL(nSurfElem)
m  = 1
DO n=1,nSurfNode,orderM
   t          =-REAL(m-1)*dt
   IF (perturb == 1) t = t+2.*dt/REAL(orderM)*(RAND(0)-.5)
   surfX(1,n) = .5*COS(t)
   surfX(2,n) = .5*SIN(t)
   m          = m+1
   !write(*,*)'knots ',n,surfX(1,n),surfX(2,n)
END DO
DO n=1,nSurfNode,orderM
   n0         = n
   n1         = n+orderM
   IF (n1 > nSurfNode) n1 = 1
   t0         = ATAN2(surfX(2,n0),surfX(1,n0))
   t1         = ATAN2(surfX(2,n1),surfX(1,n1))
   IF (t0 <= 0. .AND. t1 > 0.) t0 = t0+2.*pi
   dt         =(t1-t0)/REAL(orderM)
   DO m=1,orderM-1
      t          = t0+REAL(m)*dt
      surfX(1,n+m) = .5*COS(t)
      surfX(2,n+m) = .5*SIN(t)
      !write(*,*)'inter ',n,n+m,t0,t1,dt!surfX(1,n+m),surfX(2,n+m)!,
   END DO
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


END PROGRAM circle
