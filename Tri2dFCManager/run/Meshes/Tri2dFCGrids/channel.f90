PROGRAM channel

IMPLICIT NONE

INTEGER :: n,m,k,nx,narg,i,j,iu,perturb
INTEGER :: nNode,nTri,nEdgeBd,nCompBd
INTEGER,ALLOCATABLE,DIMENSION(:  ) :: order
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: num,edgeBd,tri
REAL :: dx,dy,x0,y0,dxp,dyp,ds,al,pi,m1,m2
REAL,PARAMETER :: fact = .5
REAL,ALLOCATABLE,DIMENSION(:,:) :: x
CHARACTER(80) :: meshfile,arg


! read inputs
narg = iargc()
m    = 0
DO WHILE (m < narg)
   m    = m+1
   CALL GETARG(m,arg)
   IF (arg(1:2) == '-n') THEN
      m = m+1
      CALL GETARG(m,arg)
      READ(arg,*)nx
   ELSE IF (arg(1:2) == '-p') THEN
      m = m+1
      CALL GETARG(m,arg)
      READ(arg,*)perturb
   ELSE IF (arg(1:2) == '-o') THEN
      m = m+1
      CALL GETARG(m,arg)
      meshfile = arg
   END IF
END DO


! determine grid dimensions
nTri    = 2*nx*nx
nNode   =(nx+1)*(nx+1)
nCompBd = 4
nEdgeBd = 4*nx


! allocate data
ALLOCATE(order(nTri))
ALLOCATE(tri(3,nTri))
ALLOCATE(x(2,nNode))
ALLOCATE(edgeBd(3,nEdgeBd))


! fill in grid data

order = 1 !all linear elements

ALLOCATE(num(nx+1,nx+1))
k = 0
DO j=1,nx+1
DO i=1,nx+1
   num(i,j) = k
   k        = k+1
END DO
END DO

k = 0
DO j=1,nx
DO i=1,nx
   k        = k+1
   tri(1,k) = num(i  ,j  )
   tri(2,k) = num(i+1,j  )
   tri(3,k) = num(i+1,j+1)
   k        = k+1
   tri(1,k) = num(i  ,j  )
   tri(2,k) = num(i+1,j+1)
   tri(3,k) = num(i  ,j+1)
END DO
END DO

k = 0
j = 1
DO i=1,nx
   k           = k+1
   edgeBd(1,k) = num(i  ,j  )
   edgeBd(2,k) = num(i+1,j  )
   edgeBd(3,k) = 0
END DO
j = nx+1
DO i=1,nx
   k           = k+1
   edgeBd(1,k) = num(i+1,j  )
   edgeBd(2,k) = num(i  ,j  )
   edgeBd(3,k) = 1
END DO
i = 1
DO j=1,nx
   k           = k+1
   edgeBd(1,k) = num(i  ,j+1)
   edgeBd(2,k) = num(i  ,j  )
   edgeBd(3,k) = 2
END DO
i = nx+1
DO j=1,nx
   k           = k+1
   edgeBd(1,k) = num(i  ,j  )
   edgeBd(2,k) = num(i  ,j+1)
   edgeBd(3,k) = 3
END DO

k  = 0
x0 = 0.
y0 = 0.
dx = 1./REAL(nx)
dy = dx
pi = 4.*ATAN(1.)
DO j=1,nx+1
DO i=1,nx+1
   k      = k+1
   x(1,k) = x0+REAL(i-1)*dx
   x(2,k) = y0+REAL(j-1)*dy
END DO
END DO
IF (perturb > 0) THEN
   DO j=2,nx
   DO i=2,nx
      CALL RANDOM_NUMBER(m1)
      CALL RANDOM_NUMBER(m2)
      k      = num(i,j)+1
      ds     = .5*fact*dx*m1
      al     = 2.*pi*m2
      dxp    = ds*COS(al)
      dyp    = ds*SIN(al)
      x(1,k) = x(1,k)+dxp
      x(2,k) = x(2,k)+dyp
   END DO
   END DO
END IF


! write to mesh file
OPEN(iu,FILE=meshfile,STATUS='replace',FORM='formatted')
WRITE(iu,*)nTri,nNode,nCompBd,nEdgeBd
DO n=1,nTri
   WRITE(iu,*)order(n),tri(:,n)
END DO
DO n=1,nNode
   WRITE(iu,*)x(:,n)
END DO
DO n=1,nEdgeBd
   WRITE(iu,*)edgeBd(:,n)
END DO
CLOSE(iu)


! deallocate data
DEALLOCATE(order)
DEALLOCATE(tri)
DEALLOCATE(x)
DEALLOCATE(edgeBd)

WRITE(*,*)
WRITE(*,*)'*** successfully wrote ',TRIM(meshfile),' ***'
WRITE(*,*)


END PROGRAM channel
