PROGRAM plate

IMPLICIT NONE

INTEGER :: n,m,k,nx,ny,narg,i,j,iu,perturb,nNode,nTri,nEdgeBd,nCompBd,npl,nen, &
           isf,nmax
INTEGER,ALLOCATABLE,DIMENSION(:  ) :: order
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: num,edgeBd,tri
REAL :: dx,dy,x0,y0,dxp,dyp,ds,al,pi,m1,m2
REAL :: x1,dx0,dx1,rmax,dx0a,dx1a
REAL,PARAMETER :: fact = .5
REAL,ALLOCATABLE,DIMENSION(:,:) :: x
REAL,ALLOCATABLE,DIMENSION(:  ) :: xsb,xsp,xsn
CHARACTER(80) :: meshFile,arg



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
   ELSE IF (arg(1:2) == '-o') THEN
      m = m+1
      CALL GETARG(m,arg)
      meshFile = arg
   END IF
END DO


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
dx0  = .01
dx1  = .48
n    = npl+1
rmax = 1.1
nmax = 200
dx0a = 1.
ALLOCATE(xsp(nmax))
CALL sfuns(isf,x0,x1,dx0,dx1,n,rmax,nmax,dx0a,dx1a,xsp)


! Compute 1-d stretching function normal to plate
isf  = 2
x0   = 0.
x1   = 1.
dx0  = .001!.002 initial wall spacing
!~ dx1  = 1.!.48
!~ dx1  = .61
dx1  = 0. ! final spacing from outer wall
!~ n    = nen+npl+1
!~ n    =(nen+npl)/4+1
!~ n    =(nen+npl)/5+3 !Specified by Shaun Feb. 2015
n=33 !Specified by Shaun Feb.-Mar. 2015, specific for a 100 npl user input
!~ n=63 !Specified by Shaun Mar. 2015, specific for a 180 npl user input
ny   = n-1
rmax = 1.1
nmax = 200
dx0a = 1.
ALLOCATE(xsn(nmax))
CALL sfuns(isf,x0,x1,dx0,dx1,n,rmax,nmax,dx0a,dx1a,xsn)


! determine grid dimensions
nx     = nen+npl
nTri   = 2*nx*ny
nNode  =(nx+1)*(ny+1)
nCompBd = 5
nEdgeBd = 2*(nx+ny)


! allocate data
ALLOCATE(order(nTri))
ALLOCATE(tri(3,nTri))
ALLOCATE(x(2,nNode))
ALLOCATE(edgeBd(3,nEdgeBd))


! fill in grid data

order = 1 !all linear elements

ALLOCATE(num(nx+1,ny+1))
k = 0
DO j=1,ny+1
DO i=1,nx+1
   num(i,j) = k
   k        = k+1
END DO
END DO

k = 0
DO j=1,ny
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
   IF (i <= nen) THEN
      edgeBd(3,k) = 0
   ELSE
      edgeBd(3,k) = 1
   END IF
END DO
j = ny+1
DO i=1,nx
   k           = k+1
   edgeBd(1,k) = num(i+1,j  )
   edgeBd(2,k) = num(i  ,j  )
   edgeBd(3,k) = 2
END DO
i = 1
DO j=1,ny
   k           = k+1
   edgeBd(1,k) = num(i  ,j+1)
   edgeBd(2,k) = num(i  ,j  )
   edgeBd(3,k) = 3
END DO
i = nx+1
DO j=1,ny
   k           = k+1
   edgeBd(1,k) = num(i  ,j  )
   edgeBd(2,k) = num(i  ,j+1)
   edgeBd(3,k) = 4
END DO

k = 0
DO j=1,ny+1
DO i=1,nx+1
   k = k+1
   IF (i <= nen) THEN
      x(1,k) = xsb(i)
   ELSE
      x(1,k) = xsp(i-nen)
   END IF
   x(2,k) = xsn(j)
END DO
END DO


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
WRITE(*,*)'*** successfully wrote ',TRIM(meshFile),' ***'
WRITE(*,*)


END PROGRAM plate
