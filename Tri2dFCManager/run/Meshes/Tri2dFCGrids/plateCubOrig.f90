PROGRAM plate

IMPLICIT NONE

INTEGER :: n,m,k,nx,ny,narg,i,j,iu,perturb,nNode,nTri,nEdgeBd,nCompBd,npl,nen, &
           isf,nmax,n1,n2,n3,ord,nne,i1,j1
INTEGER,ALLOCATABLE,DIMENSION(:  ) :: order
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: num,edgeBd,tri
REAL :: dx,dy,x0,y0,dxp,dyp,ds,al,pi,m1,m2
REAL :: x1,dx0,dx1,rmax,dx0a,dx1a,ri,si,r,lA,lB,gamma
REAL,DIMENSION(2) :: xl1,xl2,xl3,xInt
REAL,ALLOCATABLE,DIMENSION(:,:) :: x,rs,xel
REAL,ALLOCATABLE,DIMENSION(:  ) :: xsb,xsp,xsn,dsn
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
x0   =-2.
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
x1   = 4.
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
x1   = 4.
dx0  = .002
dx1  = .48
n    = nen+npl+1
!n    =(nen+npl)/4+1
ny   = n-1
rmax = 1.1
nmax = 200
dx0a = 1.
ALLOCATE(xsn(nmax))
CALL sfuns(isf,x0,x1,dx0,dx1,n,rmax,nmax,dx0a,dx1a,xsn)


! determine grid dimensions
ord     = 3 ! cubic elements
nne     =(ord+2)*(ord+1)/2
nx      = nen+npl
nTri    = 2*nx*ny
nNode   =(ord*nx+1)*(ord*ny+1)
nCompBd = 5
nEdgeBd = 2*(nx+ny)


! allocate data
ALLOCATE(order(nTri))
ALLOCATE(tri(nne,nTri))
ALLOCATE(x(2,nNode))
ALLOCATE(edgeBd(3,nEdgeBd))


! fill in grid data

order = ord !all linear elements

ALLOCATE(num(ord*nx+1,ord*ny+1))
k = 0
DO j=1,ord*ny+1
DO i=1,ord*nx+1
   num(i,j) = k
   k        = k+1
END DO
END DO

k = 0
DO j=1,ny
DO i=1,nx
   i1        = ord*(i-1)+1
   j1        = ord*(j-1)+1
   k         = k+1
   tri(1 ,k) = num(i1  ,j1  )
   tri(2 ,k) = num(i1+3,j1  )
   tri(3 ,k) = num(i1+3,j1+3)
   tri(4 ,k) = num(i1+1,j1  )
   tri(5 ,k) = num(i1+2,j1  )
   tri(6 ,k) = num(i1+3,j1+1)
   tri(7 ,k) = num(i1+3,j1+2)
   tri(8 ,k) = num(i1+2,j1+2)
   tri(9 ,k) = num(i1+1,j1+1)
   tri(10,k) = num(i1+2,j1+1)
   k         = k+1
   tri(1 ,k) = num(i1  ,j1  )
   tri(2 ,k) = num(i1+3,j1+3)
   tri(3 ,k) = num(i1  ,j1+3)
   tri(4 ,k) = num(i1+1,j1+1)
   tri(5 ,k) = num(i1+2,j1+2)
   tri(6 ,k) = num(i1+2,j1+3)
   tri(7 ,k) = num(i1+1,j1+3)
   tri(8 ,k) = num(i1  ,j1+2)
   tri(9 ,k) = num(i1  ,j1+1)
   tri(10,k) = num(i1+1,j1+2)
END DO
END DO

k  = 0
j  = 1
j1 = ord*(j-1)+1
DO i=1,nx
   i1          = ord*(i-1)+1
   k           = k+1
   edgeBd(1,k) = num(i1  ,j1  )
   edgeBd(2,k) = num(i1+3,j1  )
   IF (i <= nen) THEN
      edgeBd(3,k) = 0
   ELSE
      edgeBd(3,k) = 1
   END IF
END DO
j  = ny+1
j1 = ord*(j-1)+1
DO i=1,nx
   i1          = ord*(i-1)+1
   k           = k+1
   edgeBd(1,k) = num(i1+3,j1  )
   edgeBd(2,k) = num(i1  ,j1  )
   edgeBd(3,k) = 2
END DO
i  = 1
i1 = ord*(i-1)+1
DO j=1,ny
   j1          = ord*(j-1)+1
   k           = k+1
   edgeBd(1,k) = num(i1  ,j1+3)
   edgeBd(2,k) = num(i1  ,j1  )
   edgeBd(3,k) = 3
END DO
i  = nx+1
i1 = ord*(i-1)+1
DO j=1,ny
   j1          = ord*(j-1)+1
   k           = k+1
   edgeBd(1,k) = num(i1  ,j1  )
   edgeBd(2,k) = num(i1  ,j1+3)
   edgeBd(3,k) = 4
END DO

! compute the verticies of the triangles
DO j=1,ny+1
DO i=1,nx+1
   i1 = ord*(i-1)+1
   j1 = ord*(j-1)+1
   k  = num(i1,j1)+1
   IF (i <= nen) THEN
      x(1,k) = xsb(i)
   ELSE
      x(1,k) = xsp(i-nen)
   END IF
   x(2,k) = xsn(j)
END DO
END DO


! estimate the grid spacing at each vertex
ALLOCATE(dsn(nNode))
dsn = 0.
DO n=1,nTri
   n1      = tri(1,n)+1
   n2      = tri(2,n)+1
   n3      = tri(3,n)+1
   xl1     = x(:,n1)
   xl2     = x(:,n2)
   xl3     = x(:,n3)
   ds      = .5*((xl2(1)-xl1(1))*(xl3(2)-xl1(2)) &
                -(xl2(2)-xl1(2))*(xl3(1)-xl1(1)));
   dsn(n1) = dsn(n1)+ds
   dsn(n2) = dsn(n2)+ds
   dsn(n3) = dsn(n3)+ds
END DO
DO n=1,nNode
   dsn(n) = SQRT(dsn(n))
END DO


! compute interior vertices of each triangle
ALLOCATE(rs(2,nne))
ALLOCATE(xel(2,nne))
CALL equiTriPts(ord,nne,rs)
DO n=1,nTri
   n1  = tri(1,n)+1
   n2  = tri(2,n)+1
   n3  = tri(3,n)+1
   xl1 = x(:,n1)
   xl2 = x(:,n2)
   xl3 = x(:,n3)
   
   ! straight portion
   DO i=1,nne
      ri = rs(1,i)
      si = rs(2,i)
      xel(:,i) =((-3.*ri+2.-SQRT(3.)*si)*xl1 &
               + ( 3.*ri+2.-SQRT(3.)*si)*xl2 &
               + (2.+2.*SQRT(3.)*si)*xl3)/6.
   END DO

   ! stretched portion
   DO i=1,nne
      ri       = rs(1,i)
      si       = rs(2,i)

      ! edge 1
      r        = dsn(n2)/dsn(n1)
      gamma    =(r-1)/(r+1)
      lA       =(-3.*ri+2.-SQRT(3.)*si)/6.
      lB       =( 3.*ri+2.-SQRT(3.)*si)/6.
      xInt     = 2.*gamma*(xl1-xl2)
      xel(:,i) = xel(:,i)+xInt*lA*lB

      ! edge 2
      r        = dsn(n3)/dsn(n2)
      gamma    =(r-1)/(r+1)
      lA       =( 3.*ri+2.-SQRT(3.)*si)/6.
      lB       =(2.+2.*SQRT(3.)*si)/6.
      xInt     = 2.*gamma*(xl2-xl3)
      xel(:,i) = xel(:,i)+xInt*lA*lB

      ! edge 3
      r        = dsn(n1)/dsn(n3)
      gamma    =(r-1)/(r+1)
      lA       =(2.+2.*SQRT(3.)*si)/6.
      lB       =(-3.*ri+2.-SQRT(3.)*si)/6.
      xInt     = 2.*gamma*(xl3-xl1)
      xel(:,i) = xel(:,i)+xInt*lA*lB
   END DO

   DO i=4,nne
      k      = tri(i,n)+1
      x(:,k) = xel(:,i)
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
