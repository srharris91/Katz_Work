PROGRAM channel

IMPLICIT NONE

INTEGER :: n,m,k,nx,narg,i,j,iu,perturb,i1,j1,ord,nne,n1,n2,n3
INTEGER :: nNode,nTri,nEdgeBd,nCompBd
INTEGER,ALLOCATABLE,DIMENSION(:  ) :: order
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: num,edgeBd,tri
REAL,DIMENSION(2) :: xl1,xl2,xl3,xInt
REAL :: dx,dy,x0,y0,dxp,dyp,ds,al,pi,m1,m2,ri,si,r,lA,lB,gamma
REAL,PARAMETER :: fact = .5
REAL,ALLOCATABLE,DIMENSION(:) :: dsn
REAL,ALLOCATABLE,DIMENSION(:,:) :: x,rs,xsp
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
ord     = 3 ! cubic elements
nne     =(ord+2)*(ord+1)/2
nTri    = 2*nx*nx
nNode   =(ord*nx+1)*(ord*nx+1)
nCompBd = 4
nEdgeBd = 4*nx


! allocate data
ALLOCATE(order(nTri))
ALLOCATE(tri(nne,nTri))
ALLOCATE(x(2,nNode))
ALLOCATE(edgeBd(3,nEdgeBd))


! fill in grid data

order = ord !all cubic elements

ALLOCATE(num(ord*nx+1,ord*nx+1))
k = 0
DO j=1,ord*nx+1
DO i=1,ord*nx+1
   num(i,j) = k
   k        = k+1
END DO
END DO

k = 0
DO j=1,nx
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
   edgeBd(3,k) = 0
END DO
j  = nx+1
j1 = ord*(j-1)+1
DO i=1,nx
   i1          = ord*(i-1)+1
   k           = k+1
   edgeBd(1,k) = num(i1+3,j1  )
   edgeBd(2,k) = num(i1  ,j1  )
   edgeBd(3,k) = 1
END DO
i  = 1
i1 = ord*(i-1)+1
DO j=1,nx
   j1          = ord*(j-1)+1
   k           = k+1
   edgeBd(1,k) = num(i1  ,j1+3)
   edgeBd(2,k) = num(i1  ,j1  )
   edgeBd(3,k) = 2
END DO
i  = nx+1
i1 = ord*(i-1)+1
DO j=1,nx
   j1          = ord*(j-1)+1
   k           = k+1
   edgeBd(1,k) = num(i1  ,j1  )
   edgeBd(2,k) = num(i1  ,j1+3)
   edgeBd(3,k) = 3
END DO


! compute the verticies of the triangles
x0 = 0.
y0 = 0.
dx = 1./REAL(nx)
dy = dx
pi = 4.*ATAN(1.)
DO j=1,nx+1
DO i=1,nx+1
   i1     = ord*(i-1)+1
   j1     = ord*(j-1)+1
   k      = num(i1,j1)+1
   x(1,k) = x0+REAL(i-1)*dx
   x(2,k) = y0+REAL(j-1)*dy
END DO
END DO
IF (perturb > 0) THEN ! perturb the verticies
   DO j=2,nx
   DO i=2,nx
      CALL RANDOM_NUMBER(m1)
      CALL RANDOM_NUMBER(m2)
      i1     = ord*(i-1)+1
      j1     = ord*(j-1)+1
      k      = num(i1,j1)+1
      ds     = .5*fact*dx*m1
      al     = 2.*pi*m2
      dxp    = ds*COS(al)
      dyp    = ds*SIN(al)
      x(1,k) = x(1,k)+dxp
      x(2,k) = x(2,k)+dyp
   END DO
   END DO
END IF


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
ALLOCATE(xsp(2,nne))
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
      xsp(:,i) =((-3.*ri+2.-SQRT(3.)*si)*xl1 &
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
      xsp(:,i) = xsp(:,i)+xInt*lA*lB

      ! edge 2
      r        = dsn(n3)/dsn(n2)
      gamma    =(r-1)/(r+1)
      lA       =( 3.*ri+2.-SQRT(3.)*si)/6.
      lB       =(2.+2.*SQRT(3.)*si)/6.
      xInt     = 2.*gamma*(xl2-xl3)
      xsp(:,i) = xsp(:,i)+xInt*lA*lB

      ! edge 3
      r        = dsn(n1)/dsn(n3)
      gamma    =(r-1)/(r+1)
      lA       =(2.+2.*SQRT(3.)*si)/6.
      lB       =(-3.*ri+2.-SQRT(3.)*si)/6.
      xInt     = 2.*gamma*(xl3-xl1)
      xsp(:,i) = xsp(:,i)+xInt*lA*lB
   END DO

   DO i=4,nne
      k      = tri(i,n)+1
      x(:,k) = xsp(:,i)
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
DEALLOCATE(rs)
DEALLOCATE(xsp)

WRITE(*,*)
WRITE(*,*)'*** successfully wrote ',TRIM(meshfile),' ***'
WRITE(*,*)


END PROGRAM channel
