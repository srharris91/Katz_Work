PROGRAM circle

IMPLICIT NONE

INTEGER :: n,m,k,narg,i,j,iu,level,na,npCircFine,npRadFine,skip,npCirc,npRad
INTEGER :: nNode,nTri,nEdgeBd,nCompBd,isf,nmax
INTEGER,ALLOCATABLE,DIMENSION(:  ) :: order
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: map,edgeBd,tri
REAL :: dang,ang,pi,x0,x1,dx0,dx1,rmax,dx0a,dx1a
REAL,ALLOCATABLE,DIMENSION(:  ) :: xn,xnT
REAL,ALLOCATABLE,DIMENSION(:,:) :: x
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
pi         = 4.*ATAN(1.)
npCircFine = 512
npRadFine  = 129
skip       = 2**(level-1)
npCirc     = npCircFine/skip
npRad      =(npRadFine-1)/skip+1
nTri       = npCirc*(npRad-1)*2
nNode      = npCirc* npRad
nCompBd    = 2
nEdgeBd    = npCirc*2


! determine stretching function normal to cylinder
isf  = 2
x0   = .5
x1   = 20.
!dx0  = .0064
dx0  = .00064
!dx0  = .000064
dx1  = 1.
rmax = 1.1
nmax = npCircFine*2
dx0a = 1.
ALLOCATE(xnT(nmax))
CALL sfuns(isf,x0,x1,dx0,dx1,npRadFine,rmax,nmax,dx0a,dx1a,xnT)
k = 0
ALLOCATE(xn(npRad))
DO j=1,npRadFine,skip
   k     = k+1
   xn(k) = xnT(j)
END DO
DEALLOCATE(xnT)


! allocate data
ALLOCATE(order(nTri))
ALLOCATE(tri(3,nTri))
ALLOCATE(x(2,nNode))
ALLOCATE(edgeBd(3,nEdgeBd))


! compute finest level coordinates and structured mappings
ALLOCATE(map(npCirc,npRad))

dang = 2.*pi/REAL(npCirc)
k    = 0
DO j=1,npRad
DO n=1,npCirc
   k        = k+1
   ang      = REAL(n-1)*dang
   x(1,k)   = xn(j)*COS(ang)
   x(2,k)   = xn(j)*SIN(ang)
   map(n,j) = k-1
END DO
END DO


! fill in connectivity data
order = 1
k     = 0
DO j=1,npRad-1
DO n=1,npCirc/2
   na       = n+1
   k        = k+1
   tri(1,k) = map(na,j  ) !surface edge @j=1
   tri(2,k) = map(n ,j  ) !surface edge @j=1
   tri(3,k) = map(n ,j+1)
   k        = k+1
   tri(1,k) = map(na,j  )
   tri(2,k) = map(n ,j+1)
   tri(3,k) = map(na,j+1)
END DO
DO n=npCirc/2+1,npCirc
   na       = n+1
   IF (n==npCirc) na = 1
   k        = k+1
   tri(1,k) = map(n ,j+1) !surface edge @j=npRad
   tri(2,k) = map(na,j+1) !surface edge @j=npRad
   tri(3,k) = map(n ,j  )
   k        = k+1
   tri(1,k) = map(n ,j  )
   tri(2,k) = map(na,j+1)
   tri(3,k) = map(na,j  )
END DO
END DO


! create boundary edge array
k = 0
DO n=1,npCirc
   k            = k+1
   na           = n+1
   IF (n==npCirc) na = 1
   edgeBd(1,k)  = map(na,1)
   edgeBd(2,k)  = map(n ,1)
   edgeBd(3,k)  = 0
   k            = k+1
   edgeBd(1,k)  = map(n ,npRad)
   edgeBd(2,k)  = map(na,npRad)
   edgeBd(3,k)  = 1
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
DEALLOCATE(map)
DEALLOCATE(xn)

WRITE(*,*)
WRITE(*,*)'*** successfully wrote ',TRIM(meshfile),' ***'
WRITE(*,*)


END PROGRAM circle
