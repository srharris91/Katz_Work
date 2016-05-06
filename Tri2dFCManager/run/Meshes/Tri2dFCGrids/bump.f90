PROGRAM bump

IMPLICIT NONE

INTEGER :: n,m,k,narg,i,j,iu,level,npBumpFine,npNormFine,skip,npBump,npNorm
INTEGER :: nNode,nTri,nEdgeBd,nCompBd,isf,nmax
INTEGER,ALLOCATABLE,DIMENSION(:  ) :: order
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: map,edgeBd,tri
REAL :: dang,ang,pi,x0,x1,dx0,dx1,rmax,dx0a,dx1a,dy,yb
REAL,ALLOCATABLE,DIMENSION(:  ) :: xb,xn,xnT
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
   IF (arg(1:2) == '-n') THEN !"plate -n 2" this gives the signified level listed above
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
npBumpFine = 513
npNormFine = 129
skip       = 2**(level-1)
npBump     =(npBumpFine-1)/skip+1
npNorm     =(npNormFine-1)/skip+1
nTri       =(npBump-1)*(npNorm-1)*2
nNode      = npBump*    npNorm
nCompBd    = 5
nEdgeBd    =(npBump-1)*2+(npNorm-1)*2


! determine stretching function along the bump
isf  = 2
x0   = 0.
x1   = 1.
dx0  = .064
dx1  = .0002
rmax = 1.1
nmax = npBumpFine*2
dx0a = 1.
n    =(npBumpFine-1)/2+1
ALLOCATE(xnT(nmax))
CALL sfuns(isf,x0,x1,dx0,dx1,n,rmax,nmax,dx0a,dx1a,xnT)
k = 0
ALLOCATE(xb(npBump))
DO j=1,(npBumpFine-1)/2+1,skip
   k     = k+1
   xb(k) = 5.*(xnT(j)-1.)
END DO
k = k-1
DO j=(npBumpFine-1)/2+1,1,-skip
   k     = k+1
   xb(k) =-5.*(xnT(j)-1.)
END DO
DEALLOCATE(xnT)


! determine stretching function normal to bump
isf  = 2
x0   = 0.
x1   = 1.
dx0  = .0008 ! wall Spacing
!dx0  = .00008
!dx0  = .000008
!~ dx1  = .064
dx1  = 0. !Shauns
rmax = 1.1
nmax = npNormFine*2
dx0a = 1.
ALLOCATE(xnT(nmax))
CALL sfuns(isf,x0,x1,dx0,dx1,npNormFine,rmax,nmax,dx0a,dx1a,xnT)
k = 0
ALLOCATE(xn(npNorm))
DO j=1,npNormFine,skip
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
ALLOCATE(map(npBump,npNorm))

k = 0
DO j=1,npNorm
DO n=1,npBump
   k        = k+1
   x(1,k)   = xb(n)
   IF (x(1,k) < -.5 .OR. x(1,k) > .5) THEN
      yb    = 0.
   ELSE
      yb    = 0.05*(SIN(pi*(.5+x(1,k))))**4
   END IF
   dy = 5.-yb
   x(2,k)   = yb+dy*xn(j)
   map(n,j) = k-1
END DO
END DO


! fill in connectivity data
order = 1
k     = 0
DO j=1,npNorm-1
DO n=1,(npBump-1)/2
   k        = k+1
   tri(1,k) = map(n  ,j  ) !surface edge @j=1
   tri(2,k) = map(n+1,j  ) !surface edge @j=1
   tri(3,k) = map(n+1,j+1)
   k        = k+1
   tri(1,k) = map(n  ,j  )
   tri(2,k) = map(n+1,j+1)
   tri(3,k) = map(n  ,j+1)
END DO
DO n=(npBump-1)/2+1,npBump-1
   k        = k+1
   tri(1,k) = map(n  ,j  ) !surface edge @j=npNorm
   tri(2,k) = map(n+1,j  ) !surface edge @j=npNorm
   tri(3,k) = map(n  ,j+1)
   k        = k+1
   tri(1,k) = map(n+1,j  )
   tri(2,k) = map(n+1,j+1)
   tri(3,k) = map(n  ,j+1)
END DO
END DO


! create boundary edge array
k = 0
DO n=1,npBump-1
   k            = k+1
   edgeBd(1,k)  = map(n  ,1)
   edgeBd(2,k)  = map(n+1,1)
   IF (x(1,k) < -.5) THEN !Shauns change, set inviscid before bump
      edgeBd(3,k)  = 0
   ELSE
      edgeBd(3,k)  = 1 !viscous on and after bump
   END IF
!~    edgeBd(3,k)  = 0 !Shauns change, uncomment this line and comment if statement above to undo Shauns change
END DO
DO n=1,npBump-1
   k            = k+1
   edgeBd(1,k)  = map(n+1,npNorm)
   edgeBd(2,k)  = map(n  ,npNorm)
!~    edgeBd(3,k)  = 1
   edgeBd(3,k)  = 2 !Shauns change, switch with commented line above
END DO
DO n=1,npNorm-1
   k            = k+1
   edgeBd(1,k)  = map(1,n+1)
   edgeBd(2,k)  = map(1,n  )
!~    edgeBd(3,k)  = 2
   edgeBd(3,k)  = 3 !Shauns change, switch with commented line above
END DO
DO n=1,npNorm-1
   k            = k+1
   edgeBd(1,k)  = map(npBump,n  )
   edgeBd(2,k)  = map(npBump,n+1)
!~    edgeBd(3,k)  = 3
   edgeBd(3,k)  = 4 !Shauns change, switch with commented line above
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
DEALLOCATE(xb)
DEALLOCATE(xn)

WRITE(*,*)
WRITE(*,*)'*** successfully wrote ',TRIM(meshfile),' ***'
WRITE(*,*)


END PROGRAM bump
