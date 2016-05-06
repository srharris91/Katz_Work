PROGRAM vtk2fc

IMPLICIT NONE

INTEGER :: n,m,narg,iu,nNode,nTri,nEdgeBd,nCompBd
INTEGER,ALLOCATABLE,DIMENSION(:  ) :: order
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: edgeBd,tri
REAL :: al
REAL,ALLOCATABLE,DIMENSION(:,:) :: x
CHARACTER(80) :: vtkFile,fcFile,arg



! read inputs
narg = iargc()
m    = 0
DO WHILE (m < narg)
   m = m+1
   CALL GETARG(m,arg)
   IF (arg(1:4) == '-vtk') THEN
      m = m+1
      CALL GETARG(m,arg)
      vtkFile = arg
   ELSE IF (arg(1:3) == '-fc') THEN
      m = m+1
      CALL GETARG(m,arg)
      fcFile = arg
   END IF
END DO


iu = 11
OPEN(iu,FILE=vtkFile,STATUS='old',FORM='formatted')
READ(iu,*)
READ(iu,*)nNode,nTri,nEdgeBd,nCompBd
ALLOCATE(order(nTri))
ALLOCATE(tri(3,nTri))
ALLOCATE(x(2,nNode))
ALLOCATE(edgeBd(3,nEdgeBd))
READ(iu,*)
READ(iu,*)
READ(iu,*)
DO n=1,nNode
   READ(iu,*)x(:,n),al
END DO
READ(iu,*)
DO n=1,nTri
   READ(iu,*)m,tri(:,n)
END DO
DO n=1,nEdgeBd
   READ(iu,*)m,edgeBd(1:2,n) !delaundo gives backward
END DO
READ(iu,*)
DO n=1,nTri
   READ(iu,*)m
END DO
DO n=1,nEdgeBd
   READ(iu,*)m
END DO
READ(iu,*)
READ(iu,*)
READ(iu,*)
DO n=1,nTri
   READ(iu,*)m
END DO
DO n=1,nEdgeBd
   READ(iu,*)edgeBd(3,n)
END DO
CLOSE(iu)


order = 1 !all linear elements


! write to mesh file
OPEN(iu,FILE=fcFile,STATUS='replace',FORM='formatted')
WRITE(iu,*)nTri,nNode,nCompBd,nEdgeBd
DO n=1,nTri
   WRITE(iu,*)order(n),tri(:,n)
END DO
DO n=1,nNode
   WRITE(iu,*)x(:,n)
END DO
DO n=1,nEdgeBd
   WRITE(iu,*)edgeBd(1,n),edgeBd(2,n),edgeBd(3,n)-1
END DO
CLOSE(iu)


! deallocate data
DEALLOCATE(order)
DEALLOCATE(tri)
DEALLOCATE(x)
DEALLOCATE(edgeBd)


WRITE(*,*)
WRITE(*,*)'*** successfully wrote ',TRIM(fcFile),' ***'
WRITE(*,*)


END PROGRAM vtk2fc
