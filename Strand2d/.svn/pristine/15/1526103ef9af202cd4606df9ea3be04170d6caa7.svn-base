!> \brief
!! This subroutine forms the cells surrounding points structures.
!! \param pid Process ID.
!! \param nFaces Number of faces.
!! \param nGfaces Number of ghost faces on the partition.
!! \param nNodes Number of Nodes.
!! \param nGnodes Number of ghost nodes on the partition.
!! \param ncsp1 Dimension of csp1.
!! \param face List of faces.
!! \param csp1 Cells surrounding points index array.
!! \param csp2 Cells surrounding points data array.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-9-20
!! \par Further Documentation:
!! \par Source code:
!! \include src/Numerics/formCsp.F90


SUBROUTINE formcsp( &
     pid, &
     nFaces, &
     nGfaces, &
     nNodes, &
     nGnodes, &
     ncsp1, &
     face, &
     csp1, &
     csp2)


IMPLICIT NONE

INTEGER,INTENT(IN   ) :: &
     pid, &
     nFaces, &
     nGfaces, &
     nNodes, &
     nGnodes, &
     ncsp1
INTEGER,INTENT(INOUT),DIMENSION(2,nFaces  ) :: face
INTEGER,INTENT(  OUT),DIMENSION(  ncsp1   ) :: csp1
INTEGER,INTENT(  OUT),DIMENSION(  nNodes+1) :: csp2
INTEGER :: n,k,j,i,ii,jj,deg,m,mm
INTEGER,ALLOCATABLE,DIMENSION(:) :: flag,csp1E,csp2E,csp1EX,csp2EX
INTEGER :: ncsp1E,ncsp1EX


! make arrays 1-based
face = face+1

! Form initial arrays, which may have single cell stencils
ALLOCATE(csp2E(nNodes+1))
csp2E = 0
DO n=1,nFaces-nGfaces
DO k=1,2
   j = face(k,n)+1
   csp2E(j) = csp2E(j)+1
END DO
END DO

IF (MINVAL(csp2E(2:nNodes+1)) == 0) THEN !test for hanging nodes
   WRITE(*,*)'*** hanging node detected in grid in edgeExtract ***'
   STOP
END IF

DO n=2,nNodes+1 !add csp2E to previous values
   csp2E(n) = csp2E(n)+csp2E(n-1)
END DO
ncsp1E = csp2E(nNodes+1)
IF (ncsp1E > ncsp1) THEN
   WRITE(*,*)'*** nCsp1E exceeds dimensions in gridGeomFormCsp ***'
   STOP
END IF
ALLOCATE(csp1E(ncsp1E))

DO n=1,nFaces-nGfaces
DO k=1,2
   j        = face(k,n)
   i        = csp2E(j)+1
   csp2E(j) = i
   csp1E(i) = n
END DO
END DO

DO n=nNodes+1,2,-1 ! finalize pointer array
   csp2E(n) = csp2E(n-1)
END DO
csp2E(1) = 0


! check for single cell stencil and add more cells if necessary

ALLOCATE(flag(nFaces), &
         csp2EX(nNodes+1))
flag   = 0
csp2EX = 0
DO n=1,nNodes-nGnodes
   deg = csp2E(n+1)-csp2E(n)
IF (deg <= 1) THEN
DO mm=csp2E(n)+1,csp2E(n+1)
   j       = csp1E(mm)
   flag(j) = n
END DO
DO mm=csp2E(n)+1,csp2E(n+1)
   j       = csp1E(mm)
   DO k=1,2
      jj = face(k,j)
   DO m=csp2E(jj)+1,csp2E(jj+1)
      i  = csp1E(m)
      IF (flag(i) /= n) THEN
         flag(i)     = n
         csp2EX(n+1) = csp2EX(n+1)+1
      END IF
   END DO
   END DO
END DO
END IF
END DO

DO n=2,nNodes+1
   csp2EX(n) = csp2EX(n)+csp2EX(n-1)
END DO
nCsp1EX = csp2EX(nNodes+1)
i       = nCsp1E+nCsp1EX
IF (i > ncsp1) THEN
   WRITE(*,*)'*** ncsp1 exceeds dimensions in gridGeomFormCsp ***'
   STOP
END IF
ALLOCATE(csp1EX(nCsp1EX))

flag = 0
DO n=1,nNodes-nGnodes
   deg = csp2E(n+1)-csp2E(n)
IF (deg <= 1) THEN
DO mm=csp2E(n)+1,csp2E(n+1)
   j       = csp1E(mm)
   flag(j) = n
END DO
DO mm=csp2E(n)+1,csp2E(n+1)
   j       = csp1E(mm)
   DO k=1,2
      jj = face(k,j)
   DO m=csp2E(jj)+1,csp2E(jj+1)
      i  = csp1E(m)
      IF (flag(i) /= n) THEN
         flag(i)    = n
         ii         = csp2EX(n)+1
         csp2EX(n)  = ii
         csp1EX(ii) = i
      END IF
   END DO
   END DO
END DO
END IF
END DO
DO n=nNodes+1,2,-1
   csp2EX(n) = csp2EX(n-1)
END DO
csp2EX(1) = 0


! Combine csp and cspe arrays

k = 0
DO n=1,nNodes
   DO m=csp2E(n)+1,csp2E(n+1)
      k = k+1
      csp1(k) = csp1E(m)
   END DO
   DO m=csp2EX(n)+1,csp2EX(n+1)
      k = k+1
      csp1(k) = csp1EX(m)
   END DO
END DO
IF (k > nCsp1) THEN
   WRITE(*,*)'*** nCsp1 added incorrectly in gridGeomFormCsp ***'
   STOP
END IF

DO n=1,nNodes+1
   csp2(n) = csp2E(n)+csp2EX(n)
END DO

IF (csp2(nNodes+1) > ncsp1) THEN
   WRITE(*,*)'*** ncsp1 too small in formCsp ***'
   STOP
END IF

DEALLOCATE( &
     flag, &
     csp1E, &
     csp2E, &
     csp1EX, &
     csp2EX)

!DO n=1,nNodes
!   write(*,*)n,(csp1(m),m=csp2(n)+1,csp2(n+1))
!END DO

! make arrays 1-based
face = face-1
csp1 = csp1-1


END SUBROUTINE formcsp
