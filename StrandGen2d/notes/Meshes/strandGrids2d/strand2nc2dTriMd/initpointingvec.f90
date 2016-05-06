!> \brief
!! This subroutine initializes the strand pointing vectors.
!!
!! Comments:
!!
!! Versions:
!!   - 1.0 Katz 07/23/2010
!!
!! Additional notes:\par
!!   none
!!
!! Source code:
!!   \include initPointingVec.f90


SUBROUTINE initpointingvec(iLo,iHi,iLoN,iHiN,ig1,ig2,ig,ig1N,igN,nBNodes, &
                           ncsp1,pv,pv0,ang,csp1,csp2,xSurf,face,  &
                           bNode,bNorm)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: iLo,iHi,iLoN,iHiN,ig1,ig2,ig,ig1N,igN,nBNodes, &
                         ncsp1
INTEGER,INTENT(IN   ),DIMENSION(ncsp1    ) :: csp1
INTEGER,INTENT(IN   ),DIMENSION(  iLoN:iHiN+igN+1) :: csp2
REAL   ,INTENT(  OUT),DIMENSION(  iLoN:iHiN+igN  ) :: ang
REAL   ,INTENT(  OUT),DIMENSION(2,iLoN:iHiN+igN  ) :: pv
REAL   ,INTENT(  OUT),DIMENSION(2,iLoN:iHiN+igN  ) :: pv0
INTEGER,INTENT(IN   ),DIMENSION(2,iLo:iHi+ig) :: face
INTEGER,INTENT(IN   ),DIMENSION(  nBNodes) :: bNode
REAL,   INTENT(IN   ),DIMENSION(2,nBNodes) :: bNorm
REAL,   INTENT(IN   ),DIMENSION(2,iLoN:iHiN+igN  ) :: xSurf
INTEGER :: n1,n2,n,m,i,j
REAL :: pi,xa,ya,xb,yb,ax,ay,ds,a1,nx,ny,nx1,ny1
REAL,ALLOCATABLE,DIMENSION(:,:) :: pvc
REAL :: ang0 = 90.


pi   = 4.*ATAN(1.)
ang0 = ang0*pi/180.


! find surface cell normals

ALLOCATE(pvc(2,iLo:iHi+ig1+ig2))
DO n=iLo,iHi+ig1+ig2
   n1       = face(1,n)
   n2       = face(2,n)
   xa       = xSurf(1,n1)
   ya       = xSurf(2,n1)
   xb       = xSurf(1,n2)
   yb       = xSurf(2,n2)
   ax       = xb-xa
   ay       = yb-ya
   ds       = 1./SQRT(ax*ax+ay*ay)
   ax       = ax*ds
   ay       = ay*ds
   pvc(1,n) =-ay
   pvc(2,n) = ax
END DO


! compute node normals

DO n=iLoN,iHiN+ig1N
   ax = 0.
   ay = 0.
   DO m=csp2(n)+1,csp2(n+1)
      j  = csp1(m)
      ax = ax+pvc(1,j)
      ay = ay+pvc(2,j)
   END DO
   ds       = 1./SQRT(ax*ax+ay*ay)
   pv0(1,n) = ax*ds
   pv0(2,n) = ay*ds
END DO


! Constrain boundary vectors to lie in their planes

DO n=1,nBNodes
   i        = bNode(n)
   ax       = 0.
   ay       = 0.
   DO m=csp2(i)+1,csp2(i+1)
      j     = csp1(m)
      ax    = ax+pvc(1,j)
      ay    = ay+pvc(2,j)
   END DO
   nx1      = bNorm(1,n)
   ny1      = bNorm(2,n)
   a1       = ax*nx1+ay*ny1
   ax       = ax-a1*nx1
   ay       = ay-a1*ny1
   ds       = 1./SQRT(ax*ax+ay*ay)
   pv0(1,i) = ax*ds
   pv0(2,i) = ay*ds
END DO

pv = pv0


! find visibility cone half angle

DO n=iLoN,iHiN+ig1N
   nx = pv0(1,n)
   ny = pv0(2,n)
   a1 = ang0
   DO m=csp2(n)+1,csp2(n+1)
      j  = csp1(m)
      ax = pvc(1,j)
      ay = pvc(2,j)
      ds = nx*ax+ny*ay
      ds = .5*pi-ACOS(ds)
      a1 = MIN(a1,ds)
   END DO
   ang(n) = COS(a1)
   IF (a1 < 0.) THEN
      WRITE(*,*)
      WRITE(*,*)'*** visibility condition not met at node: ***'
      WRITE(*,'(I10,3ES14.4)')n,xSurf(:,n)
      WRITE(*,*)
   END IF
END DO

DEALLOCATE(pvc)


END SUBROUTINE initpointingvec
