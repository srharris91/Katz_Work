!> \brief
!! This subroutine computes the projected face area of a triangle.
!!
!! Comments:
!!
!! Versions:
!!   - 1.0 Katz 04/20/2010
!!
!! Additional notes:\par
!!  none
!!
!! Source code:
!!   \include triProjArea.F90


SUBROUTINE triProjArea(xa,ya,za,xb,yb,zb,xc,yc,zc,ax,ay,az)

IMPLICIT NONE

REAL,INTENT(IN   ) :: xa,ya,za,xb,yb,zb,xc,yc,zc
REAL,INTENT(  OUT) :: ax,ay,az
REAL :: x1,y1,z1,x2,y2,z2


x1 = xa-xb
y1 = ya-yb
z1 = za-zb
x2 = xc-xb
y2 = yc-yb
z2 = zc-zb

ax = .5*(y1*z2-y2*z1)
ay = .5*(x2*z1-x1*z2)
az = .5*(x1*y2-x2*y1)


END SUBROUTINE triProjArea
