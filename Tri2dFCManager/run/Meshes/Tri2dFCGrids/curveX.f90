SUBROUTINE curveX( &
     type, &
     nInfo, &
     info, &
     z, &
     xInt)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: type,nInfo
REAL   ,INTENT(IN   ),DIMENSION(nInfo) :: info
REAL   ,INTENT(IN   ) :: z
REAL   ,INTENT(  OUT),DIMENSION(2) :: xInt
REAL   :: zm,zp,xm,ym,xp,yp,a,R,tm,tp
REAL,PARAMETER :: eps=1.D-14
REAL,DIMENSION(2) :: xe0


zm =.5*(1.-z)
zp =.5*(1.+z)
a  = 1./max(zm*zp,eps)

IF      (type == 0) THEN !linear
   xe0(1) = 0.
   xe0(2) = 0.
   xm     = 0.
   ym     = 0.
   xp     = 0.
   yp     = 0.
ELSE IF (type == 1) THEN !circle
   R      = info(1)
   tm     = info(2)
   tp     = info(3)
   xe0(1) = R*COS(zm*tm+zp*tp)
   xe0(2) = R*SIN(zm*tm+zp*tp)
   xm     = R*COS(tm)
   ym     = R*SIN(tm)
   xp     = R*COS(tp)
   yp     = R*SIN(tp)
END IF

xe0(1)  = xe0(1)-xm*zm-xp*zp
xe0(2)  = xe0(2)-ym*zm-yp*zp
xInt(1) = a*xe0(1)
xInt(2) = a*xe0(2)


END SUBROUTINE curveX
