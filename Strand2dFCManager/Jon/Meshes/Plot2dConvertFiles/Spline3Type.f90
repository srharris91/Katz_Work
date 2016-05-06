MODULE Spline3Type
USE Pre
IMPLICIT NONE

REAL(wp),DIMENSION(2)::S,Y,Ypp, X, Xpp, YT, YppT
REAL(wp)::h

PRIVATE::h, S, Y, Ypp, X, Xpp, YT, YppT

CONTAINS

SUBROUTINE SetLim(So,Yo,Yppo,Xo,Xppo)
REAL(wp),DIMENSION(2),INTENT(IN)::So,Yo,Yppo
REAL(wp),DIMENSION(2),INTENT(IN),OPTIONAL::Xo,Xppo

IF (PRESENT(Xo) .and. PRESENT(Xppo)) THEN
	X   = Xo
	Xpp = Xppo
END IF

	S   = So
	Y   = Yo
	Ypp = Yppo
	h   = So(2)-So(1)
	
END SUBROUTINE SetLim

FUNCTION SplinePts( si )
! Finds location of points along the splines
! s   - peramiterized value
! Active%pts - Active Point Set (Must be filled)!
REAL(wp):: SplinePts
REAL(wp), INTENT(IN)::si
REAL(wp):: hs6 
REAL(wp):: Sp 
REAL(wp):: Sm 
hs6 = h**2/6._wp
Sp  = (S(2)-si)/h
Sm  = (si-S(1))/h
!Write(*,*)"P ",Y(1),Y(2),"Pr ",Ypp(1),Ypp(2),"S ",S(1),S(2)
SplinePts = Y(1)*Sp + Y(2)*Sm - hs6*Ypp(1)*(Sp-Sp**3) - hs6*Ypp(2)*(Sm-Sm**3)

END FUNCTION SplinePts

FUNCTION SplinePass( si )
REAL(wp)::SplinePass
REAL(wp),INTENT(IN)::si

REAL(wp)::dsx,dsy
! Doesn't Matter Which one is X and Y, as long as both of them are there !
YT=Y
YppT=Ypp
dsy = SPlineDrv( si )
YT=X
YppT=Xpp
dsx = SplineDrv( si )

SplinePass = sqrt(dsx**2+dsy**2)
!WRITE(*,*)"SplinePass ",SplinePass
END FUNCTION

FUNCTION SplineDrv( si )
REAL(wp)::SplineDrv
REAL(wp), INTENT(IN)::si
REAL(wp)::a,b,c,d,e,f,g,x,si2,sip
REAL(wp)::Yppi,Yppip,Yi,Yip
Yppi=YppT(1)
Yppip=YppT(2)
Yi = YT(1)
Yip= YT(2)
x = si
si2=S(1)
sip=S(2)
a=3._wp*(Yppi-Yppip)*x**2
b=a+6._wp*(si2*Yppip-sip*Yppi)*x
c=b-h**2*(Yppi-Yppip)
d=c-3._wp*(si2**2*Yppip-sip**2*yppi-2._wp*(Yi-Yip))

SplineDrv = -d/(6._wp*h) 

END FUNCTION

END MODULE Spline3Type
