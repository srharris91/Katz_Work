FUNCTION sinc(x)
!purpose: calc sinc=sin/x

IMPLICIT NONE
REAL,INTENT(IN)::x
REAL::sinc

REAL,PARAMETER::EPSILON=1.0E-30

!Check to see abs(x)>epsilon
IF (ABS(x)>EPSILON) THEN
    sinc=SIN(x)/x
ELSE
    sinc=1.
END IF

END FUNCTION sinc
