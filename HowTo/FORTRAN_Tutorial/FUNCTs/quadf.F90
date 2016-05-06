REAL FUNCTION quadf(x,a,b,c)
!purpose, evaluate the quadratic polynomial of the form 
!quadf = a*x^2 + b*x + c

IMPLICIT NONE
REAL,INTENT(IN)::x,a,b,c

!evaluation expression
quadf=a*x**2 + b*x + c
END FUNCTION quadf
