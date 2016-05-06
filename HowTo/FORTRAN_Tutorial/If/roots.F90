PROGRAM roots
!purpose: find the roots of a quadratic equation shown in this form
! a*x^2 + b*x + c = 0

IMPLICIT NONE

REAL::a,b,c,discriminant,imag_part,real_part,x1,x2

!Prompt user for coefficients a,b,c
WRITE(*,*) "Please input the coefficients of the quadratic equation a*x^2 + b*x + c = 0 in order of a,b, and c"
READ(*,*) a,b,c
WRITE(*,*)"The coefficients a,b and c are = ",a,b,c

!calculate the discriminant
discriminant = b**2 -4*a*c

!solve for the roots
IF(discriminant > 0.) THEN !there are two real roots
    x1 = (-b+sqrt(discriminant))/(2.*a)
    x2 = (-b-sqrt(discriminant))/(2.*a)
    WRITE(*,*) "x1 = ",x1
    WRITE(*,*) "x2 = ",x2
ELSE IF (discriminant < 0.) THEN !there are complex roots
    real_part = (-b)/(2.*a)
    imag_part = sqrt(abs(discriminant))/(2.*a)
    WRITE(*,*) "This equation has complex roots:"
    WRITE(*,*) "x1 = ",real_part,"+i",imag_part
    WRITE(*,*) "x2 = ",real_part,"-i",imag_part
ELSE IF (discriminant == 0.) THEN !there is one repeated root
    x1 = (-b)/(2.*a)
    WRITE(*,*) "This equation has two identical real roots:"
    WRITE(*,*) "x1 = x2 = ",x1
END IF


END PROGRAM roots
