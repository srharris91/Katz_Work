PROGRAM test_quadf
!purpose: test quadf function
IMPLICIT NONE
REAL::quadf,a,b,c,x

!get input data
WRITE(*,*)'Enter quadratic coefficients a,b,c: '
READ(*,*)a,b,c
WRITE(*,*)'Enter x location: '
READ(*,*) x

!write result
WRITE(*,100) 'quadf(',x,') = ',quadf(x,a,b,c)
100 FORMAT (A,F10.4,A,F12.4)
END PROGRAM test_quadf
