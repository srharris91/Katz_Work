PROGRAM test_sinc
!purpose: test sinc function
IMPLICIT NONE
!data dictionary: function types
REAL::sinc !function declaration
!data dictionary: variables
REAL::x

!get value to evaluate
WRITE(*,*)'Enter x: '
READ(*,*)x
!output answer
WRITE(*,'(1x,A,F8.5)') 'sinc(x) = ',sinc(x)
END PROGRAM test_sinc


