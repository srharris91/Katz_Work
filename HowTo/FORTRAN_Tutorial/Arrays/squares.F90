PROGRAM squares

!purpose arrays introduction

IMPLICIT NONE

INTEGER::i
INTEGER, DIMENSION(10):: number, square

! initialize number and calculate square
DO i=1,10
    number(i)=i
    square(i)=number(i)**2
END DO

! write out stuff
DO i=1,10
    WRITE(*,100) number(i), square(i)
       100 FORMAT (1X,'Number = ',I6, ' Square = ',I6)
END DO


END PROGRAM squares
