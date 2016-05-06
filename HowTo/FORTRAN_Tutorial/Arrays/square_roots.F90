PROGRAM square_roots
! initialize arrays using implied do loops

IMPLICIT NONE

INTEGER::i
REAL,DIMENSION(10)::value=(/(i,i=1,10)/)
REAL,DIMENSION(10)::square_root

!calculate the square roots of the numbers
DO i=1,10
    square_root(i)=sqrt(value(i))
END DO

!write out the values and each number
DO i=1,10
    WRITE(*,100) value(i), square_root(i)
    100 FORMAT (1x, 'Value= ',F5.1,' Square Root = ',F10.4)
END DO

END PROGRAM
