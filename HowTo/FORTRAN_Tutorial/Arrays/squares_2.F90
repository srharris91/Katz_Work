PROGRAM squares_2

IMPLICIT NONE
! use subscript alterations

INTEGER::i
INTEGER,DIMENSION(-5:5)::number,square

!initial9ize number and calculate square
DO i=-5,5
    number(i)=i
    square(i)=number(i)**2
END DO

DO i=-5,5
    WRITE(*,100) number(i),square(i)
    100 FORMAT (1x,'Number= ',I6,' Square = ',I6)
END DO




END PROGRAM squares_2

