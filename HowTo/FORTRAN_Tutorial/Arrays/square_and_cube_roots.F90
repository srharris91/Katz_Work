PROGRAM square_and_cube_roots
!purpose to calc the table of numbers, square roots, and cube roots using implied do loops

IMPLICIT NONE
!data dictionary constants
INTEGER, PARAMETER::MAX_SIZE=10

!data dictionary variables
INTEGER::j,i
REAL,DIMENSION(MAX_SIZE)::value
REAL,DIMENSION(MAX_SIZE)::square_root
REAL,DIMENSION(MAX_SIZE)::cube_root

!calc square roots and cube roots
DO j=1,MAX_SIZE
    value(j)=real(j)
    square_root(j)=sqrt(value(j))
    cube_root(j)=value(j)**(1.0/3.0)
END DO

!write out each number, square root, and cube root
WRITE(*,100)
100 FORMAT ('0',20x,'Table of Square and Cube Roots',/&
                4x, '  Number      Square Root   Cube Root',&
                3x, '  Number      Square Root   Cube Root',/,&
                4x, '=========    ============= ===========', &
                3x, '=========    ============= ===========')
WRITE(*,110) (value(j),square_root(j),cube_root(j),j=1,MAX_SIZE)
110 FORMAT (2(4x,F6.0,9x,F6.4,6x,F6.4))


!output nested implied do loop
WRITE(*,120) ((i,j,j=1,3),i=1,2)
120 FORMAT (1x,I5,1x,I5)


END PROGRAM square_and_cube_roots
