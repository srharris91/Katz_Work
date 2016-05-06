PROGRAM bounds
! check if compiler catched out of bounds stuff

IMPLICIT NONE

INTEGER::i
REAL,DIMENSION(5)::a=(/1.,2.,3.,4.,5./)
REAL,DIMENSION(5)::b=(/10.,20.,30.,40.,50./)

DO i=1,6
    WRITE(*,100)i,a(i)
    100 FORMAT (1x,'a(',I1,') = ',F6.2)
END DO


END PROGRAM bounds
