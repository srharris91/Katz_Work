PROGRAM add_arrays

!illustrate element-by-element addition and whole array addition

IMPLICIT NONE

INTEGER::i
REAL,DIMENSION(4)::a=(/1.,2.,3.,4./)
REAL,DIMENSION(4)::b=(/5.,6.,7.,8./)
REAL,DIMENSION(4)::c,d

!element-by-element addition
DO i=1,4
    c(i)=a(i)+b(i)
END DO

! whole array addition
d=a+b

!write out each of them
WRITE(*,100)'c',c
WRITE(*,100)'d',d
100 FORMAT (' ',A,' = ',5(F6.1,1x))



END PROGRAM add_arrays
