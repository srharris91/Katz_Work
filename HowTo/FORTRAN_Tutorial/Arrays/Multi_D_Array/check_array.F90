PROGRAM check_array
!purpose: illustrate array inquiry functions
IMPLICIT NONE
!list of variables:
REAL,DIMENSION(-5:5,0:3) :: a=0

!get shape and size of bounds
WRITE(*,100) SHAPE(a)
100 FORMAT (1x,'The shape of the array is           :',7I6)
WRITE(*,110) SIZE(a)
110 FORMAT (1x,'The size of the array is            :',I6)
WRITE(*,120) LBOUND(a)
120 FORMAT (1x,'The Lower bounds of the array is    :',7I6)
WRITE(*,130) UBOUND(a)
130 FORMAT (1x,'The Upper bounds of the array is    :',7I6)
END PROGRAM check_array
