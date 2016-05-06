PROGRAM test_hypotenuse
! PURPOSE: TEST SUBROUTINE CALC_HYPOTENUSE

IMPLICIT NONE

!data dictionary: variables
REAL::s1,s2,hypot

!get lengths of the two sides
WRITE(*,*)'Enter length of side 1: '
READ(*,*)s1
WRITE(*,*)'Enter length of side 2: '
READ(*,*)s2

!call calc_hypotenuse
CALL calc_hypotenuse(s1,s2,hypot)

!output
WRITE(*,1000) hypot
1000 FORMAT (1x,'The length of the hypotenuse is: ',F10.4)


END PROGRAM test_hypotenuse
