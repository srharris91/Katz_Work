SUBROUTINE calc_hypotenuse (side_1, side_2, hypotenuse)
! purpose: to calc the hypotenuse of a right triangle from the two sides

IMPLICIT NONE

!Data dictionary parameters
REAL, INTENT(IN):: side_1
REAL, INTENT(IN):: side_2
REAL, INTENT(OUT):: hypotenuse

!data dictionary: variables
REAL::temp

!calc hypotenuse
temp=side_1**2 + side_2**2
hypotenuse=sqrt(temp)

END SUBROUTINE calc_hypotenuse
