SUBROUTINE subs_as_arguments(x,y,sub,result)
!purpose: test passing names as arguments
IMPLICIT NONE

!data dictionary: parameters
EXTERNAL::sub
REAL,INTENT(IN)::x,y
REAL,INTENT(OUT)::result

CALL sub(x,y,result)

END SUBROUTINE subs_as_arguments
