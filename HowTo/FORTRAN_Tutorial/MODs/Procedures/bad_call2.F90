PROGRAM bad_call2
USE my_subs
IMPLICIT NONE
REAL::x=1.
CALL bad_argument(x)
END PROGRAM bad_call2
