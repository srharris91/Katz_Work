Module ran001
!purpose: declare data shared between subs, random0 and seed
IMPLICIT NONE
SAVE
INTEGER::n=9876
END MODULE ran001

!*******************************************************************
!*******************************************************************

SUBROUTINE random0 (ran)
!purpose: generate a pseudorandom number between 0 and 1.0

USE ran001  ! shared seed
IMPLICIT NONE
!data dictionary: parameters
REAL,INTENT(OUT)::ran

!Calc next number
n=MOD(8121*n + 28411, 134456)

!generate random value from this number
ran = REAL(n) / 134456

END SUBROUTINE random0

!*******************************************************************
!*******************************************************************

SUBROUTINE seed (iseed)
!purpose: set the seed for random number generator random0

USE ran001

IMPLICIT NONE
INTEGER,INTENT(IN)::iseed

!set seed
n=ABS(iseed)

END SUBROUTINE seed

!*******************************************************************
!*******************************************************************

