MODULE my_subs
CONTAINS
    SUBROUTINE bad_argument (i)
    IMPLICIT NONE
    INTEGER, INTENT(IN)::i
    WRITE(*,*)' i= ',i
    END SUBROUTINE bad_argument
END MODULE my_subs
