PROGRAM test_random0
!purpose: subroutine to test the random number generator random0
IMPLICIT NONE

!data dictionary: variables
REAL::ave
INTEGER::i,iseed,iseq
REAL::ran,sum

!get seed
WRITE(*,*)'Enter seed: '
READ(*,*) iseed

!set seed
CALL SEED (iseed)

!print out 10 random #'s
WRITE(*,*)'10 Random numbers: '
DO i=1,10
    CALL random0(ran)
    WRITE(*,'(3x,F16.6)') ran
END DO

! average 5 consectuve 1000-value sequences
WRITE(*,*) 'Averages of 5 consecutive 1000-sample sequences:'
DO iseq=1,5
    sum=0.
    DO i=1,1000
        CALL random0(ran)
        sum=sum+ran
    END DO
    ave=sum/1000.
    WRITE(*,'(3x,F16.6)')ave
END DO

END PROGRAM test_random0
