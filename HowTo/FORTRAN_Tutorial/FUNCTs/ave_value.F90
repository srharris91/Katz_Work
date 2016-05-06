REAL FUNCTION ave_value (func, first_value, last_value, n)
!purpose; evaluate func of function over the range by taking n evenly spaced samples
IMPLICIT NONE
!data dictionary: parameters
REAL,EXTERNAL::func     ! external function (input parameter value)
REAL,INTENT(IN)::first_value, last_value
INTEGER,INTENT(IN)::n

!data dictionary: variables
REAL::delta
INTEGER::i
REAL::sum

!get step size
delta = (last_value - first_value) / (REAL(n-1))

!accumulate sum
sum=0.
DO i=1,n
    sum=sum+func(REAL(i-1) * delta)
END DO

! get average
ave_value = sum/REAL(n)

END FUNCTION ave_value
