SUBROUTINE running_average (x,ave,std_dev,nvals,reset)
    !purpose: calculate a running average, std_dev, and number of data points as data values are received

    IMPLICIT NONE
    !data dictionary: parameters
    REAL,INTENT(IN)::x
    REAL,INTENT(OUT)::ave,std_dev
    INTEGER,INTENT(OUT)::nvals
    LOGICAL,INTENT(IN)::reset
    !data dictionary: variables
    INTEGER,SAVE::n
    REAL,SAVE::sum_x
    REAL,SAVE::sum_x2
    !if reset flag, then reset to zero
    calc_sums: IF (reset) THEN
        n       = 0
        sum_x   = 0.
        sum_x2  = 0.
        ave     = 0.
        std_dev = 0.
        nvals   = 0
    ELSE
        !accumulate sums if not
        n=n+1
        sum_x=sum_x+x
        sum_x2=sum_x2+x**2
        !calc average
        ave=sum_x/REAL(n)
        !std_dev
        IF (n>=2) THEN
            std_dev=SQRT((REAL(n) * sum_x2 - sum_x**2) &
                / (REAL(n) * REAL(n-1)) )
        ELSE 
            std_dev = 0.
        END IF 
        ! number of data points
        nvals=n
    END IF calc_sums
END SUBROUTINE running_average
