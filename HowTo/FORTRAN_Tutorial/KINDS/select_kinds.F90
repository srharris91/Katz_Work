PROGRAM select_kinds
    !purpose: illustrate the used of SELECTED_REAL_KIND, KIND(), PRECISION, and RANGE()
    IMPLICIT NONE
    !declare parameters
    INTEGER,PARAMETER::SGL=SELECTED_REAL_KIND(p=6,r=37)!(p is # of decimal digits of precision, r is range of the exponent required in powers of 10
    INTEGER,PARAMETER::DBL=SELECTED_REAL_KIND(p=13,r=200)


    !Declare variables for each type:
    REAL(kind=SGL):: var1=0.
    REAL(kind=DBL)::var2=0._DBL

    !WRITE characteristics of selected variables
    WRITE(*,100) 'var1',KIND(var1), PRECISION(var1), RANGE(var1)
    WRITE(*,100) 'var2',KIND(var2), PRECISION(var2), RANGE(var2)
    100 FORMAT (1x,A,': kind = ',I2,', Precision = ',I2,', Range = ',I3)
END PROGRAM select_kinds
