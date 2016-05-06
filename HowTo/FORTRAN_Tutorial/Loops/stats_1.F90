PROGRAM stats_1

! purpose to calculate the mean and standard deviation of an input data set containing an arbitrary number of input values

IMPLICIT NONE

INTEGER::n=0    ! number of input samples
REAL::std_dev=0.!standard deviation of input samples
REAL::sum_x=0.  !sum of input samples
REAL::sum_x2=0. !sum of squares
REAL::x=0.      !input data value
REAL::x_bar     !the average of the input samples

!while look to read input values
DO
    !read next value
    WRITE(*,*)'Enter number:'
    READ(*,*)x
    WRITE(*,*)'The number is = ',x

    !test for loop exit
    IF (x<0.) EXIT

    !otherwise accumulate sums
    n = n+1
    sum_x = sum_x+x
    sum_x2= sum_x2 + x**2
END DO

! Calculate the mean and standard deviation
x_bar = sum_x / real(n)
std_dev = sqrt((real(n) * sum_x2 - sum_x**2) / (real(n) * real(n-1)))

!output to user
WRITE(*,*)'The mean is = ',x_bar
WRITE(*,*)'Standard Dev.=',std_dev
WRITE(*,*)'Number of pts=',n

END PROGRAM stats_1

