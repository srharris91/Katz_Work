PROGRAM test_ave_value
!purpose: test function ave_value by calling it with a user-defined function my_func

IMPLICIT NONE
!data dictionary: function types
REAL::ave_value
REAL,External::my_function

!data dictionary: variables
REAL::ave

!call function
ave=ave_value(my_function,0.,1.,101)
WRITE(*,1000)'my_function',ave
1000 FORMAT (1x,'The average value of ',A,' between 0. and 1. is ', F16.6,' . ')

END PROGRAM test_ave_value

REAL FUNCTION my_function(x)
IMPLICIT NONE
REAL,INTENT(IN)::x
my_function = 3.*x
END FUNCTION my_function

