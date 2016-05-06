PROGRAM my_first_program

IMPLICIT NONE

!purpose: to illustrate some of the basic features of a fortran program

!Declare Variables
INTEGER :: i,j,k        ! all variables are integers

!Get Input from user
WRITE(*,*) 'Enter the numbers to multiply: '
READ(*,*)   i, j

! Multiply numbers together
k=i*j

!Write out the results
WRITE(*,*)  'Result = ',k

!Finish up
STOP
END PROGRAM
