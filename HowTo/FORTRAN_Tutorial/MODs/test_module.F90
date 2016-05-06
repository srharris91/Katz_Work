PROGRAM test_module
!purpose: illustrate sharing data via a module
USE shared_data
IMPLICIT NONE
REAL,PARAMETER::PI=3.141592
values=PI*(/1.,2.,3.,4.,5./)
CALL sub1

END PROGRAM test_module

!*********************************************************************
!*********************************************************************
SUBROUTINE sub1
!purpose: illustrate sharing data via a module
USE shared_data
IMPLICIT NONE

WRITE(*,*) values

END SUBROUTINE sub1
