PROGRAM test_subs_as_arguments
!purpose: test passing subroutine names as arguments

IMPLICIT NONE
!Data dictionary: parameters
EXTERNAL::sum,prod
REAL::x,y,result

!get x,y
WRITE(*,*)' ENTER x and y: '
READ(*,*) x,y

!calculate the product
CALL subs_as_arguments(x,y,prod,result)
WRITE(*,*) 'The product is ',result
!calculate the sum
CALL subs_as_arguments(x,y,sum,result)
WRITE(*,*) 'The sum is ',result


END PROGRAM test_subs_as_arguments

!************************************************************
!************************************************************

SUBROUTINE prod(x,y,result)
! calc product of two real numbers
IMPLICIT NONE
REAL,INTENT(IN)::x,y
REAL,INTENT(OUT)::result

!calc value
result = x*y

END SUBROUTINE prod

!************************************************************
!************************************************************

SUBROUTINE sum(x,y,result)
! calc sum of two real numbers
IMPLICIT NONE
REAL,INTENT(IN)::x,y
REAL,INTENT(OUT)::result

!calc value
result = x+y

END SUBROUTINE sum
