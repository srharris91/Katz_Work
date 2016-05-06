!> \brief
!! These subroutines perform LU factorization and solution of a
!! block tridiagonal system of equations.
!!
!! Comments:
!!
!! Versions:
!!   - 1.0 Katz 04/2/2012
!!
!! Additional notes:\par
!!  none
!!
!! Source code:
!!   \include nbtrluA.F90

SUBROUTINE nbtrluA(il,ih,i,nq,a,b,c,d,e)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: il,ih,i,nq
REAL   ,INTENT(IN   ),DIMENSION(nq,nq) :: a,b,c
REAL   ,INTENT(INOUT),DIMENSION(nq,il:ih) :: d
REAL   ,INTENT(INOUT),DIMENSION(nq,nq,il:ih) :: e
REAL,DIMENSION(nq) :: dd
REAL,DIMENSION(nq,nq) :: aa,bb


IF (i == il) THEN
   CALL matinv(nq,b,aa)
   d(:,i)   = MATMUL(aa,d(:,i))
   e(:,:,i) = MATMUL(aa,c)
ELSE
   aa       = b-MATMUL(a,e(:,:,i-1))
   CALL matinv(nq,aa,bb)
   dd       = d(:,i)-MATMUL(a,d(:,i-1))
   d(:,i)   = MATMUL(bb,dd)
   e(:,:,i) = MATMUL(bb,c)
END IF


END SUBROUTINE nbtrluA




!==========================================================================
SUBROUTINE nbtrbkA(il,ih,nq,d,e)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: il,ih,nq
REAL   ,INTENT(INOUT),DIMENSION(nq,il:ih) :: d
REAL   ,INTENT(IN   ),DIMENSION(nq,nq,il:ih) :: e
INTEGER :: i


DO i=ih-1,il,-1
   d(:,i) = d(:,i)-MATMUL(e(:,:,i),d(:,i+1))
END DO


END SUBROUTINE nbtrbkA



!==========================================================================
SUBROUTINE matinv(n,a,b)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: n
REAL   ,INTENT(IN   ),DIMENSION(n,n) :: a
REAL   ,INTENT(  OUT),DIMENSION(n,n) :: b
INTEGER :: info
INTEGER,DIMENSION(n) :: ipiv
REAL,DIMENSION(n) :: work


b = a
CALL sgetrf(n,n,b,n,ipiv,info)
CALL sgetri(n,b,n,ipiv,work,n,info)


END SUBROUTINE matinv
