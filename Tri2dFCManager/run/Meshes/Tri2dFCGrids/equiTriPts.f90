SUBROUTINE equiTriPts(order,nsp,rs)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: order,nsp
REAL   ,INTENT(  OUT),DIMENSION(2,nsp) :: rs
INTEGER :: i,j,m
REAL :: L0,L1,dL
REAL,DIMENSION(2,nsp) :: rsT
INTEGER,DIMENSION(nsp) :: key


! equally spaced points in the unit equilateral triangle
! with numbering consistent with gmesh format

IF (order == 0) THEN
   rs(1,1) = 0.
   rs(2,1) = 0.
ELSE IF (order > 5) THEN
   WRITE(*,*)'*** Please choose order < 6 in equiTriPts.f90 ***'
   STOP
ELSE
   dL = 1./REAL(order)
   m  = 0
   DO i=0,order
      L0 = REAL(i)*dL
      DO j=0,order-i
         L1       = REAL(j)*dL
         m        = m+1
         rsT(1,m) = L1-L0
         rsT(2,m) =(2.-3.*(L0+L1))/SQRT(3.)
      END DO
   END DO

   ! order consistent with gmesh format
   IF      (order == 1) THEN
      key(1)  = 3
      key(2)  = 2
      key(3)  = 1
   ELSE IF (order  == 2) THEN
      key(1)  = 6
      key(2)  = 3
      key(3)  = 1
      key(4)  = 5
      key(5)  = 2
      key(6)  = 4
   ELSE IF (order == 3) THEN
      key(1)  = 10
      key(2)  = 4
      key(3)  = 1
      key(4)  = 9
      key(5)  = 7
      key(6)  = 3
      key(7)  = 2
      key(8)  = 5
      key(9)  = 8
      key(10) = 6
   ELSE IF (order == 4) THEN
      key(1)  = 15
      key(2)  = 5
      key(3)  = 1
      key(4)  = 14
      key(5)  = 12
      key(6)  = 9
      key(7)  = 4
      key(8)  = 3
      key(9)  = 2
      key(10) = 6
      key(11) = 10
      key(12) = 13
      key(13) = 11
      key(14) = 8
      key(15) = 7
   ELSE IF (order == 5) THEN
      key(1)  = 21
      key(2)  = 6
      key(3)  = 1
      key(4)  = 20
      key(5)  = 18
      key(6)  = 15
      key(7)  = 11
      key(8)  = 5
      key(9)  = 4
      key(10) = 3
      key(11) = 2
      key(12) = 7
      key(13) = 12
      key(14) = 16
      key(15) = 19
      key(16) = 17
      key(17) = 10
      key(18) = 8
      key(19) = 14
      key(20) = 9
      key(21) = 13
   END IF

   ! change numbering
   DO i=1,nsp
      rs(1,i) = rsT(1,key(i))
      rs(2,i) = rsT(2,key(i))
   END DO
END IF


END SUBROUTINE equiTriPts
