!> \brief
!! This subroutine fills the entries of the csp1 array.
!!
!! Comments:
!!
!! Versions:
!!   - 1.0 Katz 08/16/2010
!!
!! Additional notes:\par
!!  none
!!
!! Source code:
!!   \include fillCsp.f90


SUBROUTINE fillcsp(iHi,iLo,iLoN,iHiN,ig1,ig2,ig,igN,ncsp1,face,csp1,csp2)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: iHi,iLo,iLoN,iHiN,ig1,ig2,ig,igN,ncsp1
INTEGER,INTENT(IN   ),DIMENSION(2,iLo:iHi+ig) :: face
INTEGER,INTENT(INOUT),DIMENSION(ncsp1) :: csp1
INTEGER,INTENT(INOUT),DIMENSION(iLoN:iHiN+igN+1) :: csp2
INTEGER :: n,i,j,k


DO n=iLo,iHi+ig1+ig2
DO k=1,2
   j       = face(k,n)
   i       = csp2(j)+1
   csp2(j) = i
   csp1(i) = n
END DO
END DO
DO n=iHiN+igN+1,iLoN+1,-1
   csp2(n) = csp2(n-1)
END DO
csp2(iLoN) = 0

!do n=iLoN,iHiN
!   write(*,'(i12,10i8)')n,(csp1(j),j=csp2(n)+1,csp2(n+1))
!end do


END SUBROUTINE fillcsp
