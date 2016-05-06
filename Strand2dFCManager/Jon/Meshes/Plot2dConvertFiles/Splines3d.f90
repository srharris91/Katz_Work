MODULE Splines_3
USE Pre
USE Spline3Type
IMPLICIT NONE

contains

SUBROUTINE NewPts(n,ptsx,ptsy, nnpts, ptsxN, ptsyN)
! Brings in points that the splines will be used to connect.
! pts_ - Input points
! nnpts- number of new points between known points
! pts_N- Output points
! n    - Number of points
INTEGER,INTENT(IN)::n,nnpts
REAL(wp),INTENT(IN),DIMENSION(n)::ptsx,ptsy
REAL(wp),INTENT(OUT),DIMENSION(nnpts*(n-1))::ptsxN,ptsyN
REAL(wp),DIMENSION(n)::S, ppx, ppy, So
REAL(wp),DIMENSION(n-1)::ds
REAL(wp)::Er, Ert, nr,si
INTEGER::i,ii,d
! Initialize the S array
CALL CreateS(n, ptsx, ptsy, ppx, ppy, ds, S, 1)
WRITE(*,*)'Building New Points'
Er = 10._wp
Ert= 0._wp
nr = real(n,wp)
ppx= 0._wp
ppy= 0._wp

! Iterate to find Correct S values
DO WHILE(Er>.00000001_wp)
	! Set So as old S value
	So = S
	! Build Matrix For Solving and Solve
	CALL CalcRHS(n,ptsx,ds,ppx)
	CALL CalcRHS(n,ptsy,ds,ppy)
	CALL BuildMat(n, S, ds, ppx)
	CALL BuildMat(n, S, ds, ppy)
	! Remake S Array
	CALL CreateS( n, ptsx, ptsy, ppx, ppy, ds, S )
	
	Ert= 0._wp
	DO i=1,n
		Ert = Ert + (So(i)-S(i))**2/nr
	END DO
	Er = sqrt(Ert)
	WRITE(*,*) "Error : ", Er
END DO

! Solve for new pts
ii = 1
DO i=1,n-1
	d=1
	CALL SetLim(S(i:i+1),ptsx(i:i+1),ppx(i:i+1))
	DO ii = (i-1)*nnpts+1, nnpts*i
		si = (S(i+1)-S(i))/(nnpts+1)
		ptsxN(ii) = SplinePts(real(d,wp)*si+S(i))
		d=d+1
		WRITE(*,*)ptsxN(ii), ptsx(i:i+1)
	END DO
	d=1
	CALL SetLim(S(i:i+1),ptsy(i:i+1),ppy(i:i+1))
	DO ii = (i-1)*nnpts+1, nnpts*i
		si = (S(i+1)-S(i))/(nnpts+1)
		ptsyN(ii) = SplinePts(real(d,wp)*si+S(i))
		d=d+1
	END DO
END DO
END SUBROUTINE NewPts

SUBROUTINE BuildMat(n, S, ds, pp)
! Builds the needed matrix for solving for the curvature of the splines
!	Then solves them and returns curvature.
! n  - Number of points
! ds - ds values
! pp - array of curvatures

INTEGER, INTENT(IN):: n
REAL(wp), DIMENSION( n ), INTENT(INOUT):: pp
REAL(wp), DIMENSION(n-1), INTENT(IN):: ds
REAL(wp), DIMENSION( n ), INTENT(IN):: S

REAL(wp), DIMENSION(n-2,n-2):: A
INTEGER::i
A = 0._wp
A(1,1) = 2._wp*(ds(1)+ds(2))
A(1,2) = ds(2)

Do i=2,n-3
	A(i,i-1)= ds(i)
	A( i,i )= 2._wp*(ds(i)+ds(i+1))
	A(i,i+1)= ds(i+1)
	!pp(i)   = RHSCalc( S(i+1), S(i), S(i-1), ds(i+1), ds(i) )
END Do

A(n-2,n-3) =ds(n-2)
A(n-2,n-2) =2._wp*(ds(n-2)+ds(n-1))

CALL linear_system(n-2,A,pp(2:n-1))
pp(1) = 0._wp
pp(n) = 0._wp
END SUBROUTINE BuildMat

SUBROUTINE CreateS(n, pts1, pts2, pppts1, pppts2, ds, S, ptsy)
! Creates the ds and S that will be used to test for convergence
! n   - the number of points
! pts - the points
! pppts-double derivitive of points (empty array for first iteration)
! ds  - arc length between each point 
! S   - Total Arc Length
! Flag- optional flag for first iteration initialization

INTEGER, INTENT(IN):: n
REAL(wp),DIMENSION(n),INTENT(IN)::pts1, pppts1, pts2, pppts2
REAL(wp),DIMENSION(n-1),INTENT(INOUT)::ds
REAL(wp),INTENT(INOUT),DIMENSION(n)::S
INTEGER,OPTIONAL::ptsy

INTEGER::i
REAL(wp),DIMENSION(n)::So
So = S
IF (present(ptsy)) THEN
	S(1) = 0._wp
	DO i = 1, n-1
		ds(i)  = sqrt((pts1(i+1)-pts1(i))**2+(pts2(i+1)-pts2(i))**2)
		S(i+1) = ds(i) + S(i)
	END DO

ELSE
	
	DO i = 1, n-1
		CALL SetLim(So(i:i+1),pts1(i:i+1),pppts1(i:i+1),pts2(i:i+1), pppts2(i:i+1))
		ds(i) = abs(gauss1d(SplinePass,So(i+1),So(i)))
		S(i+1)= ds(i) + S(i)
		Call Errors( S(i+1), 132 )
	END DO

END IF

END SUBROUTINE CreateS

SUBROUTINE Errors( Num, i )
! This will watch for NaN's and Inf
! Num - number to be tested
! i   - number for error flag
REAL(wp),INTENT(IN)::Num
INTEGER,INTENT(IN) ::i

IF (Num /= Num) THEN
	WRITE(*,*) "NaN Faul Error Line ", i
	STOP
END IF

IF (Num > 9999999999999999999999999999999999999999999999999999._wp) THEN
	WRITE(*,*) "Inf Faul Error Line ", i
	STOP
END IF
END SUBROUTINE Errors

SUBROUTINE CalcRHS( n, pts, ds, pppts )
INTEGER,INTENT(IN)::n
REAL(wp),DIMENSION(n),INTENT(in)::pts
REAL(wp),DIMENSION(n-1),INTENT(in)::ds
REAL(wp),DIMENSION(n),INTENT(out)::pppts

INTEGER::i

DO i = 2,n-1
	pppts(i) = RHSCalc( pts(i+1), pts(i), pts(i-1), ds(i), ds(i-1) )
END DO

END SUBROUTINE

FUNCTION RHSCalc( sp1, s, sm1, h, hm1 )
! Calculates the RHS of the RHS of the equation
! sp1 - y(i+1)
! s   - y(i)
! sm1 - y(i-1)
! h   - ds(i)
! hm1 - ds(i-1)
REAL(wp)::RHSCalc
REAL(wp),INTENT(in)::sp1, s, sm1, h, hm1

RHSCalc = 6._wp*((sp1-s)/h-(s-sm1)/hm1)
CALL Errors( RHSCalc, 63 )
END FUNCTION RHSCalc

! Solves a linearSystem
subroutine linear_system(N,A,b)
	integer,intent(in) :: N
    real(wp),dimension(N,N),intent(in) :: A
    real(wp),dimension(N),intent(inout) :: b

    integer,dimension(N) :: ipiv
    integer :: info
    call dgesv(N,1,A,N,ipiv,b,N,info)
    if (info /= 0) stop "Error inverting matrix"
end subroutine linear_system

! Gauss integration (External 1st point, 2nd point)
function gauss1d(f,a,b)
	real(wp), external :: f
    real(wp),intent(in) :: a,b
    real(wp) :: gauss1d

    integer :: i
    real(wp),dimension(5):: zi,wi
    real(wp) :: x

    zi(1) = 0.0_wp
    zi(2) =  (1.0_wp/3.0_wp)*sqrt(5.0_wp-2.0_wp*sqrt(10.0_wp/7.0_wp)) 
    zi(3) = -(1.0_wp/3.0_wp)*sqrt(5.0_wp-2.0_wp*sqrt(10.0_wp/7.0_wp))
    zi(4) =  (1.0_wp/3.0_wp)*sqrt(5.0_wp+2.0_wp*sqrt(10.0_wp/7.0_wp)) 
    zi(5) = -(1.0_wp/3.0_wp)*sqrt(5.0_wp+2.0_wp*sqrt(10.0_wp/7.0_wp))

    wi(1) = 128.0_wp/225.0_wp
    wi(2) = (322.0_wp+13.0_wp*sqrt(70.0_wp))/900.0_wp
    wi(3) = (322.0_wp+13.0_wp*sqrt(70.0_wp))/900.0_wp
    wi(4) = (322.0_wp-13.0_wp*sqrt(70.0_wp))/900.0_wp
    wi(5) = (322.0_wp-13.0_wp*sqrt(70.0_wp))/900.0_wp

    gauss1d = 0.0_wp

    do i=1,5
         x = ((b-a)*0.5_wp)*zi(i) + ((a+b)*0.5_wp) 
         gauss1d = gauss1d + f(x)*wi(i)
    end do
    
    gauss1d = gauss1d*(b-a)*0.5_wp

end function gauss1d

END MODULE Splines_3
