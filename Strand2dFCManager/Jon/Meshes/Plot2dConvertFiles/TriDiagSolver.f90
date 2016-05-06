MODULE StrandSolv
INTEGER, PARAMETER :: p = SELECTED_REAL_KIND(15)
INTEGER::Tag
REAL(8),ALLOCATABLE, DIMENSION(:,:)::XX
REAL(8),PARAMETER:: Infinity = 999999999999999999999999999999999._p
CONTAINS

SUBROUTINE InverseSolver(M,n,Mout)
IMPLICIT NONE
INTEGER,INTENT(IN)::n
REAL(8),DIMENSION(n,n)::M
REAL(8),DIMENSION(n,n)::Mout

INTEGER:: i, ii
REAL(8):: A
CHARACTER(30)::Formating



	!WRITE(Formating,*) "(",n,"(F5.3, X))" 
!Do i=1,n
	!WRITE(*,Formating) M(i,:)
!End DO

WRITE(*,*)
WRITE(*,*)"Solving Inverse Matrix"
WRITE(*,*)
! Make Mout Unity
Mout = 0;
DO i=1,n
	Mout(i,i)=1.
END DO


! Gauss Elimination !
DO i=1,n
	! Make 1 !
	A = M(i,i)
	M(i,:) = M(i,:)/A
	Mout(i,:)=Mout(i,:)/A
	! Subtract From other Areas !
	DO ii =1,n
		IF (ii /= i .AND. M(ii,i) /= 0.) THEN
			A = M(ii,i)
			M(ii,:)=M(ii,:)-A*M(i,:)
			Mout(ii,:)=Mout(ii,:)-A*Mout(i,:)
		END IF
	END DO
	CALL ProgressBar(i,n)
END DO



END SUBROUTINE InverseSolver

SUBROUTINE TriDiagSolver(n,Mo)
INTEGER::n
REAL(8),DIMENSION(n,n)::M,Mo
WRITE(*,*)"Building Tri-Diagnal Matrix ... "
WRITE(*,*)

M = 0
M(1,1)=4.
M(1,2)=1.
M(n,n)=4.
M(n,n-1)=1.
DO i=2,n-1
	M(i,i)=4.
	M(i,i-1)=1.
	M(i,i+1)=1.
	CALL ProgressBar(i,n-1)
END DO
CALL InverseSolver(M,n,Mo)
END SUBROUTINE

SUBROUTINE ProgressBar(i,n)
IMPLICIT NONE
INTEGER,INTENT(IN)::i,n

INTEGER::ii,nn,k
REAL::Prct
CHARACTER(30)::Bar

prct=0.
k   =0
ii  =0
nn  =0

ii=i
nn=n
prct=REAL(ii)/REAL(nn)*100
Bar="???% |                    |"
WRITE(unit=bar(1:3),fmt="(i3)")int(prct)
DO k=1,int(prct*.2)
	bar(6+k:6+k)="="
END DO 

WRITE(unit=6,fmt="(a1,a1,a30)",advance='no') '+',char(13),bar
IF (prct==100.) WRITE(*,*) "done"
RETURN
END SUBROUTINE

SUBROUTINE Splines_3(ptsx,ptsy,n,Nptsx,Nptsy)
IMPLICIT NONE
INTEGER::n
REAL(8),DIMENSION(n),INTENT(in)::ptsx,ptsy
REAL(8),DIMENSION(2*n),INTENT(out)::Nptsx, Nptsy

REAL(8),DIMENSION(n)::yppx,yppy,d
REAL(8)::h,hp,hm,xpt
REAL(8)::sxo,sxn, syo, syn, ds
REAL(8),DIMENSION(n)::sx, sy

INTEGER::i,ii, iii, check
REAL(8),DIMENSION(n-2,n-2)::Mo

WRITE(*,*)"Computing Splines..."
yppx(n)=0.
yppx(1)=0.
yppy(n)=0.
yppy(1)=0.
Tag   =0
ds    =100.
syo   =0.
sxo   =0.
Do WHILE(ds>.001)
! Find Distances s
CALL S_Get(n,ptsx,ptsy,sx,sy,sxn,syn,yppx,yppy)
ds = sqrt((syn-syo)**2)
WRITE(*,'(f10.3)',advance='no') ds
! Find the D's (SEE Numerical Methods for Scientists and Engineers - Hamming pg 351)
DO i=2,N-1
	yppy(i) = 6._p*((ptsy(i+1)-ptsy(i))/sx(i)-(ptsy(i)-ptsy(i-1))/sx(i-1))
	yppx(i) = 6._p*((ptsx(i+1)-ptsx(i))/sy(i)-(ptsx(i)-ptsx(i-1))/sy(i-1))
END DO
! Build inverse Matrix !
Call Build_Mat(ptsy,sx,yppy,n)
Call Build_Mat(ptsx,sy,yppx,n)


END DO

WRITE(*,*)
WRITE(*,*)"Calculating Points"
ii = 1
DO i=1,n-1
	h = ptsx(i+1)-ptsx(i)
	Do iii=1,2
		Nptsx(ii)=ptsx(i)+real(iii)/3.*h
		hp = (1.-(real(iii))/3.)
		hm = (real(iii)/3.)
		Nptsy(ii)=ptsy(i)*hp+ptsy(i+1)*hm-h*h/6.*yppy(i)*(hp-hp**3)-h*h/6.*yppy(i+1)*(hm-hm**3)
		IF (Nptsy(ii)/=Nptsy(ii)) THEN
		WRITE(*,*) "NAN Fault 153"
		WRITE(*,*) hp, hm
		DO check=1,n
			write(*,*) yppy(check)
		END DO
		STOP
		END IF
		ii = ii+1
	END DO
	CALL ProgressBar(i,n-1)
END DO
END SUBROUTINE Splines_3

subroutine linear_system(N,A,b)
	integer,intent(in) :: N
    real(p),dimension(N,N),intent(in) :: A
    real(p),dimension(N),intent(inout) :: b

    integer,dimension(N) :: ipiv
    integer :: info
    call dgesv(N,1,A,N,ipiv,b,N,info)
    if (info /= 0) stop "Error inverting matrix"
end subroutine linear_system

SUBROUTINE S_Get(n,ptsx,ptsy,sx,sy,stx,sty,yppx,yppy)
INTEGER,INTENT(in)::n
REAL(8),DIMENSION(n),INTENT(IN)::ptsx,ptsy,yppx,yppy
REAL(8),DIMENSION(n),INTENT(OUT)::sx,sy
REAL(8),INTENT(OUT)::stx,sty

INTEGER::i
REAL(8)::s
REAL(8),DIMENSION(5,n)::X

stx=0._p
sty=0._p
sx =0._p
sy =0._p

IF (Tag==0) THEN
WRITE(*,*)"INITIALIZED"
DO i=1,n-1
		s    = sqrt((ptsx(i)-ptsy(i))**2+(ptsy(i)-ptsy(i+1))**2)
		stx  = s+stx
		sty  = s+sty
		sx(i)= s
		sy(i)= s
END DO
Tag = 1

ELSE
	ALLOCATE(XX(5,n-1))
	DO i=1,n-1
		X(:,i)=[ptsy(i),ptsy(i+1),ptsx(i+1)-ptsx(i+1),yppy(i),yppy(i+1)]
	END DO
	stx = 0._p
	sty = 0._p
DO i=1,n-1
		Tag  = i
		s    = gauss1d(Passfx,ptsx(i),ptsx(i+1))
		CALL TESTANS(s,214)
		stx  = s+stx
		sx(i)= s
END DO

	DO i=1,n-1
		X(:,i)=[ptsx(i),ptsx(i+1),ptsy(i+1)-ptsy(i+1),yppx(i),yppx(i+1)]
	END DO

DO i=1,n-1
		Tag  = i
		s    = gauss1d(Passfx,ptsy(i),ptsy(i+1))
		CALL TESTANS(s,226)
		sty  = s+sty
		sy(i)= s
END DO
	DEALLOCATE(XX)
END IF

END SUBROUTINE S_Get

function gauss1d(f,a,b)
	real(p), external :: f
    real(p),intent(in) :: a,b
    real(p) :: gauss1d

    integer :: i
    real(p),dimension(5):: zi,wi
    real(p) :: x

    zi(1) = 0.0_p
    zi(2) =  (1.0_p/3.0_p)*sqrt(5.0_p-2.0_p*sqrt(10.0_p/7.0_p)) 
    zi(3) = -(1.0_p/3.0_p)*sqrt(5.0_p-2.0_p*sqrt(10.0_p/7.0_p))
    zi(4) =  (1.0_p/3.0_p)*sqrt(5.0_p+2.0_p*sqrt(10.0_p/7.0_p)) 
    zi(5) = -(1.0_p/3.0_p)*sqrt(5.0_p+2.0_p*sqrt(10.0_p/7.0_p))

    wi(1) = 128.0_p/225.0_p
    wi(2) = (322.0_p+13.0_p*sqrt(70.0_p))/900.0_p
    wi(3) = (322.0_p+13.0_p*sqrt(70.0_p))/900.0_p
    wi(4) = (322.0_p-13.0_p*sqrt(70.0_p))/900.0_p
    wi(5) = (322.0_p-13.0_p*sqrt(70.0_p))/900.0_p

    gauss1d = 0.0_p

    do i=1,5
         x = ((b-a)*0.5_p)*zi(i) + ((a+b)*0.5_p) 
         gauss1d = gauss1d + f(x)*wi(i)
    end do
    
    gauss1d = gauss1d*(b-a)*0.5_p

end function gauss1d

function fprime(x)
REAL(8)::fprime
REAL(8),INTENT(in) ::x


fprime = sqrt(1.+(x)**2)

end function fprime

function Findxy(x)
REAL(8)::x, xhi, xlo,h
REAL(8)::Findxy

xhi = (1._p-x)/2._p
xlo = (x-1._p)/2._p
h   = XX(3,Tag)**2

Findxy = XX(1,Tag)*(xhi)+XX(2,Tag)*(xlo)- h/6._p*XX(4,Tag)*(xhi+xhi**3)-h/6._p*XX(5,Tag)*(xlo+xlo**3)

end function

function TOrderD(func,x)
real(p)::TOrderD
real(p),external::func
real(p)::x
real(p),Dimension(4)::C
C(1)=1._p/12._p
C(2)=-2._p/3._p
C(3)=-C(2)
C(4)=-C(1)
TOrderD = (C(1)*func(x-.002_p)+C(2)*func(x-.001_p)+C(3)*func(x+.001_p)+C(4)*func(x+.002_p))/.001_p
Call TESTANS(TOrderD,298)
end function

function Passfx(x)
REAL(p)::Passfx
REAL(p)::x

Passfx = fprime(TOrderD(Findxy,x))
Call TESTANS(Passfx,306)
end function Passfx

SUBROUTINE Build_Mat(pts,h,ypp,n)
INTEGER::n, Check
REAL(p),DIMENSION(n),INTENT(in)::pts,h
REAL(p),DIMENSION(n),INTENT(OUT)::ypp

REAL(p),DIMENSION(n-2,n-2)::M
INTEGER::i
M = 0._p
M(1,1)=2._p*(h(1)+h(2))
M(1,2)=h(2)
M(n-2,n-2)=2._p*(h(n-2)+h(n-1))
M(n-2,n-3)=h(n-2)

DO i=2,n-3

	M(i,i-1)= h(i-1)
	M(i,i)  = 2._p*(h(i)+h(i-1))
	M(i,i+1)= h(i)

END DO

CALL linear_system(n-2,M,ypp(2:n-1))
IF (ypp(3)/=ypp(3)) THEN
	WRITE(*,*)
	WRITE(*,*) "NAN ERROR 329"
		DO check=1,n
			write(*,*) ypp(check)
		END DO
		DO check=1,n-2
			write(*,*) M(i,:)
		END DO
	STOP
END IF
END SUBROUTINE Build_Mat

SUBROUTINE TESTANS(x,i)
REAL(p)::x
integer::i

if (x /= x) WRITE(*,*)"NAN ERROR ", i
if (x > Infinity) WRITE(*,*)"Inf ERROR ", i

END SUBROUTINE TESTANS

END MODULE StrandSolv
