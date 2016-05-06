PROGRAM LED_FirstDWave
IMPLICIT None
! Note:  Function is of the form du/dt+a*du/dx=f(x)
!~  this takes a reconstructed with shockwave form giving a 2nd order accurate values
!Declarations
REAL, ALLOCATABLE, DIMENSION(:,:):: u !wave equation function u
REAL:: a,dphalf,dmhalf,ull,urr,url,ulr,uphalf,umhalf,dx,dt,l,time,x
REAL:: xmh,xph,trapS
INTEGER:: i,j,dxi,dti
OPEN(unit=43,file='LEDOutput.txt')
OPEN(unit=40,file='LEDInput.csv')
!User Input
READ(40,*) ! don't read the first input line
READ(40,*)l,dx,time,dt,a
dxi=int(l/dx+1) ! spatial nodes
dti=int(time/dt+1) ! time nodes
WRITE(*,'(A,I3,A,I3)') 'spatial nodes =',dxi, '   time nodes =',dti
ALLOCATE(u(-1:(dxi+2),dti)) ! add 4 to spatial nodes to get boundary values
READ(40,*) ! don't read this charater line either
!~ READ(40,*) u(1:dxi,1) ! initial values in Input.csv file
READ(40,*) ! Read initial values in Input.csv file as empty line
DO i=1,dxi,1 ! if have function for initial values
	x=real(i-1)*dx ! to calculate the x value on the bar
	u(i,1)=IV(x)
END DO
u(-1:0,:)=0. ! outside boundary values or ghost nodes, left side
u((dxi+1):(dxi+2),:)=0. ! outside boundary values or ghost nodes, right side
WRITE(43,*) u(1:dxi,1)
!Calc
DO j=1,dti-1,1 ! time steps
	DO i=1,dxi,1 ! each dx node on bar (off bar values were noted above)
!LED reconstructed method
! sign(1.,a) to handle negative a values
			urr=u(i+1,j)-sign(1.,a)*dx/2.*(u(i+2,j)-u(i,j))/(2.*dx)
			url=u(i,j)+sign(1.,a)*dx/2.*(u(i+1,j)-u(i-1,j))/(2.*dx)
			ull=u(i-1,j)+sign(1.,a)*dx/2.*(u(i-2,j)-u(i,j))/(2.*dx)
			ulr=u(i,j)-sign(1.,a)*dx/2.*(u(i+1,j)-u(i-1,j))/(2.*dx)
			!no absolute value in dphalf or dmhalf to handle negative a
			dphalf=.5*(a)*(urr-url)
			dmhalf=.5*(a)*(ulr-ull)
			umhalf=.5*abs(a)*(u(i,j)+u(i-1,j))-dmhalf
			uphalf=.5*abs(a)*(u(i,j)+u(i+1,j))-dphalf
! single trapezoidal of source term from x(i-1/2) to x(i+1/2)
			xmh=(real(i-1)-.5)*dx
			xph=(real(i-1)+.5)*dx
			trapS=.5*(xph-xmh)*(source(xph)+source(xmh))
			u(i,j+1)=-a*dt*(uphalf-umhalf)/dx+u(i,j)+trapS*dt/dx
	END DO
	WRITE(43,*) u(1:dxi,j+1)!output
END DO
CLOSE(43)
CLOSE(40)
DEALLOCATE(u)

CONTAINS
FUNCTION source(x) ! this is for the source function and term
IMPLICIT NONE
REAL:: source,x
source=-1./5.*(x-5.)
END FUNCTION
FUNCTION IV(x) ! this is for the initial values if it is a function
IMPLICIT NONE
REAL::IV,x
IV=cos(2.*x)
END FUNCTION


END PROGRAM 
