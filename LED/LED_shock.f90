PROGRAM LED_FirstDWave
IMPLICIT None
! Note:  Function is of the form du/dt+a*du/dx=f(x)
!~  this takes a reconstructed with shockwave form giving a 2nd order accurate values
!Declarations
REAL, ALLOCATABLE, DIMENSION(:,:):: u !wave equation function u
REAL:: a,dphalf,dmhalf,ull,urr,url,ulr,uphalf,umhalf,dx,dt,l,time,x,f,s
INTEGER:: i,j,dxi,dti
OPEN(unit=43,file='LEDOutput.txt')
OPEN(unit=40,file='LEDInput.csv')
!User Input
READ(40,*) ! don't read the first input line
READ(40,*)l,dx,time,dt,a
dxi=int(l/dx+1) ! spatial nodes
dti=int(time/dt) ! time steps
WRITE(*,'(A,I3,A,I3)') 'spatial nodes =', dxi, '   time steps =',dti
ALLOCATE(u(1:dxi,1:dti))
READ(40,*) ! don't read this charater line either
READ(40,*) u(1:dxi,1)
WRITE(43,*) u(1:dxi,1)
!Calc
DO j=1,dti-1,1 ! time steps
	DO i=1,dxi,1 ! spatial steps, consider left end being fixed
!first order forward diff time and backward diff space to get equation
		IF (i==1 .OR. i==2) THEN ! left boundary
			u(i,j+1)=u(i,j)!fixed left boundary (1st and 2nd nodes)
		ELSE ! for every other point
			IF (u(i,j)>u(i-1,j) .AND. u(i,j)>u(i+1,j)) THEN
				s=0 ! first order accurate for maximum
			ELSEIF (u(i-1,j)>u(i-2,j) .AND. u(i-1,j)>u(i,j)) THEN
				s=0 ! first order accurate for left of maximum
			ELSEIF (u(i-2,j)>u(i-3,j) .AND. u(i-2,j)>u(i-1,j)) THEN
				s=0 ! first order accurate for 2 left of maximum
			ELSEIF (u(i+1,j)>u(i,j) .AND. u(i+1,j)>u(i+2,j)) THEN
				s=0 ! first order accurate for right of maximum
			ELSEIF (u(i+2,j)>u(i+1,j) .AND. u(i+2,j)>u(i+3,j)) THEN
				s=0 ! first order accurate for 2 right of maximum
			ELSEIF (u(i,j)<u(i-1,j) .AND. u(i,j)<u(i+1,j)) THEN
				s=0 ! first order accurate for minimum
			ELSEIF (u(i-1,j)<u(i-2,j) .AND. u(i-1,j)<u(i,j)) THEN
				s=0 ! first order accurate for left of minimum
			ELSEIF (u(i-2,j)<u(i-3,j) .AND. u(i-2,j)<u(i-1,j)) THEN
				s=0 ! first order accurate for 2 left of minimum
			ELSEIF (u(i+1,j)<u(i,j) .AND. u(i+1,j)<u(i+2,j)) THEN
				s=0 ! first order accurate for right of minimum
			ELSEIF (u(i+2,j)<u(i+1,j) .AND. u(i+2,j)<u(i+3,j)) THEN
				s=0 ! first order accurate for 2 right of minimum
			ELSE
				s=1 ! reconstructed LED for everythihng else
			END IF
			urr=u(i+1,j)-s*dx/2*(u(i+2,j)-u(i,j))/(2*dx)
			url=u(i,j)+s*dx/2*(u(i+1,j)-u(i-1,j))/(2*dx)
			ull=u(i-1,j)+s*dx/2*(u(i-2,j)-u(i,j))/(2*dx)
			ulr=u(i,j)-s*dx/2*(u(i+1,j)-u(i-1,j))/(2*dx)
			dphalf=.5*(a)*(urr-url) !no absolute value to handle negative a
			dmhalf=.5*(a)*(ulr-ull)
			umhalf=.5*abs(a)*(u(i,j)+u(i-1,j))-dmhalf
			uphalf=.5*abs(a)*(u(i,j)+u(i+1,j))-dphalf
			x=(i-1)*dx ! to calculate the source x value
			f=source(x) ! to calculate the source value
			u(i,j+1)=-a*dt*(uphalf-umhalf)/dx+u(i,j)+f*dt
		END IF
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
source=(5.*x-x**2)/64.
END FUNCTION

END PROGRAM 
