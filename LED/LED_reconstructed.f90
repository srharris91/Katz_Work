PROGRAM LED_FirstDWave
IMPLICIT None
! Note:  Function is of the form du/dt+a*du/dx=0
!~  this takes a reconstructed form giving a 2nd order accurate values
!Declarations
REAL, ALLOCATABLE, DIMENSION(:,:):: u !wave equation function u
REAL:: a,dphalf,dmhalf,ull,urr,url,ulr,uphalf,umhalf,dx,dt,l,time
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
			urr=u(i+1,j)-dx/2*(u(i+2,j)-u(i,j))/(2*dx)
			url=u(i,j)+dx/2*(u(i+1,j)-u(i-1,j))/(2*dx)
			ull=u(i-1,j)+dx/2*(u(i-2,j)-u(i,j))/(2*dx)
			ulr=u(i,j)-dx/2*(u(i+1,j)-u(i-1,j))/(2*dx)
			dphalf=.5*(a)*(urr-url) !no absolute value to handle negative a
			dmhalf=.5*(a)*(ulr-ull)
			umhalf=.5*abs(a)*(u(i,j)+u(i-1,j))-dmhalf
			uphalf=.5*abs(a)*(u(i,j)+u(i+1,j))-dphalf
			u(i,j+1)=-a*dt*(uphalf-umhalf)/dx+u(i,j)
		END IF
	END DO
	WRITE(43,*) u(1:dxi,j+1)!output
END DO
CLOSE(43)
CLOSE(40)
DEALLOCATE(u)
END PROGRAM 