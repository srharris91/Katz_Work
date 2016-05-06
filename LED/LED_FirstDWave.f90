PROGRAM LED_FirstDWave
IMPLICIT None
! Note:  Function is of the form du/dt+a*du/dx=0
!Declarations
REAL, ALLOCATABLE, DIMENSION(:,:):: u !wave equation function u
REAL:: a,dphalf,dmhalf,ul,ur,dx,dt,l,time
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
		IF (i==1) THEN ! left boundary
			u(i,j+1)=a*dt/dx*(0.-u(i,j))+u(i,j)!imaginary node to be 0
		ELSEIF (i==dxi) THEN ! right boundary
			u(i,j+1)=a*dt/dx*(u(i-1,j)-u(i,j))+u(i,j)
		ELSE
			u(i,j+1)=a*dt/dx*(u(i-1,j)-u(i,j))+u(i,j) 
		ENDIF
	END DO
	WRITE(43,*) u(1:dxi,j+1)!output
END DO
CLOSE(43)
CLOSE(40)
DEALLOCATE(u)
END PROGRAM
