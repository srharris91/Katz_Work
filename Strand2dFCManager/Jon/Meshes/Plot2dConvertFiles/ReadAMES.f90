! READS UIUC Airfoil Data and organizes it into a usefull X,Y array
! Jon Thorne July 14. 2014

MODULE Airfoil
USE Pre
IMPLICIT NONE


CONTAINS

SUBROUTINE ReadData_UIUC(inputfile,Array,n)
CHARACTER(40),INTENT(in)::inputfile
REAL(wp),DIMENSION(:,:),ALLOCATABLE::Array
CHARACTER(40)::newFile
INTEGER,INTENT(out)::n

INTEGER::n1,n2,i,c

WRITE(*,*)'READING ',TRIM(ADJUSTL(inputfile))

OPEN(UNIT=17,FILE=TRIM(ADJUSTL(inputfile)),ACTION='READ')
READ(17,*)newFile
WRITE(*,*)'New File Name: ',TRIM(ADJUSTL(newFile)),'.mesh'
WRITE(newFile,*) TRIM(ADJUSTL(newFile)),'.mesh'
READ(17,*)n1,n2
n = n1+n2
Write(*,*)n1-1,n2-1,n
ALLOCATE(Array(n-1,2))
c = 3
DO i=n1,1,-1
	Write(*,*)c
	READ(17,*) Array(i,:)
	c = c+1
END DO

DO i=n1, n-1
	Write(*,*)c
	READ(17,*) Array(i,:)
	c = c+1
END DO

DO i=1,n-1
	WRITE(*,*) Array(i,:)
END DO
CLOSE(17)
n= n-1
END SUBROUTINE ReadData_UIUC

END MODULE Airfoil
