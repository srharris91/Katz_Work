PROGRAM extremes

! purpose: find largest and smallest values in a data set

IMPLICIT NONE

!data dictionary constants
INTEGER, PARAMETER::MAX_SIZE=10

!data dictionary: variables
INTEGER,DIMENSION(MAX_SIZE)::input
INTEGER::ilarge,ismall,j,nvals,temp

!Get number of values in data set
WRITE(*,*)'Enter number of values in data set:'
READ(*,*) nvals

size: IF(nvals<=MAX_Size) THEN
    !get input values
    in: DO j=1,nvals
        WRITE(*,100) 'Enter value ',j
        100 FORMAT (' ',A,I3,': ')
        READ(*,*)input(j)
    END DO in

    !find largest value
    temp=input(1)
    ilarge=1
    large: DO j=2,nvals
        IF (input(j)>temp) THEN
            temp=input(j)
            ilarge=j
        END IF
    END DO large

    !find smallest value
    temp=input(1)
    ismall =1
    small: DO j=2,nvals
        IF (input(j) < temp) THEN
            temp=input(j)
            ismall=j
        END IF
    END DO small

    !write out list
    WRITE(*,110)
    110 FORMAT ('0','The values are:')
    out: DO j=1,nvals
        IF (j==ilarge) THEN
            WRITE(*,'(1x,I6,2x,A)') input(j), 'LARGEST'
        ELSE IF (j==ismall) THEN
            WRITE(*,'(1x,I6,2x,A)') input(j), 'SMALLEST'
        ELSE
            WRITE(*,'(1x,I6)') input(j)
        END IF
    END DO out
ELSE size
    WRITE(*,120) nvals, MAX_SIZE
    120 FORMAT (1x, 'Too many input values: ',I6,' > ',I6)
END IF size


            

END PROGRAM extremes
