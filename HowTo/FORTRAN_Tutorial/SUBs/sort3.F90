PROGRAM sort3
!purpose: read real input data set, and sort it into ascending order using the selection sort algorithm.  Then write data to output, this calls the sort subroutine

IMPLICIT NONE

!data dictionary: constants
INTEGER,PARAMETER::MAX_SIZE=10

!data dictionary: variables
REAL,DIMENSION(MAX_SIZE)::a
LOGICAL::exceed=.FALSE.
CHARACTER(len=20)::filename
INTEGER::i,nvals=0,status
REAL::temp

!get name of the file containing input data
WRITE(*,*)'Enter the file name with the data to be sorted: '
READ(*,1000)filename
1000 FORMAT (A20)

!open input file, 
OPEN (UNIT=9, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=status)

!was open successful?
fileopen: IF(status==0) THEN
    !file was open, so read it, sort it, and write out
    DO
        READ(9,*,IOSTAT=status) temp
        IF (status/=0) EXIT     !exit on end of data
        nvals=nvals+1
        size: IF (nvals<=MAX_SIZE) THEN
            a(nvals)=temp
        ELSE
            exceed = .TRUE.
        END IF size
    END DO

    !was the array size exceeded?
    toobig: IF(exceed) THEN !yes too big
        WRITE(*,1010) nvals,MAX_SIZE
        1010 FORMAT (' Maximum array size exceeded: ',I6,' > ',I6)
    ELSE                    !no go ahead and sort the data
        CALL sort(a,nvals)
        WRITE(*,1020)' The sorted output data values are: '
        1020 FORMAT (A)
        WRITE(*,1030) (a(i), i=1,nvals)
        1030 FORMAT (4x, F10.4)
    END IF toobig
ELSE 
    !file failed to open
    WRITE(*,1040) status
    1040 FORMAT (1x, 'File open failed -- status = ',I6)
END IF fileopen

END PROGRAM sort3



!****************************************************************************
!****************************************************************************
!****************************************************************************


SUBROUTINE sort (arr,n)

!purpose: sort array "arr" into ascending order using a selection sort

IMPLICIT NONE

!data dictionary: parameters
INTEGER,INTENT(IN)::n
REAL,DIMENSION(n),INTENT(INOUT)::arr

!data dictionary: variables
INTEGER::i,iptr,j
REAL::temp

!sort array
outer: DO i=1,n-1
    !find min value in arr(i) through arr(n)
    iptr=i
    inner: DO j=i+1,n
        minval: IF (arr(j)<arr(iptr)) THEN
            iptr=j
        END IF minval
    END DO inner

    ! iptr now points to the minimum value
    swap: IF (i/=iptr) THEN
        temp = arr(i)
        arr(i)=arr(iptr)
        arr(iptr)=temp
    END IF swap

END DO outer

END SUBROUTINE sort
