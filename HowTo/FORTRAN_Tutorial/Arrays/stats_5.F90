PROGRAM stats_5
!purpose: calc mean, median, std of an input data set read from a file.  also shows allocatable arrays
IMPLICIT NONE
!data dictionary: variables
REAL,ALLOCATABLE,DIMENSION(:)::a
CHARACTER(len=20)::filename
INTEGER ::i,iptr,j
REAL    ::median
INTEGER ::nvals=0
INTEGER ::status
REAL    ::std_dev
REAL    ::sum_x=0.
REAL    ::sum_x2=0.
REAL    ::temp
REAL    ::x_bar

!get the name of the file containing the input data
WRITE(*,1000)
1000 FORMAT (1x,'Enter the file name with the data to be sorted:')
READ(*,'(A20)') filename

!open input data
OPEN (UNIT=9, FILE=filename, STATUS='OLD',ACTION='READ', IOSTAT=status)
fileopen: IF (status ==0) THEN ! successful open
    DO 
        READ(9,*,IOSTAT=status) temp
        IF (status/=0) EXIT
        nvals=nvals+1
    END DO

    !allocate memory
    WRITE(*,*)' Allocating a: size = ',nvals
    ALLOCATE(a(nvals),STAT=status)
    allocate_ok: IF (status==0) THEN
        REWIND(UNIT=9) !rewind to top of file
        !now read the file to get the array
        READ(9,*) a
        !sort data
        outer: DO i=1,nvals-1
            iptr=i
            inner: DO j=i+1,nvals
                minval: IF (a(j)<a(iptr)) THEN
                    iptr=j
                END IF minval
            END DO inner

            !iptr now points to the min value, so swap a(iptr) with a(i) if i/=iptr
            swap: IF (i/=iptr) THEN
                temp=a(i)
                a(i)=a(iptr)
                a(iptr)=temp
            END IF swap
        END DO outer

        !data is now sorted. get sums and stats
        sums: DO i=1,nvals
            sum_x=sum_x+a(i)
            sum_x2=sum_x2+a(i)**2
        END DO sums

        ! check to see if we have enough input data
        enough: IF (nvals<2) THEN
            WRITE(*,*) ' At least 2 values must be entered.'

        ELSE
            !calc stats
            x_bar=sum_x/REAL(nvals)
            std_dev = sqrt(((REAL(nvals) * sum_x2) - sum_x**2) &
                        / (REAL(nvals) * REAL(nvals-1)))
            even: IF (mod(nvals,2)==0) THEN
                median = (a(nvals/2) + a(nvals/2+1)) /2.
            ELSE 
                median = a(nvals/2+1)
            END IF even

            !output
            WRITE(*,*) ' The mean of this data set is:      ',x_bar
            WRITE(*,*) ' The median of this data set is:    ',median
            WRITE(*,*) ' The standard deviation is:         ',std_dev
            WRITE(*,*) ' The number of the data points is:  ',nvals
        END IF enough

        !deallocate
        DEALLOCATE(a,STAT=status)
    END IF allocate_ok

    ELSE !file failed to open   
        WRITE(*,1050) status
        1050 FORMAT (1x,'File open failed--status = ',I6)

    END IF fileopen

END PROGRAM stats_5
