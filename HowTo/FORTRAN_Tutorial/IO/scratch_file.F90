PROGRAM scratch_file

!purpose show scratch file use and file positioning.  Will read in arbitrary number, save to scratch file, and stop when a negative is read.  Ask user for a record number to display, rewinds the file and get that value and display it


IMPLICIT NONE

!Data dictionary, constants
INTEGER, PARAMETER::LU=8! io unit for scratch file

!Data dictionary, variables
REAL::data          ! data value stored in a disk file
INTEGER::icount=0   ! number of input data records
INTEGER::irec       ! record number to recover and display
INTEGER::j          ! loop index

!Open scratch file
OPEN (UNIT=LU, STATUS='SCRATCH')

!Prompt user for input data
WRITE(*,100)
100 FORMAT (1x,'Enter Postivie or zero input values. ',/, &
            1x,'A negative value terminates input.')

! get the input values, and write them to the scratch file
DO 
    WRITE(*,110) icount+1 !prompt for next value
    110 FORMAT (1x,'Enter sample ',I4,':')
    READ(*,*) data      ! read data
    IF(data<0.) EXIT    !if negative then exit
    icount=icount+1     ! valid value: bump count
    WRITE(LU,120) data  ! write data to scratch file
    120 FORMAT (1x, ES16.6)
END DO

! now we have all the records, ask which record to see
WRITE(*,130) icount
130 FORMAT (1x,'Which record do you want to see (1 to ',I4,' )?')
READ (*,*) irec

! do we have a legal record number?  if so, get the record. if not, tell the user and stop
IF((irec>=1) .AND. (irec<=icount) ) THEN
    REWIND(UNIT=LU) ! rewind the scratch file (start from top)
    DO j=1,irec !read forward to the desirec record
        READ(LU,*) data !read all values until the desired record
    END DO
    WRITE(*,140) irec, data
    140 FORMAT (1x, 'The value of the record ',I4, ' is ',ES14.5)
ELSE !illegal record
    WRITE(*,150) irec
        150 FORMAT (1x,'Illegal record number entered: ',I8)
END IF

! close file
CLOSE(LU)

END PROGRAM scratch_file
