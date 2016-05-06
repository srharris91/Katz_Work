PROGRAM read_file
! purpose: illustrate how to read an unknown number of values from an input data file.  detecting both any formatting errors and the end fo the file

IMPLICIT NONE

!data dictionary
CHARACTER(len=20)::filename ! name of file to open
INTEGER::nvals=0            !number of values to read in
INTEGER::status             !I/O status
REAL::value                 !The value to read in

!get the file name and echo it back to the user
WRITE(*,*)'Please enter input file name:'
READ(*,*)filename
WRITE(*,1000)filename
1000 FORMAT (' ','The input file name is: ',A)

!open the file, and check for error on open.
OPEN(UNIT=3, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=status)
openif: IF(status==0) THEN
    !open was ok, let's read the values
    readloop: DO
        READ(3,*,IOSTAT=status) value   ! get the next value
        IF (status/=0) EXIT             ! exit if not valid
        nvals = nvals+1                 ! valid: increase count
        WRITE(*,1010) nvals, value      ! Echo to screen
        1010 FORMAT (' ','Line ', I6, ': Value = ',F10.4)
    END DO readloop

    ! the readloop has now terminated, was it because the end of file?  or an error?
    readif: IF(status > 0) THEN
            WRITE(*,1020) nvals + 1
            1020 FORMAT ('0', 'An error occurred reading line ',I6)
        ELSE
            WRITE(*,1030) nvals
            1030 FORMAT ('0','End of file reached.  There were ',I6, ' values in the file.')
    END IF readif

    ELSE openif
        WRITE(*,1040) status
        1040 FORMAT (' ','Error opening file: IOSTAT = ',I6)
END IF openif

!close file
CLOSE (unit=3)

END PROGRAM read_file
