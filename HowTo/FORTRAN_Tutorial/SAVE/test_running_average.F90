PROGRAM test_running_average
    !purpose: test running average subroutine
    IMPLICIT NONE
    !declare variables:
    INTEGER::istat
    REAL::ave,std_dev,x
    INTEGER::nvals
    CHARACTER(len=20) file_name
    !clear running sums
    CALL running_average(0.,ave,std_dev,nvals,.TRUE.)

    ! get name of file containing the input data
    WRITE(*,*)' Enter the file name containing the data: '
    READ(*,'(A20)')file_name
    OPEN (UNIT=21, FILE=file_name, STATUS='OLD', ACTION='READ',IOSTAT=istat)
    openok: IF (istat==0) THEN
        calc: DO
            READ(21,*,IOSTAT=istat) x
            IF (istat/=0) EXIT
            CALL running_average(x,ave,std_dev,nvals,.FALSE.)
            WRITE(*,1020) ' Value = ',x,' Ave = ',ave, &
                          ' Std_dev = ',std_dev, &
                          ' Nvals = ',nvals
            1020 FORMAT (1x,3(A,F10.4),A,I6)
        END DO calc
    ELSE openok
        WRITE(*,1030)istat
        1030 FORMAT (1x,'File open failed -- status = ',I6)
    END IF openok

END PROGRAM test_running_average
