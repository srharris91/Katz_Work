PROGRAM generate
!purpose: calc total instantaneous power supplied by station. Show use of 2-d array

IMPLICIT NONE
!data dictionary: constants
INTEGER,PARAMETER::MAX_GEN=4,MAX_TIME=6

!data dictionary: variables
CHARACTER(len=20)::filename
INTEGER::igen,itime
REAL,DIMENSION(MAX_TIME,MAX_GEN)::power
REAL,DIMENSION(MAX_GEN)::power_ave
REAL,DIMENSION(MAX_TIME)::power_sum
INTEGER::status

!initialize to zero
power_ave=0.
power_sum=0.

!get name of file containing the input data
WRITE(*,1000)
1000 FORMAT (' Enter the file name containing the input data: ')
READ (*,'(A20)') filename 

! open input data file
OPEN (UNIT=9, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=status)

fileopen: IF (status ==0) THEN !status was opened successfully
!    READ(9,*,IOSTAT=status)power!this tranposes the read in file
    READ(9,*,IOSTAT=status) ((power(itime,igen),igen=1,MAX_GEN),itime=1,MAX_TIME) !intrinsic do loops to read in file without transposing it

    !calc power output
    sum1: DO itime=1,MAX_TIME
        sum2: DO igen=1,MAX_GEN
            power_sum(itime)=power(itime,igen)+power_sum(itime)
        END DO sum2
    END DO sum1

    !calc ave output
    ave1: DO igen=1, MAX_GEN
        ave2: DO itime=1,MAX_TIME
            power_ave(igen)=power(itime,igen)+power_ave(igen)
        END DO ave2
        power_ave(igen)=power_ave(igen)/REAL(MAX_TIME)
    END DO ave1

    ! tell user
    out1: DO itime=1,MAX_TIME
        WRITE(*,1010) itime, power_sum(itime)
        1010 FORMAT (' THe instantaneous power at time ', I1, ' is ',F7.2,' MW.')
    END DO out1

    out2: DO igen =1,MAX_GEN
        WRITE(*,1020) igen, power_ave(igen)
        1020 FORMAT (' The average power of generator ', I1, ' is ',F7.2, ' MW.')
    END DO out2
ELSE fileopen !file failed to open
    WRITE(*,1030) status
    1030 FORMAT (1x,'File open failed--status = ',I6)
END IF fileopen

END PROGRAM generate
