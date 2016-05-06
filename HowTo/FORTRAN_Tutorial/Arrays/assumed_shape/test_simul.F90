PROGRAM test_simul
!purpose: test subroutine simul, to solve N linear equations in N unknowns
    !data dictionary: constants
    INTEGER,PARAMETER::MAX_SIZE = 10
    !data dictionary: variables
    REAL,DIMENSION(MAX_SIZE,MAX_SIZE) :: a
    REAL,DIMENSION(MAX_SIZE)::b
    INTEGER::error
    CHARACTER(len=20)::file_name
    INTEGER::i,j,n,istat
    !get name of the disk file
    WRITE(*,"(' Enter the file name containing the eqns: ')")
    READ(*,'(A20)')file_name
    !open and see if it exists
    OPEN(UNIT=1,FILE=file_name,STATUS='OLD',ACTION='READ',IOSTAT=istat)
    fileopen: IF (istat==0) THEN
        READ(1,*)n !number of equations in system
        ! if not too many then
        size_ok: IF (n<=MAX_SIZE ) THEN
            DO i=1,n
                READ(1,*) (a(i,j),j=1,n),b(i)
            END DO

            ! display coeffs.
            WRITE(*,"(/,1x,'Coefficients before call:')")
            DO i=1,n
                WRITE(*,"(1x,7F11.4)") (a(i,j),j=1,n),b(i)
            END DO

            ! solve equations
            CALL simul(a,b,MAX_SIZE,n,error)
            !check error
            error_check: IF (error/=0) THEN
                WRITE(*,1010)
                1010 FORMAT (/1x,'Zero pivot encountered!',&
                //1x,'There is no unique solution to this system.')
            ELSE error_check
                ! no errors. display coefficients
                WRITE(*,"(/,1x,'Coefficients after call:')")
                DO i=1,n
                    WRITE(*,"(1x,7F11.4)") (a(i,j),j=1,n),b(i)
                END DO

                ! write final answer
                WRITE(*,"(/,1x,'The solutions are:')")
                DO i=1,n
                    WRITE(*,"(3x,'X(',I2,') = ',F16.6)") i,b(i)
                END DO

            END IF error_check
        END IF size_ok
    ELSE fileopen
        !failed to open
        WRITE(*,1020) istat
        1020 FORMAT (1x,'File open failed --status = ',I6)

    END IF fileopen
END PROGRAM test_simul
