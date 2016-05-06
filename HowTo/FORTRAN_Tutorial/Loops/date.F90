PROGRAM date

!purpose: calculates the day of the year corresponding to a specified date. Illustrates counting loops and select case construct

IMPLICIT NONE

!Data dictionary
INTEGER::day,day_of_year,i,leap_day,month,year !dd,mm,yyyy

! get day,month, and year to convert
WRITE(*,*)'Enter the current month (1-12), day(1-31), and year (yyyy) in that order'
READ(*,*)month,day,year

!check for leap year and add extra day if necessary
IF (MOD(year,400) == 0) THEN
    leap_day = 1    !years divisible by 400 are leap years
    WRITE(*,*) 'It is a leap year!'
ELSE IF (MOD(year,100) == 0) THEN
    leap_day = 0    !other centuries are not leap years
ELSE IF (MOD(year,4) == 0) THEN
    leap_day = 1    !otherwise every 4th year is a leap year
    WRITE(*,*) 'It is a leap year!'
ELSE
    leap_day = 0    !other years are not leap years
END IF

!calc day of year
day_of_year = day
DO i=1,month-1
    !Add days in months from January to last month
    SELECT CASE (i)
    CASE (1,3,5,7,8,10,12)
        day_of_year = day_of_year + 31
    CASE (4,6,9,11)
        day_of_year = day_of_year + 30
    CASE (2)
        day_of_year = day_of_year + 28 + leap_day
    END SELECT
END DO

!Tell user
WRITE(*,*)'Day      =   ',day
WRITE(*,*)'Month    =   ',month
WRITE(*,*)'Year     =   ',year
WRITE(*,*)'day of year= ',day_of_year


END PROGRAM date
