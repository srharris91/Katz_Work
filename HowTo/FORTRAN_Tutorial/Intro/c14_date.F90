PROGRAM c14_date

! Purpose:
!       Calculate the age of an organic sample from the percentage of the origincal carbon 14 remaining in the sample
IMPLICIT NONE

REAL,PARAMETER::LAMDA = 0.00012097  ! The radioactive decay constant of carbon 14 in units of 1/years
!Data dictionary
REAL::age   !age of the sample in years
REAL::percent   !percent of carbon 14 remaining at the time of the measurement in percent
REAL::ratio     !ratio of the carbon 14 remaining at the time of the measurement to the original amount of carbon 14 (no units)

!Prompt the user for percentage of c-14 remaining
WRITE(*,*) "Enter the percentage of carbon 14 remaining"
READ(*,*) percent
WRITE(*,*) "The remaining carbon 14 = ",percent, "%."

!perform calculations
ratio = percent/100.
age = (-1.0/LAMDA) * log(ratio) ! get age in years

!output
WRITE(*,*)"The age of the sample is = ",age, "years"

!Finish up
END PROGRAM c14_date
