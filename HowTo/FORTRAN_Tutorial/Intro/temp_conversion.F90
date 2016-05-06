PROGRAM temp_conversion
! purpose: convert from Fahrenheit to kelvin

IMPLICIT NONE   !Force explicit declaration of variables

! Data dictionary: declare variables types, definitions, & units
REAL :: temp_f  ! temp in fahrenheit
REAL :: temp_k  ! temp in kelvin

! prompt user for input temp
WRITE(*,*) "Enter the tempreature in degrees Fahrenheit: "
READ (*,*) temp_f

! calc by converting F to K
temp_k = (5./9.) * (temp_f - 32.) + 273.15

! write output
WRITE(*,*) temp_f,' degrees Fahrenheit = ',temp_k,' kelvins'

!Finish up
END PROGRAM temp_conversion
