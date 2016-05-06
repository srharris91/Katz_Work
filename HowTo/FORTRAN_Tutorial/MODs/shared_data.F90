MODULE shared_data

!purpose: to show how to share data between two routines

IMPLICIT NONE
SAVE
INTEGER, PARAMETER::num_vals=5
REAL,DIMENSION(num_vals)::values

END MODULE shared_data
