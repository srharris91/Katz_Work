!> \brief
!!   This subroutine reads the perturb solution input file.
!! \param fileLen Length of character string of input file.
!! \param fileName Name of input file.
!! \param nq Dimension of q (number of equations).
!! \param pert Amount by which to perturb solution.
!! \param rValue Solution initialization values.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-9-13
!! \par Further Documentation:
!! \par Source code:
!! \include src/Fortran/inputSOLUTIONReadPerturb.F90


SUBROUTINE inputsolutionreadperturb(&
     fileLen, &
     fileName, &
     nq, &
     pert, &
     rValue0)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: fileLen
CHARACTER(fileLen),INTENT(IN   ) :: fileName
INTEGER,INTENT(IN   ) :: nq
REAL   ,INTENT(  OUT) :: pert
REAL   ,INTENT(  OUT),DIMENSION(nq) :: rValue0
REAL   ,DIMENSION(100) :: rValue
INTEGER,PARAMETER :: iUnit = 11
NAMELIST/perturb/&
     rValue, &
     pert


OPEN(iUnit,FILE=fileName,STATUS='old')
READ(iUnit,perturb)
CLOSE(iUnit)

rValue0(1:nq) = rValue(1:nq)


END SUBROUTINE inputsolutionreadperturb
