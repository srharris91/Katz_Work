!> \brief
!!   This subroutine reads the ideal gas input file.
!! \param fileLen Length of character string of input file.
!! \param fileName Name of input file.
!! \param comp Component number for this ideal gas instance.
!! \param rGasN Ideal gas constant for component comp.
!! \param gammaN Ratio of specific heats for component comp.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-9-13
!! \par Further Documentation:
!! \par Source code:
!! \include src/Fortran/inputSTATEReadIdealGas.F90


SUBROUTINE inputstatereadidealgas(&
     fileLen, &
     fileName, &
     comp, &
     rGasN, &
     gammaN)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: fileLen
CHARACTER(fileLen),INTENT(IN   ) :: fileName
INTEGER,INTENT(IN   ) :: comp
REAL   ,INTENT(  OUT) :: rGasN
REAL   ,INTENT(  OUT) :: gammaN
REAL   ,DIMENSION(500) :: rGas
REAL   ,DIMENSION(500) :: gamma
INTEGER,PARAMETER :: iUnit = 11
NAMELIST/idealGas/&
     rGas, &
     gamma


OPEN(iUnit,FILE=fileName,STATUS='old')
READ(iUnit,idealGas)
CLOSE(iUnit)

rGasN  = rGas(comp)
gammaN = gamma(comp)


END SUBROUTINE inputstatereadidealgas
