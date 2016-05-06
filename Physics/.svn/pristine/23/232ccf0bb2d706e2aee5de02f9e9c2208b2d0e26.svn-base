!> \brief
!!   This subroutine reads the incompressible input file.
!! \param fileLen Length of character string of input file.
!! \param fileName Name of input file.
!! \param comp Component number for this fluid instance.
!! \param rhoN Constant density for component comp.
!! \param CpN Constant pressure specific heat for component comp.
!! \param CvN Constant volume specific heat for component comp.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-9-13
!! \par Further Documentation:
!! \par Source code:
!! \include src/Fortran/inputSTATEReadIncompressible.F90


SUBROUTINE inputstatereadincompressible(&
     fileLen, &
     fileName, &
     comp, &
     rhoN, &
     CpN, &
     CvN)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: fileLen
CHARACTER(fileLen),INTENT(IN   ) :: fileName
INTEGER,INTENT(IN   ) :: comp
REAL   ,INTENT(  OUT) :: rhoN
REAL   ,INTENT(  OUT) :: CpN
REAL   ,INTENT(  OUT) :: CvN
REAL   ,DIMENSION(500) :: rho
REAL   ,DIMENSION(500) :: Cp
REAL   ,DIMENSION(500) :: Cv
INTEGER,PARAMETER :: iUnit = 11
NAMELIST/incompressible/&
     rho, &
     Cp, &
     Cv


OPEN(iUnit,FILE=fileName,STATUS='old')
READ(iUnit,incompressible)
CLOSE(iUnit)

rhoN = rho(comp)
CpN  = Cp (comp)
CvN  = Cv (comp)


END SUBROUTINE inputstatereadincompressible
