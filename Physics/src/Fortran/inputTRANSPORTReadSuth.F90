!> \brief
!!   This subroutine reads the Sutherland Law input file.
!! \param fileLen Length of character string of input file.
!! \param fileName Name of input file.
!! \param comp Component number for this Sutherland Law instance.
!! \param prnN Prandlt number for component comp.
!! \param prnNT Turbulent Prandlt number for component comp..
!! \param t0N Sutherland reference temperature for component comp.
!! \param mu0N Sutherland reference viscosity for component comp.
!! \param sN Sutherland constant for component comp.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-9-13
!! \par Further Documentation:
!! \par Source code:
!! \include src/Fortran/inputTRANSPORTReadSuth.F90


SUBROUTINE inputtransportreadsuth(&
     fileLen, &
     fileName, &
     comp, &
     prnN, &
     prnNT, &
     t0N, &
     mu0N, &
     sN)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: fileLen
CHARACTER(fileLen),INTENT(IN   ) :: fileName
INTEGER,INTENT(IN   ) :: comp
REAL   ,INTENT(  OUT) :: prnN
REAL   ,INTENT(  OUT) :: prnNT
REAL   ,INTENT(  OUT) :: t0N
REAL   ,INTENT(  OUT) :: mu0N
REAL   ,INTENT(  OUT) :: sN
REAL   ,DIMENSION(500) :: prn
REAL   ,DIMENSION(500) :: prnT
REAL   ,DIMENSION(500) :: t0
REAL   ,DIMENSION(500) :: mu0
REAL   ,DIMENSION(500) :: s
INTEGER,PARAMETER :: iUnit = 11
NAMELIST/sutherland/&
     prn, &
     prnT, &
     t0, &
     mu0, &
     s


OPEN(iUnit,FILE=fileName,STATUS='old')
READ(iUnit,sutherland)
CLOSE(iUnit)

prnN  = prn(comp)
prnNT = prnT(comp)
t0N   = t0 (comp)
mu0N  = mu0(comp)
sN    = s  (comp)


END SUBROUTINE inputtransportreadsuth
