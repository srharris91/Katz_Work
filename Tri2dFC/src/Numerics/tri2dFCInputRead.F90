!> \brief
!!   This subroutine reads the tri2dFC numerics input file.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2013-1-1
!! \par Further Documentation:
!! \par Source code:
!! \include src/Numerics/tri2dFCInputRead.F90


SUBROUTINE tri2dfcinputread(&
     fileLen, &
     fileName, &
     len, &
     iPrint, &
     iTest, &
     iDebug, &
     iConvFile, &
     iSolnFile, &
     iResdFile, &
     iErrFile, &
     iSurfFile, &
     standAlone, &
     restartStep, &
     nRestart, &
     nOutput, &
     nSteps, &
     nPseudoSteps, &
     nPseudoSteps0, &
     nLinearSteps, &
     nRKStages, &
     implicit, &
     nLevels, &
     ordersT, &
     spacing, &
     gradMethod, &
     mgCycle, &
     limiter, &
     timeAcc, &
     dtUnsteady, &
     cfl, &
     vnn, &
     smooth, &
     convLimit, &
     relax, &
     systemType)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: fileLen,len
CHARACTER(fileLen),INTENT(IN   ) :: fileName
CHARACTER(    len),INTENT(  OUT) :: systemType
INTEGER,INTENT(  OUT),DIMENSION(len) :: ordersT
INTEGER,DIMENSION(3) :: orders
INTEGER,INTENT(  OUT) :: &
     iPrint, &
     iTest, &
     iDebug, &
     iConvFile, &
     iSolnFile, &
     iResdFile, &
     iErrFile, &
     iSurfFile, &
     standAlone, &
     restartStep, &
     nRestart, &
     nOutput, &
     nSteps, &
     nPseudoSteps, &
     nPseudoSteps0, &
     nLinearSteps, &
     nRKStages, &
     implicit, &
     nLevels, &
     spacing, &
     gradMethod, &
     mgCycle, &
     limiter, &
     timeAcc
REAL,INTENT(  OUT) :: &
     dtUnsteady, &
     cfl, &
     vnn, &
     smooth, &
     convLimit, &
     relax
INTEGER,PARAMETER :: iUnit = 9
NAMELIST/numerics/&
     iPrint, &
     iTest, &
     iDebug, &
     iConvFile, &
     iSolnFile, &
     iResdFile, &
     iErrFile, &
     iSurfFile, &
     standAlone, &
     restartStep, &
     nRestart, &
     nOutput, &
     nSteps, &
     nPseudoSteps, &
     nPseudoSteps0, &
     nLinearSteps, &
     nRKStages, &
     implicit, &
     nLevels, &
     orders, &
     spacing, &
     gradMethod, &
     mgCycle, &
     limiter, &
     timeAcc, &
     dtUnsteady, &
     cfl, &
     vnn, &
     smooth, &
     convLimit, &
     relax, &
     systemType


OPEN(iUnit,FILE=TRIM(filename),STATUS='old')
READ(iUnit,numerics)
CLOSE(iUnit)


WRITE(*,*)
WRITE(*,numerics)
WRITE(*,*)


ordersT(1:nLevels) = orders(1:nLevels)


!null terminate the characters for passing back to c++
systemType = TRIM(systemType)//CHAR(0)


END SUBROUTINE tri2dfcinputread
