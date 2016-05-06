!> \brief
!!   This subroutine reads the strand2dFC numerics input file.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2013-1-1
!! \par Further Documentation:
!! \par Source code:
!! \include src/Numerics/strand2dFCInputRead.F90


SUBROUTINE strand2dfcinputread(&
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
     surfOrdersT, &
     strandOrdersT, &
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
INTEGER,INTENT(  OUT),DIMENSION(len) :: surfOrdersT
INTEGER,INTENT(  OUT),DIMENSION(len) :: strandOrdersT
INTEGER,DIMENSION(10) :: surfOrders
INTEGER,DIMENSION(10) :: strandOrders
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
     surfOrders, &
     strandOrders, &
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


surfOrdersT(1:10)   = surfOrders(1:10)
strandOrdersT(1:10) = strandOrders(1:10)


!null terminate the characters for passing back to C++
systemType = TRIM(systemType)//CHAR(0)


END SUBROUTINE strand2dfcinputread
