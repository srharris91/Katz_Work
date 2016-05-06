!> \brief
!!   This subroutine reads the strand numerics input file.
!! \param fileLen Length of the character string containing fileName.
!! \param fileName Name of input file.
!! \param len Length of systemType name
!! \param systemType type of System layer to use for this block.
!! \param iPrint Flag to indicate a verbose mode.
!! \param iTest Flag to indicate unit test mode.
!! \param iDebug Flag to indicate debug mode.
!! \param iConvFile Flag to include temporal convergence file.
!! \param iSolnFile Flag to include solution field output file.
!! \param iResdFile Flag to include residual field output file.
!! \param iErrFile Flag to include error field output file (used for exact
!! or manufactured solutions).
!! \param iSurfFile Flag to include surface quantity output field file.
!! \param standAlone Flag indicating if solver is running in stand alone mode.
!! \param restartStep Physical time restart step.
!! \param nRestart Interval of restart step output files.
!! \param nOutput Interval to write output files.
!! \param nSteps Number of physical time steps.
!! \param nPseudoSteps Number of pseudo-time steps.
!! \param nPseudoSteps0 Number of pseudo-time steps at the first physical
!! time step.
!! \param nLinearSteps Number of linear-time steps.
!! \param nLevels Number of multigrid levels.
!! \param mgCycle Flag to determine multigrid cycle topology.
!! \param gradient Flag to indicate method of gradient computation.
!! \param nodeVal Flag to indicate method of nodal value computation.
!! \param perturb Flag to perturb nodes (used for accuracy testing).
!! \param limiter Flag to turn limiter on or off.
!! \param dtUnsteady Unsteady time step.
!! \param cflLinear CFL number for linear-time step definition.
!! \param cfl0 CFL number at beginning of CFL ramp.
!! \param vnn0 Von Neumann number at beginning of CFL ramp.
!! \param brelax Boundary relaxation, 0.=fully 1st order BC, 1.=fully 2nd order BC
!! \param gradClip Gradient clipping, 1.=liberal clipping, 1.D14=no clipping
!! \param nRamp Number of steps over which to perform CFL ramp.
!! \param cfl CFL number for pseudo-time step definition.
!! \param vnn Von Neumann number for pseudo-time step definition.
!! \param convLimit Convergence level to exit pseudo-time iterations.
!! \param relax Multigrid relaxation parameter.
!! \param coarseDis Dissipation constant used on coarse multigrid levels.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-9-11
!! \par Further Documentation:
!! \par Source code:
!! \include src/Numerics/strandInputRead.F90


SUBROUTINE strandinputread(&
     fileLen, &
     fileName, &
     len, &
     systemType, &
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
     nLevels, &
     mgCycle, &
     gradient, &
     nodeVal, &
     perturb, &
     limiter, &
     dtUnsteady, &
     cflLinear, &
     cfl0, &
     vnn0, &
     nRamp, &
     brelax, &
     gradClip, &
     cfl, &
     vnn, &
     convLimit, &
     relax, &
     coarseDis)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: fileLen,len
CHARACTER(fileLen),INTENT(IN   ) :: fileName
CHARACTER(    len),INTENT(  OUT) :: systemType
INTEGER,INTENT(INOUT) :: &
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
     nLevels, &
     mgCycle, &
     gradient, &
     nodeVal, &
     perturb, &
     limiter, &
     nRamp
REAL,INTENT(INOUT) :: &
     dtUnsteady, &
     cflLinear, &
     brelax, &
     gradClip, &
     cfl0, &
     vnn0, &
     cfl, &
     vnn, &
     convLimit, &
     relax, &
     coarseDis
INTEGER,PARAMETER :: iUnit = 9
NAMELIST/numerics/&
     iprint,&
     itest,&
     idebug,&
     iConvFile,&
     iSolnFile,&
     iResdFile,&
     iErrFile,&
     iSurfFile,&
     standAlone,&
     restartStep,&
     nRestart,&
     nOutput,&
     nSteps,&
     nPseudoSteps,&
     nPseudoSteps0,&
     nLinearSteps,&
     nLevels,&
     mgCycle,&
     gradient,&
     nodeVal,&
     perturb,&
     limiter,&
     dtUnsteady,&
     cflLinear,&
     cflLinear, &
     cfl0, &
     vnn0, &
     nRamp, &
     brelax, &
     gradClip, &
     cfl,&
     vnn,&
     convLimit,&
     relax,&
     coarseDis,&
     systemType


OPEN(iUnit,FILE=TRIM(filename),STATUS='old')
READ(iUnit,numerics)
CLOSE(iUnit)

!null terminate the characters for passing back to c++
systemType = TRIM(systemType)//CHAR(0)

!fix nonsense values of ramping parameters
nRamp = MAX(1,nRamp)
cfl0  = MIN(cfl0,cfl)
vnn0  = MIN(vnn0,vnn)


END SUBROUTINE strandinputread
