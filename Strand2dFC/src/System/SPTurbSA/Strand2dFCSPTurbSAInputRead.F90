!> \brief
!!   This subroutine reads the System input file.
!! \param fileLen Length of the character string containing fileName.
!! \param fileName Name of input file.
!! \param nq number of equations.
!! \param nComp Number of fluid components.
!! \param nBpatches Number of boundary patches.
!! \param nBpatchesT Overdimensioned nBpatches.
!! \param inviscid Flag to include inviscid terms.
!! \param dissipation Flag to include dissipation terms.
!! \param viscous Flag to include viscous terms.
!! \param sourceMMS Flag to include MMS source terms.
!! \param isolution Flag indicating which exact solution is used.
!! \param istate Integer key to the state type.
!! \param itransport Integer key to the transport type.
!! \param bType Integer key to the boundary type.
!! \param bValue Reference boundary values for each boundary patch.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2013-1-1
!! \par Further Documentation:
!! \par Source code:
!! \include src/System/SPTurbSA/Strand2dFCSPTurbSAInputRead.F90


SUBROUTINE strand2dfrspturbsainputread(&
     fileLen, &
     fileName, &
     nq, &
     nComp, &
     nBpatches, &
     nBpatchesT, &
     inviscid, &
     dissipation, &
     viscous, &
     sourceMMS, &
     solution0, &
     state0, &
     transport0, &
     bType0, &
     bValue0 &
)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: fileLen
CHARACTER(fileLen),INTENT(IN   ) :: fileName
INTEGER,INTENT(IN   ) :: nq
INTEGER,INTENT(IN   ) :: nComp
INTEGER,INTENT(  OUT) :: nBpatches
INTEGER,INTENT(IN   ) :: nBpatchesT
INTEGER,INTENT(  OUT) :: inviscid
INTEGER,INTENT(  OUT) :: dissipation
INTEGER,INTENT(  OUT) :: viscous
INTEGER,INTENT(  OUT) :: sourceMMS
INTEGER,INTENT(  OUT) :: solution0
INTEGER,INTENT(  OUT),DIMENSION(   nComp     ) :: state0
INTEGER,INTENT(  OUT),DIMENSION(   nComp     ) :: transport0
INTEGER,INTENT(  OUT),DIMENSION(   nBpatchesT) :: bType0
REAL   ,INTENT(  OUT),DIMENSION(nq,nBpatchesT) :: bValue0
INTEGER :: n
INTEGER,PARAMETER :: iUnit   = 9
INTEGER,PARAMETER :: ibcUnit = 7
CHARACTER(80) :: solution
CHARACTER(80),DIMENSION(  1   ) :: state
CHARACTER(80),DIMENSION(  1   ) :: transport
CHARACTER(80),DIMENSION(  1000) :: bType
REAL,         DIMENSION(5,1000) :: bValue
NAMELIST/SPTurbSA/&
     inviscid, &
     dissipation, &
     viscous, &
     solution, &
     state, &
     transport, &
     nBpatches, &
     bType, &
     bValue


OPEN(iUnit,FILE=fileName,STATUS='old')
READ(iUnit,SPTurbSA)
CLOSE(iUnit)


! ordered such that at corners the higher number takes precedence
DO n=1,nBPatches
   IF      (TRIM(bType(n)) == 'inviscidWall') THEN
      bType0(n) = 1
   ELSE IF (TRIM(bType(n)) == 'viscousWall' ) THEN
      bType0(n) = 2
   ELSE IF (TRIM(bType(n)) == 'farField'    ) THEN
      bType0(n) = 3
   ELSE IF (TRIM(bType(n)) == 'outflow'     ) THEN
      bType0(n) = 4
   ELSE
      WRITE(*,*)'*** boundary type not recognized in SPTurbSAinputread.F90 ***'
      STOP
   END IF
   bValue0(1:5,n) = bValue(1:5,n)
END DO

DO n=1,nComp
   IF      (TRIM(state(n)) == 'idealGas') THEN
      state0(n) = 1
   ELSE
      WRITE(*,*)'*** state type not recognized in SPTurbSAinputread.F90 ***'
      STOP
   END IF
END DO

DO n=1,nComp
   IF      (TRIM(transport(n)) == 'sutherland') THEN
      transport0(n) = 1
   ELSE
      WRITE(*,*)'*** transport type not recognized in SPTurbSAinputread.F90 ***'
      STOP
   END IF
END DO

IF      (TRIM(solution) == 'freeStream') THEN
   solution0 = 1
   sourceMMS = 0
ELSE IF (TRIM(solution) == 'perturb'   ) THEN
   solution0 = 2
   sourceMMS = 0
ELSE IF (TRIM(solution) == 'mms2d'     ) THEN
   solution0 = 3
   sourceMMS = 1
ELSE IF (TRIM(solution) == 'ringleb'   ) THEN
   solution0 = 4
   sourceMMS = 0
ELSE
   WRITE(*,*)'*** solution type not recognized in SPTurbSAinputread.F90 ***'
   STOP
END IF


END SUBROUTINE strand2dfrspturbsainputread
