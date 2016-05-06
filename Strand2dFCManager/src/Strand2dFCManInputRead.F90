!> \brief
!! This subroutine reads the input file.
!! \param fileLen Length of the character string containing fileName.
!! \param fileName Name of input file.
!! \param len Length of systemType name
!! \param meshFile name of mesh file.
!! \param nLevels Number of mesh levels to create.
!! \param meshOrdersT Desired order of the mesh on each level.
!! \param surfDist Method of spacing surface element DOFS.
!! \param strandDist Method of spacing strand nodes.
!! \param perturb Percentage by which strand distribution is perturbed.
!! \param nStrandNode Number of nodes on strands.
!! \param wallSpacing Spacing of first cell off the wall.
!! \param strandLength Length of strands.
!! \param stretchRatio Strand stretching ratio between adjacent cells.
!! \param deltaSmooth Strand smoothing threshold.
!! \param iplotmesh Flag to write Paraview mesh file.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   03-15-2013
!! \par Further Documentation:
!! \par Source code:
!! \include src/Strand2dFCManInputRead.F90


SUBROUTINE strand2dfcmaninputread(&
     fileLen, &
     fileName, &
     len, &
     meshFile, &
     nLevels, &
     meshOrdersT, &
     surfDist, &
     strandDist, &
     perturb, &
     nStrandNode, &
     wallSpacing, &
     strandLength, &
     stretchRatio, &
     deltaSmooth, &
     iplotmesh)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: fileLen,len
CHARACTER(fileLen),INTENT(IN   ) :: fileName
CHARACTER(    len),INTENT(  OUT) :: meshFile
INTEGER,INTENT(  OUT) :: nLevels
INTEGER,INTENT(  OUT),DIMENSION(len) :: meshOrdersT
INTEGER,INTENT(  OUT) :: surfDist
INTEGER,INTENT(  OUT) :: strandDist
REAL,   INTENT(  OUT) :: perturb
INTEGER,INTENT(  OUT) :: nStrandNode
REAL,   INTENT(  OUT) :: wallSpacing
REAL,   INTENT(  OUT) :: strandLength
REAL,   INTENT(  OUT) :: stretchRatio
REAL,   INTENT(  OUT) :: deltaSmooth
INTEGER,INTENT(  OUT) :: iplotmesh
INTEGER,PARAMETER :: iUnit = 9
INTEGER,DIMENSION(10) :: meshOrders
NAMELIST/Strand2dFCManager/&
     meshFile, &
     nLevels, &
     meshOrders, &
     surfDist, &
     strandDist, &
     perturb, &
     nStrandNode, &
     wallSpacing, &
     strandLength, &
     stretchRatio, &
     deltaSmooth, &
     iplotmesh


OPEN(iUnit,FILE=TRIM(fileName),STATUS='old')
READ(iUnit,Strand2dFCManager)
CLOSE(iUnit)

WRITE(*,*)
WRITE(*,Strand2dFCManager)
WRITE(*,*)

meshOrdersT(1:nLevels) = meshOrders(1:nLevels)


!null terminate the characters for passing back to c++
meshFile = TRIM(meshFile)//CHAR(0)


END SUBROUTINE strand2dfcmaninputread
