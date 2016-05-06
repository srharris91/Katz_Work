!> \brief
!! This subroutine reads the input file.
!! \param fileLen Length of the character string containing fileName.
!! \param fileName Name of input file.
!! \param iplotmesh Flag to write Paraview mesh file.
!! \param len Length of systemType name
!! \param meshFile name of mesh file.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   03-15-2013
!! \par Further Documentation:
!! \par Source code:
!! \include src/tri2dFCManInputRead.F90


SUBROUTINE tri2dfcmaninputread(&
     fileLen, &
     fileName, &
     iplotmesh, &
     len, &
     meshFile)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: fileLen,len
INTEGER,INTENT(OUT  ) :: iplotmesh
CHARACTER(fileLen),INTENT(IN   ) :: fileName
CHARACTER(    len),INTENT(  OUT) :: meshFile
INTEGER,PARAMETER :: iUnit = 9
NAMELIST/Tri2dFCManager/&
     meshFile, &
     iplotmesh


OPEN(iUnit,FILE=TRIM(fileName),STATUS='old')
READ(iUnit,Tri2dFCManager)
CLOSE(iUnit)


!null terminate the characters for passing back to c++
meshFile = TRIM(meshFile)//CHAR(0)


END SUBROUTINE tri2dfcmaninputread
