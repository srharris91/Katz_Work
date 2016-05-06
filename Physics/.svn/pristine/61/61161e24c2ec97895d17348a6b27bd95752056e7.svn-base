!> \brief
!!   This subroutine reads the freestream solution input file.
!! \param fileLen Length of character string of input file.
!! \param fileName Name of input file.
!! \param nq Dimension of q (number of equations).
!! \param rValue Solution initialization values.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-9-13
!! \par Further Documentation:
!! \par Source code:
!! \include src/Fortran/inputSOLUTIONReadFreeStream.F90


SUBROUTINE inputsolutionreadfreestream(&
     fileLen, &
     fileName, &
     nq, &
     rValue0)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: fileLen
CHARACTER(fileLen),INTENT(IN   ) :: fileName
INTEGER,INTENT(IN   ) :: nq
REAL   ,INTENT(  OUT),DIMENSION(nq) :: rValue0
REAL   ,DIMENSION(100) :: rValue
INTEGER,PARAMETER :: iUnit = 11
NAMELIST/freeStream/&
     rValue


OPEN(iUnit,FILE=fileName,STATUS='old')
READ(iUnit,freeStream)
CLOSE(iUnit)

rValue0(1:nq) = rValue(1:nq)


END SUBROUTINE inputsolutionreadfreestream
