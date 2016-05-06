!> \brief
!!   This subroutine reads the mms2d solution input file.
!! \param fileLen Length of character string of input file.
!! \param fileName Name of input file.
!! \param nq Dimension of q (number of equations).
!! \param period Period of MMS solution.
!! \param amplitude Amplitude of variation of MMS solution.
!! \param rValue Solution initialization values.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-9-13
!! \par Further Documentation:
!! \par Source code:
!! \include src/Fortran/inputSOLUTIONReadMms2d.F90


SUBROUTINE inputsolutionreadmms2d(&
     fileLen, &
     fileName, &
     nq, &
     period, &
     amplitude, &
     rValue0)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: fileLen
CHARACTER(fileLen),INTENT(IN   ) :: fileName
INTEGER,INTENT(IN   ) :: nq
REAL   ,INTENT(  OUT) :: period
REAL   ,INTENT(  OUT) :: amplitude
REAL   ,INTENT(  OUT),DIMENSION(nq) :: rValue0
REAL   ,DIMENSION(100) :: rValue
INTEGER,PARAMETER :: iUnit = 11
NAMELIST/mms2d/&
     rValue, &
     period, &
     amplitude


OPEN(iUnit,FILE=fileName,STATUS='old')
READ(iUnit,mms2d)
CLOSE(iUnit)

rValue0(1:nq) = rValue(1:nq)


END SUBROUTINE inputsolutionreadmms2d
