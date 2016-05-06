!> \brief
!! This subroutine reads general inputs for the meshing layer.
!!
!! Comments:
!!
!! Versions:
!!   - 1.0 Katz 08/03/2010
!!
!! Additional notes:\par
!!  none
!!
!! Source code:
!!   \include readinput.f90


SUBROUTINE readinput(inputfile,length,iwrite,strandDist,nPtsPerStrand, &
                     nfringe, &
                     stretchRatio,wallSpacing,strandLength,deltaSmooth, &
                     trimType)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: length
CHARACTER(length),INTENT(IN   ) :: inputfile
INTEGER,INTENT(  OUT) :: iwrite         !< verbosity flag
INTEGER,INTENT(  OUT) :: strandDist     !< type of strand point distribution
INTEGER,INTENT(  OUT) :: nPtsPerStrand  !< number of strand points
INTEGER,INTENT(  OUT) :: nfringe        !< number of fringe cells (fine level)
REAL   ,INTENT(  OUT) :: stretchRatio   !< max. allowable stretch ratio
REAL   ,INTENT(  OUT) :: wallspacing    !< spacing of first strand point
REAL   ,INTENT(  OUT) :: strandLength   !< length of strand
REAL   ,INTENT(  OUT) :: deltaSmooth    !< smoothing threshold
INTEGER,INTENT(  OUT) :: trimType       !< type of trimming to perform
INTEGER :: iinputunit
NAMELIST/strand2d/iwrite,strandDist,nPtsPerStrand,nfringe,stretchRatio, &
                  wallSpacing,strandLength,deltaSmooth,trimType


iinputunit = 7
OPEN(iinputunit,FILE=inputfile,STATUS='old',CONVERT='big_endian')

READ(iinputunit,strand2d)
WRITE(*,*)
WRITE(*,strand2d)
WRITE(*,*)

CLOSE(iinputunit)



END SUBROUTINE readinput
