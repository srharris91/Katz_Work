!> \brief
!! This subroutine creates the header file for the multiblock plot.
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
!!   \include plotheader.f90


SUBROUTINE plotheader(nStrandBlocks)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: nStrandBlocks
CHARACTER(48) :: filename
INTEGER :: n,iu


iu = 11
OPEN(iu,FILE='StrandParts.pvtu',STATUS='replace')
WRITE(iu,'(A)')'<?xml version="1.0"?>'
WRITE(iu,'(A)')'<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
WRITE(iu,'(A)')'<PUnstructuredGrid GhostLevel="0">'
WRITE(iu,'(A)')'<PCellData Scalars="ID">'
WRITE(iu,'(A)')'<PDataArray type="Int32" Name="ID"/>'
WRITE(iu,'(A)')'</PCellData>'
WRITE(iu,'(A)')'<PPoints>'
WRITE(iu,'(A)')'<PDataArray type="Float32" NumberOfComponents="3"/>'
WRITE(iu,'(A)')'</PPoints>'
DO n=0,nStrandBlocks-1
   IF (                  n < 10     ) WRITE(filename,'(I1)')n
   IF (n >= 10     .AND. n < 100    ) WRITE(filename,'(I2)')n
   IF (n >= 100    .AND. n < 1000   ) WRITE(filename,'(I3)')n
   IF (n >= 1000   .AND. n < 10000  ) WRITE(filename,'(I4)')n
   IF (n >= 10000  .AND. n < 100000 ) WRITE(filename,'(I5)')n
   IF (n >= 100000 .AND. n < 1000000) WRITE(filename,'(I6)')n
   filename = '<Piece Source="StrandPart'//TRIM(filename)//'.vtu"/>'
   WRITE(iu,'(A)')TRIM(filename)
END DO
WRITE(iu,'(A)')'</PUnstructuredGrid>'
WRITE(iu,'(A)')'</VTKFile>'
CLOSE(iu)


END SUBROUTINE plotheader
