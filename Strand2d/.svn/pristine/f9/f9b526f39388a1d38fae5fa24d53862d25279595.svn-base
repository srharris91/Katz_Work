!> \brief
!! This subroutine creates the header file for the multiblock strand plot.
!! \param step Physical time step.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-11-27
!! \par Further Documentation:
!! \par Source code:
!! \include src/System/StrandSPLam/StrandSPLamStepHeader.F90


SUBROUTINE strandsplamstepheader( &
     step)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: step
CHARACTER(48) :: cblock,cstep
INTEGER :: n,iu,nStrandBlocks


nStrandBlocks = 1


! create output directory for this step
n = step
IF (                  n < 10     ) WRITE(cstep,'(I1)')n
IF (n >= 10     .AND. n < 100    ) WRITE(cstep,'(I2)')n
IF (n >= 100    .AND. n < 1000   ) WRITE(cstep,'(I3)')n
IF (n >= 1000   .AND. n < 10000  ) WRITE(cstep,'(I4)')n
IF (n >= 10000  .AND. n < 100000 ) WRITE(cstep,'(I5)')n
IF (n >= 100000 .AND. n < 1000000) WRITE(cstep,'(I6)')n
CALL SYSTEM('mkdir -p output.'//TRIM(cstep))
CALL SYSTEM('rm -f output.'//TRIM(cstep)//'/*part*')


! need to figure out convergence file output format


iu = 11
OPEN(iu,FILE='solution'//TRIM(cstep)//'.pvtu',STATUS='replace')
WRITE(iu,'(A)')'<?xml version="1.0"?>'
WRITE(iu,'(A)')'<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
WRITE(iu,'(A)')'<PUnstructuredGrid GhostLevel="0">'
WRITE(iu,'(A)')'<PPointData Scalars="P T rho entropy" Vectors="U">'
WRITE(iu,'(A)')'<PDataArray type="Float32" Name="P"/>'
WRITE(iu,'(A)')'<PDataArray type="Float32" Name="T"/>'
WRITE(iu,'(A)')'<PDataArray type="Float32" Name="rho"/>'
WRITE(iu,'(A)')'<PDataArray type="Float32" Name="entropy"/>'
WRITE(iu,'(A)')'<PDataArray type="Float32" Name="U" NumberOfComponents="3"/>'
WRITE(iu,'(A)')'</PPointData>'
WRITE(iu,'(A)')'<PPoints>'
WRITE(iu,'(A)')'<PDataArray type="Float32" NumberOfComponents="3"/>'
WRITE(iu,'(A)')'</PPoints>'
DO n=0,nStrandBlocks-1
   IF (                  n < 10     ) WRITE(cblock,'(I1)')n
   IF (n >= 10     .AND. n < 100    ) WRITE(cblock,'(I2)')n
   IF (n >= 100    .AND. n < 1000   ) WRITE(cblock,'(I3)')n
   IF (n >= 1000   .AND. n < 10000  ) WRITE(cblock,'(I4)')n
   IF (n >= 10000  .AND. n < 100000 ) WRITE(cblock,'(I5)')n
   IF (n >= 100000 .AND. n < 1000000) WRITE(cblock,'(I6)')n
   cblock = '<Piece Source="output.'//TRIM(cstep)//'/'//'soln_part'//TRIM(cblock)//'.vtu"/>'
!   cblock = '<Piece Source="soln_part'//TRIM(cblock)//'.vtu"/>'
   WRITE(iu,'(A)')TRIM(cblock)
END DO
WRITE(iu,'(A)')'</PUnstructuredGrid>'
WRITE(iu,'(A)')'</VTKFile>'
CLOSE(iu)
!CALL SYSTEM('mv solution'//TRIM(cstep)//'.pvtu output.'//TRIM(cstep))


iu = 11
OPEN(iu,FILE='residual'//TRIM(cstep)//'.pvtu',STATUS='replace')
WRITE(iu,'(A)')'<?xml version="1.0"?>'
WRITE(iu,'(A)')'<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
WRITE(iu,'(A)')'<PUnstructuredGrid GhostLevel="0">'
WRITE(iu,'(A)')'<PCellData Scalars="r1 r2 r3 r4">'
WRITE(iu,'(A)')'<PDataArray type="Float32" Name="r1"/>'
WRITE(iu,'(A)')'<PDataArray type="Float32" Name="r2"/>'
WRITE(iu,'(A)')'<PDataArray type="Float32" Name="r3"/>'
WRITE(iu,'(A)')'<PDataArray type="Float32" Name="r4"/>'
WRITE(iu,'(A)')'</PCellData>'
WRITE(iu,'(A)')'<PPoints>'
WRITE(iu,'(A)')'<PDataArray type="Float32" NumberOfComponents="3"/>'
WRITE(iu,'(A)')'</PPoints>'
DO n=0,nStrandBlocks-1
   IF (                  n < 10     ) WRITE(cblock,'(I1)')n
   IF (n >= 10     .AND. n < 100    ) WRITE(cblock,'(I2)')n
   IF (n >= 100    .AND. n < 1000   ) WRITE(cblock,'(I3)')n
   IF (n >= 1000   .AND. n < 10000  ) WRITE(cblock,'(I4)')n
   IF (n >= 10000  .AND. n < 100000 ) WRITE(cblock,'(I5)')n
   IF (n >= 100000 .AND. n < 1000000) WRITE(cblock,'(I6)')n
   cblock = '<Piece Source="output.'//TRIM(cstep)//'/'//'resd_part'//TRIM(cblock)//'.vtu"/>'
!   cblock = '<Piece Source="resd_part'//TRIM(cblock)//'.vtu"/>'
   WRITE(iu,'(A)')TRIM(cblock)
END DO
WRITE(iu,'(A)')'</PUnstructuredGrid>'
WRITE(iu,'(A)')'</VTKFile>'
CLOSE(iu)
!CALL SYSTEM('mv residual'//TRIM(cstep)//'.pvtu output.'//TRIM(cstep))


iu = 11
OPEN(iu,FILE='error'//TRIM(cstep)//'.pvtu',STATUS='replace')
WRITE(iu,'(A)')'<?xml version="1.0"?>'
WRITE(iu,'(A)')'<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
WRITE(iu,'(A)')'<PUnstructuredGrid GhostLevel="0">'
WRITE(iu,'(A)')'<PCellData Scalars="e1 e2 e3 e4">'
WRITE(iu,'(A)')'<PDataArray type="Float32" Name="e1"/>'
WRITE(iu,'(A)')'<PDataArray type="Float32" Name="e2"/>'
WRITE(iu,'(A)')'<PDataArray type="Float32" Name="e3"/>'
WRITE(iu,'(A)')'<PDataArray type="Float32" Name="e4"/>'
WRITE(iu,'(A)')'</PCellData>'
WRITE(iu,'(A)')'<PPoints>'
WRITE(iu,'(A)')'<PDataArray type="Float32" NumberOfComponents="3"/>'
WRITE(iu,'(A)')'</PPoints>'
DO n=0,nStrandBlocks-1
   IF (                  n < 10     ) WRITE(cblock,'(I1)')n
   IF (n >= 10     .AND. n < 100    ) WRITE(cblock,'(I2)')n
   IF (n >= 100    .AND. n < 1000   ) WRITE(cblock,'(I3)')n
   IF (n >= 1000   .AND. n < 10000  ) WRITE(cblock,'(I4)')n
   IF (n >= 10000  .AND. n < 100000 ) WRITE(cblock,'(I5)')n
   IF (n >= 100000 .AND. n < 1000000) WRITE(cblock,'(I6)')n
   cblock = '<Piece Source="output.'//TRIM(cstep)//'/'//'err_part'//TRIM(cblock)//'.vtu"/>'
!   cblock = '<Piece Source="err_part'//TRIM(cblock)//'.vtu"/>'
   WRITE(iu,'(A)')TRIM(cblock)
END DO
WRITE(iu,'(A)')'</PUnstructuredGrid>'
WRITE(iu,'(A)')'</VTKFile>'
CLOSE(iu)
!CALL SYSTEM('mv error'//TRIM(cstep)//'.pvtu output.'//TRIM(cstep))


END SUBROUTINE strandsplamstepheader
