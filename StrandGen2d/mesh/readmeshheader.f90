!> \brief
!! This subroutine reads the mesh file header.
!!
!! Comments:
!!  nSurfNodes    - number of surface nodes
!!  nSurfFaces    - number of surface faces
!!  nBndNodes     - number of boundary nodes. Boundary nodes are needed if the
!!                  geometry is open as for a plate (as opposed to closed as
!!                  for an airfoil). In the case of a flat plate, nBndNodes=2,
!!                  the first being the node at the entrance of the domain, the
!!                  second at the exit of the domain. For 2d strand grids this
!!                  number is either 0 or 2 (closed or open geometry).
!!  nPtsPerStrand - number of points along a strand. Note the total number of
!!                  points along a strand is actually nPtsPerStrand+1, since
!!                  the first strand node at the surface is indexed with 0.
!!  nSurfPatches  - number of surface patches for boundary condition assignment
!!  nNodePatches  - number node patches for assignment of boundary nodes
!!
!! Versions:
!!   - 1.0 Katz 08/03/2010
!!
!! Additional notes:\par
!!  none
!!
!! Source code:
!!   \include readmeshheader.f90


SUBROUTINE readmeshheader(meshfile, &
                          length, &
                          iunit, &
                          surfaceOnly0, &
                          nMesh0, &
                          nSurfNodes0, &
                          nSurfFaces0, &
                          nBndEdges0, &
                          nSurfPatches0, &
                          nEdgePatches0, &
                          nPtsPerStrand0)

USE avdefs

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: length, &
                         iunit
CHARACTER(length),INTENT(IN   ) :: meshfile
INTEGER,INTENT(  OUT) :: surfaceOnly0, &
                         nMesh0, &
                         nSurfNodes0, &
                         nSurfFaces0, &
                         nBndEdges0, &
                         nSurfPatches0, &
                         nEdgePatches0, &
                         nPtsPerStrand0
TYPE(file_id_prefix) :: fileIdPrefix
TYPE(file_header)    :: fileHeader
CHARACTER(128)       :: fileDescription
TYPE(mesh_header), &
     ALLOCATABLE,DIMENSION(:  ) :: meshHeader
CHARACTER(128), &
     ALLOCATABLE,DIMENSION(:  ) :: meshDescription
TYPE(strand2d_header), &
     ALLOCATABLE,DIMENSION(:  ) :: strandHeader


nMesh0 = 1 !for now

OPEN(iunit,FILE=meshfile,STATUS='old',FORM='unformatted')


! read file prefix, header, and description

READ(iunit)fileIdPrefix%magic_string, &
           fileIdPrefix%magic_number, &
           fileIdPrefix%format_version

READ(iunit)fileHeader%mesh_revision, &
           fileHeader%mesh_count, &      
           fileHeader%contact_name, &  
           fileHeader%contact_org, &
           fileHeader%contact_email, &
           fileHeader%precision, &                   
           fileHeader%dimensions, &                  
           fileHeader%reference_point, &             
           fileHeader%reference_point_description, & 
           fileHeader%coordinate_system, &           
           fileHeader%model_scale, &                 
           fileHeader%grid_units, &
           fileHeader%description_size

READ(iunit)fileDescription

ALLOCATE(meshHeader(nMesh0),      &
         meshDescription(nMesh0), &
         strandHeader(nMesh0)     )


! read generic mesh header and description

READ(iunit)meshHeader(1)%mesh_name,                &
        meshHeader(1)%mesh_type,                   &
        meshHeader(1)%mesh_generator,              &
        meshHeader(1)%changed_date,                &
        meshHeader(1)%coordinate_system,           &
        meshHeader(1)%model_scale,                 &
        meshHeader(1)%grid_units,                  &
        meshHeader(1)%reynolds_number,             &
        meshHeader(1)%reference_length,            &
        meshHeader(1)%wall_distance,               & ! 0. if surface only
        meshHeader(1)%reference_point,             &
        meshHeader(1)%reference_point_description, &
        meshHeader(1)%translation_vector,          &
        meshHeader(1)%rotation_matrix,             &
        meshHeader(1)%periodicity,                 &
        meshHeader(1)%periodicity_description,     &
        meshHeader(1)%refinement_level,            &
        meshHeader(1)%iblank,                      &
        meshHeader(1)%description_size

READ(iunit)meshDescription(1)


! read strand specific mesh header

READ(iunit)strandHeader(1)%surfaceOrVolume   ,&
        strandHeader(1)%surfaceOnly          ,&
        strandHeader(1)%nSurfNodes           ,&
        strandHeader(1)%nSurfFaces           ,&
        strandHeader(1)%nBndEdges            ,&
        strandHeader(1)%nPtsPerStrand        ,& !0  if surfaceOnly
        strandHeader(1)%nSurfPatches         ,&
        strandHeader(1)%nEdgePatches         ,&
        strandHeader(1)%strandLength         ,& !0. if surfaceOnly
        strandHeader(1)%stretchRatio         ,& !0. if surfaceOnly
        strandHeader(1)%smoothingThreshold   ,& !0. if surfaceOnly
        strandHeader(1)%strandDistribution      !0. if surfaceOnly

surfaceOnly0 = 0
IF (strandHeader(1)%surfaceOnly) surfaceOnly0 = 1
nSurfNodes0    = strandHeader(1)%nSurfNodes
nSurfFaces0    = strandHeader(1)%nSurfFaces
nBndEdges0     = strandHeader(1)%nBndEdges
nSurfPatches0  = strandHeader(1)%nSurfPatches
nEdgePatches0  = strandHeader(1)%nEdgePatches
nPtsPerStrand0 = strandHeader(1)%nPtsPerStrand

DEALLOCATE(meshHeader)
DEALLOCATE(meshDescription)
DEALLOCATE(strandHeader)


! report mesh statistics

WRITE(*,*)
WRITE(*,*)'*** strand mesh statistics ***'
WRITE(*,'(A39,I8)')'number of surface nodes: ',nSurfNodes0
WRITE(*,'(A39,I8)')'number of triangular surface faces: ',nSurfFaces0
WRITE(*,'(A39,I8)')'number of boundary edges: ',nBndEdges0
WRITE(*,'(A39,I8)')'number of surface patches: ',nSurfPatches0
WRITE(*,'(A39,I8)')'number of boundary patches: ',nEdgePatches0
WRITE(*,'(A39,I8)')'number of points per strand: ',nPtsPerStrand0
WRITE(*,*)


END SUBROUTINE readmeshheader
