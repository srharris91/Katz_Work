PROGRAM airfoil

USE avdefs

IMPLICIT NONE

INTEGER :: nSurfNodes,nSurfFaces,nBndEdges,nPtsPerStrand, &
           nSurfPatches,nEdgePatches,nMesh
TYPE(file_id_prefix) :: fileIdPrefix
TYPE(file_header)    :: fileHeader
CHARACTER(128)       :: fileDescription
TYPE(mesh_header), &
     ALLOCATABLE,DIMENSION(:  ) :: meshHeader
CHARACTER(128), &
     ALLOCATABLE,DIMENSION(:  ) :: meshDescription
TYPE(strand2d_header), &
     ALLOCATABLE,DIMENSION(:  ) :: strandHeader
TYPE(strand2d_surf_patch), &
     ALLOCATABLE,DIMENSION(:,:) :: strandSurfPatch
TYPE(strand2d_edge_patch), &
     ALLOCATABLE,DIMENSION(:,:) :: strandEdgePatch
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: conn
INTEGER,ALLOCATABLE,DIMENSION(  :) :: bndEdges
INTEGER,ALLOCATABLE,DIMENSION(  :) :: nodeClip
INTEGER,ALLOCATABLE,DIMENSION(  :) :: faceClip
INTEGER,ALLOCATABLE,DIMENSION(  :) :: faceTag
INTEGER,ALLOCATABLE,DIMENSION(  :) :: edgeTag
REAL*8 ,ALLOCATABLE,DIMENSION(:,:) :: xSurf
REAL*8 ,ALLOCATABLE,DIMENSION(:,:) :: pointingVec
REAL*8 ,ALLOCATABLE,DIMENSION(  :) :: xStrand
INTEGER :: n,m,k,narg,i,j,iu,n1,n2,npc,ip
REAL :: y0,dy
REAL,ALLOCATABLE,DIMENSION(:,:) :: x
CHARACTER(80) :: fileName,arg,ptsFile


! read inputs

narg = iargc()
m    = 0
DO WHILE (m < narg)
   m    = m+1
   CALL GETARG(m,arg)
   IF (arg(1:4) == '-pts') THEN
      m = m+1
      CALL GETARG(m,arg)
      ptsFile = arg
   ELSE IF (arg(1:2) == '-m') THEN
      m = m+1
      CALL GETARG(m,arg)
      fileName = arg
   END IF
END DO


! read .pts file

OPEN(UNIT=iu,FILE=ptsFile,FORM='formatted')
READ(iu,*)
READ(iu,*)
READ(iu,*)
READ(iu,*)
READ(iu,*)npc
npc = npc-1
ALLOCATE(x(2,npc))
READ(iu,*)
READ(iu,*)
READ(iu,*)
READ(iu,*)
READ(iu,*)
READ(iu,*)
READ(iu,*)
DO n=1,npc
   READ(iu,*)x(1,n),x(2,n)
END DO
CLOSE(iu)


! compute dimensions of mesh and allocate data

nMesh         = 1
nSurfNodes    = npc
nSurfFaces    = npc
nBndEdges     = 0
nPtsPerStrand = 0
nSurfPatches  = 1
nEdgePatches  = 0

ALLOCATE(meshHeader(nMesh),      &
         meshDescription(nMesh), &
         strandHeader(nMesh)     )
ALLOCATE(strandSurfPatch(nSurfPatches,nMesh), &
         strandEdgePatch(nEdgePatches,nMesh)  )
ALLOCATE(conn       (2,nSurfFaces   ), &
         bndEdges   (  nBndEdges    ), &
         nodeClip   (  nSurfNodes   ), &
         faceClip   (  nSurfFaces   ), &
         faceTag    (  nSurfFaces   ), &
         edgeTag    (  nBndEdges    ), &
         xSurf      (2,nSurfNodes   ), &
         pointingVec(2,nSurfNodes   ), &
         xStrand    (0:nPtsPerStrand)  )

! generate file prefix, header and description

fileIdPrefix%magic_string   = 'AVMESH'
fileIdPrefix%magic_number   = 1
fileIdPrefix%format_version = '1.0'

fileHeader%mesh_revision               = '1.0'
fileHeader%mesh_count                  = nMesh
fileHeader%contact_name                = 'Aaron Katz'
fileHeader%contact_org                 = 'Utah State University'
fileHeader%contact_email               = 'aaron.katz@usu.edu'
fileHeader%precision                   = 2
fileHeader%dimensions                  = 2
fileHeader%reference_point(1)          = 0.0
fileHeader%reference_point(2)          = 0.0
fileHeader%reference_point(3)          = 0.0
fileHeader%reference_point_description = 'origin'
fileHeader%coordinate_system           = 'xByUzL'
fileHeader%model_scale                 = 1.0
fileHeader%grid_units                  = 'm'

fileDescription = 'airfoil'
fileHeader%description_size = 128

WRITE(meshHeader(1)%mesh_name, '(A,I1)') 'strand mesh # ', n
meshHeader(1)%mesh_type                   = 'strand'
meshHeader(1)%mesh_generator              = 'airfoil.f90'
meshHeader(1)%changed_date                = 'September 3, 2011'
meshHeader(1)%coordinate_system           = 'xByUzL'
meshHeader(1)%model_scale                 = 1.0
meshHeader(1)%grid_units                  = 'm'
meshHeader(1)%reynolds_number             = 1.0
meshHeader(1)%reference_length            = 1.0
meshHeader(1)%wall_distance               = 1.0
meshHeader(1)%reference_point(1)          = 0.0
meshHeader(1)%reference_point(2)          = 0.0
meshHeader(1)%reference_point(3)          = 0.0
meshHeader(1)%reference_point_description = 'origin'
meshHeader(1)%translation_vector(1)       = 0.0
meshHeader(1)%translation_vector(2)       = 0.0
meshHeader(1)%translation_vector(3)       = 0.0
meshHeader(1)%rotation_matrix(1)          = 1.0
meshHeader(1)%rotation_matrix(2)          = 0.0
meshHeader(1)%rotation_matrix(3)          = 0.0
meshHeader(1)%rotation_matrix(4)          = 0.0
meshHeader(1)%rotation_matrix(5)          = 1.0
meshHeader(1)%rotation_matrix(6)          = 0.0
meshHeader(1)%rotation_matrix(7)          = 0.0
meshHeader(1)%rotation_matrix(8)          = 0.0
meshHeader(1)%rotation_matrix(9)          = 1.0
meshHeader(1)%periodicity                 = 1.0
meshHeader(1)%periodicity_description     = '2d'
meshHeader(1)%refinement_level            = 1
meshHeader(1)%iblank                      = 0

meshDescription(1) = 'airfoil in 2d'
meshHeader(1)%description_size = 128

strandHeader(1)%surfaceOrVolume    = 'surface'
strandHeader(1)%surfaceOnly        = .true.
strandHeader(1)%nSurfNodes         = nSurfNodes
strandHeader(1)%nSurfFaces         = nSurfFaces
strandHeader(1)%nBndEdges          = nBndEdges
strandHeader(1)%nPtsPerStrand      = nPtsPerStrand
strandHeader(1)%nSurfPatches       = nSurfPatches
strandHeader(1)%nEdgePatches       = nEdgePatches
strandHeader(1)%strandLength       = 1.
strandHeader(1)%stretchRatio       = 1.
strandHeader(1)%smoothingThreshold = 1.
strandHeader(1)%strandDistribution = 'geometric'

DO n=1,nSurfPatches
   strandSurfPatch(n,1)%surfPatchID     = n
   strandSurfPatch(n,1)%surfPatchBody   = 'airfoil'
   strandSurfPatch(n,1)%surfPatchComp   = 'surface'
   strandSurfPatch(n,1)%surfPatchBCType = 'dirichlet'
END DO
do n=1,nEdgePatches
   !there are no edge patches for this grid
end do


! generate grid data

DO n=1,nSurfFaces
   conn(1,n) = n
   conn(2,n) = n+1
   faceTag(n) = 1
END DO
conn(2,nSurfFaces) = 1
nodeClip = 0
faceClip = 0
pointingVec = 0.
xSurf = x



! write the mesh to file

iu = 7
OPEN(UNIT=iu,FILE=filename,FORM='unformatted')

WRITE(iu)fileIdPrefix%magic_string, &
         fileIdPrefix%magic_number, &
         fileIdPrefix%format_version

WRITE(iu)fileHeader%mesh_revision, &
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

WRITE(iu)fileDescription

CALL write_mesh_header(iu,meshHeader(1),meshDescription(1))
CALL write_strand2d_header(iu,strandHeader(1),strandSurfPatch(:,1), &
                           strandEdgePatch(:,1))


WRITE(iu)((conn(k,n),k=1,2),n=1,nSurfFaces)
WRITE(iu)(bndEdges(n),n=1,nBndEdges)
WRITE(iu)(nodeClip(n),n=1,nSurfNodes)
WRITE(iu)(faceClip(n),n=1,nSurfFaces)
WRITE(iu)(faceTag(n),n=1,nSurfFaces)
WRITE(iu)(edgeTag(n),n=1,nBndEdges)
WRITE(iu)((xSurf(k,n),k=1,2),n=1,nSurfNodes)
WRITE(iu)((pointingVec(k,n),k=1,2),n=1,nSurfNodes)
WRITE(iu)(xStrand(n),n=0,nPtsPerStrand)

CLOSE(iu)

DEALLOCATE(meshHeader,      &
           meshDescription, &
           strandHeader)
DEALLOCATE(strandSurfPatch, &
           strandEdgePatch)
DEALLOCATE(conn       , &
           bndEdges   , &
           nodeClip   , &
           faceClip   , &
           faceTag    , &
           edgeTag    , &
           xSurf      , &
           pointingVec, &
           xStrand      )
DEALLOCATE(x)

WRITE(*,'(A,A)')  'Successfully wrote: ', filename


END PROGRAM airfoil
