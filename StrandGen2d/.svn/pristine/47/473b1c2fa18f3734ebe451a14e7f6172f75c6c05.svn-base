module avdefs

   type file_id_prefix
      character(6) :: magic_string     
      integer(4)   :: magic_number
      character(8) :: format_version   
   end type file_id_prefix

   type file_header
      character(128) :: mesh_revision
      integer(4)     :: mesh_count
      character(128) :: contact_name
      character(128) :: contact_org
      character(128) :: contact_email
      integer(4)     :: precision
      integer(4)     :: dimensions
      real(8)        :: reference_point(3)
      character(128) :: reference_point_description
      character(128) :: coordinate_system
      real(8)        :: model_scale
      character(128) :: grid_units
      integer(4)     :: description_size
   end type file_header

   !this follows file_header:
   ! character( file_header % description_size ) :: file_description

   type mesh_header
      character(128) :: mesh_name
      character(128) :: mesh_type
      character(128) :: mesh_generator
      character(128) :: changed_date
      character(128) :: coordinate_system
      real(8)        :: model_scale
      character(128) :: grid_units
      real(8)        :: reynolds_number
      real(8)        :: reference_length
      real(8)        :: wall_distance
      real(8)        :: reference_point(3)
      character(128) :: reference_point_description
      real(8)        :: translation_vector(3)
      real(8)        :: rotation_matrix(9)
      real(8)        :: periodicity
      character(128) :: periodicity_description
      integer(4)     :: refinement_level
      integer(4)     :: iblank
      integer(4)     :: description_size
   end type mesh_header

   !this follows each mesh_header:
   ! character( mesh_header % description_size ) :: mesh_description

   type bfstruc_header
      integer(4):: jmax                         
      integer(4):: kmax                   
      integer(4):: lmax                   
      integer(4):: n_patches
      integer(4):: n_patch_int_params
      integer(4):: n_patch_r8_params
      integer(4):: n_patch_c80_params
   end type bfstruc_header

   !n_patches copies of the following type follows each bfstruc_header
   type bfstruc_patch
      integer(4):: bctype
      integer(4):: direction
      integer(4):: jbegin
      integer(4):: jend
      integer(4):: kbegin
      integer(4):: kend
      integer(4):: lbegin
      integer(4):: lend
      character(128) :: description
   end type bfstruc_patch

   !these follow bfstruc_patch's:
   ! integer(4)    :: patch_int_params (n_patches, n_patch_int_params)
   ! real(8)       :: patch_r8_params  (n_patches, n_patch_r8_params)
   ! character(80) :: patch_c80_params (n_patches, n_patch_c80_params)

   type unstruc_header
      integer(4) :: full_viscous_layer_count
      integer(4) :: partial_viscous_layer_count
      integer(4) :: nNodes
      integer(4) :: nEdges
      integer(4) :: nFaces
      integer(4) :: nMaxNodesPerFace
      integer(4) :: nSumNodesPerFace
      integer(4) :: nTriFaces
      integer(4) :: nQuadFaces
      integer(4) :: nPolyFaces
      integer(4) :: nCells
      integer(4) :: nMaxNodesPerCell
      integer(4) :: nSumNodesPerCell
      integer(4) :: nMaxFacesPerCell
      integer(4) :: nSumFacesPerCell
      integer(4) :: nHexCells
      integer(4) :: nTetCells
      integer(4) :: nPriCells
      integer(4) :: nPyrCells
      integer(4) :: nPolyCells
      integer(4) :: nPatches
   end type unstruc_header

   !nPatches copies of the following type follows each unstruc_header
   type unstruc_patchheader
      character(32) :: patch_label
      character(16) :: patch_type
      integer(4) :: patch_id
      integer(4) :: nNodes
      integer(4) :: nFaces   
   end type unstruc_patchheader

   type strand_header
      character(32) :: surfaceOrVolume
      logical       :: surfaceOnly
      integer(4)    :: nSurfNodes
      integer(4)    :: nTriFaces
      integer(4)    :: nQuadFaces
      integer(4)    :: nBndEdges
      integer(4)    :: nPtsPerStrand
      integer(4)    :: nSurfPatches
      integer(4)    :: nEdgePatches
      real(8)       :: strandLength
      real(8)       :: stretchRatio
      real(8)       :: smoothingThreshold
      character(32) :: strandDistribution
   end type strand_header

   !nSurfPatches copies of the following type follows each strand_header
   type strand_surf_patch
      integer(4)    :: surfPatchID
      character(32) :: surfPatchBody
      character(32) :: surfPatchComp
      character(32) :: surfPatchBCType
   end type strand_surf_patch

   !nEdgePatches copies of the following type follows each strand_header
   type strand_edge_patch
      integer(4)    :: edgePatchID
      character(32) :: edgePatchBody
      character(32) :: edgePatchComp
      character(32) :: edgePatchBCType
      real(8)       :: nx,ny,nz
   end type strand_edge_patch

   type strand2d_header
      character(32) :: surfaceOrVolume
      logical       :: surfaceOnly
      integer(4)    :: nSurfNodes
      integer(4)    :: nSurfFaces
      integer(4)    :: nBndEdges
      integer(4)    :: nPtsPerStrand
      integer(4)    :: nSurfPatches
      integer(4)    :: nEdgePatches
      real(8)       :: strandLength
      real(8)       :: stretchRatio
      real(8)       :: smoothingThreshold
      character(32) :: strandDistribution
   end type strand2d_header

   !nSurfPatches copies of the following type follows each strand_header
   type strand2d_surf_patch
      integer(4)    :: surfPatchID
      character(32) :: surfPatchBody
      character(32) :: surfPatchComp
      character(32) :: surfPatchBCType
   end type strand2d_surf_patch

   !nEdgePatches copies of the following type follows each strand_header
   type strand2d_edge_patch
      integer(4)    :: edgePatchID
      character(32) :: edgePatchBody
      character(32) :: edgePatchComp
      character(32) :: edgePatchBCType
      real(8)       :: nx,ny
   end type strand2d_edge_patch

   type cart_header
      integer(4)    :: nLevels
      integer(4)    :: nBlocks
      integer(4)    :: nFringe
      integer(4)    :: nxc
      integer(4)    :: nyc
      integer(4)    :: nzc
      real(8)       :: domXLo
      real(8)       :: domYLo
      real(8)       :: domZLo
      real(8)       :: domXHi
      real(8)       :: domYHi
      real(8)       :: domZHi
   end type cart_header

   !nBlocks copies of the following type follows each cart_header
   type cart_block
      integer(4),allocatable,dimension(:)   :: iRatio
      integer(4),allocatable,dimension(:,:) :: ilo
      integer(4),allocatable,dimension(:,:) :: ihi
      integer(4),allocatable,dimension(:)   :: levNum
      integer(4),allocatable,dimension(:)   :: iblFLag
      integer(4),allocatable,dimension(:)   :: bdryFLag
   end type cart_block

   !if iblFlag == 2 in cart_block, indicating the block contains both 
   !blanked and unblanked patches, one copy of the following type follows
   type cart_block_iblank
      integer(4)                              :: jd
      integer(4)                              :: kd
      integer(4)                              :: ld
      integer(4),allocatable,dimension(:,:,:) :: iblank
   end type cart_block_iblank

   !if bdryFlag == 1 in cart_block, indicating the block resides on 
   !one or more farfield boundaries, one copy of the following type follows 
   type cart_block_bdry
      integer(4) :: bXLo
      integer(4) :: bYLo
      integer(4) :: bZLo
      integer(4) :: bXHi
      integer(4) :: bYHi
      integer(4) :: bZHi
   end type cart_block_bdry

end module
