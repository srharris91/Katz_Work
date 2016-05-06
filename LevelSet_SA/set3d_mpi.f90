PROGRAM set3d

!*************************************************************************************!
!
! Parallel Min/Max Flow Level Set Solver
!
! Version: 3.0  
!
! Date: 20/4/2016
!
! Type: H-J WENO 5 Serial Explicit TVD R-K Level Set Solver
!
!*************************************************************************************!


!*************************************************************************************!
! Call Modules and Libraries
!*************************************************************************************!

USE set3d_mod
USE set3d_SUBs
USE mpi
USE Lib_VTK_IO
!INCLUDE 'mpif.h'

IMPLICIT NONE

!*************************************************************************************!
! Set Data
!*************************************************************************************!
TYPE (SET_smaller),ALLOCATABLE,DIMENSION(:,:,:) :: D_s   ! Data set of phi and gradPhi
TYPE (SET),ALLOCATABLE,DIMENSION(:,:,:) :: D_all   ! Data set of phi and gradPhi and other stuff
!INTEGER :: j,sUnit,nbytePhi,offset,n1,n2,n3,fN,im,jm,km,ip,jp,kp,dd
INTEGER :: j,n1,n2,n3,fN,im,jm,km,ip,jp,kp,dd
INTEGER :: nx,ny,nz,nn
INTEGER,DIMENSION(3) :: iter
REAL :: a1,a2,a3,xLo(3),xHi(3),x1
REAL :: pX1,pX2,pX3,pY1,pY2,pY3,pZ1,pZ2,pZ3,gX,gY,gZ,dis,minD,k1,k2,k3
REAL :: B1,B2,B3,C1,C2,C3,pSx,pSy,pSz,pS,pX,pY,pZ,sgn,ddx,ddy,ddz,gM
REAL :: y1,z1,maxX,maxY,maxZ,minX,minY,minZ,c,f,phiErr,phiErr_all
REAL :: dx,t1,t2,t3,t4,t5,phiX,phiY,phiZ,phiXX,phiYY,phiZZ,phiXZ,phiXY,phiYZ
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: phi,phi1,phi2,phiS,phiN,gradPhiMag,phiO!,phiE!,phi1,phi2
INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: phiSB,phiNB
REAL,ALLOCATABLE,DIMENSION(:,:) :: centroid,surfX
!CHARACTER(LEN=1) :: lf=char(10)
!CHARACTER(LEN=1024) :: extent,origin,spacing,coffset
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: surfElem
CHARACTER header*80,filename*80
INTEGER*2 padding
!INTEGER*4 ntri,iunit,nSurfNode,k,i,n,nnn,p,kk,share,nSurfElem
INTEGER*4 ntri,iunit,nSurfNode,k,i,n,p,kk,share,nSurfElem
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: normals,triangles,nodesT
REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: gridX,grad2Phi,gradMixPhi
INTEGER,DIMENSION(3) :: SIZE_X
REAL :: a,b,pi,aa,bb,period,ax,ay,az,axy,ayz,axz,bxy,by,bz,byz,bxz,cxy,cyz,cxz,pp
REAL :: x,y,z,bx,cx,cy,cz,dMinus,dPlus,gMM
INTEGER ,DIMENSION(3):: order1,orderUp
INTEGER :: solutionType,order2
REAL :: phiErrS
REAL,DIMENSION(3) :: h,convergenceLimit
REAL :: tol,pln
CHARACTER(1024) :: input_fname,MeshFile
REAL :: dxx,CFL
!REAL :: aaa,bbb,ccc,ax1,ax2,ay1,ay2,az1,az2,Dijk


!*************************************************************************************!
! MPI Struct Variables
!*************************************************************************************!
INTEGER :: mpi_set_smaller
INTEGER :: num = 1
INTEGER, DIMENSION(1) :: blens = (/ 4 /), indx = (/ 0 /), old_types =(/ MPI_DOUBLE_PRECISION /)

!*************************************************************************************!
! MPI Group Variables
!*************************************************************************************!
INTEGER :: old_group, new_group
INTEGER :: size_grid
INTEGER,ALLOCATABLE,DIMENSION(:) :: size_grid_ranks
INTEGER :: new_comm
INTEGER :: new_size

!*************************************************************************************!
! MPI Cart Variables
!*************************************************************************************!
INTEGER :: GRID_COMM
INTEGER,DIMENSION(3) :: dims = (/ 0,0,0 /)
LOGICAL,DIMENSION(3) :: pdc = (/ .TRUE.,.TRUE.,.TRUE. /)
LOGICAL :: reorder = .TRUE.

INTEGER :: rank, numtasks
INTEGER,DIMENSION(3) :: coords
INTEGER,ALLOCATABLE,DIMENSION(:)   :: coords_dumb
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: coords_all
INTEGER,DIMENSION(6) :: nbrs
INTEGER,PARAMETER :: DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4, BACK = 5, FRONT = 6       !Positive x=DOWN, y=RIGHT, z=BACK
REAL :: test1, test2            !cpu_time

!*************************************************************************************!
! MPI Passing Variables
!*************************************************************************************!
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: req
INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: stat
INTEGER :: nPass

!*************************************************************************************!
! MPI Stuff
!*************************************************************************************!
INTEGER*4 :: error,rank_old,num_proc_old

CALL MPI_Init(error)
CALL MPI_Comm_size(MPI_COMM_WORLD, num_proc_old, error)
CALL MPI_Comm_rank(MPI_COMM_WORLD, rank_old, error)

!*************************************************************************************!
! Creation of MPI Struct
!*************************************************************************************!

CALL MPI_Type_struct(num, blens, indx, old_types, mpi_set_smaller, error)
CALL MPI_Type_commit(mpi_set_smaller, error)


!*************************************************************************************!
! Creation of MPI Group 
!*************************************************************************************!
CALL MPI_Comm_group(MPI_COMM_WORLD, old_group, error)

size_grid = FLOOR((REAL(num_proc_old))**(1.0/3.0))**3
ALLOCATE(size_grid_ranks(size_grid))

ALLOCATE(coords_all(size_grid,3))
ALLOCATE(coords_dumb(size_grid*3))

DO n = 1, size_grid
   size_grid_ranks(n) = n-1
END DO

CALL MPI_Group_incl(old_group, size_grid, size_grid_ranks, new_group, error)
CALL MPI_Comm_create(MPI_COMM_WORLD, new_group, new_comm, error)

IF (rank_old >= size_grid) THEN
   CALL MPI_Finalize(error)
END IF

CALL MPI_Comm_size(new_comm, new_size, error)

!*************************************************************************************!
! MPI Cartesian Grid construction
!*************************************************************************************!

CALL MPI_Dims_create(new_size, 3, dims, error)
CALL MPI_Cart_create(new_comm, 3, dims, pdc, reorder, GRID_COMM, error)

CALL MPI_Comm_free(new_comm, error)

CALL MPI_Comm_size(GRID_COMM, numtasks, error)
CALL MPI_Comm_rank(GRID_COMM, rank, error)

CALL MPI_Cart_coords(GRID_COMM, rank, 3, coords, error)
CALL MPI_Cart_shift(GRID_COMM, 0, 1, nbrs(LEFT),  nbrs(RIGHT), error)
CALL MPI_Cart_shift(GRID_COMM, 1, 1, nbrs(UP),    nbrs(DOWN),  error)
CALL MPI_Cart_shift(GRID_COMM, 2, 1, nbrs(FRONT), nbrs(BACK),  error)


CALL MPI_Gather(coords,3,MPI_INTEGER,coords_dumb,3,MPI_INTEGER,0,GRID_COMM,error)
IF (rank == 0) THEN
    DO i = 1,size_grid
        DO j = 1, 3
            coords_all(i,j) = coords_dumb((i-1)*3 + j)
        END DO
    END DO
END IF

DEALLOCATE(coords_dumb)
DEALLOCATE(size_grid_ranks)

!*************************************************************************************!
! Input.namelist reading
!*************************************************************************************!

CALL getarg(1,input_fname)
CALL InputRead(rank,TRIM(input_fname),MeshFile,dx,CFL,dd,tol,h,nPass,iter,order1,order2,orderUp,convergenceLimit,solutionType)


!*************************************************************************************!
! Solution Type
!*************************************************************************************!

!solutionType = 1

IF (solutionType == 1) THEN ! Read STL


!*************************************************************************************!
! Import STL Data (every processor)
!*************************************************************************************!


! start system time
CALL cpu_time(t1)

! import .stl data
!CALL getarg(1,filename)

IF (rank == 0) THEN
   WRITE(*,*)
   PRINT*, " Reading in .stl Mesh "
   WRITE(*,*)
END IF


iunit=13
OPEN(unit=iunit,file=meshfile,status='old',access='stream',form='unformatted')

! read .stl header info 
READ(iunit) header
READ(iunit) ntri
   
ALLOCATE(normals(3,ntri))
ALLOCATE(triangles(3,ntri*3))
ALLOCATE(surfElem(ntri,3))
 
! read .stl data
k=1
DO i = 1,ntri
   READ(iunit) normals(1,i),normals(2,i),normals(3,i)
   READ(iunit) triangles(1,k),triangles(2,k),triangles(3,k)
   READ(iunit) triangles(1,k+1),triangles(2,k+1),triangles(3,k+1)
   READ(iunit) triangles(1,k+2),triangles(2,k+2),triangles(3,k+2)
   READ(iunit) padding
  k=k+3
END DO
  
CLOSE(iunit)

ALLOCATE(nodesT(3,ntri*5))
nSurfElem = ntri

! search through data and put into surfX and surfElem style arrays
DO k = 1,ntri
  nodesT(1,k) = 1000000. 
  nodesT(2,k) = 1000000. 
  nodesT(3,k) = 1000000. 
END DO

! eliminate repeated nodes and clean up arrays
i = 1
nSurfNode = 3
k = 0;
DO n = 1,ntri
   DO p = 1,3
      share = 0 
      DO kk = 1,nSurfNode
         IF ((abs(nodesT(1,kk) - triangles(1,i)) < 1.e-13) .AND. &
             (abs(nodesT(2,kk) - triangles(2,i)) < 1.e-13) .AND. &
             (abs(nodesT(3,kk) - triangles(3,i)) < 1.e-13)) THEN
            share = kk
            EXIT
         END IF
      END DO
      IF (share > 0) THEN
         surfElem(n,p) = share
      ELSE
         k             = k+1 
         nodesT(:,k)   = triangles(:,i)
         surfElem(n,p) = k !1-based
      END IF
      i = i+1
   END DO
   nSurfNode = k 
END DO

CALL MPI_Barrier(GRID_COMM, error)
IF (rank == 0) THEN
!   WRITE(*,*) "Made it here"
   WRITE(*,*) nSurfNode
END IF

! allocate surfX
ALLOCATE(surfX(nSurfNode,3))

! fill in surface node data 
DO k = 1,nSurfNode
   surfX(k,1) = nodesT(1,k)
   surfX(k,2) = nodesT(2,k)
   surfX(k,3) = nodesT(3,k)
END DO

! deallocate unnecessary data
DEALLOCATE(nodesT)
DEALLOCATE(triangles)
DEALLOCATE(normals)


!*************************************************************************************!
! Determine xLo and xHi
!*************************************************************************************!

! initialize
x1 = surfX(1,1);
y1 = surfX(1,2);
z1 = surfX(1,3);

maxX = x1;
maxY = y1;
maxZ = z1;

minX = x1;
minY = y1;
minZ = z1;

! find the max and min
DO n = 2,nSurfNode 
   x1 = surfX(n,1)
   y1 = surfX(n,2)
   z1 = surfX(n,3)

   IF (x1 > maxX) THEN
      maxX = x1
   END IF
   IF (y1 > maxY) THEN
     maxY = y1
   END IF
   IF (z1 > maxZ) THEN
      maxZ = z1
   END IF

   IF (x1 < minX) THEN
      minX = x1
   END IF
   IF (y1 < minY) THEN
     minY = y1
   END IF
   IF (z1 < minZ) THEN
      minZ = z1
   END IF
   
END DO

!*************************************************************************************!
! Define Cartesiang Grid Size and Allocate Phi (every processor)
!*************************************************************************************!

! find the characteristic size of the object to add around
ddx = maxX-minX
ddy = maxY-minY
ddz = maxZ-minZ

! set dx
!dx = 1.;

! define Cartesian grid
nx = ceiling((maxX-minX)/dx)+1;
ny = ceiling((maxY-minY)/dx)+1;
nz = ceiling((maxZ-minZ)/dx)+1;

!nx = MAXVAL((/nx,ny,nz/))
!ny = nx
!nz = nx

! number of cells you want to add
!dd = 10

! adding more cells edge
nx = nx+2*dd
ny = ny+2*dd
nz = nz+2*dd

DO WHILE (MOD(nx+1,dims(2)) /=0)
    nx = nx+1
END DO

DO WHILE (MOD(ny+1,dims(1)) /=0)
    ny = ny+1
END DO

DO WHILE (MOD(nz+1,dims(3)) /=0)
    nz = nz+1
END DO

! set xLo and xHi
xLo =(/minX-dd*dx,minY-dd*dx,minZ-dd*dx/)
xHi =(/maxX+dd*dx,maxY+dd*dx,maxZ+dd*dx/)

 SIZE_X(1) = ((nx+1)/dims(1))-1
 SIZE_X(2) = ((ny+1)/dims(2))-1
 SIZE_X(3) = ((nz+1)/dims(3))-1
IF (rank == 0) THEN
   WRITE(*,*) "SIZE_X(1) = ",SIZE_X(1)
   WRITE(*,*) "SIZE_X(2) = ",SIZE_X(2)
   WRITE(*,*) "SIZE_X(3) = ",SIZE_X(3)
   WRITE(*,*) 
END IF
 
!*************************************************************************************!
! Allocation of Passing Variables
!*************************************************************************************!

ALLOCATE(req(0:3, 1:nPass)) !8*((SIZE_X(2)+1)*(SIZE_X(3)+1))-1))
!ALLOCATE(req_v(0:3, 1:nPass)) !8*((SIZE_X(1)+1)*(SIZE_X(3)+1))-1))
!ALLOCATE(req_d(0:3, 1:nPass)) !8*((SIZE_X(1)+1)*(SIZE_X(2)+1))-1))

IF (rank == 0) THEN
   WRITE(*,*) "Made it here"
END IF
ALLOCATE(stat(MPI_STATUS_SIZE, 0:3, 1:nPass)) !8*((SIZE_X(2)+1)*(SIZE_X(3)+1))-1))
!ALLOCATE(stat_v(MPI_STATUS_SIZE, -1:SIZE_X(2), -1:SIZE_X(3), 0:3, 1:nPass)) !8*((SIZE_X(1)+1)*(SIZE_X(3)+1))-1))
!ALLOCATE(stat_d(MPI_STATUS_SIZE, -1:SIZE_X(1), -1:SIZE_X(2), 0:3, 1:nPass)) !8*((SIZE_X(1)+1)*(SIZE_X(2)+1))-1))

! allocate phi
ALLOCATE(D_s (-nPass:SIZE_X(1)+nPass,-nPass:SIZE_X(2)+nPass,-nPass:SIZE_X(3)+nPass))
ALLOCATE(phi (-nPass:SIZE_X(1)+nPass,-nPass:SIZE_X(2)+nPass,-nPass:SIZE_X(3)+nPass))
ALLOCATE(phi1(-nPass:SIZE_X(1)+nPass,-nPass:SIZE_X(2)+nPass,-nPass:SIZE_X(3)+nPass))
ALLOCATE(phi2(-nPass:SIZE_X(1)+nPass,-nPass:SIZE_X(2)+nPass,-nPass:SIZE_X(3)+nPass))
ALLOCATE(D_all(-1:SIZE_X(1),-1:SIZE_X(2)  ,-1:SIZE_X(3)))
D_s%phi = 1.
D_s%gradPhiX = 0.
D_s%gradPhiY = 0.
D_s%gradPhiZ = 0.

ALLOCATE(gridX(-1:SIZE_X(1),-1:SIZE_X(2),-1:SIZE_X(3),3))

DO i = -1,SIZE_X(1)
   DO j = -1,SIZE_X(2)
      DO k = -1,SIZE_X(3)
         gridX(i,j,k,1) = xLo(1) + ((coords(1)*(SIZE_X(1)+1))+i)*dx;
         gridX(i,j,k,2) = xLo(2) + ((coords(2)*(SIZE_X(2)+1))+j)*dx;
         gridX(i,j,k,3) = xLo(3) + ((coords(3)*(SIZE_X(3)+1))+k)*dx;
      END DO
   END DO
END DO

!------------------------------------------------------------------------------------------------------------------------
! output to vtk files
D_all%x     = gridX(-1:SIZE_X(1),-1:SIZE_X(2),-1:SIZE_X(3),1)
D_all%y     = gridX(-1:SIZE_X(1),-1:SIZE_X(2),-1:SIZE_X(3),2)
D_all%z     = gridX(-1:SIZE_X(1),-1:SIZE_X(2),-1:SIZE_X(3),3)
CALL SET_smaller_to_SET(D_all,D_s,nPass)

CALL output(SIZE_X(1),SIZE_X(2),SIZE_X(3),D_all,'xml_output.pvts','xml_output',rank,dims,coords,coords_all,size_grid)
!------------------------------------------------------------------------------------------------------------------------

!*************************************************************************************!
! Determine Inside and Outside of Surface (every processor)
!*************************************************************************************!

! cut down on excess and use minimum nodes
im = floor((minX-xLo(1))/dx)-3
jm = floor((minY-xLo(2))/dx)-3
km = floor((minZ-xLo(3))/dx)-3
! cut down on excess and use maximum nodes
ip = floor((maxX-xLo(1))/dx)+3
jp = floor((maxY-xLo(2))/dx)+3
kp = floor((maxZ-xLo(3))/dx)+3

! print out grid spacing
IF (rank == 0) THEN
   PRINT*, " Setting Grid Size "
   PRINT*, " Grid Size: nx =",nx,", ny =",ny,",nz =",nz 
   PRINT*, " Grid Spacing: dx =", dx
   WRITE(*,*)
   PRINT*, " Determining Inside and Outside of Geometry "
   WRITE(*,*) 
END IF

! allocate centroid
ALLOCATE(centroid(nSurfElem,4))

DO n = 1,nSurfElem
   n1 = surfElem(n,1)
   n2 = surfElem(n,2)
   n3 = surfElem(n,3)
   pX1 = surfX(n1,1)
   pY1 = surfX(n1,2)
   pZ1 = surfX(n1,3)
   pX2 = surfX(n2,1)
   pY2 = surfX(n2,2)
   pZ2 = surfX(n2,3)
   pX3 = surfX(n3,1)
   pY3 = surfX(n3,2)
   pZ3 = surfX(n3,3)
   centroid(n,1) = (pX1+pX2+pX3)/3.
   centroid(n,2) = (pY1+pY2+pY3)/3.
   centroid(n,3) = (pZ1+pZ2+pZ3)/3.
END DO

!tol = 2.
tol = tol*dx

! find which nodes are inside and outside
fN = 0 !initialize fN
DO i = -1,SIZE_X(1) 
   DO j = -1,SIZE_X(2)  
      DO k = -1,SIZE_X(3)
        pln = gridX(i,j,k,1)
         ! search through all surface elements to find closest element
         minD = 100000.;
         DO n = 1,nSurfElem
          IF(ABS(centroid(n,1) - pln) < tol) THEN
            pX = centroid(n,1)
            pY = centroid(n,2)
            pZ = centroid(n,3)
            gX = gridX(i,j,k,1)
            gY = gridX(i,j,k,2)
            gZ = gridX(i,j,k,3)
            dis = sqrt((pX-gX)*(pX-gX) + (pY-gY)*(pY-gY) + (pZ-gZ)*(pZ-gZ))
            IF (dis < minD) THEN
               minD = dis;
               fN   = n;
            END IF
          END IF
         END DO
    IF (minD < 99999.) THEN
      ! create three vectors from our point to the surfElem points
      n1 = surfElem(fN,1)
      n2 = surfElem(fN,2)
      n3 = surfElem(fN,3)
      A1 = surfX(n1,1) - gridX(i,j,k,1)
      A2 = surfX(n1,2) - gridX(i,j,k,2)
      A3 = surfX(n1,3) - gridX(i,j,k,3)
      B1 = surfX(n2,1) - gridX(i,j,k,1)
      B2 = surfX(n2,2) - gridX(i,j,k,2)
      B3 = surfX(n2,3) - gridX(i,j,k,3)
      C1 = surfX(n3,1) - gridX(i,j,k,1)
      C2 = surfX(n3,2) - gridX(i,j,k,2)
      C3 = surfX(n3,3) - gridX(i,j,k,3)
    
      ! cross product two of the vectors
      pSx = A2*B3-A3*B2;
      pSy = -(A1*B3-B1*A3);
      pSz = A1*B2-B1*A2;

      ! dot product last vector with your cross product result
      pS =-(pSx*C1+pSy*C2+pSz*C3);
      
      gM = 1. ! gradiant magnitude
      ! return sign of the dot product
      CALL phiSign(pS,sgn,dx,gM)
              
      D_s(i,j,k)%phi = sgn
    ELSE
      D_s(i,j,k)%phi = 1.
    END IF
      

      END DO
   END DO
END DO


! deallocate surface mesh data
DEALLOCATE(surfX)
DEALLOCATE(surfElem)
DEALLOCATE(centroid)

IF (rank == 0) THEN
   CALL cpu_time(t2)
   PRINT*, " Search Run Time: ",t2-t1," Seconds"
   WRITE(*,*)
END IF

CALL SET_smaller_to_SET(D_all,D_s,nPass)
CALL MPI_Barrier(GRID_COMM, error)

CALL output(SIZE_X(1),SIZE_X(2),SIZE_X(3),D_all,'in_out.pvts','in_out',rank,dims,coords,coords_all,size_grid)

CALL MPI_Barrier(GRID_COMM, error)
!STOP
!*************************************************************************************!
! Fast Marching Method (parallel processor)
!*************************************************************************************!

CALL MPI_Barrier(GRID_COMM, error)

! allocate arrays used for FMM
ALLOCATE(phiS(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3)))
ALLOCATE(phiO(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3)))
ALLOCATE(phiN(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3)))
ALLOCATE(gradPhiMag(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3)))
!ALLOCATE(phiE(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3)))
!ALLOCATE(phi1(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3)))
!ALLOCATE(phi2(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3)))
!ALLOCATE(phi3(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3)))

IF (rank == 0) THEN
   PRINT*, " Level Set Time Integration "
   WRITE(*,*) 
END IF

OPEN(UNIT=55,FILE='it_time.txt',STATUS='replace')


! set the phi sign array
phiS = D_s(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3))%phi

! number of iterations
!iter = 1000
!iter1 =10000 

!orderUp = 1
!order1 = 4

! time step
!h1 = 0.1
dxx = dx/sqrt(ddx*ddx+ddy*ddy+ddz*ddz) 
h(1) = CFL*dxx*h(1)

! iterate
DO n=0,iter(1)


    phi=D_s%phi

      !********************* Explicit Forward Euler Scheme **************************!

        ! DO i = 1,SIZE_X(1)-1
        !     DO j = 1,SIZE_X(2)-1
        !         DO k = 1,SIZE_X(3)-1
         DO i = 0,SIZE_X(1)
             DO j = 0,SIZE_X(2)
                 DO k = 0,SIZE_X(3)
                     CALL weno(gM,i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi,dPlus,dMinus,nPass,coords,dims)
                     CALL phiSign(phiS(i,j,k),sgn,dx,gM)
                     k1  = sgn*(1.-gM) 
                     phiN(i,j,k) = phi(i,j,k)+h(1)*k1
                 END DO
             END DO
         END DO


    !****************** Explicit Third Order TVD RK Scheme ************************!
!    DO i = 0,SIZE_X(1)
!        DO j = 0,SIZE_X(2)
!            DO k = 0,SIZE_X(3)
!                !CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order1,gMM)
!                CALL upwind(gM,i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi,dPlus,dMinus,orderUp(1),nPass)
!                !CALL weno(gM,i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi,dPlus,dMinus)
!
!                CALL phiSign(phiS(i,j,k),sgn,dx,gM)
!                k1  = sgn*(1.-gM) 
!                D_s(i,j,k)%phi1 = phi(i,j,k)+h(1)*k1
!            END DO
!        END DO
!    END DO
!
!    phi1=D_s%phi1
!
!    DO i = 0,SIZE_X(1)
!        DO j = 0,SIZE_X(2)
!            DO k = 0,SIZE_X(3)
!                !CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order1,gMM)
!                CALL upwind(gM,i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi1,dPlus,dMinus,orderUp(1),nPass)
!                CALL phiSign(phiS(i,j,k),sgn,dx,gM)
!                k2  = sgn*(1.-gM) !sgn*(1.-gM)
!                D_s(i,j,k)%phi2 = 3./4.*phi(i,j,k)+1./4.*phi1(i,j,k)+1./4.*h(1)*k2
!            END DO
!        END DO
!    END DO
!
!    phi2=D_s%phi2
!
!    DO i = 0,SIZE_X(1)
!        DO j = 0,SIZE_X(2)
!            DO k = 0,SIZE_X(3)
!                !CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order1,gMM)
!                CALL upwind(gM,i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi2,dPlus,dMinus,orderUp(1),nPass)
!                CALL phiSign(phiS(i,j,k),sgn,dx,gM)
!                k3  = sgn*(1.-gM) !sgn*(1.-gM)
!                phiN(i,j,k) = 1./3.*phi(i,j,k)+2./3.*phi2(i,j,k)+2./3.*h(1)*k3
!            END DO
!        END DO
!    END DO

   IF (rank == 0) THEN
      CALL cpu_time(test1)
   END IF
    ! extrapolation boundary condition
!   DO i = 0,SIZE_X(1)
!       DO j = 0,SIZE_X(2)
!           DO k = 0,SIZE_X(3)
!               IF (i == 0 .AND. coords(1) == 0) THEN
!                   phiN(i,j,k) = D_s(i+1,j+1,k+1)%phi + dx/sqrt(2.) 
!               ELSE IF (j == 0 .AND. coords(2) == 0) THEN
!                   phiN(i,j,k) = D_s(i+1,j+1,k+1)%phi + dx/sqrt(2.) 
!               ELSE IF (k == 0 .AND. coords(3) == 0) THEN
!                   phiN(i,j,k) = D_s(i+1,j+1,k+1)%phi + dx/sqrt(2.) 
!               ELSE IF (i == SIZE_X(1) .AND. coords(1) == dims(1)-1) THEN
!                   phiN(i,j,k) = D_s(i-1,j-1,k-1)%phi + dx/sqrt(2.) 
!               ELSE IF (j == SIZE_X(2) .AND. coords(2) == dims(2)-1) THEN
!                   phiN(i,j,k) = D_s(i-1,j-1,k-1)%phi + dx/sqrt(2.) 
!               ELSE IF (k == SIZE_X(3) .AND. coords(3) == dims(3)-1) THEN
!                   phiN(i,j,k) = D_s(i-1,j-1,k-1)%phi + dx/sqrt(2.) 
!               END IF
!           END DO
!       END DO
!   END DO


   phiErr = 0.
   phiErr_all = 0.
  
   !phiErr = sum((phi-phiN)
 
   ! calculate RMS
   DO i = 0,SIZE_X(1)
      DO j = 0,SIZE_X(2)
         DO k = 0,SIZE_X(3)
            phiErr = phiErr + (phi(i,j,k)-phiN(i,j,k))**2
         END DO
      END DO
   END DO

   CALL MPI_Barrier(GRID_COMM, error)
   CALL MPI_Allreduce(phiErr, phiErr_all, 1, MPI_DOUBLE_PRECISION,&
       MPI_SUM, GRID_COMM, error)


   ! check error
   phiErr_all = sqrt(phiErr_all/(nx*ny*nz))

      IF (phiErr_all < convergenceLimit(1)) THEN
       IF (rank == 0) THEN
         PRINT*, " Distance function time integration has reached steady state "
       END IF
         EXIT
      END IF
   
   D_s(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3))%phi = phiN

   IF (rank == 0) THEN
      PRINT*, " Iteration: ",n," ", " RMS Error: ",phiErr_all
   END IF

   ! check for NAN
   IF (isnan(phiErr_all)) STOP 
       
   !***************************** Passing Phi ************************************!

   CALL MPI_Barrier(GRID_COMM, error)


CALL Pass(D_s,nbrs,SIZE_X,dims,req,stat,mpi_set_smaller,GRID_COMM,error,nPass)
  
   IF (rank == 0) THEN
      CALL cpu_time(test2)
      WRITE(55,*) n, test2-test1
   END IF
END DO
   

CLOSE(55)
!WRITE(*,*)

! correct distances
!phi = phi*dx

phiO = D_s(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3))%phi
! print out run time
IF (rank == 0) THEN
   CALL cpu_time(t3)
   PRINT*, " Initialization Run Time: ",t3-t1," Seconds"
   WRITE(*,*)
END IF

ELSEIF (solutionType == 2) THEN ! MMS

    ! set dx
    !dx = 1.;

    ! define Cartesian grid
    nx = ceiling((2.)/dx)+1;
    ny = ceiling((2.)/dx)+1;
    nz = ceiling((2.)/dx)+1;

    ! set xLo and xHi
    xLo =(/-1.,-1.,-1./)
    xHi =(/1.,1.,1./)

    ! allocate phi
    ALLOCATE(phi(0:nx,0:ny,0:nz))
    phi = 1.

    ! allocate a grid of x,y,z points 
    !ALLOCATE(gridX(0:nx,0:ny,0:nz,3))
    DO i = 0,SIZE_X(1)
        DO j = 0,SIZE_X(2)
            DO k = 0,SIZE_X(3)
                gridX(i,j,k,1) = xLo(1) + i*dx
                gridX(i,j,k,2) = xLo(2) + j*dx
                gridX(i,j,k,3) = xLo(3) + k*dx
            END DO
        END DO
    END DO



    period = 10.
    pi = 4.*atan(1.);
    a = .1
    aa = a/3.
    bb = 2.*pi/period;


    ! set phi
    DO i = 0,SIZE_X(1)
        DO j = 0,SIZE_X(2)
            DO k = 0,SIZE_X(3)

                CALL RANDOM_NUMBER (b)      
                ax     = 2.*a*(b-.5); ! random number between -a and a
                CALL RANDOM_NUMBER (b)
                ay     = 2.*a*(b-.5); ! random number between -a and a
                CALL RANDOM_NUMBER (b)
                az     = 2.*a*(b-.5); ! random number between -a and a
                CALL RANDOM_NUMBER (b)
                axy    = 2.*a*(b-.5); ! random number between -a and a    
                CALL RANDOM_NUMBER (b)
                axz    = 2.*a*(b-.5); ! random number between -a and a
                CALL RANDOM_NUMBER (b)
                ayz    = 2.*a*(b-.5); ! random number between -a and a

                CALL RANDOM_NUMBER (b)
                c      = 0.;
                IF (b < .5) c = 1.; ! 0. or 1. randomly
                bx     = c*.5*pi;
                CALL RANDOM_NUMBER (b)
                c      = 0.;
                IF (b < .5) c = 1.; ! 0. or 1. randomly
                by     = c*.5*pi;
                CALL RANDOM_NUMBER (b)
                c      = 0.;
                IF (b < .5) c = 1.; ! 0. or 1. randomly
                bz     = c*.5*pi;
                CALL RANDOM_NUMBER (b)
                c      = 0.;
                IF (b < .5) c = 1.; ! 0. or 1. randomly
                bxy    = c*.5*pi;
                CALL RANDOM_NUMBER (b)
                c      = 0.;
                IF (b < .5) c = 1.; ! 0. or 1. randomly
                bxz    = c*.5*pi;
                CALL RANDOM_NUMBER (b)
                c      = 0.;
                IF (b < .5) c = 1.; ! 0. or 1. randomly
                byz    = c*.5*pi;

                CALL RANDOM_NUMBER (b)
                b      = 1.+2.*A*(b-.5); ! random number between 1.-A and 1.+A
                cx     = b*bb;
                CALL RANDOM_NUMBER (b)
                b      = 1.+2.*A*(b-.5); ! random number between 1.-A and 1.+A
                cy     = b*bb;
                CALL RANDOM_NUMBER (b)
                b      = 1.+2.*A*(b-.5); ! random number between 1.-A and 1.+A
                cz     = b*bb;
                CALL RANDOM_NUMBER (b)
                b      = 1.+2.*A*(b-.5); ! random number between 1.-A and 1.+A
                cxy    = b*bb/period;
                CALL RANDOM_NUMBER (b)
                b      = 1.+2.*A*(b-.5); ! random number between 1.-A and 1.+A
                cxz    = b*bb/period;
                CALL RANDOM_NUMBER (b)
                b      = 1.+2.*A*(b-.5); ! random number between 1.-A and 1.+A
                cyz    = b*bb/period;

                x = gridX(i,j,k,1)
                y = gridX(i,j,k,2)
                z = gridX(i,j,k,3)


                phi(i,j,k) =  sin(x*y/0.1) !+ay*sin(cy*y)+axy*sin(bxy*x*y)                   
            END DO
        END DO
    END DO



END IF


!*************************************************************************************!
! Paraview Output (parallel process)
!*************************************************************************************!
IF(rank == 0) THEN
    PRINT*, " Writing Out Cartesian Grid to Paraview Format "
    WRITE(*,*) 
END IF


CALL SET_smaller_to_SET(D_all,D_s,nPass)
CALL output(SIZE_X(1),SIZE_X(2),SIZE_X(3),D_all,'blockBefore.pvts','blockBefore',rank,dims,coords,coords_all,size_grid)

!STOP

!*************************************************************************************!
! Determine Narrow Band
!*************************************************************************************!

ALLOCATE(phiSB(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3)))
ALLOCATE(phiNB(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3)))

CALL narrowBand(SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,D_s,phiNB,phiSB,nPass)

!*************************************************************************************!
! Initialize Gradients
!*************************************************************************************!


!ALLOCATE(gradPhi(0:nx,0:ny,0:nz,3))
ALLOCATE(grad2Phi(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3),3))
ALLOCATE(gradMixPhi(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3),3))

D_all%gradPhiXX =0.
D_all%gradPhiYY =0.
D_all%gradPhiZZ =0.
D_all%gradPhiXY =0.
D_all%gradPhiXZ =0.
D_all%gradPhiYZ =0.
D_all%gradPhiMag=0.
grad2Phi = 0.
gradMixPhi = 0.
gradPhiMag = 0.

phiN = D_s(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3))%phi


!orderUp = 1
!order1 = 2
!order2 = 2

!*************************************************************************************!
! Min/Max Flow (parallel process)
!*************************************************************************************!

!iter2 = 10000
!h1 = 0.1
h(2) = h(2)*dx*dx

!phi(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3))     = phiN

DO n = 1,iter(2)
    !****************** Explicit Third Order TVD RK Stage 1 ***********************!

    ! Calculate first derivative if it falls within stencil band
    phi=D_s%phi
    DO i = 0,SIZE_X(1)
        DO j = 0,SIZE_X(2)
            DO k = 0,SIZE_X(3)
                IF (phiSB(i,j,k) == 1) THEN
                    CALL firstDeriv(i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi,phiX,phiY,phiZ,order1(2),gMM,nPass,coords,dims)
                    D_s(i,j,k)%gradPhiX = phiX
                    D_s(i,j,k)%gradPhiY = phiY
                    D_s(i,j,k)%gradPhiZ = phiZ   
                END IF
            END DO
        END DO
    END DO

    ! Calculate second derivative flow if it is in the narrow band
    DO i = 0,SIZE_X(1)
        DO j = 0,SIZE_X(2)
            DO k = 0,SIZE_X(3)
                IF (phiNB(i,j,k) == 1) THEN
                    CALL secondDeriv(i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi &
                        ,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,order2,nPass)
                    D_all(i,j,k)%gradPhiXX = phiXX
                    D_all(i,j,k)%gradPhiYY = phiYY
                    D_all(i,j,k)%gradPhiZZ = phiZZ
                    D_all(i,j,k)%gradPhiXY = phiXY
                    D_all(i,j,k)%gradPhiXZ = phiXZ
                    D_all(i,j,k)%gradPhiYZ = phiYZ          
                END IF

            END DO
        END DO
    END DO

    ! Caluclate the min/max flow
    DO i = 0,SIZE_X(1)
        DO j = 0,SIZE_X(2)
            DO k = 0,SIZE_X(3)
                IF (phiNB(i,j,k) == 1) THEN
                    CALL firstDeriv(i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi,phiX,phiY,phiZ,order1(2),gMM,nPass,coords,dims)
                    CALL minMax(i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),phi,D_s,D_all,F,nPass)
                    k1 = F !F*gMM !max(F,0.)*dPlus + min(F,0.)*dMinus ! F*gM  
                    phi1(i,j,k) = D_s(i,j,k)%phi + h(2)*k1
                END IF

            END DO
        END DO
    END DO

    !phi1 = D_s%phi1

    !****************** Explicit Third Order TVD RK Stage 2 ***********************!

    ! Calculate first derivative if it falls within stencil band
    DO i = 0,SIZE_X(1)
        DO j = 0,SIZE_X(2)
            DO k = 0,SIZE_X(3)
                IF (phiSB(i,j,k) == 1) THEN
                    CALL firstDeriv(i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi1,phiX,phiY,phiZ,order1(2),gMM,nPass,coords,dims)
                    D_s(i,j,k)%gradPhiX = phiX
                    D_s(i,j,k)%gradPhiY = phiY
                    D_s(i,j,k)%gradPhiZ = phiZ   
                END IF
            END DO
        END DO
    END DO

    ! Calculate second derivative flow if it is in the narrow band
    DO i = 0,SIZE_X(1)
        DO j = 0,SIZE_X(2)
            DO k = 0,SIZE_X(3)
                IF (phiNB(i,j,k) == 1) THEN
                    CALL secondDeriv(i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi1 &
                        ,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,order2,nPass)
                    D_all(i,j,k)%gradPhiXX = phiXX
                    D_all(i,j,k)%gradPhiYY = phiYY
                    D_all(i,j,k)%gradPhiZZ = phiZZ
                    D_all(i,j,k)%gradPhiXY  = phiXY
                    D_all(i,j,k)%gradPhiXZ  = phiXZ
                    D_all(i,j,k)%gradPhiYZ  = phiYZ          
                END IF

            END DO
        END DO
    END DO

    ! Caluclate the min/max flow
    DO i = 0,SIZE_X(1)
        DO j = 0,SIZE_X(2)
            DO k = 0,SIZE_X(3)
                IF (phiNB(i,j,k) == 1) THEN
                    CALL firstDeriv(i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi1,phiX,phiY,phiZ,order1(2),gMM,nPass,coords,dims)
                    CALL minMax(i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),phi1,D_s,D_all,F,nPass)
                    k2 = F !F*gMM !max(F,0.)*dPlus + min(F,0.)*dMinus ! F*gM  
                    phi2(i,j,k) = 3./4.*D_s(i,j,k)%phi + 1./4.*phi1(i,j,k) + 1./4.*h(2)*k2
                END IF

            END DO
        END DO
    END DO

    !phi2 = D_s%phi2

    !****************** Explicit Third Order TVD RK Stage 3 ***********************!

    ! Calculate first derivative if it falls within stencil band
    DO i = 0,SIZE_X(1)
        DO j = 0,SIZE_X(2)
            DO k = 0,SIZE_X(3)
                IF (phiSB(i,j,k) == 1) THEN
                    CALL firstDeriv(i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi2,phiX,phiY,phiZ,order1(2),gMM,nPass,coords,dims)
                    D_s(i,j,k)%gradPhiX = phiX
                    D_s(i,j,k)%gradPhiY = phiY
                    D_s(i,j,k)%gradPhiZ = phiZ   
                END IF
            END DO
        END DO
    END DO

    ! Calculate second derivative flow if it is in the narrow band
    DO i = 0,SIZE_X(1)
        DO j = 0,SIZE_X(2)
            DO k = 0,SIZE_X(3)
                IF (phiNB(i,j,k) == 1) THEN
                    CALL secondDeriv(i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi2 &
                        ,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,order2,nPass)
                    D_all(i,j,k)%gradPhiXX = phiXX
                    D_all(i,j,k)%gradPhiYY = phiYY
                    D_all(i,j,k)%gradPhiZZ = phiZZ
                    D_all(i,j,k)%gradPhiXY = phiXY
                    D_all(i,j,k)%gradPhiXZ = phiXZ
                    D_all(i,j,k)%gradPhiYZ = phiYZ          
                END IF

            END DO
        END DO
    END DO

    ! Caluclate the min/max flow
    DO i = 0,SIZE_X(1)
        DO j = 0,SIZE_X(2)
            DO k = 0,SIZE_X(3)
                IF (phiNB(i,j,k) == 1) THEN
                    CALL firstDeriv(i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi2,phiX,phiY,phiZ,order1(2),gMM,nPass,coords,dims)
                    CALL minMax(i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),phi2,D_s,D_all,F,nPass)
                    k3 = F !F*gMM !max(F,0.)*dPlus + min(F,0.)*dMinus ! F*gM  
                    phiN(i,j,k) = 1./3.*D_s(i,j,k)%phi + 2./3.*phi2(i,j,k) + 2./3.*h(2)*k3
                END IF

            END DO
        END DO
    END DO


   phiErr = 0.
   phiErr_all = 0.
  
   !phiErr = sum((phi-phiN)
 
   ! calculate RMS
   DO i = 0,SIZE_X(1)
      DO j = 0,SIZE_X(2)
         DO k = 0,SIZE_X(3)
            phiErr = phiErr + (phi(i,j,k)-phiN(i,j,k))**2
         END DO
      END DO
   END DO


CALL MPI_Barrier(GRID_COMM, error)
CALL MPI_Allreduce(phiErr, phiErr_all, 1, MPI_DOUBLE_PRECISION,&
    MPI_SUM, GRID_COMM, error)


   ! check error
   phiErr_all = sqrt(phiErr_all/(nx*ny*nz))

      IF (phiErr_all < convergenceLimit(2)) THEN
       IF (rank == 0) THEN
         PRINT*, " Distance function time integration has reached steady state "
       END IF
         EXIT
      END IF

   D_s(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3))%phi = phiN
   !D_s%phi = phiN

   IF (rank == 0) THEN
      PRINT*, " Iteration: ",n," ", " RMS Error: ",phiErr_all
   END IF
  
   !***************************** Passing Phi ************************************!

CALL Pass(D_s,nbrs,SIZE_X,dims,req,stat,mpi_set_smaller,GRID_COMM,error,nPass)
!   CALL MPI_Barrier(GRID_COMM, error)
  

   ! set phi value for sign
   phiS = D_s(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3))%phi





   !orderUp = 1
   !h2 = .08 !.0001
   !h(2) = h(2)*dx*dx
   ! reinitialize phi for a number of time steps
   IF (MOD(n,100)==10000000) THEN
       DO nn=1,20000

           !****************** Explicit Third Order TVD RK Scheme ************************!

           DO i = 0,SIZE_X(1)
               DO j = 0,SIZE_X(2)
                   DO k = 0,SIZE_X(3)
                       CALL upwind(gM,i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,D_s%phi,dPlus,dMinus,orderUp(2),nPass)
                       CALL phiSign(phiS(i,j,k),sgn,dx,gM)
                       k1  = sgn*(1.-gM) 
                       phi1(i,j,k) = D_s(i,j,k)%phi+h(2)*k1
                   END DO
               END DO
           END DO

           DO i = 0,SIZE_X(1)
               DO j = 0,SIZE_X(2)
                   DO k = 0,SIZE_X(3)
                       CALL upwind(gM,i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi1,dPlus,dMinus,orderUp(2),nPass)
                       CALL phiSign(phiS(i,j,k),sgn,dx,gM)
                       k2  = sgn*(1.-gM) 
                       phi2(i,j,k) = 3./4.*D_s(i,j,k)%phi+1./4.*phi1(i,j,k)+1./4.*h(2)*k2
                   END DO
               END DO
           END DO

           DO i = 0,SIZE_X(1)
               DO j = 0,SIZE_X(2)
                   DO k = 0,SIZE_X(3)
                       CALL upwind(gM,i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi2,dPlus,dMinus,orderUp(2),nPass)
                       CALL phiSign(phiS(i,j,k),sgn,dx,gM)
                       k3  = sgn*(1.-gM) 
                       phiN(i,j,k) = 1./3.*D_s(i,j,k)%phi+2./3.*phi2(i,j,k)+2./3.*h(2)*k3
                   END DO
               END DO
           END DO

           phiErrS = 0.
           phiErr_all=0.
           pp = 0.

           DO i = 0,SIZE_X(1)
               DO j = 0,SIZE_X(2)
                   DO k = 0,SIZE_X(3)
                       !IF (phiSB(i,j,k) == 1) THEN         
                       phiErrS = phiErrS + (D_s(i,j,k)%phi-phiN(i,j,k))*(D_s(i,j,k)%phi-phiN(i,j,k))
                       pp = pp+1.
                       !END IF
                   END DO
               END DO
           END DO
           CALL MPI_Barrier(GRID_COMM, error)
           CALL MPI_Allreduce(phiErrS, phiErr_all, 1, MPI_DOUBLE_PRECISION,&
               MPI_SUM, GRID_COMM, error)

           ! check error
           phiErr_all = sqrt(phiErr_all/(pp))
           IF (phiErr_all < convergenceLimit(3)) THEN
               PRINT*, " Distance function time integration has reached steady state "
               EXIT
           END IF

           D_s%phi = phiN

           PRINT*, " #:",n," ","RMS:",phiErr,"Sub #:",nn," ","Sub-RMS:",phiErrS

       END DO
   END IF

   !****************************** Narrow Band **********************************!

   CALL narrowBand(SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,D_s,phiNB,phiSB,nPass)

!   h(2) = 0.0001
!   ! reinitialize phi for 10 time steps
!   DO nn=1,0
!      ! call upwind routine
!      
!   DO i = 0,SIZE_X(1)
!      DO j = 0,SIZE_X(2)
!         DO k = 0,SIZE_X(3)
!               IF (phiNB(i,j,k) == 1) THEN
!                  CALL upwind(gM,i,j,k,nx,ny,nz,dx,D_s)
!                  CALL phiSign(phiS(i,j,k),sgn)
!                  k1  = sgn*(1.-gM)
!                  phiN(i,j,k) = D_s(i,j,k)%phi+h(2)*k1
!               END IF
!            END DO
!         END DO
!      END DO
!      D_s%phi = phiN
!   END DO


END DO
!WRITE(*,*)


! print out run time
IF(rank == 0) THEN
    CALL cpu_time(t4)
    PRINT*, " Total Run Time: ",t4-t1," Seconds"
    WRITE(*,*)
END IF

!*************************************************************************************!
! Fast Marching Method
!*************************************************************************************!

IF (rank == 0) THEN
    PRINT*, " Level Set Reinitialization "
    PRINT*, 
END IF

! set the phi sign array
phiS = D_s(0:SIZE_X(1),0:SIZE_X(2),0:SIZE_X(3))%phi

! number of iterations
!iter(3) = 0

! time step
!h3 = .01
h(3) = h(3)*dx*dx

! order of accuracy for the upwinding
!orderUp = 1
!order1 = 4

! iterate
!DO n=0,iter(3)
!
!    !****************** Explicit Third Order TVD RK Scheme ************************!
!
!    DO i = 0,SIZE_X(1)
!        DO j = 0,SIZE_X(2)
!            DO k = 0,SIZE_X(3)
!                !CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order1(3),gMM)
!                CALL upwind(gM,i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,D_s%phi,dPlus,dMinus,orderUp(3),nPass)
!                CALL phiSign(phiS(i,j,k),sgn)
!                k1  = sgn*(1.-gM) 
!                D_s(i,j,k)%phi1 = D_s(i,j,k)%phi+h(3)*k1
!            END DO
!        END DO
!    END DO
!
!    DO i = 0,SIZE_X(1)
!        DO j = 0,SIZE_X(2)
!            DO k = 0,SIZE_X(3)
!                !CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order1(3),gMM)
!                CALL upwind(gM,i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,D_s%phi1,dPlus,dMinus,orderUp(3),nPass)
!                CALL phiSign(phiS(i,j,k),sgn)
!                k2  = sgn*(1.-gM) !sgn*(1.-gM)
!                D_s(i,j,k)%phi2 = 3./4.*D_s(i,j,k)%phi+1./4.*D_s(i,j,k)%phi1+1./4.*h(3)*k2
!            END DO
!        END DO
!    END DO
!
!
!    DO i = 0,SIZE_X(1)
!        DO j = 0,SIZE_X(2)
!            DO k = 0,SIZE_X(3)
!                !CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order1(3),gMM)
!                CALL upwind(gM,i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,D_s%phi2,dPlus,dMinus,orderUp(3),nPass)
!                CALL phiSign(phiS(i,j,k),sgn)
!                k3  = sgn*(1.-gM) !sgn*(1.-gM)
!                phiN(i,j,k) = 1./3.*D_s(i,j,k)%phi+2./3.*D_s(i,j,k)%phi2+2./3.*h(3)*k3
!            END DO
!        END DO
!    END DO
!
!    ! extrapolation boundary condition
!    DO i = 0,SIZE_X(1)
!        DO j = 0,SIZE_X(2)
!            DO k = 0,SIZE_X(3)
!                IF (i == 0) THEN
!                    phiN(i,j,k) = D_s(1,j,k)%phi + dx/sqrt(2.) 
!                ELSE IF (j == 0) THEN
!                    phiN(i,j,k) = D_s(i,1,k)%phi + dx/sqrt(2.) 
!                ELSE IF (k == 0) THEN
!                    phiN(i,j,k) = D_s(i,j,1)%phi + dx/sqrt(2.) 
!                ELSE IF (i == nx) THEN
!                    phiN(i,j,k) = D_s(nx-1,j,k)%phi + dx/sqrt(2.) 
!                ELSE IF (j == ny) THEN
!                    phiN(i,j,k) = D_s(i,ny-1,k)%phi + dx/sqrt(2.) 
!                ELSE IF (k == nz) THEN
!                    phiN(i,j,k) = D_s(i,j,nz-1)%phi + dx/sqrt(2.) 
!                END IF
!            END DO
!        END DO
!    END DO
!
!    !********************************* RMS ***************************************!
!
!    phiErr = 0.
!
!    ! calculate RMS
!    DO i = 0,SIZE_X(1)
!        DO j = 0,SIZE_X(2)
!            DO k = 0,SIZE_X(3)
!                phiErr = phiErr + (D_s(i,j,k)%phi-phiN(i,j,k))*(D_s(i,j,k)%phi-phiN(i,j,k))
!            END DO
!        END DO
!    END DO
!
!    ! check error
!    phiErr = sqrt(phiErr/(nx*ny*nz))
!    IF (phiErr < 1.E-10) THEN
!        PRINT*, " Distance function time integration has reached steady state "
!        EXIT
!    END IF
!
!    ! set new phi 
!    phi = phiN
!
!    PRINT*, " Iteration: ",n," ", " RMS Error: ",phiErr
!
!    ! check for NAN
!    IF (isnan(phiErr)) STOP 
!
!
!END DO
!PRINT*,
!
!CALL MPI_BARRIER(GRID_COMM, error)
!
!! print out run time
!CALL cpu_time(t5)
!PRINT*, " Reinitialization Run Time: ",t5-t1," Seconds"
!PRINT*,


!*************************************************************************************!
! Asymptotic Error
!*************************************************************************************!

phiErr = 0.

! calculate RMS
DO i = 0,SIZE_X(1)
    DO j = 0,SIZE_X(2)
        DO k = 0,SIZE_X(3)
            phiErr = phiErr + (D_s(i,j,k)%phi-phiO(i,j,k))*(D_s(i,j,k)%phi-phiO(i,j,k))
        END DO
    END DO
END DO

phiErr = sqrt(phiErr/(nx*ny*nz))
IF (rank == 0) THEN
    PRINT*, " Asymptotic Error: ",phiErr
END IF

!*************************************************************************************!
! Output Grad Phi Mag
!*************************************************************************************!


! grad phi
!order1 = 4
DO i = 0,SIZE_X(1)
    DO j = 0,SIZE_X(2)
        DO k = 0,SIZE_X(3)
            CALL firstDeriv(i,j,k,SIZE_X(1),SIZE_X(2),SIZE_X(3),dx,phi,phiX,phiY,phiZ,order1(3),gMM,nPass,coords,dims)
            D_all(i,j,k)%gradPhiMag = gMM
        END DO
    END DO
END DO




!*************************************************************************************!
! Paraview Output
!*************************************************************************************!

IF(rank == 0) THEN
    PRINT*, " Writing Out Cartesian Grid to Paraview Format "
    WRITE(*,*) 
END IF

CALL SET_smaller_to_SET(D_all,D_s,nPass)
CALL output(SIZE_X(1),SIZE_X(2),SIZE_X(3),D_all,'blockAfter.pvts','blockAfter',rank,dims,coords,coords_all,size_grid)

!*************************************************************************************!
! Deallocate arrays
!*************************************************************************************!

DEALLOCATE(D_s)
DEALLOCATE(phiO)
DEALLOCATE(D_all)
DEALLOCATE(phiS)
DEALLOCATE(phiSB)
DEALLOCATE(phiNB)
DEALLOCATE(grad2Phi)
DEALLOCATE(gradPhiMag)
DEALLOCATE(gradMixPhi)
DEALLOCATE(coords_all)
DEALLOCATE(req)
!DEALLOCATE(req_v)
!DEALLOCATE(req_d)
DEALLOCATE(stat)
!DEALLOCATE(stat_v)
!DEALLOCATE(stat_d)
DEALLOCATE(gridX)
CALL MPI_Finalize(error)

!*************************************************************************************!
! Program End
!*************************************************************************************!

END PROGRAM set3d
