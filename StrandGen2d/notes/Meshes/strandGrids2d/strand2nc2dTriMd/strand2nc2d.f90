PROGRAM strand2nc2d

IMPLICIT NONE

CHARACTER(12) :: inputFile
CHARACTER(64) :: strandFile,nc2dFile,arg
INTEGER :: inputFileSize,iwrite,strandDist,nPtsPerStrand,nFringe,trimType
INTEGER :: istrandunit,surfaceOnly,nSurfNodes,nSurfFaces,nBndNodes, &
           nSurfPatches,nNodePatches
INTEGER :: k,n,j,nmax,n1,n2,p1,p2,p3,p4,m,narg
INTEGER :: inc2dunit,nverts,ncells,nedgeb,npoibf,ncomps
INTEGER :: ilo,ihi,ilon,ihin,ig1,ig2,ig,ig1n,ign,nNodes,nBnodes,csp1d,psp1d
REAL :: stretchRatio,wallSpacing,strandLength,deltaSmooth
REAL :: rSmoothL2,xn,yn,nx,ny,ds
INTEGER,ALLOCATABLE,DIMENSION(:) :: bndNodes,nodeTag,faceTag,nodeClip, &
                                    faceClip,csp1,csp2,psp1,psp2
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: surfFaces,face,bNode,cell,edgb,bpts,flag
REAL,ALLOCATABLE,DIMENSION(:) :: xStrand,xs,ang
REAL,ALLOCATABLE,DIMENSION(:,:) :: bndNormal,xSurf,pointingVec,bNorm,pv0,x


narg = iargc()
m    = 0
DO WHILE (m < narg)
   m    = m+1
   CALL GETARG(m,arg)
   IF (arg(1:2) == '-i') THEN
      m = m+1
      CALL GETARG(m,arg)
      strandFile = arg
   ELSE IF (arg(1:2) == '-o') THEN
      m = m+1
      CALL GETARG(m,arg)
      nc2dFile = arg
   END IF
END DO


! read strand grid input file
inputFile = 'strand.input'
inputFileSize = 12
CALL readinput(inputFile, &
               inputFileSize, &
               iwrite, &
               strandDist, &
               nPtsPerStrand, &
               nFringe, &
               stretchRatio, &
               wallSpacing, &
               strandLength, &
               deltaSmooth, &
               trimType);


! read strand grid mesh header
OPEN(istrandunit,FILE=strandFile,STATUS='old',FORM='unformatted', &
     CONVERT='big_endian')
READ(istrandunit)surfaceOnly
READ(istrandunit)nSurfNodes
READ(istrandunit)nSurfFaces
READ(istrandunit)nBndNodes
READ(istrandunit)nPtsPerStrand
READ(istrandunit)nSurfPatches
READ(istrandunit)nNodePatches


! allocate strand grid data
ALLOCATE(surfFaces(2,nSurfFaces))
ALLOCATE(bndNodes(nBndNodes))
ALLOCATE(bndNormal(2,nBndNodes))
ALLOCATE(nodeTag(nBndNodes))
ALLOCATE(faceTag(nSurfFaces))
ALLOCATE(xSurf(2,nSurfNodes))
ALLOCATE(nodeClip(nSurfNodes))
ALLOCATE(faceClip(nSurfFaces))
ALLOCATE(pointingVec(2,nSurfNodes))
ALLOCATE(xStrand(0:nPtsPerStrand))

! read strand grid data
READ(istrandunit)((surfFaces  (k,n),k=1,2),n=1,nSurfFaces   )
READ(istrandunit) (bndNodes   (  n)       ,n=1,nBndNodes    )
READ(istrandunit)((bndNormal  (k,n),k=1,2),n=1,nBndNodes    )
READ(istrandunit) (faceTag    (  n)       ,n=1,nSurfFaces   )
READ(istrandunit) (nodeTag    (  n)       ,n=1,nBndNodes    )
READ(istrandunit)((xSurf      (k,n),k=1,2),n=1,nSurfNodes   )
READ(istrandunit) (nodeClip   (  n)       ,n=1,nSurfNodes   )
READ(istrandunit) (faceClip   (  n)       ,n=1,nSurfFaces   )
READ(istrandunit)((pointingVec(k,n),k=1,2),n=1,nSurfNodes   )
IF (surfaceOnly == 0) READ(istrandunit) (xStrand(n),n=0,nPtsPerStrand)
CLOSE(istrandunit)


! generate 1d stretching function if necessary
IF (surfaceOnly == 1) THEN
   nmax = 10000;
   ALLOCATE(xs(0:nmax))
   CALL strand1ddist(nmax, &
                     xs, &
                     nPtsPerStrand, &
                     stretchRatio, &
                     strandDist, &
                     strandLength, &
                     wallSpacing)
   DEALLOCATE(xStrand)
   ALLOCATE(xStrand(0:nPtsPerStrand))
   DO n=0,nPtsPerStrand
      xStrand(n) = xs(n)
   END DO
   DEALLOCATE(xs)
END IF


! smooth normals
ilo = 1
ihi = nSurfFaces
ilon = 1
ihin = nSurfNodes
ig1 = 0
ig2 = 0
ig = 0
ig1n = 0
ign = 0
nNodes = nSurfNodes
nBnodes = 0
ALLOCATE(face(2,nSurfFaces))
ALLOCATE(bNode(2,nBnodes))
ALLOCATE(bNorm(2,nBnodes))
face = surfFaces

ALLOCATE(csp2(nNodes+1))
ALLOCATE(psp2(nNodes+1))
CALL findcspdim(ihi, &
                ilo, &
                ilon, &
                ihin, &
                ig1, &
                ig2, &
                ig, &
                ign, &
                csp2, &
                face, &
                csp1d)
ALLOCATE(csp1(csp1d))
CALL fillcsp(ihi, &
	     ilo, &
             ilon, &
             ihin, &
             ig1, &
             ig2, &
             ig, &
             ign, &
             csp1d, &
             face, &
             csp1, &
             csp2)
CALL findpspdim(ihi, &
		ilo, &
                ilon, &
                ihin, &
                ig, &
                ign, &
                csp1d, &
                psp2, &
                face, &
                csp1, &
                csp2, &
                psp1d)
ALLOCATE(psp1(psp1d))
CALL fillpsp(ihi, &
	     ilo, &
             ilon, &
             ihin, &
             ig, &
             ign, &
             csp1d, &
             psp1d, &
             face, &
             csp1, &
             csp2, &
             psp1, &
             psp2)
ALLOCATE(pV0(2,nNodes))
ALLOCATE(ang(nNodes))
CALL initpointingvec(ilo, &
		     ihi, &
                     ilon, &
                     ihin, &
                     ig1, &
                     ig2, &
                     ig, &
                     ig1n, &
                     ign, &
                     nBnodes, &
                     csp1d, &
                     pointingVec, &
                     pV0, &
                     ang, &
                     csp1, &
                     csp2, &
                     xSurf, &
                     face, &
                     bNode, &
                     bNorm)

n = 0
smooth: DO
   n = n+1
   CALL smoothingiter(ilo, &
                      ihi, &
                      ilon, &
                      ihin, &
                      ig, &
                      ign, &
                      nNodes, &
                      nBnodes, &
                      psp1d, &
                      bNode, &
                      bNorm, &
                      pointingVec, &
                      pV0, &
                      ang, &
                      psp1, &
                      psp2, &
                      rSmoothL2)
   WRITE(*,*)'smoothing iter: ',n,rSmoothL2
   IF (rSmoothL2 <= deltaSmooth) EXIT smooth
END DO smooth


! compute nc2d format data
nverts = nSurfNodes*(nPtsPerStrand+1)
ncells = nSurfFaces*nPtsPerStrand*2
nedgeb = 2*nSurfFaces
npoibf = 2*nSurfNodes
ncomps = 2 !just set to 2 for now

ALLOCATE(x(2,nverts))
ALLOCATE(cell(3,ncells))
ALLOCATE(edgb(3,nedgeb))
ALLOCATE(bpts(2,npoibf))
ALLOCATE(flag(0:nPtsPerStrand,nSurfNodes))

k = 0
DO n=1,nSurfNodes
   xn = xSurf(1,n)
   yn = xSurf(2,n)
   nx = pointingVec(1,n)
   ny = pointingVec(2,n)
DO j=0,nPtsPerStrand
   k = k+1
   flag(j,n) = k
   ds = xStrand(j)
   x(1,k) = xn+ds*nx
   x(2,k) = yn+ds*ny
END DO
END DO

k = 0
DO n=1,nSurfFaces
   n1 = surfFaces(1,n)
   n2 = surfFaces(2,n)
DO j=1,nPtsPerStrand
   p1 = flag(j-1,n1)
   p2 = flag(j-1,n2)
   p3 = flag(j  ,n1)
   p4 = flag(j  ,n2)
   k = k+1
   cell(1,k) = p1
   cell(2,k) = p2
   cell(3,k) = p4
   k = k+1
   cell(1,k) = p1
   cell(2,k) = p4
   cell(3,k) = p3
END DO
END DO

k = 0
DO n=1,nSurfFaces
   n1 = surfFaces(1,n)
   n2 = surfFaces(2,n)
   j = 0
   k = k+1
   edgb(1,k) = flag(j,n1)
   edgb(2,k) = flag(j,n2)
   edgb(3,k) = 1
   j = nPtsPerStrand
   k = k+1
   edgb(1,k) = flag(j,n1)
   edgb(2,k) = flag(j,n2)
   edgb(3,k) = 2
END DO

k = 0
DO n=1,nSurfNodes
   j = 0
   k = k+1
   bpts(1,k) = flag(j,n)
   bpts(2,k) = 1
   j = nPtsPerStrand
   k = k+1
   bpts(1,k) = flag(j,n)
   bpts(2,k) = 2
END DO


! write nc2d format to file
OPEN(inc2dunit,FILE=nc2dFile,STATUS='replace',FORM='formatted')
WRITE(inc2dunit,'(A1,I8)')'#',nverts
WRITE(inc2dunit,'(A1,I8)')'#',ncells
WRITE(inc2dunit,'(A1,I8)')'#',nedgeb
WRITE(inc2dunit,'(A1,I8)')'#',npoibf
WRITE(inc2dunit,'(A1,I8)')'#',ncomps
WRITE(inc2dunit,'(A20)')'VARIABLES = "X", "Y"'
WRITE(inc2dunit,'(A8,I8,A4,I8,A24)')'ZONE N= ',nverts,' E= ',ncells, &
                                    ' ,F=FEPOINT, ET=TRIANGLE'
WRITE(inc2dunit,*)
DO n=1,nverts
   WRITE(inc2dunit,*)x(1,n),x(2,n)
END DO
WRITE(inc2dunit,*)
DO n=1,ncells
   WRITE(inc2dunit,*)cell(1,n),cell(2,n),cell(3,n)
END DO
WRITE(inc2dunit,*)
DO n=1,nedgeb
   WRITE(inc2dunit,'(A1,3I8)')'#',edgb(1,n),edgb(2,n),edgb(3,n)
END DO
IF (npoibf > 0) THEN
   WRITE(inc2dunit,*)
   DO n=1,npoibf
      WRITE(inc2dunit,'(A1,2I8)')'#',bpts(1,n),bpts(2,n)
   END DO
END IF
CLOSE(inc2dunit)


! deallocate data


END PROGRAM strand2nc2d
