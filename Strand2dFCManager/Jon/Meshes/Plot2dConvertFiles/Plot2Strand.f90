SUBROUTINE ReadPlot2Strand(inputfile)
!USE StrandSolv
USE Splines_3
USE Airfoil
IMPLICIT NONE

!/ Inputs /
CHARACTER(20)::Inputfile

!/ Subroutine Values /
REAL(8), DIMENSION(:,:), ALLOCATABLE :: Bounds,BNormal
INTEGER, DIMENSION(:),ALLOCATABLE :: Boundc,BNode
INTEGER:: NBound, Method, NSideBound, MeshOrder, FType
CHARACTER(40):: FName, PlotName
NAMELIST/Manager/FType,FName, PlotName
NAMELIST/Plot2Str/NBound, Method, NSideBound, MeshOrder
NAMELIST/Bound/Bounds,BNode,BNormal

! / Values from Read NASA /
! / Reading Varibles /
INTEGER::nbl, n, i, j
INTEGER, DIMENSION(:),ALLOCATABLE:: idim, jdim
INTEGER,PARAMETER::iFile=20,oFile=21
REAL(8),DIMENSION(:,:,:),ALLOCATABLE:: x,y

! / Output Varibles /
INTEGER::TotalCells, NNodes, NInterior, Nodenum, SWNode, NWNode, SENode, NENode
INTEGER, DIMENSION(:,:), ALLOCATABLE:: Cells

! / Spline Varibles / !
REAL(8),ALLOCATABLE,DIMENSION(:)::Xpts,Ypts
INTEGER::ptsEnd,ii
REAL(8),ALLOCATABLE,DIMENSION(:,:)::xyArray
!/ END HEADER /!

OPEN(UNIT=13, FILE=TRIM(ADJUSTL(Inputfile)),STATUS="OLD",ACTION="READ")
WRITE(*,*)"FILE OPENED"
READ(13,Manager)
READ(13,Plot2Str)
ALLOCATE(Bounds(2,0:NBound-NSideBound),BNode(1:NSideBound),BNormal(1:NSideBound,1:2))
READ(13,Bound)
CLOSE(13)

WRITE(*,*)"Inputfile Read Successfull..."

! / Read NASA /

WRITE(*,*) "Read input file: ", TRIM(ADJUSTL(PlotName))
IF (FType==1) THEN
	OPEN(UNIT=iFile, FILE=Trim(ADJUSTL(PlotName)), ACTION='read',STATUS='Old')
	READ(iFile,*)nbl
	WRITE(*,*) "nbl: ", nbl
	ALLOCATE(idim(1:nbl),jdim(1:nbl))
	READ(iFile,*)(idim(n),jdim(n),n=1,nbl)
	WRITE(*,*) "idim :", idim(1)
	WRITE(*,*) "jdim :", jdim(1)
	ALLOCATE(x(1:MAXVAL(idim),1:MAXVAL(jdim),1:nbl),y(1:MAXVAL(idim),1:MAXVAL(jdim),1:nbl))

	DO n=1,nbl

	READ(iFile,*)((x(i,j,n),i=1,idim(n)),j=1,jdim(n)),&
					((y(i,j,n),i=1,idim(n)),j=1,jdim(n))

	END DO
	CLOSE(iFile)
ELSEIF (FType==2) THEN
	CALL ReadData_UIUC(PlotName,xyArray,nbl)
	ALLOCATE(x(nbl,1,1),y(nbl,1,1))
	x(:,1,1)=xyArray(:,1)
	y(:,1,1)=xyArray(:,2)
	ALLOCATE(idim(1:nbl),jdim(1:nbl))
	idim = 0.
	idim(1)=nbl
ELSE
WRITE(*,*)'Missunderstood Type'
END IF

WRITE(*,*)"READ COMPLETED"
WRITE(*,*)"PROCESSING...."

TotalCells = (MAXVAL(idim)-1)
WRITE(*,*) "Number of Cells:", TotalCells
NNodes = (MAXVAL(idim))
WRITE(*,*) "Number of Nodes:", NNodes

ALLOCATE(BoundC(0:idim(1)-1))
!/Set up Boundaries/
IF (Method==1) THEN
	DO i=0,idim(1)-1
		DO j=0,NBound-1
			IF(x(i+1,1,1)<=Bounds(2,j) .AND. x(i+1,1,1)>=Bounds(1,j)) THEN
				BoundC(i)=j
				EXIT
			END IF
		END DO
	END DO
END IF

IF (Method==0) THEN
	DO i=0,idim(1)-1
		DO j=0,NBound-2
			IF(i<=Bounds(2,j) .AND. i>=Bounds(1,j)) THEN
				BoundC(i)=j
				EXIT
			END IF
		END DO
	END DO
END IF

OPEN(UNIT=oFile, FILE = TRIM(ADJUSTL(FName)), Action="WRITE") 

! Conversion

IF (MeshOrder==1) THEN
	ptsEnd=0
	ALLOCATE(Xpts(ptsEnd),Ypts(ptsEnd))
	WRITE(oFile,*) TotalCells, NNodes+ptsEnd, NBound, NSideBound
	DO i=0,NNodes-2
	
		WRITE(oFile,*) 1,i,i+1,BoundC(i)

	END DO
ELSEIF(MeshOrder==2) THEN
	ii=0
	ptsEnd=1*NNodes
	ALLOCATE(Xpts(ptsEnd),Ypts(ptsEnd))
	WRITE(oFile,*) TotalCells, NNodes+ptsEnd-2, NBound, NSideBound
	WRITE(*,*)X(1,1,1)
	CALL NewPts(NNodes,X(:,1,1),Y(:,1,1), 1, Xpts, Ypts)
	DO i=0,NNodes-2
		
		WRITE(oFile,*) 2,i,i+1,NNodes+ii, BoundC(i)
		ii=ii+1
	END DO
ELSEIF(MeshOrder==3) THEN
	ii=0
	ptsEnd=2*NNodes
	ALLOCATE(Xpts(ptsEnd),Ypts(ptsEnd))
	WRITE(oFile,*) TotalCells, NNodes+ptsEnd-2, NBound, NSideBound
	WRITE(*,*)X(1,1,1)
	CALL NewPts(NNodes,X(:,1,1),Y(:,1,1), 2, Xpts, Ypts)
	DO i=0,NNodes-2
		
		WRITE(oFile,*) 3,i,i+1,NNodes+ii, NNodes+ii+1, BoundC(i)
		ii=ii+2
	END DO
ELSEIF(MeshOrder==4) THEN
	ii=0
	ptsEnd=3*NNodes
	ALLOCATE(Xpts(ptsEnd),Ypts(ptsEnd))
	WRITE(oFile,*) TotalCells, NNodes+ptsEnd-2, NBound, NSideBound
	WRITE(*,*)X(1,1,1)
	CALL NewPts(NNodes,X(:,1,1),Y(:,1,1), 3, Xpts, Ypts)
	DO i=0,NNodes-2
		
		WRITE(oFile,*) 4,i,i+1,NNodes+ii, NNodes+ii+1, NNodes+ii+2, BoundC(i)
		ii=ii+3
	END DO
ELSE
	Write(*,*)"Order Not Implemented!"
END IF

DO i=1,NNodes

	WRITE(oFile,*) X(i,1,1), Y(i,1,1)

END DO

DO i=1,ptsEnd-2

	WRITE(oFile,*) Xpts(i), Ypts(i)

END DO

BNode(NSideBound)=NNodes-1
DO i=1,NSideBound

	WRITE(oFile,*) BNode(i), NBound-((NSideBound-i)), BNormal(i,:)
	
END DO

CLOSE(oFile)
DEALLOCATE(idim,jdim,x,y,Bounds,Boundc,BNode,BNormal,Xpts,Ypts)
END SUBROUTINE ReadPlot2Strand
