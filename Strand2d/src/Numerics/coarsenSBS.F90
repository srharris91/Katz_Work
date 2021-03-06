!> \brief
!! This subroutine performs multigrid coarsening in the preprocess step
!! to agglomerate fine grids into coarse grids.
!! \param mglevel Flag to determine if stand alone run.
!! \param mglevel Multigrid level.
!! \param nFacesT Number of faces on this coarse level.
!! \param nGfacesT Number of ghost faces on this coarse level.
!! \param nPstrT Number of cells per strand on this coarse level.
!! \param nBedgesT Number of boundary edges on this coarse level.
!! \param nPedgesT Number of partition edges on this coarse level.
!! \param nEdgesT Number of edges on this coarse level.
!! \param nFacesP Number of faces on parent level.
!! \param nGfacesP Number of ghost faces on parent level.
!! \param nEdgesP Number of edges on parent level.
!! \param nBedgesP Number of boundary edges on parent level.
!! \param nPedgesP Number of partition edges on parent level.
!! \param nPstrP Number of cells per strand on parent level.
!! \param nFringeP Number of fringe cells on parent level.
!! \param ndim Number of spatial dimensions.
!! \param edgeP Edge list on parent level.
!! \param bTagP Boundary tags on parent level.
!! \param fTagP Face tags on parent level.
!! \param fClipP Clipping index on parent level.
!! \param f2cc Fine to coarse array for cells.
!! \param f2ce Fine to coarse array for edges.
!! \param f2cs Fine to coarse array for f2cs.
!! \param fTag Face tags on this coarse level.
!! \param bTag Boundary tags on this coarse level.
!! \param fClip Clipping index on this coarse level.
!! \param edge Edge list on this coarse level.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-10-16
!! \par Further Documentation:
!! \par Source code:
!! \include src/Numerics/coarsenSBS.F90


SUBROUTINE coarsensbs( &
     standAlone, &
     mglevel, &
     nFacesT, &
     nGfacesT, &
     nPstrT, &
     nBedgesT, &
     nPedgesT, &
     nEdgesT, &
     nFacesP, &
     nGfacesP, &
     nEdgesP, &
     nBedgesP, &
     nPedgesP, &
     nPstrP, &
     nFringeP, &
     ndim, &
     edgeP, &
     bTagP, &
     fTagP, &
     fClipP, &
     f2cc, &
     f2ce, &
     f2cs, &
     fTag, &
     bTag, &
     fClip, &
     edge)
  
IMPLICIT NONE

INTEGER,INTENT(IN   ) :: &
     standAlone, &
     mglevel, &
     nFacesP, &
     nGfacesP, &
     nEdgesP, &
     nBedgesP, &
     nPedgesP, &
     nPstrP, &
     nFringeP, &
     ndim
INTEGER,INTENT(  OUT) :: &
     nFacesT, &
     nGfacesT, &
     nPstrT, &
     nBedgesT, &
     nPedgesT, &
     nEdgesT
INTEGER,INTENT(INOUT),DIMENSION(ndim,nEdgesP) :: edgeP
INTEGER,INTENT(INOUT),DIMENSION(nBedgesP) :: bTagP
INTEGER,INTENT(INOUT),DIMENSION(nFacesP) :: fTagP
INTEGER,INTENT(INOUT),DIMENSION(nFacesP+nBedgesP) :: fClipP
INTEGER,INTENT(  OUT),DIMENSION(nFacesP+nBedgesP) :: f2cc
INTEGER,INTENT(  OUT),DIMENSION(nEdgesP) :: f2ce
INTEGER,INTENT(  OUT),DIMENSION(0:nPstrP+1) :: f2cs
INTEGER,INTENT(  OUT),DIMENSION(nFacesP) :: fTag
INTEGER,INTENT(  OUT),DIMENSION(nBedgesP) :: bTag
INTEGER,INTENT(  OUT),DIMENSION(nFacesP+nBedgesP) :: fClip
INTEGER,INTENT(  OUT),DIMENSION(ndim,nEdgesP) :: edge
INTEGER :: n,i,j,k,m,ii,mm,c1,c2,ncscA1,cc,fc,ce,cb
INTEGER,ALLOCATABLE,DIMENSION(:) :: flag,cscA1,cscA2,cscA1t,cscA2t, &
                                    que,deg,esc1,esc2
LOGICAL :: still_have_vols,additional


! make connectivity 1-based
edgeP  = edgeP+1
bTagP  = bTagP+1
fTagP  = fTagP+1


! form cells surrounding cells
ALLOCATE(cscA2(nFacesP+1))
cscA2 = 0
DO n=1,nEdgesP-nBedgesP-nPedgesP
   c1          = edgeP(1,n)
   c2          = edgeP(2,n)
   cscA2(c1+1) = cscA2(c1+1)+1
   cscA2(c2+1) = cscA2(c2+1)+1
END DO

DO n=2,nFacesP+1
   cscA2(n) = cscA2(n)+cscA2(n-1)
END DO
ncscA1 = cscA2(nFacesP+1)
ALLOCATE(cscA1(ncscA1))

DO n=1,nEdgesP-nBedgesP-nPedgesP
   c1        = edgeP(1,n)
   c2        = edgeP(2,n)
   i         = cscA2(c1)+1
   cscA2(c1) = i
   cscA1(i)  = c2
   i         = cscA2(c2)+1
   cscA2(c2) = i
   cscA1(i)  = c1
END DO
DO n=nFacesP+1,2,-1
   cscA2(n) = cscA2(n-1)
END DO
cscA2(1) = 0


! perform agglomeration by advancing front
ALLOCATE(que(nFacesP))
ALLOCATE(deg(nFacesP))
ALLOCATE(flag(nFacesP))
cc      = 0
c1      = 1
c2      = 1
f2cc    = 0
deg     = 0
flag    = 0
que(1)  = 1
flag(1) = 1

advance: DO
   still_have_vols = .FALSE.
   next_vol: DO n=c1,c2 ! find the next volume in the front
      fc = que(n)
      IF (f2cc(fc) == 0) THEN
         c1 = n+1
         still_have_vols = .TRUE.
         EXIT next_vol
      END IF
   END DO next_vol

   IF (.NOT. still_have_vols) EXIT advance ! all done then exit

   cc       = cc+1
   f2cc(fc) = cc
   deg(fc)  = deg(fc)+1
   DO n=cscA2(fc)+1,cscA2(fc+1) ! add a coarse cell
      i = cscA1(n)
      IF (f2cc(i) == 0 .AND. fTagP(fc) == fTagP(i)) THEN
      !IF (f2cc(i) == 0) THEN
         f2cc(i) = cc
         deg(fc) = deg(fc)+1
      END IF
   END DO

   DO n=cscA2(fc)+1,cscA2(fc+1) ! set degree and advance front
      i = cscA1(n)
      IF (f2cc(i) == f2cc(fc)) THEN
         deg(i) = deg(fc)
         DO m=cscA2(i)+1,cscA2(i+1)
            k = cscA1(m)
            IF (flag(k) == 0) THEN
               c2      = c2+1
               que(c2) = k
               flag(k) = 1
            END IF
         END DO
      ELSE IF (f2cc(i) == 0 .AND. flag(i) == 0) THEN
         c2      = c2+1
         que(c2) = i
         flag(i) = 1
      END IF
   END DO

   IF (deg(fc) == 1) THEN ! fix singleton cells
      ii = 1000
      k  = fc
      DO n=cscA2(fc)+1,cscA2(fc+1)
         i = cscA1(n)
         IF (deg(i) < ii .AND. fTagP(fc) == fTagP(i)) THEN
         !IF (deg(i) < ii) THEN
            ii = deg(i)
            k  = i
         END IF
      END DO
      IF (k /= fc) THEN
         deg(k)   = deg(k)+1
         cc       = cc-1
         f2cc(fc) = f2cc(k)
      END IF
      DO n=cscA2(k)+1,cscA2(k+1)
         i = cscA1(n)
         IF (f2cc(i) == f2cc(k)) deg(i) = deg(k)
      END DO
   END IF

END DO advance

!do n=1,nFacesP+nBedgesP
!   write(*,*)mglevel,n,f2cc(n)
!END DO

DEALLOCATE(que)
DEALLOCATE(deg)
DEALLOCATE(flag)


! form edges surrounding the coarse cells
ALLOCATE(esc2(cc+1))
esc2 = 0
DO n=1,nEdgesP
   c1 = f2cc(edgeP(1,n))
   c2 = f2cc(edgeP(2,n))
   IF (c1 /= c2) THEN
      IF (c1 > 0) esc2(c1+1) = esc2(c1+1)+1
      IF (c2 > 0) esc2(c2+1) = esc2(c2+1)+1
   END IF
END DO

DO n=2,cc+1
   esc2(n) = esc2(n)+esc2(n-1)
END DO
ALLOCATE(esc1(esc2(cc+1)))

DO n=1,nEdgesP
   c1 = f2cc(edgeP(1,n))
   c2 = f2cc(edgeP(2,n))
   IF (c1 /= c2) THEN
   IF (c1 > 0) THEN
      i        = esc2(c1)+1
      esc2(c1) = i
      esc1(i)  = n
   END IF
   IF (c2 > 0) THEN
      i        = esc2(c2)+1
      esc2(c2) = i
      esc1(i)  = n
   END IF
   END IF
END DO
DO n=cc+1,2,-1
   esc2(n) = esc2(n-1)
END DO
esc2(1) = 0

!DO n=1,cc
!   write(*,'(I12,20I6)')n,(esc1(i),i=esc2(n)+1,esc2(n+1))
!END DO


! set f2cc at boundaries to the negative fTag value to avoid agglomerating
! different edge types. Fringe cells should remain f2cc = 0.
k = 0
DO n=nEdgesP-nBedgesP+1,nEdgesP
   k        = k+1
   c2       = edgeP(2,n)
   f2cc(c2) =-bTagP(k)
END DO


! add interior coarse level edges
f2ce = 0
ce   = 0
DO n=1,cc
DO m=esc2(n)+1,esc2(n+1)
   i = esc1(m)
   IF (f2ce(i) == 0) THEN
      c1         = f2cc(edgeP(1,i))
      c2         = f2cc(edgeP(2,i))
   IF (MIN(c1,c2) > 0) THEN
      ce         = ce+1
      edge(1,ce) = c1
      edge(2,ce) = c2
      f2ce(i)    = ce
   END IF
   END IF
END DO
END DO
k = ce


! add partition coarse level edges
DO n=1,cc
DO m=esc2(n)+1,esc2(n+1)
   i = esc1(m)
   IF (f2ce(i) == 0) THEN
      c1         = f2cc(edgeP(1,i))
      c2         = f2cc(edgeP(2,i))
   IF (MIN(c1,c2) == 0) THEN
      ce         = ce+1
      edge(1,ce) = c1
      edge(2,ce) = c2
      f2ce(i)    = ce
   END IF
   END IF
END DO
END DO
nPedgesT = ce-k
k        = ce


! add boundary coarse level edges
DO n=1,cc
DO m=esc2(n)+1,esc2(n+1)
   i = esc1(m)
   IF (f2ce(i) == 0) THEN
      c1         = f2cc(edgeP(1,i))
      c2         = f2cc(edgeP(2,i))
   IF (MIN(c1,c2) < 0) THEN
      ce         = ce+1
      edge(1,ce) = c1
      edge(2,ce) = c2
      f2ce(i)    = ce
   END IF
   END IF
END DO
END DO
nBedgesT = ce-k
nEdgesT  = ce

!DO n=1,ce
!   write(*,'(I12,10I6)')n,edge(:,n),ce,nedgesP-igeP
!END DO

!DO n=1,nEdgesP-igeP
!   write(*,*)n,f2ce(n)
!END DO

DEALLOCATE(esc1)
DEALLOCATE(esc2)


! fix boundary and fringe values
nGfacesT = nPedgesT
nFacesT  = cc+nGfacesT
DO n=nEdgesT-nPedgesT-nBedgesT+1,nEdgesT
   cc        = cc+1
   edge(2,n) = cc
END DO


! set partition edge agglomeration array
DO n=nEdgesP-nPedgesP-nBedgesP+1,nEdgesP-nBedgesP
   c2         = edgeP(2,n)
   i          = f2ce(n)
   f2cc(c2)   = edge(2,i)
END DO


! set boundary edges agglomeration array
k = 0
ii = nEdgesT-nBedgesT
DO n=nEdgesP-nBedgesP+1,nEdgesP
   k          = k+1
   c2         = edgeP(2,n)
   i          = f2ce(n)
   f2cc(c2)   = edge(2,i)
   bTag(i-ii) = bTagP(k)
END DO


! set component face tags
fTag = 0
DO n=1,nFacesP
   i       = f2cc(n)
   j       = fTag(i)
   k       = fTagP(n)
   IF (j /= 0 .AND. j /= k) THEN
      WRITE(*,*)'*** agglomeration across surface types in gmgGEOMCoarsen***'
      STOP
   END IF
   fTag(i) = k
END DO

!do n=1,nFacesP+nBedgesP
!   write(*,*)n,f2cc(n),cc
!END DO
!stop

!DO n=1,nEdgesT
!   write(*,'(I12,10I6)')n,edge(:,n)
!END DO
!stop

!DO n=1,nEdgesP
!   write(*,*)n,f2ce(n)
!END DO
!stop

!do n=1,nFacesT
!   write(*,*)n,fTag(n)
!end do
!stop

!do n=1,nBedgesT
!   write(*,*)n,bTag(n)
!end do
!stop

! agglomerate the strand direction

nPstrT =(nPstrP-nFringeP)/2+1
f2cs   = 0
IF (standAlone == 1) THEN
   i = 0
   DO n=1,nPstrT
      i       = i+1
      f2cs(i) = n
      i       = i+1
      f2cs(i) = n
   END DO
   IF (nPstrP /= i) THEN !3-cell agglom. at strand end
      i       = i+1
      f2cs(i) = nPstrT
   END IF
ELSE
   IF (nPstrT == 1) THEN
      nPstrT  = 2
      f2cs(1) = 1
      f2cs(2) = 2
   ELSE
      i = 0
      DO n=1,nPstrT-1
         i       = i+1
         f2cs(i) = n
         i       = i+1
         f2cs(i) = n
      END DO
      IF (nPstrP-nFringeP /= i) THEN !3-cell agglom. at strand end
         i       = i+1
         f2cs(i) = nPstrT-1
      END IF
      f2cs(i+1:nPstrP) = nPstrT
   END IF
END IF


! set clipping on coarse level. Note: this is a conservative clip
IF (standAlone == 1) THEN
   fClip = nPstrT
ELSE
   fClip = nPstrT-1
   DO n=1,nFacesP
      i = f2cc(n)
      j = fClipP(n)
      IF (j > 0) THEN
         k = MAX(0,f2cs(j)-1)
         IF (f2cs(j+1) /= f2cs(j)) k = f2cs(j)
      ELSE
         k = 0
      END IF
      fClip(i) = MIN(fClip(i),k)
   END DO
   DO n=nEdgesT-nBedgesT+1,nEdgesT
      c1        = edge(1,n)
      c2        = edge(2,n)
      fClip(c2) = fClip(c1)
   END DO
END IF


!WRITE(*,'(A40,I6)'  )'mesh level: ',mglevel
!DO n=1,nFacesP
!   write(*,*)n,fClipP(n)
!END DO
!DO n=1,nFacesT+nBedgesT
!   write(*,*)n,fclip(n)
!END DO

!do n=1,nPstrP
!   write(*,*)n,f2cs(n)
!end do
!stop

! print out coarse block information
WRITE(*,*)
WRITE(*,'(A40,I6)'  )'mesh level: ',mglevel
WRITE(*,'(A40,I8,A2,F4.2,A1)')'number interior faces (f/c ratio): ', &
                      nFacesT-nGfacesT,' (', &
                      REAL(nFacesP-nGfacesP)/REAL(nFacesT-nGfacesT),')'
WRITE(*,'(A40,I8,A2,F4.2,A1)')'number of edges (f/c ratio): ', &
                      nEdgesT,' (', &
                      REAL(nEdgesP)/REAL(MAX(nEdgesT,1)),')'
WRITE(*,'(A40,I8,A2,F4.2,A1)')'number of boundary edges (f/c ratio): ', &
                      nBedgesT,' (', &
                      REAL(nBedgesP)/REAL(MAX(nBedgesT,1)),')'
WRITE(*,'(A40,I8,A2,F4.2,A1)')'number of cells per strand (f/c ratio): ', &
                      nPstrT,' (', &
                      REAL(nPstrP)/REAL(nPstrT),')'
WRITE(*,*)


! make connectivity 0-based
edgeP = edgeP-1
bTagP = bTagP-1
fTagP = fTagP-1
f2cc  = f2cc-1
!f2ce  = f2ce-1 !keep 1-based to allow for meaningful signs
edge  = edge-1
bTag  = bTag-1
fTag  = fTag-1


END SUBROUTINE coarsensbs
