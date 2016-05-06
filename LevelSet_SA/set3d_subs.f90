MODULE set3d_SUBs
    USE set3d_mod
    USE Lib_VTK_IO
    IMPLICIT NONE
CONTAINS

    !*************************************************************************************!
    ! Input.namelist readfile
    !*************************************************************************************!
    SUBROUTINE InputRead(&
            rank    ,&
            filename,&
            meshfile,&
            dx      ,&
            CFL     ,&
            dd      ,&
            tol     ,&
            h       ,&
            nPass   ,&
            iter    ,&
            order1  ,&
            order2  ,&
            orderUp ,&
            convergenceLimit,&
            solutionType&
            )
        INTEGER ,INTENT(IN)     :: rank
        CHARACTER(*),INTENT(IN) :: filename
        CHARACTER(*),INTENT(OUT):: meshfile

        REAL                    ,INTENT(OUT)    :: dx,tol,CFL
        REAL    ,DIMENSION(3)   ,INTENT(OUT)    :: h,convergenceLimit
        INTEGER                 ,INTENT(OUT)    :: nPass,dd,order2,solutionType
        INTEGER ,DIMENSION(3)   ,INTENT(OUT)    :: iter,order1,orderUp
        INTEGER                 ,PARAMETER      :: iUnit=9

        NAMELIST/LEVELSET_SA/&
            meshfile,&
            dx      ,&
            CFL     ,&
            dd      ,&
            tol     ,&
            h       ,&
            nPass   ,&
            iter    ,&
            order1  ,&
            order2  ,&
            orderUp ,&
            convergenceLimit ,&
            solutionType

        OPEN(iUnit,FILE=TRIM(filename),STATUS='old')
        READ(iUnit,LEVELSET_SA)
        CLOSE(iUnit)

        IF (rank == 0) THEN
            WRITE(*,*)
            WRITE(*,LEVELSET_SA)
            WRITE(*,*)
        END IF

    END SUBROUTINE InputRead

    !*************************************************************************************!
    ! Output to VTS file using the following subroutine
    !*************************************************************************************!
    SUBROUTINE Pass(D_s,nbrs,SIZE_X,dims,req,stat,mpi_set_smaller,GRID_COMM,error,nPass)

    TYPE (SET_smaller),DIMENSION(-nPass:,-nPass:,-nPass:),INTENT(INOUT) :: D_s
    INTEGER,DIMENSION(:),INTENT(IN) :: nbrs,SIZE_X,dims
    INTEGER,INTENT(IN) :: nPass
    INTEGER,DIMENSION(0:,:),INTENT(INOUT) :: req
    INTEGER,DIMENSION(:,0:,:),INTENT(INOUT) :: stat
    INTEGER,INTENT(IN) :: GRID_COMM,error
    INTEGER,INTENT(IN) :: mpi_set_smaller 
    INTEGER :: nn,nnn
    INTEGER :: i,j,k
    INTEGER,PARAMETER :: DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4, BACK = 5, FRONT = 6       !Positive x=DOWN, y=RIGHT, z=BACK

    CALL MPI_Barrier(GRID_COMM,error)

    DO nn = 1, nPass
    IF (dims(2) > 1) THEN
      DO i = -1, SIZE_X(1)
        DO k = -1, SIZE_X(3) 
          CALL MPI_Isend(D_s(i,nn-1,            k),1,mpi_set_smaller,nbrs(UP)  ,0,GRID_COMM,req(0,nn),error)
          CALL MPI_Isend(D_s(i,SIZE_X(2)-(nn-1),k),1,mpi_set_smaller,nbrs(DOWN),1,GRID_COMM,req(1,nn),error)
          
          CALL MPI_Irecv(D_s(i,SIZE_X(2)+nn,    k),1,mpi_set_smaller,nbrs(DOWN),0,GRID_COMM,req(2,nn),error)
          CALL MPI_Irecv(D_s(i,-nn,             k),1,mpi_set_smaller,nbrs(UP)  ,1,GRID_COMM,req(3,nn),error)

          CALL MPI_Barrier(GRID_COMM,error)
          DO nnn = 0,3
            CALL MPI_WAIT(req(nnn,nn),stat(:,nnn,nn), error)
          END DO

        END DO
      END DO
    END IF

    IF (dims(1) > 1) THEN
      DO j = -1, SIZE_X(2)
        DO k = -1, SIZE_X(3) 
          CALL MPI_Isend(D_s(nn-1,            j,k),1,mpi_set_smaller,nbrs(LEFT) ,0,GRID_COMM,req(0,nn),error)
          CALL MPI_Isend(D_s(SIZE_X(1)-(nn-1),j,k),1,mpi_set_smaller,nbrs(RIGHT),1,GRID_COMM,req(1,nn),error)
          
          CALL MPI_Irecv(D_s(SIZE_X(1)+nn,    j,k),1,mpi_set_smaller,nbrs(RIGHT),0,GRID_COMM,req(2,nn),error)
          CALL MPI_Irecv(D_s(-nn,             j,k),1,mpi_set_smaller,nbrs(LEFT) ,1,GRID_COMM,req(3,nn),error)

          CALL MPI_Barrier(GRID_COMM,error)
          DO nnn = 0,3
            CALL MPI_WAIT(req(nnn,nn),stat(:,nnn,nn), error)
          END DO

        END DO
      END DO
    END IF

    IF (dims(3) > 1) THEN
      DO i = -1, SIZE_X(1)
        DO j = -1, SIZE_X(2) 
          CALL MPI_Isend(D_s(i,j,nn-1            ),1,mpi_set_smaller,nbrs(FRONT),0,GRID_COMM,req(0,nn),error)
          CALL MPI_Isend(D_s(i,j,SIZE_X(3)-(nn-1)),1,mpi_set_smaller,nbrs(BACK) ,1,GRID_COMM,req(1,nn),error)
          
          CALL MPI_Irecv(D_s(i,j,SIZE_X(3)+nn    ),1,mpi_set_smaller,nbrs(BACK) ,0,GRID_COMM,req(2,nn),error)
          CALL MPI_Irecv(D_s(i,j,-nn             ),1,mpi_set_smaller,nbrs(FRONT),1,GRID_COMM,req(3,nn),error)

          CALL MPI_Barrier(GRID_COMM,error)
          DO nnn = 0,3
            CALL MPI_WAIT(req(nnn,nn),stat(:,nnn,nn), error)
             END DO

        END DO
      END DO
    END IF


    END DO

!WRITE(*,*) "Finished Pass"

    END SUBROUTINE Pass
    !*************************************************************************************!
    ! Output to VTS file using the following subroutine
    !*************************************************************************************!
    SUBROUTINE SET_smaller_to_SET(S,sm,nPass)
        INTEGER,INTENT(IN) :: nPass
        TYPE (SET),DIMENSION(-1:,-1:,-1:),INTENT(INOUT) :: S
        TYPE (SET_smaller),DIMENSION(-nPass:,-nPass:,-nPass:),INTENT(IN) :: sm
        INTEGER :: nx,ny,nz
        nx = size(sm(:,0,0))
        ny = size(sm(0,:,0))
        nz = size(sm(0,0,:))


        ! nx-5 to not output the ghost nodes... and only output from -1 to end
        S%phi       = sm(-1:nx-1-nPass*2,-1:ny-1-nPass*2,-1:nz-1-nPass*2)%phi !copy small data type to larger type
        S%gradPhiX  = sm(-1:nx-1-nPass*2,-1:ny-1-nPass*2,-1:nz-1-nPass*2)%gradPhiX !copy small data type to larger type
        S%gradPhiY  = sm(-1:nx-1-nPass*2,-1:ny-1-nPass*2,-1:nz-1-nPass*2)%gradPhiY !copy small data type to larger type
        S%gradPhiZ  = sm(-1:nx-1-nPass*2,-1:ny-1-nPass*2,-1:nz-1-nPass*2)%gradPhiZ !copy small data type to larger type

    END SUBROUTINE SET_smaller_to_SET
    SUBROUTINE output_vtk(nx,ny,nz,Dat,output_filename,rank,coords)
        !INTEGER,DIMENSION(:),INTENT(IN)::dims
        INTEGER,DIMENSION(:),INTENT(IN):: coords
        CHARACTER(LEN=*),INTENT(IN)  ::  output_filename
        INTEGER,INTENT(IN)  ::  nx,ny,nz,rank
        TYPE (SET),DIMENSION(-1:,-1:,-1:),INTENT(IN) :: Dat
        CHARACTER(LEN=1024) ::  output_filename_rank,output_folder,cmd
        INTEGER :: E_IO
        INTEGER :: nnx1,nnx2,nny1,nny2,nnz1,nnz2

        output_folder = 'output/' // trim(output_filename) // '/'
        cmd = 'mkdir -p ' // trim(output_folder)
        CALL SYSTEM(trim(cmd))
        
        
        IF (rank < 10) THEN
            WRITE(output_filename_rank,'(A,A,I1,A4)') trim(output_folder),trim(output_filename),rank,'.vts'
        ELSE IF (rank >=10 .AND. rank<100) THEN
            WRITE(output_filename_rank,'(A,A,I2,A4)') trim(output_folder),trim(output_filename),rank,'.vts'
        ELSE IF (rank >=100 .AND. rank<1000) THEN
            WRITE(output_filename_rank,'(A,A,I3,A4)') trim(output_folder),trim(output_filename),rank,'.vts'
        ELSE 
            WRITE(*,*) 'error in putting output_filename with this rank',rank
        END IF
        ! output to vtk files
            nnx1 = (coords(1))*nx + coords(1)-1
            nnx2 = (coords(1)+1)*nx + coords(1)

            nny1 = (coords(2))*ny + coords(2) - 1
            nny2 = (coords(2)+1)*ny + coords(2)

            nnz1 = (coords(3))*nz + coords(3) - 1
            nnz2 = (coords(3)+1)*nz + coords(3)
        WRITE(*,*) 'Writing file ',trim(output_filename_rank)
        E_IO = VTK_INI_XML_write(fformat='binary', filename=trim(output_filename_rank),mesh_topology='StructuredGrid',&
            nx1=nnx1        ,&
            nx2=nnx2        ,&
            ny1=nny1        ,&
            ny2=nny2        ,&
            nz1=nnz1        ,&
            nz2=nnz2        )
        E_IO = VTK_FLD_XML(fld_action='open')
        E_IO = VTK_FLD_XML(fld=0.e1,fname='TIME')
        E_IO = VTK_FLD_XML(fld=1,fname='CYCLE')
        E_IO = VTK_FLD_XML(fld_action='close')
        E_IO = VTK_GEO_XML_WRITE(& 
            nx1=nnx1        ,&
            nx2=nnx2        ,&
            ny1=nny1        ,&
            ny2=nny2        ,&
            nz1=nnz1        ,&
            nz2=nnz2        ,&
            NN=INT((nx+2)*(ny+2)*(nz+2)),&
            X=Dat%x         ,&
            Y=Dat%y         ,&
            Z=Dat%z)
        E_IO = VTK_DAT_XML(var_location='node',var_block_action='open')
        E_IO = VTK_VAR_XML(             &
            NC_NN=(nx+2)*(ny+2)*(nz+2) ,&
            varname='Phi'              ,&
            var=Dat%phi)
        E_IO = VTK_VAR_XML(             &
            NC_NN=(nx+2)*(ny+2)*(nz+2) ,&
            varname='gradPhiX'         ,&
            var=Dat%gradPhiX)
        E_IO = VTK_VAR_XML(             &
            NC_NN=(nx+2)*(ny+2)*(nz+2) ,&
            varname='gradPhiY'         ,&
            var=Dat%gradPhiY)
        E_IO = VTK_VAR_XML(             &
            NC_NN=(nx+2)*(ny+2)*(nz+2) ,&
            varname='gradPhiZ'         ,&
            var=Dat%gradPhiZ)
        E_IO = VTK_DAT_XML(var_location='node',var_block_action='close')
        E_IO = VTK_GEO_XML_WRITE()
        E_IO = VTK_END_XML()
    END SUBROUTINE output_vtk

    SUBROUTINE output_pvtk(nx,ny,nz,output_filename,input_filename,rank,dims,coords,num_proc)
        INTEGER,DIMENSION(:),INTENT(IN)::dims
        INTEGER,DIMENSION(:,:),INTENT(IN):: coords
        INTEGER         ,INTENT(IN)  ::  num_proc
        CHARACTER(LEN=*),INTENT(IN)  ::  output_filename
        CHARACTER(LEN=*),INTENT(IN)  ::  input_filename
        CHARACTER(LEN=1024),DIMENSION(num_proc)::  input_filenames
        CHARACTER(LEN=1024) :: output_file
        INTEGER,INTENT(IN)  ::  nx,ny,nz,rank
        INTEGER ::  i
        !TYPE (SET),DIMENSION(-1:,-1:,-1:),INTENT(IN) :: Dat
        !CHARACTER(LEN=1024) ::  output_filename_rank
        INTEGER :: E_IO
        INTEGER :: nnx1,nnx2,nny1,nny2,nnz1,nnz2
        ! output to vtk files
        output_file= 'output/' // trim(input_filename) // '/' // trim(output_filename)

        IF (rank == 0) THEN
                DO i=1,num_proc
                    IF (i <= 10) THEN
                        WRITE(input_filenames(i),'(A,I1,A4)') &
                            trim(input_filename),i-1,'.vts'
                    ELSE IF (i >10 .AND. rank<=100) THEN
                        WRITE(input_filenames(i),'(A,I2,A4)') &
                            trim(input_filename),i-1,'.vts'
                    ELSE IF (i >100 .AND. rank<=1000) THEN
                        WRITE(input_filenames(i),'(A,I3,A4)') &
                            trim(input_filename),i-1,'.vts'
                    ELSE 
                        WRITE(*,*) 'error in putting input_filenames in output_pvtk with this rank',rank
                    END IF
                END DO
                WRITE(*,*) 'Writing file ',trim(output_file)
                E_IO = PVTK_INI_XML(filename = trim(output_file),mesh_topology = 'PStructuredGrid',&
                    nx1=-1,             &
                    nx2=nx*dims(1)+1,   &
                    ny1=-1,             &
                    ny2=ny*dims(2)+1,   &
                    nz1=-1,             &
                    nz2=nz*dims(3)+1,   &
                    tp='Float64')
                  DO i = 1, num_proc
                      nnx1 = (coords(i,1))*nx + coords(i,1)-1
                      nnx2 = (coords(i,1)+1)*nx + coords(i,1)

                      nny1 = (coords(i,2))*ny + coords(i,2) - 1
                      nny2 = (coords(i,2)+1)*ny + coords(i,2)

                      nnz1 = (coords(i,3))*nz + coords(i,3) - 1
                      nnz2 = (coords(i,3)+1)*nz + coords(i,3)
                    E_IO = PVTK_GEO_XML(&
                        nx1 = nnx1,     &
                        nx2 = nnx2,     &
                        ny1 = nny1,     &
                        ny2 = nny2,     &
                        nz1 = nnz1,     &
                        nz2 = nnz2,     &
                        source=trim(input_filenames(i)))
                  END DO
                E_IO = PVTK_DAT_XML(var_location='node',var_block_action='open')
                E_IO = PVTK_VAR_XML(varname='Phi',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradPhiX',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradPhiY',tp='Float64')
                E_IO = PVTK_VAR_XML(varname='gradPhiZ',tp='Float64')
                E_IO = PVTK_DAT_XML(var_location='node',var_block_action='close')
                E_IO = PVTK_END_XML()
        END IF

    END SUBROUTINE output_pvtk

    SUBROUTINE output(nx,ny,nz,Dat,output_filename,input_filename,rank,dims,coords,coords_all,num_proc)
        INTEGER,DIMENSION(:),INTENT(IN)::dims,coords
        INTEGER,DIMENSION(:,:),INTENT(IN):: coords_all
        INTEGER         ,INTENT(IN)  ::  num_proc
        CHARACTER(LEN=*),INTENT(IN)  ::  output_filename
        CHARACTER(LEN=*),INTENT(IN)  ::  input_filename
        INTEGER,INTENT(IN)  ::  nx,ny,nz,rank
        TYPE (SET),DIMENSION(-1:,-1:,-1:),INTENT(IN) :: Dat


        CALL output_vtk(nx,ny,nz,Dat,input_filename,rank,coords)
        CALL output_pvtk(nx,ny,nz,output_filename,input_filename,rank,dims,coords_all,num_proc)

    END SUBROUTINE output

    !*************************************************************************************!
    ! Calculate upwind derivatives
    !*************************************************************************************!

    SUBROUTINE upwind(gM,i,j,k,nx,ny,nz,dx,phi,dPlus,dMinus,order,nPass)

        IMPLICIT NONE

        INTEGER,INTENT(IN) :: nPass
        INTEGER,INTENT(IN) :: nx,ny,nz,i,j,k,order
        REAL,INTENT(IN) :: dx
        !REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
        !TYPE (SET_smaller), DIMENSION(-2:nx+2,-2:ny+2,-2:nz+2),INTENT(IN) :: D_s
        REAL, DIMENSION(-nPass:nx+nPass,-nPass:ny+nPass,-nPass:nz+nPass),INTENT(IN) :: phi
        REAL,INTENT(OUT):: gM,dPlus,dMinus
        REAL :: gradX,gradY,gradZ,pa,pb,pc,pd,pe,pf,na,nb,nc,nd,ne,nf,a,b,c,d,e,f
        REAL :: u
        REAL :: uxp1,uxp2,uxp3,uxm1,uxm2,uxm3
        REAL :: uyp1,uyp2,uyp3,uym1,uym2,uym3
        REAL :: uzp1,uzp2,uzp3,uzm1,uzm2,uzm3
        INTEGER :: im,jm,km
        INTEGER :: ip,jp,kp
        INTEGER :: ip1,jp1,kp1
        INTEGER :: im1,jm1,km1
        INTEGER :: ip2,jp2,kp2
        INTEGER :: im2,jm2,km2
        INTEGER :: ip3,jp3,kp3
        INTEGER :: im3,jm3,km3

        ! calculate upwind derivatives

        IF (order == 1) THEN

            ! set plus and minus integers
            im = i-1
            jm = j-1
            km = k-1
            ip = i+1
            jp = j+1
            kp = k+1


            ! calculate derivatives
            a = (phi(i ,j ,k ) - phi(im,j ,k ))/dx
            b = (phi(ip,j ,k ) - phi(i ,j ,k ))/dx
            c = (phi(i ,j ,k ) - phi(i ,jm,k ))/dx
            d = (phi(i ,jp,k ) - phi(i ,j ,k ))/dx
            e = (phi(i ,j ,k ) - phi(i ,j ,km))/dx
            f = (phi(i ,j ,kp) - phi(i ,j ,k ))/dx

        ELSEIF (order == 3) THEN


            im1 = i-1
            jm1 = j-1
            km1 = k-1
            im2 = i-2
            jm2 = j-2
            km2 = k-2
            ip1 = i+1
            jp1 = j+1
            kp1 = k+1
            ip2 = i+2
            jp2 = j+2
            kp2 = k+2


            a = (2.*phi(ip1,j,  k  )+3.*phi(i,j,k)-6.*phi(im1,j  ,k  )+phi(im2,j  ,k  ))
            b =-(2.*phi(im1,j,  k  )+3.*phi(i,j,k)-6.*phi(ip1,j  ,k  )+phi(ip2,j  ,k  ))
            c = (2.*phi(i  ,jp1,k  )+3.*phi(i,j,k)-6.*phi(i  ,jm1,k  )+phi(i  ,jm2,k  ))
            d =-(2.*phi(i  ,jm1,k  )+3.*phi(i,j,k)-6.*phi(i  ,jp1,k  )+phi(i  ,jp2,k  ))
            e = (2.*phi(i  ,j  ,kp1)+3.*phi(i,j,k)-6.*phi(i  ,j  ,km1)+phi(i  ,j  ,km2))
            f =-(2.*phi(i  ,j  ,km1)+3.*phi(i,j,k)-6.*phi(i  ,j  ,kp1)+phi(i  ,j  ,kp2))

            a = a/(6.*dx)
            b = b/(6.*dx)
            c = c/(6.*dx)
            d = d/(6.*dx)
            e = e/(6.*dx)
            f = f/(6.*dx)



        ELSEIF (order == 5) THEN

            PRINT*, " Working on this ... "
            STOP

            im1 = i-1
            jm1 = j-1
            km1 = k-1
            im2 = i-2
            jm2 = j-2
            km2 = k-2
            im3 = i-3
            jm3 = j-3
            km3 = k-3

            ip1 = i+1
            jp1 = j+1
            kp1 = k+1
            ip2 = i+2
            jp2 = j+2
            kp2 = k+2
            ip3 = i+3
            jp3 = j+3
            kp3 = k+3

            u = phi(i,j,k)

            uxm1 = phi(im1,j,k)
            uxp1 = phi(ip1,j,k)
            uxm2 = phi(im2,j,k)
            uxp2 = phi(ip2,j,k)
            uxm3 = phi(im3,j,k)
            uxp3 = phi(ip3,j,k)

            uym1 = phi(i,jm1,k)
            uyp1 = phi(i,jp1,k)
            uym2 = phi(i,jm2,k)
            uyp2 = phi(i,jp2,k)
            uym3 = phi(i,jm3,k)
            uyp3 = phi(i,jp3,k)

            uzm1 = phi(i,j,km1)
            uzp1 = phi(i,j,kp1)
            uzm2 = phi(i,j,km2)
            uzp2 = phi(i,j,kp2)
            uzm3 = phi(i,j,km3)
            uzp3 = phi(i,j,kp3)

            ! calculate derivatives
            a =-(4.*uxm3 - 30.*uxm2 + 120.*uxm1 - 40.*u - 60.*uxp1 + 6.*uxp2)
            b = (4.*uxp3 - 30.*uxp2 + 120.*uxp1 - 40.*u - 60.*uxm1 + 6.*uxm2)
            c =-(4.*uym3 - 30.*uym2 + 120.*uym1 - 40.*u - 60.*uyp1 + 6.*uyp2)
            d = (4.*uyp3 - 30.*uyp2 + 120.*uyp1 - 40.*u - 60.*uym1 + 6.*uym2)
            e =-(4.*uzm3 - 30.*uzm2 + 120.*uzm1 - 40.*u - 60.*uzp1 + 6.*uzp2)
            f = (4.*uzp3 - 30.*uzp2 + 120.*uzp1 - 40.*u - 60.*uzm1 + 6.*uzm2)

            a = a/(120.*dx)
            b = b/(120.*dx)
            c = c/(120.*dx)
            d = d/(120.*dx)
            e = e/(120.*dx)
            f = f/(120.*dx)

        ELSE
            PRINT*," This order not set yet in upwind"
            STOP
        END IF


        ! check positive   
        pa = max(a,0.)
        pb = max(b,0.)
        pc = max(c,0.)
        pd = max(d,0.)
        pe = max(e,0.)
        pf = max(f,0.)

        ! check negative
        na = min(a,0.)
        nb = min(b,0.)
        nc = min(c,0.)
        nd = min(d,0.)
        ne = min(e,0.)
        nf = min(f,0.)

        ! gradient depends on which direction data is traveling
        IF (phi(i,j,k) > 0.) THEN
            gradX = max(pa*pa,nb*nb)
            gradY = max(pc*pc,nd*nd)
            gradZ = max(pe*pe,nf*nf)
        ELSE
            gradX = max(pb*pb,na*na)
            gradY = max(pd*pd,nc*nc)
            gradZ = max(pf*pf,ne*ne)
        END IF

        ! magnitude of the gradient of phi
        gM = sqrt(gradX+gradY+gradZ)

        dPlus = max(a,0.)**2+min(b,0.)**2+max(c,0.)**2+min(d,0.)**2+max(e,0.)**2+min(f,0.)**2
        dPlus = sqrt(dPlus)

        dMinus = min(a,0.)**2+max(b,0.)**2+min(c,0.)**2+max(d,0.)**2+min(e,0.)**2+max(f,0.)**2
        dMinus = sqrt(dMinus)



    END SUBROUTINE upwind 




    !*************************************************************************************!
    ! Determine Sign
    !*************************************************************************************!

    !SUBROUTINE phiSign(pS,sgn,dx)
    SUBROUTINE phiSign(pS,sgn,dxx,gM)

        REAL,INTENT(IN) :: pS,dxx,gM
        REAL,INTENT(OUT) :: sgn

        ! non smeared
       ! IF (pS < 0.) THEN
       !     sgn =-1.
       ! ELSEIF (pS > 0.) THEN
       !     sgn = 1.
       ! ELSE
       !     sgn = 0.
       ! END IF

        ! smeared
        sgn = pS/sqrt(pS*pS + dxx*dxx*gM)


    ENDSUBROUTINE phiSign



    !*************************************************************************************!
    ! Determine Narrow Band
    !*************************************************************************************!

    SUBROUTINE narrowBand(nx,ny,nz,dx,D_s,phiNB,phiSB,nPass)
        IMPLICIT NONE


        INTEGER,INTENT(IN) :: nPass
        INTEGER,INTENT(IN) :: nx,ny,nz
        INTEGER :: i,j,k
        REAL,INTENT(IN) :: dx
        !REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
        TYPE (SET_smaller),DIMENSION(-nPass:nx+nPass,-nPass:ny+nPass,-nPass:nz+nPass),INTENT(IN) :: D_s
        INTEGER,DIMENSION(0:nx,0:ny,0:nz),INTENT(INOUT) :: phiNB,phiSB

        phiNB = 0
        phiSB = 0


        DO i = 0,nx
            DO j = 0,ny
                DO k = 0,nz

                    ! narrow band
                    IF (abs(D_s(i,j,k)%phi) < 5.1*dx) THEN
                        phiNB(i,j,k) = 1
                    END IF

                    ! stencil band
                    IF (abs(D_s(i,j,k)%phi) < 10.1*dx) THEN
                        phiSB(i,j,k) = 1
                    END IF

                END DO
            END DO
        END DO



    END SUBROUTINE narrowBand

    !*************************************************************************************!
    ! Calculate First Derivative
    !*************************************************************************************!

    SUBROUTINE firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order,gMM,nPass,coords,dims)
        IMPLICIT NONE

        INTEGER,DIMENSION(3),INTENT(IN) :: coords,dims
        INTEGER,INTENT(IN) :: nPass
        INTEGER,INTENT(IN) :: i,j,k,order
        INTEGER,INTENT(IN) :: nx,ny,nz
        REAL,INTENT(IN) :: dx
        !REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
        !TYPE (SET_smaller),DIMENSION(-2:nx+2,-2:ny+2,-2:nz+2),INTENT(IN) :: D_s
        REAL              ,DIMENSION(-nPass:nx+nPass,-nPass:ny+nPass,-nPass:nz+nPass),INTENT(IN) :: phi
        REAL,INTENT(OUT) :: phiX,phiY,phiZ,gMM
        REAL :: aa1,aa2,aa3,aa4,aa5,aa6,aa7
        REAL::dx12
        INTEGER :: ip1,jp1,kp1
        INTEGER :: im1,jm1,km1
        INTEGER :: ip2,jp2,kp2
        INTEGER :: im2,jm2,km2
        INTEGER :: ip3,jp3,kp3
        INTEGER :: im3,jm3,km3

        IF (order == 2) THEN

            ! second order first derivative
            phiX = (-1./2.*phi(i-1,j,k) +  1./2.*phi(i+1,j,k))/dx;
            phiY = (-1./2.*phi(i,j-1,k) +  1./2.*phi(i,j+1,k))/dx;
            phiZ = (-1./2.*phi(i,j,k-1) +  1./2.*phi(i,j,k+1))/dx;

        ELSEIF (order == 4) THEN 

            dx12 = 1./(12.*dx)

   IF (i<2 .AND. coords(1)==0) THEN
      phiX = (-25./12.*phi(i,j,k)+4.*phi(i+1,j,k)-3.*phi(i+2,j,k)+4./3.*phi(i+3,j,k)-1./4.*phi(i+4,j,k))/dx
   ELSEIF (j<2 .AND. coords(2)==0) THEN 
      phiY = (-25./12.*phi(i,j,k)+4.*phi(i,j+1,k)-3.*phi(i,j+2,k)+4./3.*phi(i,j+3,k)-1./4.*phi(i,j+4,k))/dx
   ELSEIF (k<2 .AND. coords(3)==0) THEN 
      phiZ = (-25./12.*phi(i,j,k)+4.*phi(i,j,k+1)-3.*phi(i,j,k+2)+4./3.*phi(i,j,k+3)-1./4.*phi(i,j,k+4))/dx
   ELSEIF (i>nx-2 .AND. coords(1)==dims(1)-1) THEN 
      phiX =-(-25./12.*phi(i,j,k)+4.*phi(i-1,j,k)-3.*phi(i-2,j,k)+4./3.*phi(i-3,j,k)-1./4.*phi(i-4,j,k))/dx
   ELSEIF (j>ny-2 .AND. coords(2)==dims(2)-1) THEN
      phiY =-(-25./12.*phi(i,j,k)+4.*phi(i,j-1,k)-3.*phi(i,j-2,k)+4./3.*phi(i,j-3,k)-1./4.*phi(i,j-4,k))/dx
   ELSEIF (k>nz-2 .AND. coords(3)==dims(3)-1) THEN
      phiZ =-(-25./12.*phi(i,j,k)+4.*phi(i,j,k-1)-3.*phi(i,j,k-2)+4./3.*phi(i,j,k-3)-1./4.*phi(i,j,k-4))/dx
   ELSE
        ! fourth order first derivatives
        phiX = (-phi(i+2,j,k)+8.*phi(i+1,j,k)-8.*phi(i-1,j,k)+phi(i-2,j,k))*dx12
        phiY = (-phi(i,j+2,k)+8.*phi(i,j+1,k)-8.*phi(i,j-1,k)+phi(i,j-2,k))*dx12
        phiZ = (-phi(i,j,k+2)+8.*phi(i,j,k+1)-8.*phi(i,j,k-1)+phi(i,j,k-2))*dx12
   END IF

        ELSEIF (order == 6) THEN ! summation by parts
            PRINT*, " Working on this ... "
            STOP

            im1 = i-1
            jm1 = j-1
            km1 = k-1
            im2 = i-2
            jm2 = j-2
            km2 = k-2
            im3 = i-3
            jm3 = j-3
            km3 = k-3

            ip1 = i+1
            jp1 = j+1
            kp1 = k+1
            ip2 = i+2
            jp2 = j+2
            kp2 = k+2
            ip3 = i+3
            jp3 = j+3
            kp3 = k+3



            aa1 =-1./60.;
            aa2 = 3./20.;
            aa3 =-3./4.;
            aa4 = 0.;
            aa5 = 3./4.;
            aa6 =-3./20.;
            aa7 = 1./60.;

            ! calculate derivatives
            phiX = (phi(im3,j,k)*aa1 + phi(im2,j,k)*aa2 + phi(im1,j,k)*aa3 + &
                &         phi(i,j,k)*aa4 + phi(ip1,j,k)*aa5 + phi(ip2,j,k)*aa6 + phi(ip3,j,k)*aa7)/dx    

            phiY = (phi(i,jm3,k)*aa1 + phi(i,jm2,k)*aa2 + phi(i,jm1,k)*aa3 + &
                &         phi(i,j,k)*aa4 + phi(i,jp1,k)*aa5 + phi(i,jp2,k)*aa6 + phi(i,jp3,k)*aa7)/dx   

            phiZ = (phi(i,j,km3)*aa1 + phi(i,j,km2)*aa2 + phi(i,j,km1)*aa3 + &
                &         phi(i,j,k)*aa4 + phi(i,j,kp1)*aa5 + phi(i,j,kp2)*aa6 + phi(i,j,kp3)*aa7)/dx   

        ELSE 
            PRINT*," This order not set yet in firstDeriv"
            STOP
        END IF

        gMM = phiX*phiX + phiY*phiY + phiZ*phiZ
        gMM = sqrt(gMM)

    END SUBROUTINE firstDeriv



    !*************************************************************************************!
    ! Calculate Second Derivative 
    !*************************************************************************************!

    SUBROUTINE secondDeriv(i,j,k,nx,ny,nz,dx,phi,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,order,nPass)
        IMPLICIT NONE



        INTEGER,INTENT(IN) :: nPass
        REAL :: dxx,dx12
        INTEGER,INTENT(IN) :: i,j,k,order
        INTEGER,INTENT(IN) :: nx,ny,nz
        REAL,INTENT(IN) :: dx
        !REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
        !TYPE (SET_smaller),DIMENSION(-2:nx+2,-2:ny+2,-2:nz+2),INTENT(IN) :: D_s
        REAL              ,DIMENSION(-nPass:nx+nPass,-nPass:ny+nPass,-nPass:nz+nPass),INTENT(IN) :: phi
        !REAL,DIMENSION(0:nx,0:ny,0:nz,3),INTENT(IN) :: gradPhi
        REAL,INTENT(OUT) :: phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ
        REAL :: aa1,aa2,aa3,aa4,aa5,aa6,aa7
        REAL :: bb1,bb2,bb3,bb4,bb5,bb6,bb7
        INTEGER :: im2,ip2,jm2,jp2,km2,kp2
        INTEGER :: ip1,jp1,kp1,im1,jm1,km1,ip3,jp3,kp3,im3,jm3,km3



        IF (order == 2) THEN
            dxx = 1./(dx*dx); 
            ! second order second derivatives
            phiXX = (-2.*phi(i,j,k) + phi(i+1,j,k) + phi(i-1,j,k))*dxx;          
            phiYY = (-2.*phi(i,j,k) + phi(i,j+1,k) + phi(i,j-1,k))*dxx;
            phiZZ = (-2.*phi(i,j,k) + phi(i,j,k+1) + phi(i,j,k-1))*dxx;

            ! second order mixed derivatives
            phiXY = phi(i+1,j+1,k  )-phi(i+1,j-1,k  )-phi(i-1,j+1,k  )+phi(i-1,j-1,k  )
            phiYZ = phi(i  ,j+1,k+1)-phi(i  ,j+1,k-1)-phi(i  ,j-1,k+1)+phi(i  ,j-1,k-1)
            phiXZ = phi(i+1,j  ,k+1)-phi(i+1,j  ,k-1)-phi(i-1,j  ,k+1)+phi(i-1,j  ,k-1)

            phiXY = phiXY*dxx/4.
            phiYZ = phiYZ*dxx/4.
            phiXZ = phiXZ*dxx/4.


            !phiXY = (-1./2.*gradPhi(i-1,j,k,2) +  1./2.*gradPhi(i+1,j,k,2))/dx
            !phiYZ = (-1./2.*gradPhi(i,j-1,k,3) +  1./2.*gradPhi(i,j+1,k,3))/dx
            !phiXZ = (-1./2.*gradPhi(i,j,k-1,1) +  1./2.*gradPhi(i,j,k+1,1))/dx


!            phiXY = (-1./2.*D_s(i-1,j,k)%gradPhiY +  1./2.*D_s(i+1,j,k)%gradPhiY)/dx;
!            phiYZ = (-1./2.*D_s(i,j-1,k)%gradPhiZ +  1./2.*D_s(i,j+1,k)%gradPhiZ)/dx;
!            phiXZ = (-1./2.*D_s(i,j,k-1)%gradPhiX +  1./2.*D_s(i,j,k+1)%gradPhiX)/dx;

        ELSEIF (order == 4) THEN 

            dx12 = 1./(12.*dx)

            !fourth order second derivatives
            phiXX = (-phi(i+2,j,k)+8.*phi(i+1,j,k)-8.*phi(i-1,j,k)+phi(i-2,j,k))*dx12
            phiYY = (-phi(i,j+2,k)+8.*phi(i,j+1,k)-8.*phi(i,j-1,k)+phi(i,j-2,k))*dx12
            phiZZ = (-phi(i,j,k+2)+8.*phi(i,j,k+1)-8.*phi(i,j,k-1)+phi(i,j,k-2))*dx12
            ! fourth order mixed derivatives
            phiXY = (-phi(i+2,j,k)+8.*phi(i+1,j,k)-8.*phi(i-1,j,k)+phi(i-2,j,k))*dx12
            phiXZ = (-phi(i,j,k+2)+8.*phi(i,j,k+1)-8.*phi(i,j,k-1)+phi(i,j,k-2))*dx12
            phiYZ = (-phi(i,j+2,k)+8.*phi(i,j+1,k)-8.*phi(i,j-1,k)+phi(i,j-2,k))*dx12


        ELSEIF (order == 6) THEN 
            PRINT*, " Working on this ... "
            STOP

            im1 = i-1
            jm1 = j-1
            km1 = k-1
            im2 = i-2
            jm2 = j-2
            km2 = k-2
            im3 = i-3
            jm3 = j-3
            km3 = k-3

            ip1 = i+1
            jp1 = j+1
            kp1 = k+1
            ip2 = i+2
            jp2 = j+2
            kp2 = k+2
            ip3 = i+3
            jp3 = j+3
            kp3 = k+3

            aa1 =-1./60.
            aa2 = 3./20.
            aa3 =-3./4.
            aa4 = 0.
            aa5 = 3./4.
            aa6 =-3./20.
            aa7 = 1./60.

            bb1 = 1./90.
            bb2 =-3./20.
            bb3 = 3./2.
            bb4 =-49/18.
            bb5 = 3./2.
            bb6 =-3./20.
            bb7 = 1./90.


            !sixth order second derivatives
            phiXX = (phi(im3,j,k)*bb1 + phi(im2,j,k)*bb2 + phi(im1,j,k)*bb3 + &
                &          phi(i,j,k)*bb4 + phi(ip1,j,k)*bb5 + phi(ip2,j,k)*bb6 + phi(ip3,j,k)*bb7)/dx/dx

            phiYY = (phi(i,jm3,k)*bb1 + phi(i,jm2,k)*bb2 + phi(i,jm1,k)*bb3 + &
                &          phi(i,j,k)*aa4 + phi(i,jp1,k)*bb5 + phi(i,jp2,k)*bb6 + phi(i,jp3,k)*bb7)/dx/dx   

            phiZZ = (phi(i,j,km3)*bb1 + phi(i,j,km2)*bb2 + phi(i,j,km1)*bb3 + &
                &          phi(i,j,k)*bb4 + phi(i,j,kp1)*bb5 + phi(i,j,kp2)*bb6 + phi(i,j,kp3)*bb7)/dx/dx      


            !sixth order mixed derivatives
            phiXY = (phi(im3,j,k)*aa1 + phi(im2,j,k)*aa2 + phi(im1,j,k)*aa3 + &
                &          phi(i,j,k)*aa4 + phi(ip1,j,k)*aa5 + phi(ip2,j,k)*aa6 + &
                &          phi(ip3,j,k)*aa7)/dx   

            phiYZ = (phi(i,jm3,k)*aa1 + phi(i,jm2,k)*aa2 + phi(i,jm1,k)*aa3 + &
                &          phi(i,j,k)*aa4 + phi(i,jp1,k)*aa5 + phi(i,jp2,k)*aa6 + &
                &          phi(i,jp3,k)*aa7)/dx   

            phiXZ = (phi(i,j,km3)*aa1 + phi(i,j,km2)*aa2 + phi(i,j,km1)*aa3 + &
                &          phi(i,j,k)*aa4 + phi(i,j,kp1)*aa5 + phi(i,j,kp2)*aa6 + & 
                &          phi(i,j,kp3)*aa7)/dx   



        ELSE 

            PRINT*," This order not set yet in secondDeriv"
            STOP

        END IF


    END SUBROUTINE secondDeriv

    !*************************************************************************************!
    ! Calculate Min/Max Flow
    !*************************************************************************************!

    SUBROUTINE minMax(i,j,k,nx,ny,nz,phi,D_s,D_all,F,nPass)
        IMPLICIT NONE

        INTEGER,INTENT(IN) :: nPass
        REAL,DIMENSION(-nPass:,-nPass:,-nPass:),INTENT(IN) :: phi
        REAL :: curv,b1,b2,b3,b4,b5,b6,gM,gradMag2,pAve,phiX,phiY,phiZ,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,thresh
        INTEGER,INTENT(IN) :: i,j,k
        INTEGER,INTENT(IN) :: nx,ny,nz
        !REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
        TYPE (SET_smaller),DIMENSION(-nPass:nx+nPass,-nPass:ny+nPass,-nPass:nz+nPass),INTENT(IN) :: D_s
        !REAL              , DIMENSION(-nPass:nx+nPass,-nPass:ny+nPass,-nPass:nz+nPass) :: phi
        !REAL,DIMENSION(0:nx,0:ny,0:nz,3),INTENT(IN) :: gradPhi,grad2Phi,gradMixPhi
        !REAL,DIMENSION(-nPass:nx+nPass,-nPass:ny+nPass,-nPass:nz+nPass,3),INTENT(IN) :: grad2Phi,gradMixPhi,gridX
        TYPE (SET),DIMENSION(-1:,-1:,-1:),INTENT(INOUT) :: D_all
        REAL,INTENT(OUT) :: F
        INTEGER :: h


        phiX = D_s(i,j,k)%gradPhiX
        phiY = D_s(i,j,k)%gradPhiY
        phiZ = D_s(i,j,k)%gradPhiZ

        phiXX = D_all(i,j,k)%gradPhiXX
        phiYY = D_all(i,j,k)%gradPhiYY
        phiZZ = D_all(i,j,k)%gradPhiZZ

        phiXY = D_all(i,j,k)%gradPhiXY
        phiXZ = D_all(i,j,k)%gradPhiXZ
        phiYZ = D_all(i,j,k)%gradPhiYZ

        ! gradient magnitude
        gradMag2 = phiX*phiX+phiY*phiY+phiZ*phiZ             
        IF (gradMag2 < 1.E-13) THEN
            gM = 0.
        ELSE
            gM = sqrt(gradMag2)
        END IF  

        ! calculate curvature
        IF (gM < 1.E-13) THEN
            curv = 0.0;
        ELSE 

            b1 = phiYY + phiZZ
            b2 = phiXX + phiZZ
            b3 = phiXX + phiYY
            b4 = phiX*phiY*phiXY
            b5 = phiX*phiZ*phiXZ
            b6 = phiY*phiZ*phiYZ

            !curv = b1*phiX*phiX + b2*phiY*phiY + b3*phiZ*phiZ - 2.*b4 - 2.*b5 - 2.*b6
            !curv  = phiXX + phiYY + phiZZ
            !denom = gM
            !curv = curv/denom
            curv = phiXX+phiYY+phiZZ

        END IF
        !IF (D_all(i,j,k)%x < -50.) THEN
            h = 1
        !ELSE
            !h = 10
        !END IF

        pAve = phi(i,j,k)+phi(i-1,j,k)+phi(i+1,j,k)+phi(i,j+1,k)+phi(i,j-1,k)+phi(i,j,k+1)+phi(i,j,k-1)
        !pAve = phi(i,j,k)+phi(i-h,j,k)+phi(i+h,j,k)+phi(i,j+h,k)+phi(i,j-h,k)+phi(i,j,k+h)+phi(i,j,k-h)

        pAve = pAve/7.

        ! threshold value
        thresh = 0.      

        ! min/max switch
        IF (pAve < thresh) THEN
            F = min(curv,0.0)
        ELSE
            F = max(curv,0.0)
        END IF

    END SUBROUTINE minMax

SUBROUTINE weno(gM,i,j,k,nx,ny,nz,dx,phi,dPlus,dMinus,nPass,coords,dims)

    IMPLICIT NONE

    INTEGER,DIMENSION(3),INTENT(IN) :: coords,dims
    INTEGER,INTENT(IN) :: nx,ny,nz,i,j,k,nPass
    REAL,INTENT(IN) :: dx
    REAL,DIMENSION(-nPass:nx+nPass,-nPass:ny+nPass,-nPass:nz+nPass),INTENT(IN) :: phi
    REAL,INTENT(OUT):: gM,dPlus,dMinus
    REAL :: gradX,gradY,gradZ,pa,pb,pc,pd,pe,pf,na,nb,nc,nd,ne,nf,a,b,c,d,e,f,uzm1,uzm2,uzm3
    REAL :: u,ap,am,bp,bm,cp,cm,dp,dm,IS0p,IS0m,IS1p,IS1m,IS2p,IS2m,epsp,epsm
    REAL :: a0p,a0m,a1p,a1m,a2p,a2m,w0p,w0m,w2p,w2m,PWp,PWm,p0,p1,p2,p3,p4,p5
    INTEGER :: im,jm,km,ip,jp,kp,im2,ip2,jm2,jp2,km2,kp2
    INTEGER :: ip1,jp1,kp1,im1,jm1,km1,ip3,jp3,kp3,im3,jm3,km3

    ! calculate WENO5 derivatives

    IF ((i>3.OR.coords(1)>0).AND.(i<nx-4.OR.coords(1)<dims(1)-1).AND.(j>3.OR.coords(2)>0).AND.&
    (j<ny-4.OR.coords(2)<dims(2)-1).AND.(k>3.OR.coords(3)>0).AND.(k<nz-4.OR.coords(3)<dims(3)-1)) THEN

        ! x-direction
        ap = (phi(i+3,j,k)-2.*phi(i+2,j,k)+phi(i+1,j,k))/dx
        am = (phi(i-3,j,k)-2.*phi(i-2,j,k)+phi(i-1,j,k))/dx
        bp = (phi(i+2,j,k)-2.*phi(i+1,j,k)+phi(i  ,j,k))/dx
        bm = (phi(i-2,j,k)-2.*phi(i-1,j,k)+phi(i  ,j,k))/dx
        cp = (phi(i+1,j,k)-2.*phi(i  ,j,k)+phi(i-1,j,k))/dx
        cm = cp
        dp = bm
        dm = bp

        IS0p = 13.*(ap-bp)*(ap-bp) + 3.*(ap-3.*bp)*(ap-3.*bp)
        IS0m = 13.*(am-bm)*(am-bm) + 3.*(am-3.*bm)*(am-3.*bm)
        IS1p = 13.*(bp-cp)*(bp-cp) + 3.*(bp +  cp)*(bp +  cp)
        IS1m = 13.*(bm-cm)*(bm-cm) + 3.*(bm +  bm)*(bm +  bm)
        IS2p = 13.*(cp-dp)*(cp-dp) + 3.*(3.*cp-dp)*(3.*cp-dp)
        IS2m = 13.*(cm-dm)*(cm-dm) + 3.*(3.*cm-dm)*(3.*cm-dm)

        p0 = (phi(i-2,j,k) - phi(i-3,j,k))/dx
        p1 = (phi(i-1,j,k) - phi(i-2,j,k))/dx
        p2 = (phi(i  ,j,k) - phi(i-1,j,k))/dx
        p3 = (phi(i+1,j,k) - phi(i  ,j,k))/dx
        p4 = (phi(i+2,j,k) - phi(i+1,j,k))/dx
        p5 = (phi(i+3,j,k) - phi(i+2,j,k))/dx

        ! epsilon with scaling term
        epsp = (1.E-6)*max(p1*p1,max(p2*p2,max(p3*p3,max(p4*p4,p5*p5)))) + 1.E-99
        epsm = (1.E-6)*max(p0*p0,max(p1*p1,max(p2*p2,max(p3*p3,p4*p4)))) + 1.E-99

        a0p = 1./((epsp+IS0p)*(epsp+IS0p))
        a0m = 1./((epsm+IS0m)*(epsm+IS0m))
        a1p = 6./((epsp+IS1p)*(epsp+IS1p))
        a1m = 6./((epsm+IS1m)*(epsm+IS1m))
        a2p = 3./((epsp+IS2p)*(epsp+IS2p))
        a2m = 3./((epsm+IS2m)*(epsm+IS2m))  

        w0p = a0p/(a0p+a1p+a2p)
        w0m = a0m/(a0m+a1m+a2m)
        w2p = a2p/(a0p+a1p+a2p)
        w2m = a2m/(a0m+a1m+a2m)

        PWp = 1./3.*w0p*(ap-2.*bp+cp)+1./6.*(w2p-0.5)*(bp-2.*cp+dp)
        PWm = 1./3.*w0m*(am-2.*bm+cm)+1./6.*(w2m-0.5)*(bm-2.*cm+dm)

        a = 1./12.*(-p1+7.*p2+7.*p3-p4) - PWm
        b = 1./12.*(-p1+7.*p2+7.*p3-p4) + PWp

        ! y-direction
        ap = (phi(i,j+3,k)-2.*phi(i,j+2,k)+phi(i,j+1,k))/dx
        am = (phi(i,j-3,k)-2.*phi(i,j-2,k)+phi(i,j-1,k))/dx
        bp = (phi(i,j+2,k)-2.*phi(i,j+1,k)+phi(i,j  ,k))/dx
        bm = (phi(i,j-2,k)-2.*phi(i,j-1,k)+phi(i,j  ,k))/dx
        cp = (phi(i,j+1,k)-2.*phi(i,j  ,k)+phi(i,j-1,k))/dx
        cm = cp
        dp = bm
        dm = bp

        IS0p = 13.*(ap-bp)*(ap-bp) + 3.*(ap-3.*bp)*(ap-3.*bp)
        IS0m = 13.*(am-bm)*(am-bm) + 3.*(am-3.*bm)*(am-3.*bm)
        IS1p = 13.*(bp-cp)*(bp-cp) + 3.*(bp +  cp)*(bp +  cp)
        IS1m = 13.*(bm-cm)*(bm-cm) + 3.*(bm +  bm)*(bm +  bm)
        IS2p = 13.*(cp-dp)*(cp-dp) + 3.*(3.*cp-dp)*(3.*cp-dp)
        IS2m = 13.*(cm-dm)*(cm-dm) + 3.*(3.*cm-dm)*(3.*cm-dm)

        p0 = (phi(i,j-2,k) - phi(i,j-3,k))/dx
        p1 = (phi(i,j-1,k) - phi(i,j-2,k))/dx
        p2 = (phi(i,j  ,k) - phi(i,j-1,k))/dx
        p3 = (phi(i,j+1,k) - phi(i,j  ,k))/dx
        p4 = (phi(i,j+2,k) - phi(i,j+1,k))/dx
        p5 = (phi(i,j+3,k) - phi(i,j+3,k))/dx

        ! epsilon with scaling term
        epsp = (1.E-6)*max(p1*p1,max(p2*p2,max(p3*p3,max(p4*p4,p5*p5)))) + 1.E-99
        epsm = (1.E-6)*max(p0*p0,max(p1*p1,max(p2*p2,max(p3*p3,p4*p4)))) + 1.E-99

        a0p = 1./((epsp+IS0p)*(epsp+IS0p))
        a0m = 1./((epsm+IS0m)*(epsm+IS0m))
        a1p = 6./((epsp+IS1p)*(epsp+IS1p))
        a1m = 6./((epsm+IS1m)*(epsm+IS1m))
        a2p = 3./((epsp+IS2p)*(epsp+IS2p))
        a2m = 3./((epsm+IS2m)*(epsm+IS2m))  

        w0p = a0p/(a0p+a1p+a2p)
        w0m = a0m/(a0m+a1m+a2m)
        w2p = a2p/(a0p+a1p+a2p)
        w2m = a2m/(a0m+a1m+a2m)

        PWp = 1./3.*w0p*(ap-2.*bp+cp)+1./6.*(w2p-0.5)*(bp-2.*cp+dp)
        PWm = 1./3.*w0m*(am-2.*bm+cm)+1./6.*(w2m-0.5)*(bm-2.*cm+dm)

        c = 1./12.*(-p1+7.*p2+7.*p3-p4) - PWm
        d = 1./12.*(-p1+7.*p2+7.*p3-p4) + PWp

        ! z-direction
        ap = (phi(i,j,k+3)-2.*phi(i,j,k+2)+phi(i,j,k+1))/dx
        am = (phi(i,j,k-3)-2.*phi(i,j,k-2)+phi(i,j,k-1))/dx
        bp = (phi(i,j,k+2)-2.*phi(i,j,k+1)+phi(i,j,k  ))/dx
        bm = (phi(i,j,k-2)-2.*phi(i,j,k-1)+phi(i,j,k  ))/dx
        cp = (phi(i,j,k+1)-2.*phi(i,j,k  )+phi(i,j,k-1))/dx
        cm = cp
        dp = bm
        dm = bp

        IS0p = 13.*(ap-bp)*(ap-bp) + 3.*(ap-3.*bp)*(ap-3.*bp)
        IS0m = 13.*(am-bm)*(am-bm) + 3.*(am-3.*bm)*(am-3.*bm)
        IS1p = 13.*(bp-cp)*(bp-cp) + 3.*(bp +  cp)*(bp +  cp)
        IS1m = 13.*(bm-cm)*(bm-cm) + 3.*(bm +  bm)*(bm +  bm)
        IS2p = 13.*(cp-dp)*(cp-dp) + 3.*(3.*cp-dp)*(3.*cp-dp)
        IS2m = 13.*(cm-dm)*(cm-dm) + 3.*(3.*cm-dm)*(3.*cm-dm)

        p0 = (phi(i,j,k-2) - phi(i,j,k-3))/dx
        p1 = (phi(i,j,k-1) - phi(i,j,k-2))/dx
        p2 = (phi(i,j,k  ) - phi(i,j,k-1))/dx
        p3 = (phi(i,j,k+1) - phi(i,j,k  ))/dx
        p4 = (phi(i,j,k+2) - phi(i,j,k+1))/dx
        p5 = (phi(i,j,k+3) - phi(i,j,k+2))/dx

        ! epsilon with scaling term
        epsp = (1.E-6)*max(p1*p1,max(p2*p2,max(p3*p3,max(p4*p4,p5*p5)))) + 1.E-99
        epsm = (1.E-6)*max(p0*p0,max(p1*p1,max(p2*p2,max(p3*p3,p4*p4)))) + 1.E-99

        a0p = 1./((epsp+IS0p)*(epsp+IS0p))
        a0m = 1./((epsm+IS0m)*(epsm+IS0m))
        a1p = 6./((epsp+IS1p)*(epsp+IS1p))
        a1m = 6./((epsm+IS1m)*(epsm+IS1m))
        a2p = 3./((epsp+IS2p)*(epsp+IS2p))
        a2m = 3./((epsm+IS2m)*(epsm+IS2m))  

        w0p = a0p/(a0p+a1p+a2p)
        w0m = a0m/(a0m+a1m+a2m)
        w2p = a2p/(a0p+a1p+a2p)
        w2m = a2m/(a0m+a1m+a2m)

        PWp = 1./3.*w0p*(ap-2.*bp+cp)+1./6.*(w2p-0.5)*(bp-2.*cp+dp)
        PWm = 1./3.*w0m*(am-2.*bm+cm)+1./6.*(w2m-0.5)*(bm-2.*cm+dm)

        e = 1./12.*(-p1+7.*p2+7.*p3-p4) - PWm
        f = 1./12.*(-p1+7.*p2+7.*p3-p4) + PWp


    ELSE

        ! set plus and minus integers
        im = i-1
        jm = j-1
        km = k-1
        ip = i+1
        jp = j+1
        kp = k+1

        ! calculate derivatives
        a = (phi(i ,j ,k ) - phi(im,j ,k ))/dx
        b = (phi(ip,j ,k ) - phi(i ,j ,k ))/dx
        c = (phi(i ,j ,k ) - phi(i ,jm,k ))/dx
        d = (phi(i ,jp,k ) - phi(i ,j ,k ))/dx
        e = (phi(i ,j ,k ) - phi(i ,j ,km))/dx
        f = (phi(i ,j ,kp) - phi(i ,j ,k ))/dx

    END IF


    ! check positive   
    pa = max(a,0.)
    pb = max(b,0.)
    pc = max(c,0.)
    pd = max(d,0.)
    pe = max(e,0.)
    pf = max(f,0.)

    ! check negative
    na = min(a,0.)
    nb = min(b,0.)
    nc = min(c,0.)
    nd = min(d,0.)
    ne = min(e,0.)
    nf = min(f,0.)

    ! gradient depends on which direction data is traveling
    IF (phi(i,j,k) > 0.) THEN
        gradX = max(pa*pa,nb*nb)
        gradY = max(pc*pc,nd*nd)
        gradZ = max(pe*pe,nf*nf)
    ELSE
        gradX = max(pb*pb,na*na)
        gradY = max(pd*pd,nc*nc)
        gradZ = max(pf*pf,ne*ne)
    END IF

    ! magnitude of the gradient of phi
    gM = sqrt(gradX+gradY+gradZ)

    dPlus = max(a,0.)**2+min(b,0.)**2+max(c,0.)**2+min(d,0.)**2+max(e,0.)**2+min(f,0.)**2
    dPlus = sqrt(dPlus)

    dMinus = min(a,0.)**2+max(b,0.)**2+min(c,0.)**2+max(d,0.)**2+min(e,0.)**2+max(f,0.)**2
    dMinus = sqrt(dMinus)


END SUBROUTINE weno 










!*************************************************************************************!
! End of Code 
!*************************************************************************************!

END MODULE set3d_SUBs
