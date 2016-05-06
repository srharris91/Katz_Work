PROGRAM reduce_ave
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER,PARAMETER::NPROCS=8
    INTEGER::rank,new_rank,sendbuf,recvbuf,numtasks
    INTEGER::ranks1(4) = (/ 0,1,2,3 /),ranks2(4) = (/ 4,5,6,7 /),ierr
    INTEGER:: orig_group,new_group,new_comm
    INTEGER:: sum_local=0,sum_all=0
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,  rank,       ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,  numtasks,   ierr)
    IF (numtasks .NE. NPROCS) THEN
        WRITE(*,*) 'Must specifiy NPROCS = ',NPROCS,'TERMINATING!!'
        CALL MPI_FINALIZE(ierr)
        STOP
    ENDIF
    WRITE(*,*) 'Number of processors and rank and numtasks = ',NPROCS,rank,numtasks
    sum_local = rank
    CALL MPI_ALLREDUCE(sum_local,sum_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    WRITE(*,*) 'local and all sum = ',sum_local,sum_all
    CALL MPI_FINALIZE(ierr)
END PROGRAM reduce_ave
