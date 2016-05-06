      SUBROUTINE NBTRLU(A,B,C,E,D,IL,IU,I,JL,JU,JV,ORDER,IPS,X)  
      implicit real(a-h,o-z)

!
!     Subroutine: solverUtilNbtrlu
!
!     Block tridiagonal solver
!
!     Location:
!      file:Code/Solver/nbtrlu.F
!
!     Notes:
!      None
!
!     Arguments:
!      None
!

      DIMENSION A(1),B(1),C(1),E(1),D(1),IPS(JV,1),X(JV,1)
      INTEGER ORDER,ORDSQ,ORDJV,OSQJV
      ORDSQ=ORDER**2
      ORDJV=ORDER*JV
      OSQJV=ORDSQ*JV
      IF(I.EQ.IL) THEN
        I0MAT=1+(I-1)*OSQJV
        I0VEC=1+(I-1)*ORDJV
        CALL LUDPVT(B,JL,JU,ORDER,JV,IPS)
        CALL LUSPVT(B,D(I0VEC),D(I0VEC),X,JL,JU,ORDER,JV,IPS)
        DO 1 J=1,ORDER
          I0MATJ=I0MAT+(J-1)*ORDJV
          IVEC=1+(J-1)*ORDJV
          CALL LUSPVT(B,C(IVEC),E(I0MATJ),X,JL,JU,ORDER,JV,IPS)
    1   CONTINUE
      END IF
      IF(I.GT.IL.AND.I.LT.IU) THEN
        I0MAT=1+(I-1)*OSQJV
        I0VEC=1+(I-1)*ORDJV
        I1MAT=I0MAT-OSQJV
        I1VEC=I0VEC-ORDJV
        CALL MULPUT(A,D(I1VEC),D(I0VEC),JL,JU,ORDER,JV,X)
        DO 2 J=1,ORDER
          IVEC=1+(J-1)*ORDJV
          I1MATJ=I1MAT+(J-1)*ORDJV
          CALL MULPUT(A,E(I1MATJ),B(IVEC),JL,JU,ORDER,JV,X)
    2   CONTINUE
C,AJK 1-10-11        CALL LUDECO(B,JL,JU,ORDER,JV,X)
        IF (ORDER > 1) CALL LUDECO(B,JL,JU,ORDER,JV,X)
        CALL LUSOLV(B,D(I0VEC),D(I0VEC),JL,JU,ORDER,JV,X)
        DO 3 J=1,ORDER
          I0MATJ=I0MAT+(J-1)*ORDJV
          IVEC=1+(J-1)*ORDJV
          CALL LUSOLV(B,C(IVEC),E(I0MATJ),JL,JU,ORDER,JV,X)
    3   CONTINUE
      END IF
      IF(I.EQ.IU) THEN
        I0MAT=1+(I-1)*OSQJV
        I0VEC=1+(I-1)*ORDJV
        I1MAT=I0MAT-OSQJV
        I1VEC=I0VEC-ORDJV
        CALL MULPUT(A,D(I1VEC),D(I0VEC),JL,JU,ORDER,JV,X)
        DO 4 J=1,ORDER
          IVEC=1+(J-1)*ORDJV
          I1MATJ=I1MAT+(J-1)*ORDJV
          CALL MULPUT(A,E(I1MATJ),B(IVEC),JL,JU,ORDER,JV,X)
    4   CONTINUE
        CALL LUDPVT(B,JL,JU,ORDER,JV,IPS)
        CALL LUSPVT(B,D(I0VEC),D(I0VEC),X,JL,JU,ORDER,JV,IPS)
      END IF
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE NBTRBK(E,D,IL,IU,JL,JU,JV,ORDER,X)
      implicit real(a-h,o-z)
      DIMENSION E(1),D(1),X(JV,1)
      INTEGER ORDER,ORDSQ,ORDJV,OSQJV
      ORDSQ=ORDER**2
      ORDJV=ORDER*JV
      OSQJV=ORDSQ*JV
      I=IU
   10 CONTINUE
        I=I-1
        I0MAT=1+(I-1)*OSQJV
        I0VEC=1+(I-1)*ORDJV
        I1VEC=I0VEC+ORDJV
        CALL MULPUT(E(I0MAT),D(I1VEC),D(I0VEC),JL,JU,ORDER,JV,X)
      IF(I.GT.IL) GO TO 10
      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE LUDPVT(A,JL,JU,ORDER,JV,IPS)
      implicit real(a-h,o-z)
      INTEGER ORDER
      DIMENSION A(JV,ORDER,1),IPS(JV,ORDER)
      NM1=ORDER-1
      DO 1 I=1,ORDER
        DO 1 II=JL,JU
          IPS(II,I)=I
    1 CONTINUE
      DO 2 II=JL,JU
        DO 17 K=1,NM1
        BIG=0.0
        DO 11 I=K,ORDER
        IP=IPS(II,I)
        SIZE=ABS(A(II,IP,K))
        IF((SIZE-BIG).LE.0) GO TO 11
   10   BIG=SIZE
        IDXPIV=I
   11   CONTINUE
        IF((IDXPIV-K).EQ.0) GO TO 15
   14   J=IPS(II,K)
        IPS(II,K)=IPS(II,IDXPIV)
        IPS(II,IDXPIV)=J
   15   KP=IPS(II,K)
        PIVOT=A(II,KP,K)
        KP1=K+1
        DO 16 I=KP1,ORDER
        IP=IPS(II,I)
        EM=-A(II,IP,K)/PIVOT
        A(II,IP,K)=-EM
        DO 16 J=KP1,ORDER
        A(II,IP,J)=A(II,IP,J)+EM*A(II,KP,J)
   16   CONTINUE
   17   CONTINUE
    2 CONTINUE
      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE LUSPVT(A,B,C,X,JL,JU,ORDER,JV,IPS)
      implicit real(a-h,o-z)
      INTEGER ORDER
      DIMENSION A(JV,ORDER,1),B(JV,1),C(JV,1),X(JV,1),IPS(JV,1)
      NP1=ORDER+1
      DO 101 II=JL,JU
      IP=IPS(II,1)
      X(II,1)=B(II,IP)
      DO 2 I=2,ORDER
      IP=IPS(II,I)
      IM1=I-1
      SUM=0.0
      DO 1 J=1,IM1
    1 SUM=SUM+A(II,IP,J)*X(II,J)
    2 X(II,I)=B(II,IP)-SUM
      IP=IPS(II,ORDER)
      C(II,ORDER)=X(II,ORDER)/A(II,IP,ORDER)
      DO 4 IBACK=2,ORDER
      I=NP1-IBACK
      IP=IPS(II,I)
      IP1=I+1
      SUM=0.0
      DO 3 J=IP1,ORDER
    3 SUM=SUM+A(II,IP,J)*C(II,J)
    4 C(II,I)=(X(II,I)-SUM)/A(II,IP,I)
  101 CONTINUE
      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE LUDECO(A,JL,JU,ORDER,JV,SUM)
      implicit real(a-h,o-z)
      INTEGER ORDER
      DIMENSION A(JV,ORDER,1),SUM(JV)
      DO 1 JC=2,ORDER
        DO 2 II=JL,JU
          SUM(II)=A(II,1,JC)/A(II,1,1)
    2   CONTINUE
        DO 1 II=JL,JU
          A(II,1,JC)=SUM(II)
    1 CONTINUE
      JRJC=1
   10 CONTINUE
      JRJC=JRJC+1
      JRJCM1=JRJC-1
      JRJCP1=JRJC+1
      DO 3 JR=JRJC,ORDER
        DO 4 II=JL,JU
          SUM(II)=A(II,JR,JRJC)
    4   CONTINUE
        DO 5 JM=1,JRJCM1
          DO 5 II=JL,JU
            SUM(II)=SUM(II)-A(II,JR,JM)*A(II,JM,JRJC)
    5   CONTINUE
        DO 3 II=JL,JU
          A(II,JR,JRJC)=SUM(II)
    3 CONTINUE
      IF(JRJC.EQ.ORDER) RETURN
      DO 6 JC=JRJCP1,ORDER
        DO 7 II=JL,JU
          SUM(II)=A(II,JRJC,JC)
    7   CONTINUE
        DO 8 JM=1,JRJCM1
          DO 8 II=JL,JU
            SUM(II)=SUM(II)-A(II,JRJC,JM)*A(II,JM,JC)
    8   CONTINUE
        DO 9 II=JL,JU
          SUM(II)=SUM(II)/A(II,JRJC,JRJC)
    9   CONTINUE
        DO 6 II=JL,JU
          A(II,JRJC,JC)=SUM(II)
    6 CONTINUE
      GO TO 10
      END
C------------------------------------------------------------------
      SUBROUTINE LUSOLV(A,B,C,JL,JU,ORDER,JV,SUM)
      implicit real(a-h,o-z)
      INTEGER ORDER
      DIMENSION A(JV,ORDER,1),B(JV,1),C(JV,1),SUM(JV)
      DO 1 II=JL,JU
        C(II,1)=B(II,1)/A(II,1,1)
    1 CONTINUE
      DO 2 JR=2,ORDER
        JRM1=JR-1
        DO 3 II=JL,JU
          SUM(II)=B(II,JR)
    3   CONTINUE
        DO 4 JM=1,JRM1
          DO 4 II=JL,JU
            SUM(II)=SUM(II)-A(II,JR,JM)*C(II,JM)
    4   CONTINUE
        DO 2 II=JL,JU
          C(II,JR)=SUM(II)/A(II,JR,JR)
    2 CONTINUE
      DO 5 JRJR=2,ORDER
        JR=ORDER-JRJR+1
        JRP1=JR+1
        DO 6 II=JL,JU
          SUM(II)=C(II,JR)
    6   CONTINUE
        DO 7 JMJM=JRP1,ORDER
          JM=ORDER-JMJM+JRP1
          DO 7 II=JL,JU
            SUM(II)=SUM(II)-A(II,JR,JM)*C(II,JM)
    7   CONTINUE
        DO 5 II=JL,JU
          C(II,JR)=SUM(II)
    5 CONTINUE
      RETURN
      END
C---------------------------------------------------------------
      SUBROUTINE MULPUT(A,B,C,JL,JU,ORDER,JV,SUM)
      implicit real(a-h,o-z)
      INTEGER ORDER
      DIMENSION A(JV,1),B(JV,1),C(JV,1),SUM(JV)
      DO 1 JR=1,ORDER
        DO 2 II=JL,JU
          SUM(II)=0.0
    2   CONTINUE
        DO 3 JC=1,ORDER
          IA=JR+(JC-1)*ORDER
          DO 4 II=JL,JU
            SUM(II)=SUM(II)+A(II,IA)*B(II,JC)
    4     CONTINUE
    3   CONTINUE
        DO 5 II=JL,JU
          C(II,JR)=C(II,JR)-SUM(II)
    5   CONTINUE
    1 CONTINUE
      RETURN
      END
