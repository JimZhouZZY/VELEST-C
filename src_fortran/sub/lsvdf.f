c
c
c
c
c
c
      subroutine LSVDF  (A,IA,M,N,B,IB,NB,S,WK,IER)
c     This is an IMSL-subroutine
      implicit none
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,M,N,IB,NB,IER
      real               A(IA,N),B(IB,1),S(N),WK(N,2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,JP1,K,L,MM,NN,NNP1,NS,NSP1
      real               ZERO,ONE,T
      integer            nm
      real               f
      DATA               ZERO/0.0e0/,ONE/1.0e0/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
C                                  BEGIN SPECIAL FOR ZERO ROWS AND
C                                    COLS. PACK THE NONZERO COLS TO THE
C                                    LEFT
      NN=N
      IER=34
      IF (NN.LE.0.OR.M.LE.0) GO TO 9000
      IER=0
      J=NN
    5 CONTINUE
      DO 10 I=1,M
         IF (A(I,J).NE.ZERO) GO TO 25
   10 CONTINUE
C                                  COL J IS ZERO. EXCHANGE IT WITH COL
C                                    N
      IF (J.EQ.NN) GO TO 20
      DO 15 I=1,M
   15 A(I,J)=A(I,NN)
   20 CONTINUE
      A(1,NN)=J
      NN=NN-1
   25 CONTINUE
      J=J-1
      IF (J.GE.1) GO TO 5
C                                  IF N=0 THEN A IS ENTIRELY ZERO AND
C                                    SVD COMPUTATION CAN BE SKIPPED
      NS=0
      IF (NN.EQ.0) GO TO 120
C                                  PACK NONZERO ROWS TO THE TOP QUIT
C                                    PACKING IF FIND N NONZERO ROWS
      I=1
      MM=M
   30 IF (I.GT.N.OR.I.GE.MM) GO TO 75
      IF (A(I,I).NE.ZERO) GO TO 40
      DO 35 J=1,NN
         IF (A(I,J).NE.ZERO) GO TO 40
   35 CONTINUE
      GO TO 45
   40 I=I+1
      GO TO 30
C                                  ROW I IS ZERO EXCHANGE ROWS I AND M
   45 IF (NB.LE.0) GO TO 55
      DO 50 J=1,NB
         T=B(I,J)
         B(I,J)=B(MM,J)
         B(MM,J)=T
   50 CONTINUE
   55 DO 60 J=1,NN
   60 A(I,J)=A(MM,J)
      IF (MM.GT.NN) GO TO 70
      DO 65 J=1,NN
   65 A(MM,J)=ZERO
   70 CONTINUE
C                                  EXCHANGE IS FINISHED
      MM=MM-1
      GO TO 30
C
   75 CONTINUE
C                                  END SPECIAL FOR ZERO ROWS AND
C                                    COLUMNS
C                                  BEGIN SVD ALGORITHM..
C                                  (1) REDUCE THE MATRIX TO UPPER
C                                    BIDIAGONAL FORM WITH HOUSEHOLDER
C                                    TRANSFORMATIONS.
C                                    H(N)...H(1)AQ(1)...Q(N-2) =
C                                    (D**T,0)**T WHERE D IS UPPER
C                                    BIDIAGONAL.
C                                  (2) APPLY H(N)...H(1) TO B. HERE
C                                    H(N)...H(1)*B REPLACES B IN
C                                    STORAGE.
C                                  (3) THE MATRIX PRODUCT W=
C                                    Q(1)...Q(N-2) OVERWRITES THE FIRST
C                                    N ROWS OF A IN STORAGE.
C                                  (4) AN SVD FOR D IS COMPUTED. HERE K
C                                    ROTATIONS RI AND PI ARE COMPUTED
C                                    SO THAT RK...R1*D*P1**(T)...PK**(T)
C                                    = DIAG(S1,...,SM) TO WORKING
C                                    ACCURACY. THE SI ARE NONNEGATIVE
C                                    AND NONINCREASING. HERE RK...R1*B
C                                    OVERWRITES B IN STORAGE WHILE
C                                    A*P1**(T)...PK**(T) OVERWRITES A
C                                    IN STORAGE.
C                                  (5) IT FOLLOWS THAT,WITH THE PROPER
C                                    DEFINITIONS, U**(T)*B OVERWRITES
C                                    B, WHILE V OVERWRITES THE FIRST N
C                                    ROW AND COLUMNS OF A.
      L=MIN0(MM,NN)
C                                  THE FOLLOWING LOOP REDUCES A TO
C                                    UPPER BIDIAGONAL AND ALSO APPLIES
C                                    THE PREMULTIPLYING TRANSFORMATIONS
C                                    TO B.
      DO 85 J=1,L
         IF (J.GE.MM) GO TO 80
         JP1=MIN0(J+1,NN)
         CALL VHS12 (1,J,J+1,MM,A(1,J),1,T,A(1,JP1),1,IA,NN-J)
         CALL VHS12 (2,J,J+1,MM,A(1,J),1,T,B,1,IB,NB)
   80    IF (J.GE.NN-1) GO TO 85
         CALL VHS12 (1,J+1,J+2,NN,A(J,1),IA,WK(J,2),A(J+1,1),IA,1,MM-J)
   85 CONTINUE
C                                  COPY THE BIDIAGONAL MATRIX INTO THE
C                                    ARRAY S FOR LSVDB
      IF (L.EQ.1) GO TO 95
      DO 90 J=2,L
         S(J)=A(J,J)
         WK(J,1)=A(J-1,J)
   90 CONTINUE
   95 S(1)=A(1,1)
C
      NS=NN
      IF (MM.GE.NN) GO TO 100
      NS=MM+1
      S(NS)=ZERO
      WK(NS,1)=A(MM,MM+1)
  100 CONTINUE
C                                  CONSTRUCT THE EXPLICIT N BY N
C                                    PRODUCT MATRIX, W=Q1*Q2*...*QL*I
C                                    IN THE ARRAY A
      DO 115 K=1,NN
         I=NN+1-K
         IF (I.GT.MIN0(MM,NN-2)) GO TO 105
         CALL VHS12 (2,I+1,I+2,NN,A(I,1),IA,WK(I,2),A(1,I+1),1,IA,NN-I)
  105    DO 110 J=1,NN
  110    A(I,J)=ZERO
         A(I,I)=ONE
  115 CONTINUE
C                                  COMPUTE THE SVD OF THE BIDIAGONAL
C                                    MATRIX
C
C      LEVEL=1
C      CALL UERSET(LEVEL,LEVOLD)
      CALL LSVDB (S(1),WK(1,1),NS,A,IA,NN,B,IB,NB,IER)
C                                  TEST FOR IER=33
C
      IF (IER.GT.128) GO TO 9000
C      CALL UERSET(LEVOLD,LEVOLD)
      IF (IER.NE.33) GO TO 120
cuk      WRITE(6,4500) IER
 4500 FORMAT(2X,'SUBR. LSVDF: IER=',I4,'AFTER CALL OF SUBR.LSVDB')
      T=0.0e0
      NM=MIN0(M,N)
      IF (S(1).NE.ZERO) T=S(NM)/S(1)
      F=100.0e0+T
      IF (F.EQ.100.0e0) GO TO 120
      IER=0
  120 CONTINUE
      IF (NS.GE.NN) GO TO 130
      NSP1=NS+1
      DO 125 J=NSP1,NN
  125 S(J)=ZERO
  130 CONTINUE
      IF (NN.EQ.N) GO TO 155
      NNP1=NN+1
C                                  MOVE RECORD OF PERMUTATIONS AND
C                                    STORE ZEROS
      DO 140 J=NNP1,N
         S(J)=A(1,J)
         IF (NN.LT.1) GO TO 140
         DO 135 I=1,NN
  135    A(I,J)=ZERO
  140 CONTINUE
C                                  PERMUTE ROWS AND SET ZERO SINGULAR
C                                    VALUES
      DO 150 K=NNP1,N
         I=S(K)
         S(K)=ZERO
         DO 145 J=1,N
            A(K,J)=A(I,J)
  145    A(I,J)=ZERO
         A(I,K)=ONE
  150 CONTINUE
C                                  END SPECIAL FOR ZERO ROWS AND
C                                    COLUMNS
  155 IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
C      CALL UERTST (IER,'LSVDF ')
cuk
      if(ier.ne.33)then
      WRITE(6,4501) IER
 4501 FORMAT(2X,'SUBR. LSVDF: AT END, IER=',I4,'CHECK SUBR. HEAD FOR',
     11X,'MEANING')
      endif
cuk
 9005 RETURN
      END
