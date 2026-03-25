c
c
c
c
c
c
      subroutine VHS12  (MODE,LP,L1,M,U,INCU,UP,C,INCC,ICV,NCV)
C
      implicit none
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            MODE,LP,L1,M,INCU,INCC,ICV,NCV
cuk      real               U(1),UP,C(1)
cuk      real               U(m+1),UP,C(ncv*m)  ! did NOT work...
      real               U(100),UP,C(100)   ! in this case only!!! (single-ev-l)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IJ,ILP,IL1,IM,INCR,I2,I3,I4,J
      real               SM,B
      real               ONE,CL,CLINV,SM1
C                                  FIRST EXECUTABLE STATEMENT
      ONE=1.e0
C
      IF (0.GE.LP.OR.LP.GE.L1.OR.L1.GT.M) GO TO 9005
      ILP=(LP-1)*INCU+1
      IL1=(L1-1)*INCU+1
      IM=(M-1)*INCU+1
      CL=ABS(U(ILP))
      IF (MODE.EQ.2) GO TO 15
C                                  CONSTRUCT THE TRANSFORMATION.
      DO 5 IJ=IL1,IM,INCU
    5 CL=MAX1(ABS(U(IJ)),CL)
      IF (CL.LE.0.0e0) GO TO 9005
      CLINV=ONE/CL
      SM=(U(ILP)*CLINV)**2
      DO 10 IJ=IL1,IM,INCU
   10 SM=SM+(U(IJ)*CLINV)**2
C                                  CONVERT DBLE. PREC. SM TO SNGL.
C                                    PREC. SM1
      SM1=SM
      CL=CL*SQRT(SM1)
      IF (U(ILP).GT.0.0e0) CL=-CL
      UP=U(ILP)-CL
      U(ILP)=CL
      GO TO 20
C                                  APPLY THE TRANSFORMATION
C                                    I+U*(U**T)/B TO C.
   15 IF (CL.LE.0.0e0) GO TO 9005
   20 IF (NCV.LE.0) GO TO 9005
      B=UP*U(ILP)
C                                  B MUST BE NONPOSITIVE HERE. IF B =
C                                    0., RETURN.
      IF (B.GE.0.0e0) GO TO 9005
      B=ONE/B
      I2=1-ICV+INCC*(LP-1)
      INCR=INCC*(L1-LP)
      DO 35 J=1,NCV
         I2=I2+ICV
         I3=I2+INCR
         I4=I3
         SM=C(I2)*UP
         DO 25 IJ=IL1,IM,INCU
            SM=SM+C(I3)*U(IJ)
            I3=I3+INCC
   25    CONTINUE
         IF (SM.EQ.0.0e0) GO TO 35
         SM=SM*B
         C(I2)=C(I2)+SM*UP
         DO 30 IJ=IL1,IM,INCU
            C(I4)=C(I4)+SM*U(IJ)
            I4=I4+INCC
   30    CONTINUE
   35 CONTINUE
 9005 RETURN
      END
