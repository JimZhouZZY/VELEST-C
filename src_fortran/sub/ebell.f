c
c
c
c
c
c
      subroutine EBELL(YL,XB,L,B,MY)
C
C  SCHWEIZ. PROJEKTIONSSYSTEM FORMELN VON H. ODERMATT
C  TRANSFORMATION EBENE ELLIPSOID
C  YL,XB LANDESKOORDINATEN IN METER (BERN 600/200)
C  L,B LAENGE UND BREITE  (GRAD)
C  MY MERIDIEANKONVERGENZ (SEXAG.SEK.)
C
      implicit none
      REAL YL,XB,L,B,MY
      integer i
      DOUBLE PRECISION A(8),BB(8),C(8),RZ(8),QZ(8),P,Q,X,Y,DB,DL
cccc      dimension A(8),BB(8),C(8),RZ(8),QZ(8)
cccc      real P,Q,X,Y,DB,DL
C
      X=1000.*XB-200000.
      Y=1000.*YL-600000.
C
      A(1)=1.4623614572021
      BB(1)=3.4564252673326D-2
      A(2)=1.225255821052D-7
      BB(2)=2.89600437564D-9
      A(3)=1.3687923002D-14
      BB(3)=4.651046030D-16
      A(4)=1.971224191D-21
      BB(4)=6.43850954D-23
      A(5)=2.97898051D-28
      BB(5)=9.600412D-30
      A(6)=4.650273D-35
      BB(6)=1.50512D-36
C     A(7)=7.48203D-42
      A(7)=0.D00
C     BB(7)=2.422D-43
      BB(7)=0.D00
C     A(8)=1.229D-48
      A(8)=0.D00
      C(1)= 2.2146704979846D-2
      C(2)=-1.280815253730D-9
      C(3)= 7.4775676024D-18
      C(4)= 4.691943327D-24
      C(5)=-3.6550101D-31
      C(6)= 3.71615D-39
C     C(7)= 1.6901D-45
      C(7)= 0.D00
C     C(8)= 1.96D-52
      C(8)= 0.D00
      RZ(1)= X
      QZ(1)= Y
cuk      DO 1 I=2,8   !  sufficiant to loop until 6 only
      do 1 I=2,6      !  otherwise VAX has an overflow...
      RZ(I)= X*RZ(I-1)-Y*QZ(I-1)
      QZ(I)= Y*RZ(I-1)+X*QZ(I-1)
    1 CONTINUE
      I=8
      Q=A(I)*RZ(I)
      P=A(I)*QZ(I)
      MY=0.
    2 I=I-1
      Q=Q+A(I)*RZ(I)
      P=P+A(I)*QZ(I)
      MY=MY+BB(I)*QZ(I)
      IF(I.GT.1) GOTO 2
      DL=3.2343101932327D-2 * P
      DB=Q*C(8)
      I=8
    3 I=I-1
      DB=Q*(DB+C(I))
      IF(I.GT.1) GOTO 3
C
cc      write(6,*)'vor L= und B= in subr. ebell!'
      L=(DL + 26782.5D00)/3600.
      B=(DB + 169028.66D00)/3600.
      RETURN
      END ! of subr. EBELL
