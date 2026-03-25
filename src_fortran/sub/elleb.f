c
c
c
c
c
c
      subroutine ELLEB(L,B,x,y)
      implicit none
      real b,x,y
      REAL L
C
C  SCHWEIZ. PROJEKTIONSSYSTEM  FORMELN VON H. ODERMATT
C  TRANSFORMATION ELLIPSOID - EBENE
C  L,B  UNTERSCHIEDE IN LAENGE UND BREITE ZU BERN (SEXAG.SEK.)
C  Y,X LANDESKOORDINATEN IN METER ( BERN 0/0)
C  MY  MERIDIANKONVERGENZ ( SEXAG. SEK.)
C
      integer i
      DOUBLE PRECISION D(5),E(5),F(5),RW(5),IW(5),P,Q,A,C,bb,bl
      DATA             BB,BL/169028.66,26782.5/
      A= L*3600.-BL
      C= B*3600.-BB
      D(1)= 0.68382546262761
      E(1)= 2.363591607471D-2
      D(2)=-3.91798328045D-8
      E(2)= 0.
      D(3)= 1.4965410352D-15
      E(3)= 4.527219881D-17
      D(4)=-8.039471422D-23
      E(4)=-3.89081120D-24
      D(5)= 7.0021390D-30
      E(5)= 2.3407700D-31
      F(1)= 4.515344386039D1
      F(2)= 1.17912305209D-4
      F(3)= 5.8474201864D-10
      F(4)= 2.73386187D-15
      F(5)= 1.4308547D-20
C
      P=30.91849390613 * A
      Q=C*F(5)
      I=5
  1   I=I-1
      Q=C*(Q+F(I))
      IF(I.GT.1) GOTO 1
C
      RW(1)=Q
      IW(1)=P
      DO 2 I=2,5
      RW(I)=Q*RW(I-1)-P*IW(I-1)
      IW(I)=P*RW(I-1)+Q*IW(I-1)
  2   CONTINUE
C
      y=D(5)*RW(5)
      x=D(5)*IW(5)
      I=5
  3   I=I-1
      y=y+D(I)*RW(I)
      x=x+D(I)*IW(I)
      IF(I.GT.1) GOTO 3
      y= y/1000.+200.       !  output in kilometers
      x= x/1000.+600.
C
      RETURN
      END ! of subr. ELLEB
