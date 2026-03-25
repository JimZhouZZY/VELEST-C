c
c
c
c
c
c
      subroutine REGWORLD(ityp,cord1,cord2,place,nreg)
      implicit none
      integer ityp,nreg,last,num,in,nlo,latin,mhemis,i,iuk
      integer nval,i1,i2
      real cord1,cord2,xnrlon,galat,galon,pi
      real galor,galar,gelat,gelon,xazr,xdr,q,as
      real pirad,pi2,gelac,xlatr,gelor,ck1,ck2,ac
      real alon,xlat,xlon,xreg
      integer indfe,lt50,lt25
      character irname*24300
      COMMON /REGIONcom/ INDFE(730),LT50(111),LT25(421),
     &                   irname,XNRLON(14400)
C     PROGRAM TO COMPUTE REGION NUMBER AND NAME BASED ON THE INPUT
C        ITYP <=0         LATITUDE&LONGITUDE
C        ITYP > 0         DISTANCE & AZIMUTH
C    THE REGION NAMES ARE FROM FLINN & ENGDAHL#S TABLES
      character place*32
      real BUF(40)
      character*1 cbk
      DATA GALAT/46.7706/,GALON/8.09428/
      data cbk/' '/
C
      PI=3.141593
      PIRAD=PI/180.
      PI2=PI*2.
      LAST=89
C
C     CONVERT ARRAY COORDINATES TO RADIANS
      GALOR=GALON*PIRAD
      GALAR=GALAT*PIRAD
      IF(ITYP) 850,850,860
C
  850 GELAT=CORD1
      GELON=CORD2
      NUM=1
      GO TO 212
C
  860 XAZR=CORD2*PIRAD
      XDR= CORD1*PIRAD
      NUM=1
C
      Q=SIN(GALAR)*COS(XDR)+COS(GALAR)*SIN(XDR)*COS(XAZR)
      IF(ABS(Q).GT.1.) GO TO 891
      AS=ASIN(Q)
      GO TO 890
C
  891 CONTINUE
      GELON=0.
      GELAT=0.
      NUM=NUM-1
      GO TO 212
C
  890 GELAC=AS
      XLATR=GELAC
      IF(COS(GELAC))45,46,45
   45 IF(COS(GALAR)) 47,46,47
   46 GELOR=0.0
      GO TO 34
C
   47 CONTINUE
      CK1=(COS(XDR)-SIN(GALAR)*SIN(XLATR))/
     *(COS(GALAR)*COS(XLATR))
C
      CK2=(SIN(XAZR)*SIN(XDR))/COS(XLATR)
      IF(CK1.GE.0.) GO TO 37
      IF(CK2.GE.0.) GO TO 882
C
      IN=1
      GO TO 881
C
  882 IN=2
      GO TO 883
C
   37 IF(CK2.GE.0.) GO TO 885
C
      IN=3
      GO TO 881
C
  885 IN=4
C
  883 IF(ABS(CK1).GT.1.0) GO TO 887
      AC=ACOS(CK1)
      GO TO 888
C
  881 IF(ABS(CK2).GT.1.0) GO TO 887
      AS=ASIN(CK2)
C
  888 GO TO (33,32,31,36),IN
C
   33 GELOR=GALOR-PI-AS
      GO TO 34
C
   32 GELOR=GALOR+AC
      GO TO 34
C
   31 GELOR=GALOR+AS
      GO TO 34
C
   36 GELOR=GALOR+AC
      GO TO 34
C
  887 CONTINUE
      GELOR=0.0
C
C     CONVERT LONGITUDE AND LATITUDE FROM RADIANS TO DEGREES
C     AND REDUCE LONGITUDE TO A FORM LESS THAN 360 DEGREES
C
   34 CONTINUE
      ALON=ABS(GELOR)
  506 IF(ALON.LT.PI2) GO TO 116
      ALON=ALON-PI2
      GO TO 506
C
  116 GELOR=SIGN(ALON,GELOR)
      IF(ABS(GELOR).GT.PI) GELOR=GELOR-SIGN(PI2,GELOR)
C
C  LONGITUDE IN DEGREES
C
      GELON=GELOR/PIRAD
C
C     LATITUDE IS CONVERTED FROM GEOCENTRIC TO GEODETIC COORDINATES
C
      GELAT=GELAC/PIRAD
C
  212 IF(NUM.EQ.1) GO TO 5000
      place=' '
      GO TO 900
C
 5000 CONTINUE
C
C      CLEAR BUFFER
C
C     ASSIGN REGION NUMBER AND NAME
C
      place=' '
  215 XLAT=GELAT
      XLON=GELON
      IF(XLON.LT.0.) XLON=XLON+360.
C
C     DETERMINE EARTH QUADRANT
C
      IF(XLAT.LT.0.) GO TO 4
      IF(XLON.GT.180.) GO TO 6
C
C     LATITUDE POSITIVE,LONGITUDE 0-180
C
      MHEMIS=0
      GO TO 9
C
C     LATITUDE POSITIVE,LONGITUDE 181-360
C
    6 MHEMIS=3600
      XLON=360.-XLON
      GO TO 9
C
    4 IF(XLON.GT.180.) GO TO 8
C
C     LATITUDE NEGATIVE,LONGITUDE 0-180
C
      MHEMIS=7200
      GO TO 9
C
C     LATITUDE NEGATIVE,LONGITUDE 181-360
C
    8 MHEMIS=10800
      XLON=360.-XLON
C
    9 CONTINUE
C
C     CHECK FOR ZERO LATITUDE
C
      IF(ABS(XLAT).GE.1.0) GO TO 13
      LATIN=0
      GO TO 17
C
   13 LATIN=ABS(XLAT)
      IF(LATIN.GT.LAST) LATIN=LAST
C
   17 CONTINUE
C
      NLO=ABS(XLON)
      I=MHEMIS + 40*LATIN + 1
      do iuk=1,40
         buf(iuk)=xnrlon(i+iuk-1)
      enddo
      NVAL=BUF(1)*1000.+1.1
      DO 56 I=2,NVAL
      IF(NLO.LT.int(BUF(I))) GOTO 57
56    CONTINUE
      I=NVAL+1
57    NREG=BUF(I-1)
      XREG=BUF(I-1)-NREG
      NREG=XREG*1000. + 0.1
      I1=INDFE(NREG)+1
      I2=INDFE(NREG+1)
      place=irname(i1:i2)
C
  900 continue
      return
C
      END ! of subr. regworld
c
cek    end of vel_topo.f
c
cek    begin of matrinv.f
c
      SUBROUTINE matrinv(N,A,B)
c
c MATRINV invertiert eine SYMMETRISCHE NxN-Matrix
c
c Aufruf: call MATRINV(N,A,B)
c
c         N = Anz. Zeilen/Kolonnen der Matrix
c         A = zu invertierende Matrix (Input)
c         B = inv(A)    (Output)
c
      DIMENSION A(N*N),B(N*N)
      B(1)= 1.0/A(1)
      IF(N-1) 60,60,2
2     NN= N*N
      DO 4 I=2,NN
4     B(I)= 0.0
      MM=1
      KN= 0
      DO 50 M=2,N
      K= M-1
      MM= MM+N+1
      KN= KN+N
      EK= A(MM)
      MI= M-N
      DO 10 I=1,K
      MI= MI+N
      IJ= I-N
      JM= KN
      DO 10 J= 1,K
      IJ= IJ+N
      JM= JM+1
10    EK= EK-A(MI)*B(IJ)*A(JM)
      B(MM)= 1.0/EK
      MI= M-N
      IM= KN
      DO 30 I=1,K
      IM= IM+1
      IJ= I-N
      JM= KN
      DO 20 J= 1,K
      IJ= IJ+N
      JM= JM+1
20    B(IM)= B(IM)-B(IJ)*A(JM)*B(MM)
      MI= MI+N
30    B(MI)= B(IM)
      IM= KN
      DO 40 I=1,K
      IM= IM+1
      MJ= M-N
      IJ= I-N
      DO 40 J=1,K
      MJ=MJ+N
      IJ= IJ+N
40    B(IJ)= B(IJ)+B(IM)*B(MJ)*EK
50    CONTINUE
60    continue
      return
      END
