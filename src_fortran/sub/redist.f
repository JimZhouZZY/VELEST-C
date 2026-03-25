c
c
c
c
c
c
      subroutine REDIST(Xkm,Ykm,XLAT,XLON)   ! 7. Jan. 1988
C   by Edi
C CALCULATES DISTANCE BETWEEN TWO POINTS,  INPUT IN RECTANGULAR XY-COOR.
C   OUTPUT IN LATITUDE/LONGITUDE
C
C     CONVERSION OF  KILOMETERS TO LATITUDE AND LONGITUDE
c
      implicit none
      real xkm,ykm,xlat,xlon
c
      real olat,olon,aa,bb,bc,sint,cost,rotate
      real xx,yy,y,x,q,lat,yp,bcl,p,lon
      integer icoordsystem
      common/GEO_COORSYSTEM/olat,olon,
     &                      rearth,ellip,rlatc,rad,
     &                      aa,bb,bc,sint,cost,rotate, icoordsystem
      double precision rearth,ellip,rlatc,rad
c
      DOUBLE PRECISION LAT1,LAT2,LAT3,clat1       ! neu
C
c
      xx=xkm
      yy=ykm
c*****	Following code is commented-out
c	and replaced by correct code to undo rotation
c	WLE 1/3/95

c      rphi=0.0   ! 5.4.91 u.k.
cc
c      SINT=DSIN(RPHI*RAD)
c      COST=DCOS(RPHI*RAD)
c      TANT=DTAN(RPHI*RAD)
cC---- ROTATE COORDINATES CLOCKWISE BACK BY RPHI
c      IF(RPHI.LT.0.1) GOTO 120
c      Y=YY*COST-XX*SINT
c      X=XX/SINT+Y*TANT
c      GOTO 130
c  120 X=XX
c      Y=YY

c	rotate coordinate system into original orientation

	x=xx*cost-yy*sint
	y=yy*cost+xx*sint
c	WLE 1/3/95

  130 CONTINUE
      IF(ABS(AA).LT.0.0000001) GOTO 900
      Q=Y/AA
      LAT=(Q+OLAT)/60.
      XLAT=Q+OLAT - 60.*LAT
      YP=60.*LAT+XLAT
      LAT1=DATAN(RLATC*DTAN(YP*RAD/60.0))
      LAT2=DATAN(RLATC*DTAN(OLAT*RAD/60.))
      LAT3=(LAT1+LAT2)/2.
      CLAT1=DCOS(LAT3)
      BCL=BB*CLAT1
      IF(ABS(BCL).LT.0.000001) GOTO 900
      P=X/(BB*CLAT1)
      LON=(P+OLON)/60.
      XLON=P+OLON - 60.*LON
      xlat=lat+xlat/60.
      xlon=lon+xlon/60.
      RETURN
  900 WRITE(6,1000) AA,BB,CLAT1
 1000 FORMAT(/,2X,' SUBR. REDISt: AA=',F10.5,2X,'BB=',F10.5,2X,
     1'COS(LAT1)=',F10.7,5X,'DIVISION BY ZERO, RUN STOPS HERE',/)
      stop'REDIST>>> division by zero!!'
      END ! of subr. redist
