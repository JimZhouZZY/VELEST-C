c
c
c
c
c
c
      subroutine DIST(xlat,xlon,xkm,ykm)   ! 7. Jan. 1988
c
      implicit none
      real xlat,xlon,xkm,ykm
c      
      real olat,olon,aa,bb,bc,sint,cost,rotate
      real q,yp,xx
      integer icoordsystem
      common/GEO_COORSYSTEM/olat,olon,
     &                      rearth,ellip,rlatc,rad,
     &                      aa,bb,bc,sint,cost,rotate, icoordsystem
      double precision rearth,ellip,rlatc,rad
c
      DOUBLE PRECISION LAT1,LAT2,LAT3
c     conversion of latitude and longitude to kilometers relative
c     to center of coordinates
c
      write(*,*) 'DEBUG dist inputs:'
      write(*,*) '  xlat=', xlat, ' xlon=', xlon

      write(*,*) 'DEBUG globals:'
      write(*,*) '  olat=', olat, ' olon=', olon
      write(*,*) '  rlatc=', rlatc, ' rad=', rad
      write(*,*) '  aa=', aa, ' bb=', bb
      write(*,*) '  cost=', cost, ' sint=', sint
      q=60*xlat-olat
      yp=q+olat
      lat1=datan(rlatc*dtan(RAD*yp/60.0))
      lat2=datan(rlatc*dtan(RAD*OLAT/60.0))
      LAT3=(LAT2+LAT1)/2.
      xx=60*xlon-olon
      q=q*aa
      xx=xx*bb*dcos(LAT3)
cc** rotate coordinate system clockwise
      yp=cost*q-sint*xx
      xx=cost*xx+sint*q
      q=yp
c
      xkm=xx
      ykm=q
      write(*,*) 'DEBUG intermediate variables:'
      write(*,*) '  q=', q, ' yp=', yp, ' xx=', xx
      write(*,*) '  lat1=', lat1, ' lat2=', lat2, ' lat3=', lat3
c
      return
      end ! of subr. dist
