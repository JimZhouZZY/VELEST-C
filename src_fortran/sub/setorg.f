c
c
c
c
c
c
      subroutine SETORG(orlat,orlon,rrotate,ifil) 
c
      implicit none
      real orlat,orlon,rrotate
      integer ifil
c
      real olat,olon,aa,bb,bc,sint,cost,rotate
      integer icoordsystem
      common/GEO_COORSYSTEM/olat,olon,
     &                      rearth,ellip,rlatc,rad,
     &                      aa,bb,bc,sint,cost,rotate, icoordsystem
      double precision rearth,ellip,rlatc,rad
c
c
      double precision r,lat1,lat2,dela,delb,PHI,BETA
c
c     O(r)LAT & O(r)LON : origin of cartesian coordinate system
c      north  w e s t
c
c
cEK Dez94      rotate=0.0  !......
c
      rotate=rrotate
cEK
      if(orlat.eq.0.0.and.orlon.eq.0.0)then
         olat=46.95240  ! BERN North
         olon=-7.439583  ! BERN West
      else
         olat=orlat
         olon=orlon
      endif
c
      olat=olat*60. ! minutes N
      olon=olon*60. ! minutes W
c
       rad=0.017453292d0
C
C  OLD ELLIPSOID:
C      rlatc = 0.99330647d0
C      data rearth /6378.163d0/, ellip /298.26d0/
C
C  NEW ELLIPSOID FOR WHOLE EARTH:   WGS72 == WORLD GEODETIC SYSTEM 1972
C
C  ALSO SET RLATC ACCORDING TO ORIGIN
C
      REARTH=6378.135D0
      ELLIP=298.26         ! (flattening)
C
c rlatc = tan(geocentric LAT) / TAN(GEOGRAPHIC LAT)
C
C
C CALCULATE RLATC:  CONVERSION FROM GEOGRAPHICAL LAT TO GEOCENTRICAL LAT
C
      PHI=OLAT*RAD/60.0              !  phi=geogr. lat
      BETA=PHI-DSIN(PHI*2.)/ELLIP    !  beta=geoc. lat
      RLATC=DTAN(BETA)/DTAN(PHI)
c
C
C WRITE ELLIPSOIDE CONSTANTS
C
      if(ifil.gt.0)then
         write(ifil,*)
         write(ifil,*)
         write(ifil,*)'SHORT DISTANCE CONVERSION on ELLIPSOIDE of'//
     &                ' WORLD GEODETIC SYSTEM 1972 (WGS72)'
         write(ifil,*)'=========================================='//
     &                '==================================='
         write(ifil,*)
         write(ifil,'('' Radius at equator (REARTH)= '',f10.5,
     &                ''  km'')') rearth
         write(ifil,'(''   1. / (ellipticity)      = '',f10.3)') ellip
         write(ifil,*)
         write(ifil,*)'Origin of cartesian coordinates [degrees]:'
         if(orlat.eq.0.0.and.orlon.eq.0.0)then
            write(ifil,*)' (Origin = city of BERNE, Switzerland)'
         endif
         write(ifil,'(1x,f12.7,'' N'',5x,f12.7,'' W'')')
     &               olat/60.,olon/60.
         write(ifil,*)
         write(ifil,*)' Rotation angle (in degr.) clockwise from'
         write(ifil,'(''North  rotate= '',f6.1)') rotate
         write(ifil,*)
      endif
c
c   calculate aa &  bb
c   length of one minute of lat and lon in km
c
      lat1=datan(rlatc*dtan(olat*rad/60.0))       ! geoc. lat for OLAT
      lat2=datan(rlatc*dtan((olat+1.)*rad/60.0))  ! geoc. lat for (OLAT+1min)
      dela=lat2 - lat1
      r=rearth*(1.0 - (dsin(lat1)**2)/ellip)      ! kugelradius fuer lat=OLAT
      aa=r*dela   ! aa = 1 min geogr. lat
      delb=dacos(dsin(lat1)**2 + dcos(rad/60.)*dcos(lat1)**2)
      bc=r*delb     ! bc = 1 min geogr. lon
      bb=r*delb/dcos(lat1)
      if(ifil.gt.0)then
         write(ifil,'('' Radius of sphere at OLAT = '',f10.3,'' km'')')r
         write(ifil,*)
         write(ifil,*)'Conversion of GEOGRAPHICAL LATITUDE to '//
     &                'GEOCENTRICAL LATITUDE:'
         write(ifil,*)'RLATC = TAN(GEOCENTR.LAT) / TAN(GEOGRAPH.LAT)'
         write(ifil,'(1x,''RLATC = '',f12.8)') rlatc
         write(ifil,*)
         write(ifil,4) aa, bc
 4       format (10x,'Short distance conversions:',/,
     &           10x,'one min lat ~ ', f7.4,' km ',/,
     &           10x,'one min lon ~ ', f7.4,' km ',/)
         write(ifil,*)
         write(ifil,*)
      endif
C
c***  convert coordinates with rotation cosines
      sint=sin(rotate*rad)
      cost=cos(rotate*rad)
31    continue
c
      return
      end ! of subr. setorg
