c
c
c
c
c
c
      subroutine DELAZ(alat,alon,blat,blon,del,dist,az)
c
c      computes distance and azimuth from a to b
c      a and b are in decimal degrees and n-e coordinates
c      del -- delta in degrees
c      dist -- distance in km
c      az -- azimuth from a to b clockwise from north in degrees
c
      implicit none
      real alat,alon,blat,blon,del,dist,az
c
      real aa,bb,bc,olat,olon,sint,cost,rotate
      common/corect/ aa,bb,bc,olat,olon
      common/coords/ sint,cost,rotate
      COMMON/ORIGI/ REARTH,ELLIP,RLATC,RAD          !  neu
      DOUBLE PRECISION REARTH,ELLIP,RLATC,RAD       !  neu
c
cc      double precision pi2,rad,flat
      double precision pi2,   flat
      double precision alatr,alonr,blatr,blonr
      double precision tana,geoa,acol,tanb,geob,bcol
      double precision diflon,cosdel,delr,top,den,azr,colat,radius
      double precision dtan,datan,dsin,dcos,dacos,datan2
      data pi2/1.570796d0/
cc      data rad/1.745329d-02/   ! defined in subr. SETORG
calt      data flat/.993231d0/
ctest      data flat/.99330647d0/ ! neu wie in setorg
ccc      call setorg(orlat,orlon)  ! is called separately !!
      flat=rlatc !!! from subr. SETORG
c-----convert to radians
      alatr=alat*rad
      alonr=alon*rad
      blatr=blat*rad
      blonr=blon*rad
c-----convert latitudes to geocentric colatitudes
      tana=flat*dtan(alatr)
      geoa=datan(tana)
      acol=pi2-geoa
      tanb=flat*dtan(blatr)
      geob=datan(tanb)
      bcol=pi2-geob
c-----calculate delta
      diflon=blonr-alonr
      cosdel=dsin(acol)*dsin(bcol)*dcos(diflon)+dcos(acol)*
     &dcos(bcol)
      delr=dacos(cosdel)
c-----calculate azimuth from a to b
      top=dsin(diflon)
      den=(dsin(acol)/dtan(bcol))-dcos(diflon)*dcos(acol)
      azr=datan2(top,den)
c----- convert to degrees
      del=delr/rad
      az=azr/rad
      if(az.lt.0.0) az=360.+az
c-----compute distance in kilometers
      colat=pi2-(alatr+blatr)/2.d0
cold      radius=6371.277d0*
cold     & (1.d0+(3.37853d-3)*((1.d0/3.d0)-((dcos(colat))**2)))
      radius=6378.163d0*      !  neu wie in setorg
     & (1.d0+(3.35278d-3)*((1.d0/3.d0)-((dcos(colat))**2))) ! neu wie in setorg
      dist=delr*radius
      return
c
c  1/298.26 = 0.0033527795   =  flattening   ( WGS72 )
c
      end ! of subr. delaz
