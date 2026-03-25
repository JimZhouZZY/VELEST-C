c
c
c
c
c
c
      subroutine STRPATH (xe,ye,ze,xr,yr,zr,rp,nrp)
c
c       strpath assigns raypoint coordinates for a straight path
c  between the event and the receiver.
c
c  input:  (xe,ye,ze) - event coordinates
c          (xr,yr,zr) - receiver coordinates
c
c  output:   rp(1,2,or 3,i) - coordinates of raypoint i
c                       nrp - number of raypoints = 2
c
      implicit none
      real xe,ye,ze,xr,yr,zr
      integer nrp
      real rp(3,2)
      nrp=2
      rp(1,1)=xe
      rp(2,1)=ye
      rp(3,1)=ze
      rp(1,2)=xr
      rp(2,2)=yr
      rp(3,2)=zr
      return
      end ! of subr. strpath
