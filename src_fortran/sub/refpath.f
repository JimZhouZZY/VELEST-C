c
c
c
c
c
c
      subroutine REFPATH (xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,
     &tkj,jl,kk,didjkk,rp,nrp,refrerr)
c
c       refpath assigns raypaths with coordinates rp(1,i), rp(2,i),
c  rp(3,i) along a refracted seismic raypath in a flat layered earth.
c
c  input:  xe,ye,ze - event coordinates
c          xr,yr,zr - receiver coordinates
c             delta - horizontal distance between the event
c                          and the receiver
c                nl - total number of layers in model
c                     (used only to dimension v, vsq, and thk)
c              v(l) - velocity of layer l
c            vsq(l) = v(l) ** 2
c            thk(l) - thickness of layer l
c                jl - event layer
c               tkj - depth of event from top of event layer
c                kk - refracting layer
c            didjkk - critical distance for the refracted ray
c
c  output:  rp(1,2,or 3,i) - raypoint coordinates
c                      nrp - number of raypoints
c
cek max. nr of raypoints =50
c
      implicit none
c
      real xe,ye,ze,xr,yr,zr,delta,tkj,didjkk,refrerr
      integer nl,jl,kk,nrp
c
      real thkjl,d1,d2,tng
      integer nrpd,i,m,nrp1,nrpd2
      real v(nl),thk(nl),vsq(nl),rp(3,200)
      thkjl=thk(jl)
c
c   d1 and d2 are the x and y directions of the raypath
c
      d1=(xr-xe)/delta
      d2=(yr-ye)/delta
      nrpd=kk-jl+1
c
c   nrpd= number of raypoints in down-going part of ray path
c
      nrp=nrpd+kk
c
c   down-going part of path
c
      thk(jl)=thkjl-tkj
      rp(1,1)=xe
      rp(2,1)=ye
      rp(3,1)=ze
c   m is the number of the layer above raypoint i
c
      do 23204 i=2,nrpd
      m=jl+i-2
c
c   tng is the tangent of the incident angle in layer m
c
      tng=v(m)/sqrt(vsq(kk)-vsq(m))
      rp(1,i)=rp(1,i-1)+thk(m)*tng*d1
      rp(2,i)=rp(2,i-1)+thk(m)*tng*d2
      rp(3,i)=rp(3,i-1)+thk(m)
23204 continue
23205 continue
c
c   up-going part of path
c
      rp(1,nrpd+1)=rp(1,nrpd)+(delta-didjkk)*d1
      rp(2,nrpd+1)=rp(2,nrpd)+(delta-didjkk)*d2
      rp(3,nrpd+1)=rp(3,nrpd)
      thk(jl)=thkjl
      nrp1=nrp-1
c
c   nrpd2 is the number of the raypoint at the
c   top of layer (kk-1)
c
      nrpd2=nrpd+2
cuk      do 23206i=nrpd2,nrp1
      do 23206 i=nrpd2,nrp1+1  ! calculate up to the top, then compute error !!
c
c   m is the number of the layer below raypoint i
c
      m=(kk-1)-(i-nrpd2)
c
c   tng is the tangent of the incidence angle in layer m
c
      tng=v(m)/sqrt(vsq(kk)-vsq(m))
      rp(1,i)=rp(1,i-1)+thk(m)*tng*d1
      rp(2,i)=rp(2,i-1)+thk(m)*tng*d2
      rp(3,i)=rp(3,i-1)-thk(m)
23206 continue
c
      refrerr=sqrt( (rp(1,nrp)-xr)**2 + (rp(2,nrp)-yr)**2 +
     &              (rp(3,nrp)-zr)**2 ) *1000.
      rp(1,nrp)=xr
      rp(2,nrp)=yr
      rp(3,nrp)=zr
      return
      end ! of subr. refpath
