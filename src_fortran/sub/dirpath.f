c
c
c
c
c
c
      subroutine DIRPATH (xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,
     &jl,tkj,salpha,deljl,rp,nrp,direrr)
c
c       dirpath assigns raypoints with coordinates rp(1,i),rp(2,i),
c  rp(3,i) along a direct seismic raypath in a flat layered earth.
c
c  input:  xe,ye,ze - event coordinates
c          xr,yr,zr - receiver coordinates
c             delta - horizontal distance between event and receiver
c                nl - total number of layers in model
c                     (used only to dimension v vsq and thk)
c              v(l) - velocity of layer l
c            vsq(l) = v(l) ** 2
c            thk(l) - thickness of layer l
c                jl - event layer
c               tkj - depth of event from top of event layer
c            salpha - sine of takeoff angle
c             deljl - horizontal travel distance in event layer
c
c  output:  rp(1,2,or 3,i) - raypoint coordinates
c                      nrp - number of raypoints
c
      implicit none
      real xe,ye,ze,xr,yr,zr,delta,tkj,salpha,deljl
      real direrr
      integer nl,jl,nrp,nrp1,i,m
      real d1,d2,tng
cek max. nr of raypoints =50
      real v(nl),thk(nl),vsq(nl),rp(3,200)
      d1=(xr-xe)/delta
      d2=(yr-ye)/delta
      nrp=jl+1
c
c   assign raypoint coordinates
c
      rp(1,1)=xe
      rp(2,1)=ye
      rp(3,1)=ze
      rp(1,2)=rp(1,1)+deljl*d1
      rp(2,2)=rp(2,1)+deljl*d2
      rp(3,2)=rp(3,1)-tkj
cuk      nrp1=nrp-1
      nrp1=nrp  ! calculate up to the top and then compute an error !
      do 23208i=3,nrp1
      m=nrp-(i-1)
c
c   m is the number of the layer below ray point i
c
      tng=v(m)*salpha/sqrt(vsq(jl)-vsq(m)*salpha**2)
c
c   tng is the tangent of the angle of incidence in layer m
c
      rp(1,i)=rp(1,i-1)+thk(m)*tng*d1
      rp(2,i)=rp(2,i-1)+thk(m)*tng*d2
      rp(3,i)=rp(3,i-1)-thk(m)
23208 continue
23209 continue
      direrr=sqrt( (rp(1,nrp)-xr)**2 + (rp(2,nrp)-yr)**2 +
     &             (rp(3,nrp)-zr)**2 ) *1000.
      rp(1,nrp)=xr
      rp(2,nrp)=yr
      rp(3,nrp)=zr
      return
      end ! of subr. dirpath
