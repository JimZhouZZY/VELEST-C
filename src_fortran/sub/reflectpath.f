c
c
c
c
c
c
      subroutine REFLECTPATH(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,
     &jl,tkj,ain,mll,rp,nrp,reflerr)
c
c     Urs Kradolfer, Nov. 1986
c
c       reflectpath assigns raypoints with coordinates rp(1,i),rp(2,i),
c  rp(3,i) along a reflected seismic raypath in a flat layered earth.
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
c               ain - sine of takeoff angle (with respect to downward vertical)
c               mll - reflector is bottom of layer MLL
c
c  output:  rp(1,2,or 3,i) - raypoint coordinates
c                      nrp - number of raypoints
c
c max. nr of raypoints =50
c
      implicit none
c
      real xe,ye,ze,xr,yr,zr,delta,tkj,ain,reflerr
      integer nl,jl,mll,nrp
c
      real d1,d2,salpha,tng,deljl
      integer m,ld,j,i
      real v(nl),thk(nl),vsq(nl),rp(3,200)
c
      d1=(xr-xe)/delta
      d2=(yr-ye)/delta
      nrp=jl+1+(mll-jl)*2+1
c
      salpha=ain
c
c   assign raypoint coordinates
c
      rp(1,1)=xe
      rp(2,1)=ye
      rp(3,1)=ze
c      tng=tan(asin(salpha))  ! compiler produces warning with range checking
      tng=asin(salpha)
      tng=tan(tng)       ! must be defined here; in case ld is zero
      deljl=tng*(thk(jl)-tkj) !horiz. travel distance in ev.layer
      rp(1,2)=rp(1,1)+deljl*d1
      rp(2,2)=rp(2,1)+deljl*d2
      rp(3,2)=rp(3,1)+(thk(jl)-tkj)
c---- now arrived at bottom of event-layer
      m=jl
      ld=mll-jl    ! nr of layer to go down to reflector
      j=2
      do i=1,ld
         j=j+1
         m=jl+i
         tng=v(m)*salpha/sqrt(vsq(jl)-vsq(m)*salpha**2)
         rp(1,j)=rp(1,j-1)+thk(m)*tng*d1
         rp(2,j)=rp(2,j-1)+thk(m)*tng*d2
         rp(3,j)=rp(3,j-1)+thk(m)       ! ray is going downward
      enddo
c---- now arrived down at reflector; going upward from now on
      if(m.ne.mll)stop'subr. REFLECTPATH >>> Error in ray-path...'
      j=j+1
      deljl=tng*thk(m)   ! horiz. travel-distance in lowest layer
      rp(1,j)=rp(1,j-1)+deljl*d1
      rp(2,j)=rp(2,j-1)+deljl*d2
      rp(3,j)=rp(3,j-1)-thk(m)
c---- now arrived at top of lowest layer
c---- going up to top of model:
c     do i=j+1,nrp-1   if one would raytrace only to the second-last point
c in order to calculate an error, we calculate up to the receiver:
      do i=j+1,nrp
         m=m-1
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
      enddo
c
      reflerr=sqrt( (rp(1,nrp)-xr)**2 + (rp(2,nrp)-yr)**2 +
     &              (rp(3,nrp)-zr)**2 ) *1000.
cc      write(6,*)'REFLECTPATH: raytracer-error was : ',reflerr,' meter'
c---- now assign the correct receiver-coordinates to the last raypoint:
      rp(1,nrp)=xr
      rp(2,nrp)=yr
      rp(3,nrp)=zr
      return
      end ! of subr. reflectpath
