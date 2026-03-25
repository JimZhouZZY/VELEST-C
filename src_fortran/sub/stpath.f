c
c
c
c
c
c
      subroutine STPATH(xe,ye,ze,xr,yr,zr,rp,nrp,tt,nl,jl,tkj,v,d,thk,
     &                  sterr)
      implicit none
      real xe,ye,ze,xr,yr,zr,tt,tkj,sterr
      integer nrp,nl,jl
c
      integer jl1,j,j1,jb,jt
      real v(nl),d(nl),thk(nl)
c
      real rp(3,200)
      nrp=jl+1
c  assign points along raypath
      rp(1,1)=xe
      rp(2,1)=ye
      rp(3,1)=ze
c
      if (nrp.eq.2) go to 50
c
      jl1=jl-1
cuk      do 25 j=1,jl1    ! calculate up to the top and then compute an error
      do 25 j=1,jl1+1
      j1=j+1
      rp(1,j1)=xe
      rp(2,j1)=ye
      jb=jl-j+1
      rp(3,j1)=d(jb)
   25 continue
c
   50 continue
      if(nrp.gt.2)then
         sterr=sqrt( (rp(1,nrp)-xr)**2 + (rp(2,nrp)-yr)**2 +
     &               (rp(3,nrp)-zr)**2 ) *1000.
      endif
      rp(1,nrp)=xr
      rp(2,nrp)=yr
      rp(3,nrp)=zr
c
c  compute travel time
      tt=tkj/v(jl)
      if (nrp.eq.2) return
c
      do 75 j=1,jl1
      jt=jl-j
      tt=tt+thk(jt)/v(jt)
   75 continue
      return
      end ! of subr. stpath
