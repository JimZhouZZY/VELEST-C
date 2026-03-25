c
c
c
c
c
c
      subroutine TRAVDERIV(raytype,nl,mll,v1,vsq1,
     &                     rp,nrp,x2,y2,z2,ss,r1,r2,ievent,inobs)
c
c     compute traveltime derivatives with respect to all the unknowns inverted
c     for !
c
      implicit none
      include '../inc/vel_com.inc'
c
      character*(*) raytype
      integer nl,mll,nrp,ievent,inobs,ii,j,jx,jndex,l,jb
      integer idownward,i
      real rp(3,inrpmax), f(inltot),zmax,r1,r2,dtdd
      real v1(nl),vsq1(nl),x2(nrp),y2(nrp),z2(nrp),ss(nrp)
c
      if(isingle.ne.0.and.iturbo.eq.1)then
         continue
      else
         do 2 ii=1,nl
2        dtdv(ii)=0.
      endif
      do 3 ii=1,3
3     dtdr(ii)=0.
c
      if(ifixsolution.eq.9) RETURN   ! hold LAT/LON/DEPTH fixed !!!
c
c     No trvdrv's for reseted obs:
c
      if(isingle.ne.0.and.w(inobs,ievent).eq.0.0) RETURN
c
      if(raytype.eq.'direct')    goto 5
      if(raytype.eq.'refracted') goto 6
      if(raytype.eq.'reflected') goto 8
      stop'TRAVDERIV>>> illegal raytype!!!'
c
c   direct ray
c
5     continue
      do 9 j=1,nrp-1
      jx=nrp-j
      x2(j)=rp(1,j+1)-rp(1,j)
      y2(j)=rp(2,j+1)-rp(2,j)
      z2(j)=rp(3,j+1)-rp(3,j)
      ss(j)=sqrt(x2(j)**2 + y2(j)**2 + z2(j)**2)
      if(isingle.ne.0.and.iturbo.eq.1)then
         continue
      else
         dtdv(jx)= -ss(j)/vsq1(jx)
      endif
9     continue
      dtdr(1)=-x2(1)/(v(jl)*ss(1))
      dtdr(2)=-y2(1)/(v(jl)*ss(1))
      dtdr(3)=-z2(1)/(v(jl)*ss(1))
      goto 40
c
c    refracted first path
c
6     continue
      do 171 j=1,nrp
      z2(j)=rp(3,j)
171   continue
      call MAXRI(nrp,z2,zmax,jndex) ! determine MAX (=zmax) of z2 (vertical
c                                     path-lengths in each layer)
      l=nl
12    if(h(l).le.(zmax+.01)) goto 11
      l=l-1
      goto 12
11    continue
      jb=l
c
      dtdd=1/v(jb)
      dtdr(1)=(r1/delta)*dtdd
      dtdr(2)=(r2/delta)*dtdd
      dtdr(3)=-sqrt(vsq(jb)-vsq(jl))/(v(jb)*v(jl))
c
      if(isingle.ne.0.and.iturbo.eq.1) goto 40  ! do not calc. dt/dv  !!!
c
c     thk(j) -- thickness of layer j
c     jl -- event layer
c     jb -- bottomming layer
c     tkj -- depth of event within event layer
c
      do 201 j=1,nl
      f(j)=1.
      if(j.ge.jl) f(j)=2.
201   if(j.gt.jb) f(j)=0.
      do 202 j=1,jb-1
      dtdv(jb)=dtdv(jb)+ thk(j)*v1(j)*f(j)/(vsq1(jb)*
     &sqrt(vsq1(jb)-vsq1(j)))
      dtdv(j)=-thk(j)*v1(jb)*f(j)/(vsq1(j)*
     &sqrt(vsq1(jb)-vsq1(j)))
202   continue
      dtdv(jl)=dtdv(jl)+tkj*v1(jb)/(vsq1(jl)*
     &sqrt(vsq1(jb)-vsq1(jl)))
      dtdv(jb)=dtdv(jb)-tkj*v1(jl)/(vsq1(jb)*
     &sqrt(vsq1(jb)-vsq1(jl))) -delta/vsq1(jb)
c
      goto 40
c
c     reflected phase
c
  8   continue
c  avoid hypocenter on layer-boundary:
chrm      if(rp(1,1).eq.rp(1,2).and.
chrm     &   rp(2,1).eq.rp(2,2).and.
chrm     &   rp(3,1).eq.rp(3,2))then
chrm            rp(3,1)=rp(3,1)-0.0001   ! move hypocenter 10 cm up
chrm      endif
c
ctest      dtdr(1)=DTDDrefl*(r1/delta)      ! valid for REFLECT
ctest      dtdr(2)=DTDDrefl*(r2/delta)
ctest      dtdr(3)=DTDHrefl
c
      jx=0
      do j=1,nrp-1
         x2(j)=rp(1,j+1)-rp(1,j)
         y2(j)=rp(2,j+1)-rp(2,j)
         z2(j)=rp(3,j+1)-rp(3,j)
         ss(j)=sqrt(x2(j)**2 + y2(j)**2 + z2(j)**2)
         if(j.eq.1)then
            dtdr(1)=-x2(1)/(v(jl)*ss(1))          ! valid for REFLECT1
            dtdr(2)=-y2(1)/(v(jl)*ss(1))
            dtdr(3)=-z2(1)/(v(jl)*ss(1))
         endif
         if(isingle.ne.0.and.iturbo.eq.1) goto 40 ! do not calc. dt/dv  !!!
         idownward=0
         if(z2(j).gt.0.0) idownward=1
         if(j.eq.(nrp-1).and.idownward.eq.1) idownward=0
         if(idownward.eq.1)then  !  DOWNwards
            jx=jl+j-1
c            if(vsq1(jx).le.0)then
c               write(6,*)'jx=',jx,'jl=',jl,'tkj=',tkj
c               write(6,*)'j=',j,'vsq1(jx)=',vsq1(jx),'nrp=',nrp
c            endif
            dtdv(jx)=-ss(j)/vsq1(jx)
         else   !  UPward
            if(jx.eq.mll)then
               dtdv(mll)=dtdv(mll)-ss(j)/vsq1(mll)
               jx=jx-1
            else
               dtdv(jx)=dtdv(jx)-ss(j)/vsq1(jx)
               jx=jx-1
            endif
         endif
      enddo
      if(jx.ne.0)stop'TRAVDERIV>>> did not reach the top...!'
      goto 40
c
c
40    continue
c
      if(ifixsolution.eq.1) dtdr(3)=0.0   ! hold depth fixed !!!
c
cek   test if nrp gt. inrpmax:::>
      if(nrp.gt.inrpmax) then
       if(.not.single_turbo)then
          write(16,2000) i,inobs,nrp
 2000  format(//,2x,'travderiv>>>, label 2000:  ',
     & 'nrp greater than inrpmax:',/,
     & 5x,'event nr=',i6,2x,'nr. of obs=',i4,5x,4hnrp=,i4,/,
     & 3x,'program stops here',/)
       endif
       stop'subr. TRAVDERIV >>> nrp > nrpmax'
      endif
c
      return
      end ! of subr. travderiv
