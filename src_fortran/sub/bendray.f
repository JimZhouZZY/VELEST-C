c
c
c
c
c
c
      subroutine BENDRAY(rp,nrp,staname,vtop,ttt)
c
      implicit none
      include '../inc/vel_com.inc'
c
      real rp(3,inrpmax)
      integer nrp
      real vtop,ttt
      character*4 staname
c
      logical hypocintop
      real rpn(3,10)
      real xydist,takeoff_angle,arrive_angle
      real deltat1,deltat2,dx1,dy1,dz1,xyz1
      real dx2,dy2,dz2,xyz2,zzz,ttt1new,xyz1n
      real ttt2new,xyz2n,ttt1old,ttt2old
      integer nrpm1,jnrp,jnrpm1
c
      if(icoordsystem.ne.2) RETURN   ! topo-array only for switzerland available
c
      xydist=
     &   SQRT( (rp(1,nrp)-rp(1,1))**2
     &        +(rp(2,nrp)-rp(2,1))**2 )
      if(xydist.gt.10.0) RETURN !epicentral-distance too big for 'airy' rays!(?)
c
      takeoff_angle=
     &   SQRT( (rp(1,2)-rp(1,1))**2
     &        +(rp(2,2)-rp(2,1))**2 )
     &   /
     &   SQRT( (rp(1,2)-rp(1,1))**2
     &        +(rp(2,2)-rp(2,1))**2
     &        +(rp(3,2)-rp(3,1))**2 )
      takeoff_angle=57.296*ASIN(takeoff_angle)
      if( (rp(3,2)-rp(3,1)) .lt. 0.0 )then
         takeoff_angle=180.-takeoff_angle     ! ray is going upwards
      endif
c
      arrive_angle=
     &   SQRT( (rp(1,nrp)-rp(1,(nrp-1)))**2
     &        +(rp(2,nrp)-rp(2,(nrp-1)))**2 )
     &   /
     &   SQRT( (rp(1,nrp)-rp(1,(nrp-1)))**2
     &        +(rp(2,nrp)-rp(2,(nrp-1)))**2
     &        +(rp(3,nrp)-rp(3,(nrp-1)))**2 )
      arrive_angle=57.296*ASIN(arrive_angle)
      if( (rp(3,2)-rp(3,1)) .lt. 0.0 )then
         arrive_angle=180.-arrive_angle     ! ray is coming upwards
      endif
c
c     write(6,'(1x,''takeoff_angle='',f6.1,''  arrive_angle='',
c    &          f6.1)') takeoff_angle, arrive_angle
c
      deltat1=0.0
      deltat2=0.0
      hypocintop=.false.
      if(rp(3,1).lt.0.0)then
         hypocintop=.true.
c        compute XYZ1: the length between the first 2 raypoints
         dx1=rp(1,2)-rp(1,1)
         dy1=rp(2,2)-rp(2,1)
         dz1=rp(3,2)-rp(3,1)
         xyz1=SQRT(dx1**2 + dy1**2 + dz1**2)
      endif
c
c     if takeoff_angle at surface < 45  --> ray too steep to be 'airy'
c
      if(hypocintop.and.takeoff_angle.le.45.0) hypocintop=.false.
      if(hypocintop)then
         continue   ! stay here; ray-start could be in the air...
      else
         if(arrive_angle.ge.135.0) RETURN ! ray too steep to be 'airy'
      endif
c
      nrpm1=nrp-1
c     compute XYZ2: the length between the last 2 raypoints
      dx2=rp(1,nrp)-rp(1,nrpm1)
      dy2=rp(2,nrp)-rp(2,nrpm1)
      dz2=rp(3,nrp)-rp(3,nrpm1)
      xyz2=SQRT(dx2**2 + dy2**2 + dz2**2)
c
c     average (3D-)ray-length in the top-layer is about 2 kilometers;
c     now split ray in top-layer into 4 parts of about equal length:
      if(hypocintop)then
         dx1=dx1/4.
         dy1=dy1/4.
         dz1=dz1/4.
      endif
      dx2=dx2/4.
      dy2=dy2/4.
      dz2=dz2/4.
      if(hypocintop)then
         do jnrp=1,4  ! first raypoint = first old raypoint!
            rpn(1,jnrp)=rp(1,1)+(jnrp-1)*dx1
            rpn(2,jnrp)=rp(2,1)+(jnrp-1)*dy1
            rpn(3,jnrp)=rp(3,1)+(jnrp-1)*dz1
         enddo
         jnrp=5   ! set last raypoint to second old raypoint to keep accuracy!
         rpn(1,jnrp)=rp(1,2)
         rpn(2,jnrp)=rp(2,2)
         rpn(3,jnrp)=rp(3,2)
      endif
      do jnrp=6,9  ! first raypoint = second-last old raypoint!
         rpn(1,jnrp)=rp(1,nrpm1)+(jnrp-6)*dx2
         rpn(2,jnrp)=rp(2,nrpm1)+(jnrp-6)*dy2
         rpn(3,jnrp)=rp(3,nrpm1)+(jnrp-6)*dz2
      enddo
      jnrp=10 ! set last raypoint to last old raypoint to keep accuracy!
      rpn(1,jnrp)=rp(1,nrp)
      rpn(2,jnrp)=rp(2,nrp)
      rpn(3,jnrp)=rp(3,nrp)
c
c     now set raypoints min. 100meters below surface (if necessary):
c
      if(hypocintop)then
         do jnrp=2,5
            call CHTOP(-rpn(1,jnrp),rpn(2,jnrp),zzz,
     &                 topo1file,topo2file)
            zzz=zzz+0.1  ! 100meters below surface
            if(rpn(3,jnrp).lt.zzz) rpn(3,jnrp)=zzz  ! set hypoc. below surface !
         enddo
      endif
      do jnrp=6,9
         call CHTOP(-rpn(1,jnrp),rpn(2,jnrp),zzz,
     &              topo1file,topo2file)
         zzz=zzz+0.1  ! 100meters below surface
         if(rpn(3,jnrp).lt.zzz) rpn(3,jnrp)=zzz   ! set hypoc. below surface !
      enddo
c
c     calculate ray-length-difference and correct traveltime:
c
      ttt1new=0.0
      if(hypocintop)then
         do jnrp=2,5
            jnrpm1=jnrp-1
            xyz1n=SQRT( (rpn(1,jnrp)-rpn(1,jnrpm1))**2
     &                 +(rpn(2,jnrp)-rpn(2,jnrpm1))**2
     &                 +(rpn(3,jnrp)-rpn(3,jnrpm1))**2 )
            ttt1new=ttt1new+xyz1n/vtop
         enddo
      endif
      ttt2new=0.0
      do jnrp=7,10
         jnrpm1=jnrp-1
         xyz2n=SQRT( (rpn(1,jnrp)-rpn(1,jnrpm1))**2
     &              +(rpn(2,jnrp)-rpn(2,jnrpm1))**2
     &              +(rpn(3,jnrp)-rpn(3,jnrpm1))**2 )
         ttt2new=ttt2new+xyz2n/vtop
      enddo
      if(hypocintop)then
         ttt1old=xyz1/vtop
         deltat1=(ttt1new-ttt1old)
      endif
      ttt2old=xyz2/vtop
      deltat2=(ttt2new-ttt2old)
      if(ABS(deltat1).gt.1e-5.or.ABS(deltat2).gt.1e-5)then
         write(6,'(1x,''BENDRAY>>> ray bended below surface!'',
     &                '' Station: '',a4)') staname
         write(6,'(1x,''dt1='',f6.3,''   dt2='',f6.3)') deltat1,deltat2
      endif
      ttt=ttt+deltat1+deltat2
c
      RETURN
      end ! of subr. bendray
