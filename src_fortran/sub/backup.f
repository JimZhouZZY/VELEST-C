c
c
c
c
c
c
      subroutine BACKUP
c
c     go back in direction of the previous solution by doing
c     only half of the adjustments calculated in the last iteration
c
      implicit none
      include '../inc/vel_com.inc'
      integer jjj,i,n,k,iccc,k2,nl,j1,j2,ifl,k1,ksta1,j,kk1,m
      integer ksta2
      real zzz,avelo
c
      real cc(ist)
      character*1 reflch
c
      jjj=0
      do 2 i=1,legs
c---------- hypocenter backup
      n=4
c
      if(i.gt.neqs) n=1
      do 2 k=1,n
      jjj=jjj+1
      b(jjj)=b(jjj)/2.0
      if(k.eq.4)then                      ! backup depth
c
c     concept of effdeltaz no longer in use!
c
c         if(effdeltaz(i).eq.0.0)then
            e(k,i)=e(k,i)-b(jjj)          ! depth not constrained
c         else
c            effdeltaz(i)=effdeltaz(i)/2.0
c            e(k,i)=e(k,i)-effdeltaz(i)           !  depth was constrained:
c         endif                                   ! go back only half of the
         if(itopo.gt.0.and.e(k,i).lt.0.0)then
c           depth above zero... is it also above surface?
            call CHTOP(-e(2,i),e(3,i),zzz,
     &                 topo1file,topo2file)
            if(e(k,i).lt.zzz) e(k,i)=zzz  ! set hypocenter down to surface !
            if(ifixsolution.gt.0) e(k,i)=zzz ! fix depth to min_depth allowed!!
         endif
      else                                       ! effective depth-change !
            e(k,i)=e(k,i)-b(jjj)
      endif
   2  continue
ccc      iccc=mod((nitt-1),invertratio)
      iccc=mod(nitt,invertratio)
      if(iccc.ne.0) goto 510
c---- velocity readjustments for each model
      do 26 k2=1,nmod
      nl=nplay(k2)
      j1=4*neqs+nshot+laysum(k2)
      j2=j1-1+nl
      k=0
      if(.not.single_turbo)then
         write(16,3)
3        format(1x,'Velocity readjustments:')
         write(16,27) k2
27       format(1x,'Velocity model',i4)
      endif
      do 4 jjj=j1,j2
      k=k+1
      b(jjj)=b(jjj)/2.
      vp(k2,k)=vp(k2,k)-b(jjj)
      if(lowveloclay.eq.0)then
         if(k.gt.1)then
            if(vp(k2,k).lt.vp(k2,k-1))then
               i f (.not.single_turbo)then
               write(16,*)'WARNING:'
               write(16,'(1x,''WARNING: Tried to introduce a '',
     &                    ''low-velocity-layer! (Layer '',i2,'')'')') k
               write(16,'(1x,''Setting DVP from '',f5.2,'' to 0.0 and'',
     &                    '' VP to vp(layer_above)+0.001'')') b(jjj)
               e n d i f
               write(6,'(1x,''WARNING: Tried to introduce a '',
     &                    ''low-velocity-layer! (Layer '',i2,'')'')') k
               write(6,'(1x,''Setting DVP from '',f5.2,'' to 0.0 and'',
     &                    '' VP to vp(layer_above)+0.001'')') b(jjj)
               write(6,*)
               b(jjj)=0.0
               vp(k2,k)=vp(k2,k-1)+0.001
            endif
         endif
      endif
5     format(1x,3f7.3,3x,a1)
      if(k.eq.ireflector)reflch=reflchar
      if(k.ne.ireflector)reflch=' '
      if(.not.single_turbo)then
         write(16,5) vp(k2,k),b(jjj),hp(k2,k),reflch
      endif
   4  continue
c
c    calculate and print average velocities of the model k2 :
c
      i f (.not.single_turbo) t h e n
      ifl=1
      write(16,*)
      write(16,*)'Calculation of average velocity starts at layer # ',
     &           ifl
      avelo=0
      do k=ifl+1,nplay(k2)
         avelo=avelo + ( hp(k2,k)-hp(k2,k-1) ) * vp(k2,k-1)
         write(16,'('' z ='',f6.2,'' ...'',f6.2,''    v = '',f4.2,
     &             ''    AVEL = '',f4.2,'' downto z = '',f6.2)')
     &       hp(k2,k-1),hp(k2,k),vp(k2,k-1),avelo/(hp(k2,k)-hp(k2,ifl)),
     &       hp(k2,k)
      enddo
      write(16,*)
c
      ifl=2
      write(16,*)
      write(16,*)'Calculation of average velocity starts at layer # ',
     &           ifl
      avelo=0
      do k=ifl+1,nplay(k2)
         avelo=avelo + ( hp(k2,k)-hp(k2,k-1) ) * vp(k2,k-1)
         write(16,'('' z ='',f6.2,'' ...'',f6.2,''    v = '',f4.2,
     &             ''    AVEL = '',f4.2,'' downto z = '',f6.2)')
     &       hp(k2,k-1),hp(k2,k),vp(k2,k-1),avelo/(hp(k2,k)-hp(k2,ifl)),
     &       hp(k2,k)
      enddo
      write(16,*)
      write(16,*)
      e n d i f
c
26    continue
510   if(nsinv.eq.0) goto 511
ccc      if(mod((nitt-1),invertratio).ne.0) goto 511
      if(mod(nitt,invertratio).ne.0) goto 511
      k1=4*neqs+nshot+nltot+1
      ksta1=ksta
      if(nsp.eq.2) ksta1=(ksta/2)
      k2=k1+ksta1-1
      do 60 j=k1,k2
60    b(j)=b(j)/2.
      do 7 j=1,nsta
      cc(j)=0.0
      if(map1(j).eq.0) goto 7
      if(map1(j).gt.ksta1) goto 7
      kk1=k1-1+map1(j)
      cc(j)=b(kk1)
      ptcor(j)=ptcor(j)-cc(j)
7     continue
      if(.not.single_turbo)then
         write(16,8) (stn(m),ptcor(m),cc(m),m=1,nsta)
8        format(5(2x,a4,2f7.3))
 6       format(1x,'P correction readjustments:')
         write(16,*)
         write(16,*)'Half adjustments made'
      endif
      if(nsp.ne.2) goto 29
      k1=4*neqs+nshot+nltot+ksta1+1
      ksta2=ksta-ksta1
      k2=k1+ksta2-1
      do 61 j=k1,k2
61    b(j)=b(j)/2.
      if(.not.single_turbo)then
         write(16,39)
39       format(1x,'S correction readjustments:')
      endif
      do 28 j=1,nsta
      cc(j)=0.0
      if(map1(j).eq.0) goto 28
      if(map1(j).gt.ksta2) goto 28
      kk1=k1-1+map1(j)
      cc(j)=b(kk1)
      stcor(j)=stcor(j)-cc(j)
28    continue
      if(.not.single_turbo)then
         write(16,8) (stn(m),stcor(m),cc(m),m=1,nsta)
      endif
511   continue
      if(.not.single_turbo)then
         write(16,*)
         write(16,*)'Half adjustments made'
      endif
29    continue
c
      return
      end ! of subr. backup
