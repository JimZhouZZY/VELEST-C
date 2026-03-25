c
c
c
c
c
c
      subroutine NITTOUTPUT(damp)
c
c     output the results of iteration NITT:
c
      implicit none
      real damp
      include '../inc/vel_com.inc'
c
      real avdt,avdx,avdy,avdz,aavdt,aavdx,aavdy,aavdz
      real avelo
      integer i,j2,j1,j,k,k2,j11,j22,kj,jjj,ifl,ksta1,k1,kk1
      integer ksta2
      real cc(ist)
      character*1 reflch
c
c     output hypocenters of this iteration:
c
      avdt=0.0
      avdx=0.0
      avdy=0.0
      avdz=0.0
      aavdt=0.0
      aavdx=0.0
      aavdy=0.0
      aavdz=0.0
      write(16,*)
      write(16,'(4h  eq, 7x, 2hot, 5x, 1hx, 6x, 1hy, 6x, 1hz, 6x,3hrms,
     &           4x,5havres,''   dot     dx     dy     dz'' )')
c
      do 29 i=1,legs
c
c     print constrain-info if necessary:
c
      if(isingle.ne.0)then
         if(isconstrain(1).eq.1)then
            write(16,*)' *** nitt<2 --> ',
     &                 'depth-adjustment := 0.0 for event ',i
         endif
         if(isconstrain(2).eq.1)then
            write(16,*)' *** igap>250 --> ',
     &                 'depth-adjustment := 0.0 for event ',i
         endif
         if(isconstrain(3).eq.1)then
            write(16,*)' *** depth-adjustment constrained for event ',i
         endif
      endif
      if(iconstrain(i).eq.1)then
         write(16,'('' ***** depth constrained for event '',i5)') i
      endif
c
c   output new hypocenter-results got in this iteration:
c
      j2=4*i
      j1=j2-3
      if(i.gt.neqs)then  ! for shots
         j1=3*neqs+i
         j2=j1
      endif
      if(icoordsystem.eq.2)then
         write(16,37) i,e(1,i),-e(2,i),(e(j,i),j=3,4),rms(i),avres(i),
     &             b(j1), -b(j1+1), +b(j1+2), b(j1+3)      ! dot  dx  dy  dz
      else
         write(16,37) i,(e(j,i),j=1,4),rms(i),avres(i),
     &                (b(j),j=j1,j2)              ! dot  dx  dy  dz
      endif
 37   format (1x,i4,3x,6f7.2,2x,4f7.3)
      avdt=avdt+b(j1)
      avdx=avdx+b(j1+1)
      avdy=avdy+b(j1+2)
      avdz=avdz+b(j1+3)
      aavdt=aavdt+ABS(b(j1))
      aavdx=aavdx+ABS(b(j1+1))
      aavdy=aavdy+ABS(b(j1+2))
      aavdz=aavdz+ABS(b(j1+3))
c
  29  continue  ! loop over all events (printout)
c
      if(icoordsystem.eq.2) avdx=-avdx
      avdt=avdt/float(legs)
      avdx=avdx/float(legs)
      avdy=avdy/float(legs)
      avdz=avdz/float(legs)
      aavdt=aavdt/float(legs)
      aavdx=aavdx/float(legs)
      aavdy=aavdy/float(legs)
      aavdz=aavdz/float(legs)
      write(16,*)
      write(16,'('' A V E R A G E   of ADJUSTMENTS :'',19x,4f7.3)')
     &           avdt,avdx,avdy,avdz
      write(16,'('' A V E R A G E   of ABSOLUTE ADJUSTMENTS :'',10x,
     &           4f7.3)') aavdt,aavdx,aavdy,aavdz
      write(16,*)
c
c     all new hypocenters printed out now.
c
      if(damp.ne.1.0)then
         write(16,124) damp
124      format(/,' Step length damping of ',f7.5,' was applied.',/)
      else
         write(16,'(/,''NO step length damping applied'',/)')
      endif
c
c     print velocity adjustments and velocity-model here:
c
      if(scale(6).eq.0.0) goto 44  ! no velocity-adjustments...
      j1=4*neqs+nshot+1
      j2=j1+nltot-1
      k=0
      write(16,*)
      write(16,40)
 40   format(' Velocity adjustments:')
      write(16,41)
 41   format(5x,'vp    dvp      hp   reflector')
c
c     do it for all velocity-models
c
      do i=1,nmod
         write(16,71) i
71       format(1x,'Velocity model',i4)
         k2=laysum(i)
         j11=j1+k2-1
         j22=j11+nplay(i)-1
         kj=1
c
c        do it for all layers in this model ( = i )
c
         do jjj=j11,j22
            if(kj.eq.ireflector)reflch=reflchar
            if(kj.ne.ireflector)reflch=' '
            write(16,63) vp(i,kj),b(jjj),hp(i,kj),reflch
 63         format(1x,2f7.3,2x,f7.3,5x,a1)
            kj=kj+1
         enddo
c
c        calculate and print average velocities of the model i :
c
         ifl=1
         write(16,*)
         write(16,*)'Calculation of average velocity starts at layer #',
     &              ifl
         avelo=0
         do kj=ifl+1,nplay(i)
            avelo=avelo + ( hp(i,kj)-hp(i,kj-1) ) * vp(i,kj-1)
            write(16,'('' z ='',f6.2,'' ...'',f6.2,''    v = '',f4.2,
     &                ''    AVEL = '',f4.2,'' downto z = '',f6.2)')
     &           hp(i,kj-1),hp(i,kj),vp(i,kj-1),
     &           avelo/(hp(i,kj)-hp(i,ifl)),hp(i,kj)
         enddo
         write(16,*)
c
         ifl=2
         write(16,*)
         write(16,*)'Calculation of average velocity starts at layer #',
     &              ifl
         avelo=0
         do kj=ifl+1,nplay(i)
            avelo=avelo + ( hp(i,kj)-hp(i,kj-1) ) * vp(i,kj-1)
            write(16,'('' z ='',f6.2,'' ...'',f6.2,''    v = '',f4.2,
     &                ''    AVEL = '',f4.2,'' downto z = '',f6.2)')
     &           hp(i,kj-1),hp(i,kj),vp(i,kj-1),
     &           avelo/(hp(i,kj)-hp(i,ifl)),hp(i,kj)
         enddo
         write(16,*)
         write(16,*)
      enddo  ! loop over all models
c
c     each velocity-model and its changes printed out now.
c
c     print p- & s- station-corrections here:
c
 44   if(scale(5).eq.0.0) goto 72 ! no station-correction adjustments...
      write(16,*)
      write(16,'('' Adjusted station corrections:'')')
      write(16,48)
 48   format(2x,' stn  ptcor  dpcor  ')
      k1=4*neqs+nshot+nltot+1
      ksta1=ksta
      if(nsp.eq.2) ksta1=ksta/2
      k2=k1+ksta1-1
c
c     print p-correction for all stations:
c
      do j=1,nsta
         cc(j)=0.
         if(map1(j).eq.0) goto 49
         if(map1(j).gt.ksta1) goto 49
         kk1=k1-1+map1(j)
         cc(j)=b(kk1)             !  STATION-CORRECTION-
cccc         ptcor(j)=ptcor(j)+cc(j)  !  ADJUSTMENT !!! (p-correction)
 49      continue
      enddo
      write(16,50) (stn(j),ptcor(j),cc(j),j=1,nsta)
 50   format(4(2x,a4,2f7.3))
      if(nsp.ne.2) goto 72
      write(16,'('' Adjusted station corrections:'')')
      write(16,73)
73    format(2x,' stn  stcor  dscor  ')
      k1=4*neqs+nshot+nltot+1+ksta1
      ksta2=ksta-ksta1
      k2=k1+ksta2-1
c
c     print s-correction for all stations:
c
      do j=1,nsta
         cc(j)=0.
         if(map1(j).eq.0) goto 74
         if(map1(j).gt.ksta2) goto 74
         kk1=k1-1+map1(j)
         cc(j)=b(kk1)              ! STATION-CORRECTION-
ccc         stcor(j)=stcor(j)+cc(j)   ! ADJUSTMENT !!! (s-correction)
74       continue
      enddo
      write(16,50) (stn(j),stcor(j),cc(j),j=1,nsta)
c
72    continue
      write(16,*)
c
      RETURN
c
      end ! of subr. nittoutput
