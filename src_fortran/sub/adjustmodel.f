c
c
c
c
c
c
      subroutine ADJUSTMODEL(damp)
c
c     adjust model-vector by the solution just obtained
c
      implicit none
      real damp
      include '../inc/vel_com.inc'
c
      integer jj,i,n,iminold,j,itime,j1,j2,k,k2,j11,j22
      integer juliam
      integer kj,jjj,k1,ksta1,kk1,ksta2
      real cc(ist)
      real dmin
c
c     adjust hypocenters of this iteration and output them:
c
      jj=0
c
      do 29 i=1,legs
      n=4
      if(i.gt.neqs) n=1   ! for shots adjust only the origin-time !!!
ccc      call TIMECLEAR(iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i),iminold)
      if(iyr(i).lt.100)then
         iminold=JULIAM(iyr(i)+1900,imo(i),iday(i),ihr(i),imin(i))
      else
         iminold=JULIAM(iyr(i),imo(i),iday(i),ihr(i),imin(i))
      endif
      write(*,'(A,I5,5F12.4)') 'DEBUG: i=',i,(e(j,i),j=1,4)
      do j=1,n
         jj=jj+1
         e(j,i)=e(j,i) + b(jj)
      enddo
      if(e(1,i).lt.0.0)then   ! change reference minute
        call TIMECLEAR(iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i),
     &                                                      itime)
        do j=1,knobs(i)
           pt(j,i)=pt(j,i)+(iminold-itime)*60.
        enddo
      endif
      do j=1,3
         isconstrain(j)=0
      enddo
      iconstrain(i)=0
c
c     constrain focal depth if necessary
c
      if(isingle.ne.0)then
         if(nitt.lt.2)then
            isconstrain(1)=1
c            write(16,*)' *** nitt<2 --> ',
c     &                 'depth-adjustment := 0.0 for event ',i
            e(4,i)=e(4,i)-b(4)
            b(4)=0.0
         endif     
         if(igap(1).gt.250)then
            dmin=999.9   ! minimum distance (epicenter --> receiver)
            do k=1,knobs(1)
               delta=SQRT(  (e(2,1)-d(k,1,1))**2
     &                    + (e(3,1)-d(k,2,1))**2  )
               if(delta.lt.dmin) then
	          dmin=delta
	       endif
            enddo
	    if (dmin.gt.15.) then
               isconstrain(2)=1
               if(iconstrain(i).eq.1) iconstrain(i)=3
c              write(16,*)' *** igap>250 --> ',
c     &                  'depth-adjustment := 0.0 for event ',i
               e(4,i)=e(4,i)-b(4)
               b(4)=0.0
	    endif
         endif
         if(abs(b(4)).gt.zadj)then
            isconstrain(3)=1
c            write(16,*)' *** depth-adjustment constrained for event ',i
            if(b(4).gt.0.0)then
               e(4,i)=e(4,i)-b(4)+zadj
               b(4)=zadj
            else
               e(4,i)=e(4,i)-b(4)-zadj
               b(4)=-zadj
            endif
         endif
      else    ! simultanous inversion
         if(i.le.neqs.and.abs(b(jj)).gt.zadj)then ! b(jj)=depth-adj. if event i
            iconstrain(i)=1
c            write(16,*)' *** depth-adjustment constrained for event ',i
            if(b(jj).gt.0.0)then
               e(4,i)=e(4,i)-b(jj)+zadj
               b(jj)=zadj
            else
               e(4,i)=e(4,i)-b(jj)-zadj
               b(jj)=-zadj
            endif
         endif
      endif
c      effdeltaz(i)=0.0
      if(itopo.gt.0)then
         if(e(4,i).lt.0.0)then  ! depth above sea-level...so near surface !
            call CHTOP(-e(2,i),e(3,i),zmin,
     &                 topo1file,topo2file) ! zmin:==surface at this point
         else
            zmin=0.0 ! depth below zero ... so zmin:==0.0 is sufficiant !
         endif
      endif
      if(e(4,i).lt.zmin.or.(ifixsolution.gt.0.and.e(4,i).le.0.0))then
c         effdeltaz(i)=b(jj)-(e(4,i)-zmin)
c
c        instead of calculating effdeltaz constrain adjustment-vector-element
c        b(jj) directly:
c
         b(jj)=b(jj)-(e(4,i)-zmin)
c
         e(4,i)=zmin
         iconstrain(i)=1
c         write(16,'('' ***** depth constrained for event '',i5)') i
      endif
c
  29  continue  ! loop over all events (adjust)
c
c     all new hypocenters adjusted now.
c
c
c  do velocity adjustments here:
c
      if(scale(6).eq.0.0) goto 44  ! no velocity-adjustments...
      j1=4*neqs+nshot+1
      j2=j1+nltot-1
      k=0
      if(.not.single_turbo)then
         write(16,*)
         write(16,40)
 40      format(' doing velocity adjustments now...')
      endif
c
c     do it for all velocity-models
c
      do 4242 i=1,nmod
      k2=laysum(i)
      j11=j1+k2-1
      j22=j11+nplay(i)-1
      kj=1
c
c     do it for all layers in this model ( = i )
c
      do 42 jjj=j11,j22
      vp(i,kj)=vp(i,kj)+b(jjj)   !   <<---- VELOCITY-ADJUSTMENT
c
c     If reflected phases are used, the velocity of the layer just above
c     the reflector is NOT allowed to have a lower velocity than any of
c     the above layers!
c     Because usually the reflector is almost at the bottom of the model,
c     namely the MOHO, we simplify the condition and don't allow any
c     low-velocity-layers at all in the whole model (of course  o n l y
c     if reflected phases are used in this VELEST-run !)
      if(lowveloclay.eq.0)then
         if(kj.gt.1)then
            if(vp(i,kj).lt.vp(i,kj-1))then
               i f (.not.single_turbo)then
               write(16,'(1x,''WARNING: Tried to introduce a '',
     &                    ''low-velocity-layer! (Layer '',i2,'')'')') kj
               write(16,'(1x,''Setting DVP from '',f5.2,'' to 0.0 and'',
     &                    '' VP to vp(layer_above)+0.001'')') b(jjj)
               endif
               write(6,'(1x,''WARNING: Tried to introduce a '',
     &                    ''low-velocity-layer! (Layer '',i2,'')'')') kj
               write(6,'(1x,''Setting DVP from '',f5.2,'' to 0.0 and'',
     &                    '' VP to vp(layer_above)+0.001'')') b(jjj)
               write(6,*)
               b(jjj)=0.0
               vp(i,kj)=vp(i,kj-1)+0.001
            endif
         endif
      endif
c
      kj=kj+1
42    continue
c
4242  continue
c
c     each velocity-model adjusted now.
c
c     adjust p- & s- station-corrections here:
c
 44   if(scale(5).eq.0.0) goto 72 ! no station-correction adjustments...
      if(.not.single_turbo)then
         write(16,*)
         write(16,'('' doing station-correction adjustments...'')')
         write(16,*)
      endif
      k1=4*neqs+nshot+nltot+1
      ksta1=ksta
      if(nsp.eq.2) ksta1=ksta/2
      k2=k1+ksta1-1
c
c     do p-correction for all stations:
c
      do j=1,nsta
         cc(j)=0.
         if(map1(j).eq.0) goto 49
         if(map1(j).gt.ksta1) goto 49
         kk1=k1-1+map1(j)
         cc(j)=b(kk1)             !  STATION-CORRECTION-
         ptcor(j)=ptcor(j)+cc(j)  !  ADJUSTMENT !!! (p-correction)
 49      continue
      enddo
      if(nsp.ne.2) goto 72
      k1=4*neqs+nshot+nltot+1+ksta1
      ksta2=ksta-ksta1
      k2=k1+ksta2-1
c
c     do s-correction for all stations:
c
      do j=1,nsta
         cc(j)=0.
         if(map1(j).eq.0) goto 74
         if(map1(j).gt.ksta2) goto 74
         kk1=k1-1+map1(j)
         cc(j)=b(kk1)              ! STATION-CORRECTION-
         stcor(j)=stcor(j)+cc(j)   ! ADJUSTMENT !!! (s-correction)
74       continue
      enddo
c
72    continue
c
      RETURN
c
      end ! of subr. adjustmodel
