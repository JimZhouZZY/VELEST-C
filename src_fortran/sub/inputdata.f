c
c
c
c
c
c
      subroutine INPUTDATA(i)  ! old name:  INPUT2
c
c     reads in the phaselists with the earthquake and/or shot data
c
      implicit none
      include '../inc/vel_com.inc'
c
      integer     i
      integer     nobs,ntobs,jj,nobsread,i1,i2,i3,icc,jjmin1,jk,
     &            jjmin,jshot,j,k,iunit,ll,itest,ie,iphaseteststopflag,
     &            l,nobsp,nobss
      real  wsum,xlat,xlon,alon,depth,z,ss1
      real  xxlat,aalon,xxx
      character*1 sc,ss,cphase(maxobsperevent)
      character*1 rmk1(maxobsperevent), rmk2(maxobsperevent),
     &            eventtype, cns,cew
      character*4 sta(maxobsperevent)
      integer ipwt(maxobsperevent)
      real sec(maxobsperevent)
c
      SAVE iphaseteststopflag  ! forget never that error was detected!
c
      data sc,ss/'s','S'/
c
c     input phase list and set-up initial trial hypocenter
      iunit=8
      if(i.eq.1)then
         if(neqs.gt.0) open(iunit,file=phasefile,
cVMS     &                      status='old',err=998,readonly)
     &                      status='old',err=998)
         if(nshot.gt.0) open(9,file=shotfile,
cVMS     &                       status='old',err=999,readonly)
     &                       status='old',err=999)
      endif
cek
cek      re-entry for reading next event if first event has less
cek      than three observations
10001 continue
cek
      nobs=0
      nobswithw0=0
      ntobs=0
      wsum= 0.0
c
      if(i.le.neqs)then
         iunit=8    ! EQS
      else
         iunit=9    ! SHOTS
      endif
c
      xlat=0.0
      xlon=0.0
      do jj=1,maxobsperevent
         cphase(jj)=' '
      enddo
c
cek  i1-switch is dummy
      i1=0
c
cek     EQS:   i1=ifx(i)     ! no longer in use !!!
cek     SHOTS: i1=icc
c
      if(ised.eq.0)then
         call INPUTCNV(iunit,nobsread,
     &                 sta,iyr(i),imo(i),iday(i),ihr(i),imin(i),sec,
     &                 rmk1,rmk2,cphase,ipwt,amx,prx,
     &                 xlat,alon,emag(i),depth,e(1,i),
     &                 i1,i2,i3,eventtype)
      endif
      if(ised.eq.1)then
         call INPUTARCVEL(iunit,nobsread,
     &                 sta,iyr(i),imo(i),iday(i),ihr(i),imin(i),sec,
     &                 rmk1,rmk2,cphase,ipwt,amx,prx,
     &                 xlat,alon,emag(i),depth,e(1,i),
     &                 i1,i2,i3,eventtype)
      endif
      if(ised.eq.2)then
         call INPUTSED_NEW(iunit,nobsread,
     &                 sta,iyr(i),imo(i),iday(i),ihr(i),imin(i),sec,
     &                 rmk1,rmk2,cphase,ipwt,amx,prx,
     &                 xlat,alon,emag(i),depth,e(1,i),
     &                 i1,ifixsolution,i3,eventtype,itrial)
      endif
      if(nobsread.eq.-1)then
         knobs(i)=-1       ! means: end of input-file detected!
         RETURN
      endif
      if(.NOT.(ised.eq.0.or.ised.eq.1.or.ised.eq.2))
     & stop'INPUTDATA: ised flag has no supported value!'
cek
cek check for events with less than three observations
cek
      if(isingle.ne.0) then
        if(nobsread.lt.3) then
            write(6,'(1x,''Event #'',i4)') i
            write(6,*) ' skipped because it has fewer than 3 obs.'
           goto 10001
        endif
      endif
      if(isingle.ne.0)then
         write(2,'(''1 E V E N T   N R .   '',i6,
     &             ''                 '',
     &             ''           0                    0'')')
     &             isingle
         write(6,'(1x,''Event #'',i6)')
     &            isingle
      endif
c
      if(i.le.neqs)then
cEK Dez94  next statement put in effect (ek: i1 set to zero, see above)
         ifx(i)=i1
cEK Dez94 next statement put out of use:
cEK         ifx(i)=0   ! NO ROTATION; do not set dtdr(2) to 0.0  UK87
      else
         icc=i1
      endif
c
c     if ITRIAL = 1 take first station as trial epicenter and
c                   ztrial as trial depth
c     Do the same, if no hypocenter has been read from input !!!
c
      if(itrial.gt.0.or.
     &      (xlat.eq.0.0.and.alon.eq.0.0))then
c
c        find first station:
c
         if(itrial.gt.0)then
            jjmin1=1
            do jj=1,nobsread
               if(ipwt(jj).lt.5)then
                  do jk=1,nsta
                     if(sta(jj).eq.stn(jk)) goto 22222
                  enddo
                  goto 2222   ! jj-th station not on station-list !!!
22222             jjmin1=jj
                  goto 222
               endif
2222           continue
            enddo
 222        jjmin=jjmin1  ! first station in data which is on station-list!!
            do jj=jjmin1,nobsread
               if(ipwt(jj).lt.5)then
                  if(sec(jj).lt.sec(jjmin1))then
                     jjmin=jj
                  endif
               endif
            enddo
cc          write(16,'('' first station is '',a4)') sta(jjmin)
            jjmin1=0
            do jj=1,nsta
               if(sta(jjmin).eq.stn(jj))then
                  jjmin1=jj
                  goto 11111
               endif
            enddo
         endif
11111    continue    ! jj is the first station
c        If no first station found, take 1st in station-list !! :
cuk         if(smn(jj,i).eq.' ') jj=1  ! is wrong!!
         if(jjmin1.eq.0) jj=1
         if(ifixsolution.ne.9)then
            xlat=xla(jj)+0.001
            alon=xlo(jj)+0.001
            if(ifixsolution.ne.1)then
               depth=ztrial
            else
               if(.not.single_turbo) write(16,*)'DEPTH fixed !!!'
               write(6,*) 'DEPTH fixed !'
            endif
         else
            if(.not.single_turbo) write(16,*)'HYPOCENTER fixed !!!'
            write(6,*)'HYPOCENTER fixed !'
            if(icoordsystem.eq.2.and.alon.gt.0.) alon=-alon ! fixed LAT
         endif                                 ! was given in LON E
      endif
c
c     in case the analyst has fixed the depth to less/equal 0.0, he probably
c     wanted to fix it at the surface... set depth to 3km above sea-level;
c     program VELEST will set it properly to the surface!
c
      if(ifixsolution.ne.0.and.itopo.gt.0.and.depth.le.0.0) depth=-3.0
      if(i.gt.neqs)then
         jshot=i-neqs       ! SHOTS
         map2(jshot)=icc
         z=depth
      else
         z=depth+zshift     ! EQS
      endif
c
c     now do transformation 'LAT/LON --> Xkm/Ykm' for trial hypocenter:
c      write(*,'(A,I0)') '[DEBUG] icoordsystem = ', icoordsystem
c
      if(icoordsystem.eq.2)then
         call GEOKO(e(2,i),e(3,i),xlat,-alon,-1) ! calc. cart. coord.
         e(2,i)=-e(2,i)
      else
         call SDC(e(2,i),e(3,i),xlat,alon,-1) ! calc. cart. coord.
      endif
      e(4,i)=z
      depthsofinput(i)=z
c
      do 15 j=1,nobsread
c
      if(isingle.eq.0)then
         if(ipwt(j).ge.4)then     ! do not accept phase-weights >= 4  !!!
            goto 15
         endif
      else    
         if(ipwt(j).ge.6)then     ! do not accept phase-weights >= 6  !!!
            goto 15               ! but take 4's & 5's for magnitude-calculation
         endif                    ! and 4's also for residual-computation !!!
      endif
      if(nsp.eq.1.and.cphase(j).eq.'s') goto 15
      if(nsp.eq.1.and.cphase(j).eq.'S') goto 15
      if(nsp.eq.1.and.cphase(j).eq.'-') goto 15 ! skip also s-p phases
c
c     if cphase(j).ne. 's' or 'S', a P-phase is assumed, so it's not necessary
c     to check, whether it's a real P or a reflected P ( 'm' or 'M').
c     But if no reflections should be assumed (e.g. ireflector.eq.0),
c     do NOT use any phases marked with 'm' or 'M' !
      if(cphase(j).ne.'p'.and.cphase(j).ne.'P'.and.ipwt(j).lt.5)then
         if(cphase(j).ne.'s'.and.cphase(j).ne.'S')then
            if(cphase(j).ne.'-')then
               if(cphase(j).ne.'m'.and.cphase(j).ne.'M')then	    
                  if(.not.single_turbo)then
                     write(16,*)'WARNING:'
                     write(16,*)'what phase is this ?  ',
     &		          cphase(j),' ???'
		  endif   
                  write(6,*)'what phase is this ?  ',cphase(j),' ???'
                  write(6,*)
                  if(isingle.gt.0)then
                     write(2,'('' DELETED: '',a4,
     &                      '' unknown phase is: '',a1)')
     &                      sta(j),cphase(j)
                     goto 15
               endif  ! if .ne. M
               else ! here, if phase is a reflected one!
                  if(ireflector.eq.0)then
                     i f (.not.single_turbo) t h e n
                        write(16,*)'WARNING:'
                        write(16,*)'subr. INPUTDATA >>> Phase is : ',
     &		             cphase(j)
                        write(16,*)'but ireflector is: ',ireflector
                        write(16,*)'Phase therefore neglected !!'
                     e n d i f  ! single_turbo
                     write(6,*)'subr. INPUTDATA >>> Phase is : ',
     &		          cphase(j)
                     write(6,*)'but ireflector is: ',ireflector
                     write(6,*)'Phase therefore neglected !!'
                     write(6,*)
                     goto 15
                  endif  ! if ireflector
               endif  ! else
	    endif ! if s-p   
         endif  ! if s
      endif  ! if p
c
      do 16 k=1,nsta
      if(sta(j).eq.stn(k)) goto 17
   16 continue
      if(.not.single_turbo)then
         if(isingle.eq.0) write(16,*)'WARNING:     Event # ',i
         if(isingle.gt.0) write(16,*)'WARNING:     Event # ',isingle
         write(16,'('' WARNING:  Station: >>>'',a4,
     &              ''<<< not found in stationlist!'')') sta(j)
         write(16,*)'Phase therefore skipped'
      endif
cc      write(6,*)'Event # ',i
      write(6,'('' WARNING:  Station: >>>'',a4,
     &           ''<<< not found in stationlist!'')') sta(j)
      write(6,*)'Phase therefore skipped'
      write(6,*)
      if(isingle.gt.0) write(2,'('' DELETED: '',a4,
     &                           '' not on station-list'')') sta(j)
      goto 15
17    continue
      ss1=sqrt( (x(k,1)-e(2,i))**2 + (x(k,2)-e(3,i))**2)
      if(ss1.gt.dmax)then
         if(.not.single_turbo)then
            write(16,*)'WARNING:'
            write(16,'('' epicentral distance:'',f6.1,
     &                 '' > dmax ('',f6.1,'') ==> skipping phase !'')')
     &                 ss1,dmax
         endif
         if(isingle.ne.0)then
            write(6,'('' epicentral distance:'',f6.1,
     &                 '' > dmax ('',f6.1,'') ==> skipping phase !'')')
     &                 ss1,dmax
         write(2,'('' DELETED: '',a4,
     &             '' epicentral-distance too large'')') sta(j)
         endif
         goto 15
      endif
c
c     test for only one p-reading of same station per event
c
      if(ipwt(j).eq.5) goto 888  ! dont test for "weight 5" phases
c      
      do ll=1,j-1
         itest=9
         i f (nsp.eq.1.and.(cphase(j).eq.'s'.or.cphase(j).eq.'S'))then
            continue
         e l s e
         if(k.eq.istm(ll,i).and.ipwt(ll).lt.4)then
            if(cphase(j).eq.'p')itest=0
            if(cphase(j).eq.'P')itest=0
            if(cphase(j).eq.'s')itest=1
            if(cphase(j).eq.'S')itest=1
            if(cphase(j).eq.'m')itest=-1
            if(cphase(j).eq.'M')itest=-1
            if(cphase(j).eq.'-')itest=2
            if(sphase(ll,i).eq.itest)then  ! twice the same phase in this event
               if(isingle.eq.0) ie=i
               if(isingle.gt.0) ie=isingle
               if(.not.single_turbo)then
                 write(16,*)'WARNING:'
                 write(16,*)'PHASETEST: POSSIBLE ERROR in phaselist !!'
                 write(16,'(1x,''---> '',3i2.2,1x,2i2.2)')
     &                    iyr(i),imo(i),iday(i),ihr(i),imin(i)
                 write(16,*)'Event=',ie,' Obs-nr. = ',j,' >>> Station ',
     &                      sta(j),' & Phase = ',cphase(j),
     &                      ' already occured!'
               endif
               write(6,*)'PHASETEST: POSSIBLE ERROR in phaselist !!'
               write(6,'(1x,''---> '',3i2.2,1x,2i2.2)')
     &                  iyr(i),imo(i),iday(i),ihr(i),imin(i)
                 write(6,*)'Event=',ie,' Obs-nr. = ',j,' >>> Station ',
     &                      sta(j),' & Phase = ',cphase(j),
     &                      ' already occured!'
               write(6,*)
               if(isingle.eq.0)then
                  if(.not.single_turbo)then
                     write(16,*)'subr. INPUTDATA >>> program will stop'
                     write(16,*)
                  endif
                  iphaseteststopflag=1
               else
                  if(.not.single_turbo)then
                     write(16,*)'INPUTDATA>>> nevertheless, '
     &                           //'program continues'
                     write(16,*)
                  endif
               endif
c              iphaseteststopflag is SAVED at beginning of subr. inputdata
c             ! forget never that error was detected!
c                 stop'PHASE-TEST: FATAL ERROR in phaselist !'
            endif
         endif
         e n d i f
      enddo
c
c     arrive here, if station on stationlist and phase is accepted by all tests:
c
 888  nobs=nobs+1
      amx(nobs)=amx(j)    
      prx(nobs)=prx(j)    ! added 18.9.91 / uk
      if(isingle.gt.0)then
         if(ised.eq.2)then
            prmk(nobs,1)=rmk1(j)
            prmk(nobs,2)=rmk2(j)
         else
            prmk(nobs,1)=' '
            prmk(nobs,2)=' '
         endif
      endif
      do 66 l=1,3
66    d(nobs,l,i)=x(k,l)
      pt(nobs,i)=sec(j)
      kpwt(nobs,i)=ipwt(j)
      istm(nobs,i)=k
      smn(nobs,i)=sta(j)
c***    kpwt(nobs,i) - p or s weight for observation nobs, event i
c***    istm(nobs,i) - station number for observation
c***    smn(nobs,i) - station name for observation
c**    iphase(nobs,i) gives the model number for the observation
      if(nsp.eq.2) goto 21
c
c
        iphase(nobs,i)=model(1)
c
c***     sphase(nobs,i) is 0. for a p observation
c***                       1. for an s observation
cuk*                      -1. for a reflected P observation
chrm                       2. for a S-P phase
c
      sphase(nobs,i)=0.
      if(cphase(j).eq.sc.or.cphase(j).eq.ss) sphase(nobs,i)=1.0
      if(cphase(j).eq.'m'.or.cphase(j).eq.'M') sphase(nobs,i)=-1.0
      if(cphase(j).eq.'-') sphase(nobs,i)=2.0
c
c *** w(nobs,i) is an observation weighting factor :
c
      w(nobs,i)=1.0/(2**(ipwt(j)*2))
c
      if(cphase(j).eq.sc.or.cphase(j).eq.ss.or.cphase(j).eq.'-')
     & w(nobs,i)=swtfac*w(nobs,i)
      if(ipwt(j).gt.4)then
         w(nobs,i)=0.0    ! necessary for single-event-mode !!!
         nobswithw0=nobswithw0+1
      endif
      wsum=wsum+w(nobs,i)
      goto 15
21    continue    !nsp=2
      if(cphase(j).eq.sc.or.cphase(j).eq.ss) goto 22 ! s-phase
      if(cphase(j).eq.'-') goto 22  ! s-p phase
      sphase(nobs,i)=0.
      if(cphase(j).eq.'m'.or.cphase(j).eq.'M') sphase(nobs,i)=-1.0 ! refl. P
      iphase(nobs,i)=model(1)      
c      iphase(nobs,i)=model(2*k-1)
      w(nobs,i)=1.0/(2**(ipwt(j)*2))
      if(ipwt(j).ge.4)then
         w(nobs,i)=0.0    ! necessary for single-event-mode !!!
         nobswithw0=nobswithw0+1
      endif
      wsum=wsum+w(nobs,i)
      goto 15
22    continue
      sphase(nobs,i)=1.
      if(cphase(j).eq.'-') sphase(nobs,i)=2.0 ! s-p phase
c
c     For s-p phases the model parameter is set to the s-wave model
c
c
      iphase(nobs,i)=model(2)
      w(nobs,i)=swtfac*1.0/(2**(ipwt(j)*2))
      if(ipwt(j).ge.4)then
         w(nobs,i)=0.0    ! necessary for single-event-mode !!!
         nobswithw0=nobswithw0+1
      endif
      wsum=wsum+w(nobs,i)
15    continue
c
c     one event is read-in now :
c
 14   continue
      knobs(i)=nobs
cek list P and S obs separateley:
      if(nsp.eq.2) then
         nobsp=0
         nobss=0
         do j=1,nobs
            if(sphase(j,i).eq.0.0) nobsp=nobsp+1
            if(sphase(j,i).eq.1.0) nobss=nobss+1
         enddo
      endif
      if(.not.single_turbo)then
         if(xlat.lt.0.0)then
            cns='S'
            xxlat=-xlat
         else
            cns='N'
            xxlat=xlat
         endif
         if(alon.lt.0.0)then
            cew='E'
            aalon=-alon
         else
            cew='W'
            aalon=alon
         endif
         if(icoordsystem.eq.2)then
            xxx=-e(2,i)
         else
            xxx=e(2,i)
         endif
         if(i.le.neqs)then
cek print P and S number of obs
           if(nsp.eq.2) then
            write(16,'(1x,i3,1x,3i2.2,1x,2i2.2,1x,f5.2,1x,f7.4,a1,
     &               1x,f8.4,a1,1x,f6.2,3f7.2,f5.2,i2,i4,3x,2i5)')
     &      i,iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i),
     &      xxlat,cns,aalon,cew,depth,xxx,(e(j,i),j=3,4),
     &      emag(i),ifx(i),knobs(i)-nobswithw0,nobsp,nobss
           else
            write(16,'(1x,i3,1x,3i2.2,1x,2i2.2,1x,f5.2,1x,f7.4,a1,
     &                 1x,f8.4,a1,1x,f6.2,3f7.2,f5.2,i2,i4)')
     &      i,iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i),
     &      xxlat,cns,aalon,cew,depth,xxx,(e(j,i),j=3,4),
     &      emag(i),ifx(i),knobs(i)-nobswithw0
           endif

         endif
         if(i.gt.neqs)then  !shots
            write(16,'(1x,i3,1x,3i2.2,1x,2i2.2,1x,f5.2,1x,f7.4,a1,
     &                1x,f8.4,a1,1x,f6.2,3f7.2,f5.2,i2,i4)')
     &      i,iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i),
     &      xxlat,cns,aalon,cew,depth,xxx,(e(j,i),j=3,4),
     &      emag(i),map2(i-neqs),knobs(i)-nobswithw0
         endif
      endif
c
c      No more active (format 1607)
c 1607 format(1x,i3,1x,3i2.2,1x,2i2.2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,1x,
c     &       f6.2,3f7.2,f5.2,i2,i4)
c
c     normalize weights:
c
      if(isingle.gt.0)then
         if( (knobs(i)-nobswithw0) .lt. nvar ) RETURN
         if(wsum.le.0.0) RETURN
      endif
      do j=1,nobs
         w(j,i)=w(j,i)*(nobs-nobswithw0)/wsum
      enddo
c
      if(iphaseteststopflag.eq.1.and.i.eq.(neqs+nshot))then
         stop'INPUTDATA >>> PHASE-TEST: FATAL ERROR in phaselist !'
      endif
      RETURN
c
80    continue
      if(.not.single_turbo)then
         write(16,85)
c---- read past end of data
   85    format(' ***** end of data encountered *****')
      write(16,*)'WARNING:  subr. INPUTDATA >>> end of data encountered'
      endif
      write(6,*)'WARNING:'
      stop'subr. INPUTDATA >>> error: end of data!'
c
  998 call OPENERROR('inputdata','EQ-data-input-file FOR008')
  999 call OPENERROR('inputdata','SHOT-data-input-file FOR009')
      return
c
      end ! of subr. inputdata
