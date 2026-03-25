c
c
c
c
c
c
      subroutine INPUTSED(iunit,nobs,
     &                  sta,iyr,imo,iday,ihr,imin,sec,
     &                  rmk1,rmk2,cphase,iwt,amx,prx,
     &                  xlat,xlon,xmagni,depth,origtime,
     &                  iswt,ifixsolution,ievnr,eventtype)
c
c     read inputfile (phase-list in *.SED format)
c     if NOBS=-1 ==> end of input-file detected!
c
      implicit none
c
      integer maxobsperevent
      parameter (maxobsperevent=180)
      integer nobs
      character*4 sta(maxobsperevent)
      integer*4 itime,itime1(maxobsperevent)
      real sec(maxobsperevent)
      character*1 rmk1(maxobsperevent),rmk2(maxobsperevent),
     &            cphase(maxobsperevent),clayP,clayS
      integer iunit,iwt(maxobsperevent)
      real amx(maxobsperevent), prx(maxobsperevent),
     &     xlat,xlon
      real xmagni, depth, origtime
      integer iswt,ifixsolution,ievnr
      character*1 eventtype
c
      character*80 cline
      integer iyr1(maxobsperevent),imo1(maxobsperevent),
     &        iday1(maxobsperevent),ihr1(maxobsperevent),
     &        kmin1(maxobsperevent)
      integer iyr,imo,iday,ihr,imin
      integer j,j1, jjmin,jjmin1
c
      ifixsolution=0
      nobs=0
      j=0
c
   1  read(iunit,'(a)',end=999) cline
      if(INDEX(cline,'INST').gt.0) goto 1
      if(INDEX(cline,'SED').gt.0) goto 1
      if(INDEX(cline,'BOL').gt.0) goto 1
      if(INDEX(cline,'EVENT').gt.0) goto 1
c
      j=j+1
      if(cline(1:4).eq.' ')then
         read(cline,'(17x,i1,i1,f5.2,2f8.3)') iswt,ifixsolution,
     &                                        depth,xlat,xlon
         nobs=j-1
         goto 99  ! one event finished!!!
      endif
      read(cline,'(a4,a1,a1,a1,i1,a1,5i2,f5.2,7x,
     &            f5.2,2x,a1,i1,4x,f3.0,f3.1)',err=9999)
     &            sta(j),rmk1(j),cphase(j),rmk2(j),iwt(j),clayP,
     &            iyr1(j),imo1(j),iday1(j),ihr1(j),kmin1(j),sec(j),
     &            sec(j+1),clayS,iwt(j+1),amx(j),prx(j)
      if(j.eq.1)then
         read(cline(64:64),'(a1)') eventtype
         read(cline(76:79),'(i4)') ievnr
      endif
      if(iwt(j).eq.9)then
         j=j-1
         goto 1   ! p-weight 9 not accepted !(means: only T(s-p) )
      endif
      if(clayP.eq.'m'.or.clayP.eq.'M') cphase(j)='M'
      if(cline(20:24).eq.' ') iwt(j)=5 ! NO P-arrival read !
ccc   VELEST cannot handle reflected S-phases; therefore skip them!
      if(clayS.ne.' ') cline(32:36)=' '    ! <-- set NO S-data if REFLECTED
      if(cline(32:36).ne.' ')then   ! field with observed S-time in sec
         cphase(j+1)='S'
         sta(j+1)=sta(j)
         rmk1(j+1)=' '
         rmk2(j+1)=' '
         iyr1(j+1)=iyr1(j)
         imo1(j+1)=imo1(j)
         iday1(j+1)=iday1(j)
         ihr1(j+1)=ihr1(j)
         kmin1(j+1)=kmin1(j)
         if(iwt(j).eq.5)then ! NO P arrival here; store amx&prx to S-arrival !
            amx(j+1)=amx(j)
            prx(j+1)=prx(j)
         endif
         j=j+1              ! increment J, because two phases have been added!!!
         goto 1
      else
         goto 1   ! no S-phase-data: elements J+1 will be overwritten
      endif
c
   99 continue
      do j=1,nobs
         if(iwt(j).lt.5)then
            jjmin1=j
            goto 222
         endif
      enddo
 222  j=jjmin1
      call TIMECLEAR(iyr1(j),imo1(j),iday1(j),ihr1(j),kmin1(j),sec(j),
     &               itime1(j))
      itime=itime1(j)
      jjmin=jjmin1
      do j=1,nobs
         if(iwt(j).lt.5)then
            call TIMECLEAR(iyr1(j),imo1(j),iday1(j),ihr1(j),kmin1(j),
     &                     sec(j),itime1(j))
            if(itime1(j).lt.itime)then
               itime=itime1(j)
               jjmin=j
            endif
         endif
      enddo
c
c     itime is the earliest arrival-minute;
c     adjust all seconds according to this (earliest) minute
c
  100 continue
      do j1=1,nobs
         if(iwt(j1).lt.5) sec(j1)=(itime1(j1)-itime)*60. + sec(j1)
      enddo
      origtime=sec(jjmin)
      iyr=iyr1(jjmin)
      imo=imo1(jjmin)
      iday=iday1(jjmin)
      ihr=ihr1(jjmin)
      imin=kmin1(jjmin)
c
      RETURN
c
  999 continue   ! come here, if end-of-input-file detected !
      nobs=-1
      RETURN
c
 9999 continue
      write(6,*)'INPUTSED>>> read-error! input-line is:'
      write(6,'(a)') cline
      stop'subr. INPUTSED >>> error!'
c
      end ! of subr. inputsed
