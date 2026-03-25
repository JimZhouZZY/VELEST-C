c
c
c
c
c
c
      subroutine INPUTSED_NEW(iunit,nobs,
     &                  sta,iyr,imo,iday,ihr,imin,sec,
     &                  rmk1,rmk2,cphase,iwt,amx,prx,
     &                  xlat,xlon,xmagni,depth,origtime,
     &                  iswt,ifixsolution,ievnr,eventtype,itrial)
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
      integer*4 itime,itime_o,itime1(maxobsperevent)
      real sec(maxobsperevent)
      character*1 rmk1(maxobsperevent),rmk2(maxobsperevent),
     &            cphase(maxobsperevent)
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
      integer ev_year,ev_month,ev_day,ev_hour,ev_min
      character*8 phase_id,amx_type
      integer     use_flag,amx_flag
      integer     found
      integer     itrial
c
      ifixsolution=0
      nobs=0
      j=0
c
   1  read(iunit,'(a)',end=999) cline
c
c     The following cards will be skipped
c
      if(INDEX(cline,'SED').gt.0) goto 1
      if(INDEX(cline,'BOL').gt.0) goto 1
      if(INDEX(cline,'EVENT').gt.0) goto 1
c
c     The INST card will be partially read
c
      if(INDEX(cline,'INST').gt.0) then
         read(cline,'(29x,i2)') iswt
         goto 1
      endif
c
cccccccccccccccccccccccccccc
c     Hypon input section  c 
cccccccccccccccccccccccccccc
c      read(line,*,err=900)
c     +		latr,lonr,ztr,y0tr,motr,dytr,h0tr,
c     +		m0tr,s0tr,inst
c latr:	trial latitude
c lonr:	  "   longitude
c ztr:	  "   depth
c y0tr:	  "   origin year
c motr:	  "      "   month
c dytr:	  "      "   day
c h0tr:	  "      "   hour
c m0tr:	  "      "   minute
c s0tr:	  "      "   second
c inst: 0=  free solution
c       1=  fix depth
c       2=  fix lat/lon depth
c       3=  fix origin time
c       4=  fix all
ccccccccccccccccccccccccccccccccccccccccccc
c
      if(INDEX(cline,'TRIAL').gt.0) then
         read(cline,*) xlat,xlon,depth,iyr,imo,iday,ihr,imin,
     &   origtime,ifixsolution 
	 xlon = -xlon
	 iyr = iyr - 1900
         if (ifixsolution.eq.4) then 
	    ifixsolution = 9
	 endif
         call TIMECLEAR(iyr,imo,iday,ihr,imin,origtime,itime_o)  	 
c	 
c        Reading the line after the TRIAL card
c
         read(iunit,'(a)',end=999) cline
         read(cline,'(2x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') 
     &        ev_year,ev_month,ev_day,ev_hour,ev_min 
	 read(cline(55:55),'(a1)') eventtype
	 read(cline(71:78),'(i8)') ievnr   
         goto 1
      endif
c
c     If the SKIP card is reached, the input for one event is finished
c
      if(INDEX(cline,'SKIP').gt.0) then
         nobs=j
         goto 99  ! one event finished!!!
      endif
c
c     If here, a phase will be read
c 
      j=j+1
      read(cline,'(a4,4x,a8,a1,a1,f8.3,i2,f5.2,f9.0,1x,a8,i3)',
     &     err=9999)	   
     &     sta(j),phase_id,rmk1(j),rmk2(j),sec(j),use_flag,
     &     prx(j),amx(j),amx_type,amx_flag
      call casefold(phase_id)
      call casefold(amx_type)
      call casefold(rmk1(j))
      call casefold(rmk2(j))
c
c     Treat the supported phase types
c
      found = 0
      if (phase_id.eq.'P       ') then
         found=1
	 cphase(j) = 'P'
      endif
      if (phase_id.eq.'S       ') then
         found=1
	 cphase(j) = 'S'
      endif
      if (phase_id.eq.'PMP     ') then
         found=1
	 cphase(j) = 'M'
      endif
      if (phase_id.eq.'S-P     ') then
         found=1
	 cphase(j) = '-'
      endif
c
c     If the read phase is not supported, it will be skipped
c
      if (found.eq.0) then
         j=j-1
	 goto1
      endif	 
c
c     Conversion of the weights
c
      if (rmk1(j).eq.'I') then
         iwt(j)=0
      endif	 
      if (rmk1(j).eq.'E') then
         iwt(j)=1
      endif	 
      if (rmk1(j).eq.'Q') then
         iwt(j)=2
      endif
      if (rmk1(j).eq.' ') then  ! Should never happen
         iwt(j)=5
      endif
chrm      if (use_flag.eq.0) then
chrm         j=j-1
chrm	 goto 1
chrm      endif
      if (use_flag.eq.0) then
         iwt(j) = 4
      endif
c
c     Assigning year,month,day,hour,min to each phase
c
      iyr1(j)=ev_year
      imo1(j)=ev_month
      iday1(j)=ev_day
      ihr1(j)=ev_hour
      kmin1(j)=ev_min
c
c     Read the next event
c
      goto 1
c
c     If here, all phases of one event are read in
c
c
   99 continue
c
c     Searching for the first phase with weight less than 5
c
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
      if (itrial.gt.0) then
         do j1=1,nobs
            if(iwt(j1).lt.5) sec(j1)=(itime1(j1)-itime)*60. + sec(j1)
chrm	    itime1(j1) = itime  ! Activate this, if you want to use itime1(j1)
         enddo
         origtime=sec(jjmin)
         iyr=iyr1(jjmin)
         imo=imo1(jjmin)
         iday=iday1(jjmin)
         ihr=ihr1(jjmin)
         imin=kmin1(jjmin)
      else
         do j1=1,nobs
            if(iwt(j1).lt.5)  then
	       sec(j1)=(itime1(j1)-itime_o)*60. + sec(j1)
	       if (sec(j1).lt.0) then
	          sec(j1) = 0
	       endif
	    endif
chrm	    itime1(j1) = itime  ! Activate this, if you want to use itime1(j1)
         enddo  
      endif
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
      end ! of subr. inputsed_new
