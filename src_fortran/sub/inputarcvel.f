c
c
c
c
c
c
      subroutine INPUTARCVEL(iunit,nobs,
     &                  sta,iyr,imo,iday,ihr,imin,sec,
     &                  rmk1,rmk2,cphase,iwt,amx,prx,
     &                  xlat,xlon,xmagni,depth,origtime,
     &                  i1,i2,i3,eventtype)
c
c  implemented by EK 9.11.93 to read velest-archive type input
c  data for single event location mode , replacing older
c  input routine for SED-format. called by ised=1
c
c     on Output:
c     ---------
c     EQS:   i1=ifx(i)   ! no longer in use !!!
c     SHOTS: i1=icc
c
      implicit none
c
      integer maxobsperevent
      parameter (maxobsperevent=180)
      integer nobs
      character*4 sta(maxobsperevent)
      integer*4 itime
      real sec(maxobsperevent)
      character*1 rmk1(maxobsperevent),rmk2(maxobsperevent),
     &            cphase(maxobsperevent)
      integer iunit,iwt(maxobsperevent)
      real amx(maxobsperevent), prx(maxobsperevent),
     &     xlat,xlon
      real xmagni, depth, origtime
      real ttime(maxobsperevent)
      character*1 eventtype
      integer i1,i2,i3
      character*1 cns,cew
c
      character*80 cline
      integer iyr,imo,iday,ihr,imin
      integer j,j1
c
  1   eventtype='L'
      nobs=-1
c
   2  read(iunit,'(a)',end=999) cline
      if(cline.eq.' ') goto 2
c
      read(cline,'(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,1x,f7.2,
     1      2x,f5.2)',
     &           err=9999,end=999)
     &           iyr,imo,iday,ihr,imin,origtime,
     &           xlat,cns,xlon,cew,depth,xmagni !  i1  =ifx(i)
                !=======================!   trial hypocenter !!
cek  i1-switch is dummy
      i1=0
c
cek     EQS:   i1=ifx(i)     ! no longer in use !!!
c     SHOTS: i1=icc
c
      if(cns.eq.'S') xlat=-xlat
      if(cew.eq.'E') xlon=-xlon
c
      call TIMECLEAR(iyr,imo,iday,ihr,imin,origtime,itime)
c
      j1=1
      do j=1,maxobsperevent
         read(iunit,'(2x,a4,2x,a1,3x,i1,3x,f6.2)',end=99)
     &            sta(j),cphase(j),iwt(j),ttime(j)
         if(sta(j).eq.' ')then
            j1=j1-1          ! event is completely read-in
            goto 99
         endif
         j1=j
      enddo
c
   99 continue
      nobs=j1
c
      do j=1,nobs
         sec(j)=origtime+ttime(j)
      enddo
c
      RETURN
c
  999 continue
      nobs=-1
      RETURN
c
 9999 continue
      write(6,*)'INPUTARCVEL>>> read-error! input-line is:'
      write(6,'(a)') cline
      stop'subr. INPUTARCVEL >>> error!'
c
      end ! of subr. inputarcvel
