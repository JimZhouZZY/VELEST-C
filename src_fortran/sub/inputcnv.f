c
c
c
c
c
c
      subroutine INPUTCNV(iunit,nobs,
     &                  sta,iyr,imo,iday,ihr,imin,sec,
     &                  rmk1,rmk2,cphase,iwt,amx,prx,
     &                  xlat,xlon,xmagni,depth,origtime,
     &                  i1,i2,i3,eventtype)
c
c     on Output:
c     ---------
c     EQS:   i1=ifx(i)   
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
      integer j,j1,j2
c
  1   eventtype='L'
      j2=-1
      nobs=0
c
   2  read(iunit,'(a)',end=99) cline
      if(cline.eq.' ') goto 2
cek search for end of *.cnv file, marked by 9999
      if(cline.eq.'9999') goto 999
c
cek   next line adjusted by EK 3.12.90
      read(cline,'(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,1x,f7.2,
     1      2x,f5.2)',
     &           err=9999,end=999)
     &           iyr,imo,iday,ihr,imin,origtime,
     &           xlat,cns,xlon,cew,depth,xmagni ! i1 =ifx(i)
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
      j2=0
   10 j1=j2+1
      j2=j1+5
      read(iunit,'(6(a4,a1,i1,f6.2))',end=999)
     &            (sta(j),cphase(j),iwt(j),ttime(j),j=j1,j2)
      do j=j1,j2
         if(sta(j).eq.' ')then
            if(j.eq.1.and.j1.eq.1) goto 1  ! blank-line read !
            j2=j-1   ! event is completely read-in
            goto 99
         endif
      enddo
      goto 10         ! read next input-line
c
   99 continue
      nobs=j2
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
      write(6,*)'INPUTCNV>>> read-error! input-line is:'
      write(6,'(a)') cline
      stop'subr. INPUTCNV >>> error!'
c
      end ! of subr. inputcnv
