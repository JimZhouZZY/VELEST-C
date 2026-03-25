c
c
c
c
c
c
      subroutine FINALHYPOCOUT
c
c     output all the final hypocenters, print output simultaneous mode
c     and *cnv format
c
      implicit none
      include '../inc/vel_com.inc'
c
      real zero,xlat,xlon,xxx,sec
      integer izero,i,nin,j
      real tt(ist)
      character*1 phzz(ist), cns,cew
      parameter (zero=0.0,izero=0)
      integer year19,sec10,lat1000,lon1000,xout,yout,dep10
      integer mag10,rms100,idmin,aa
      real dist,dmin,xsta,ysta
c
      if(.not.single_turbo) write(16,*)
      if(isingle.ne.0)then
         if(.not.single_turbo) write(16,1111)
         write(6,1111)
      else
         if(.not.single_turbo)then
            write(16,11)
         endif
      endif
      do 12 i=1,legs
c
c     convert hypocenter into degrees:
c
      if(icoordsystem.eq.2)then
         call GEOKO( -e(2,i), e(3,i), xlat, xlon , 1 ) ! calc. LAT/LON
         xlon=-xlon
         xxx=-e(2,i)
      else
         call SDC( e(2,i), e(3,i), xlat, xlon , 1 ) ! calc. LAT/LON
         xxx=e(2,i)
      endif
      if(xlat.lt.0.0)then
         cns='S'
         xlat=-xlat
      else
         cns='N'
      endif
      if(xlon.lt.0.0)then
         cew='E'
         xlon=-xlon
      else
         cew='W'
      endif
c
      sec=e(1,i)
      nin=imin(i)
   23 if(sec.lt.0.) goto 13
   15 if(sec.lt.60.) goto 14
      sec=sec-60.
      nin=nin+1
      goto 15
   13 sec=sec+60.
      nin=nin-1
      goto 23
   11 format(1h ,4x,' date    origin   latitude longitude ',   ! file016
     &     ' depth  mag  no  rms      x      y      z')
 1111 format('  date    origin   latitude longitude',          ! screen (unit=6)
     &     '  depth  mag  no  rms      x      y      z')
c
c     output summary information
c
   14 continue
      if(nin.lt.0)then     !  U.K. 3.Feb.87
         nin=nin+60
         ihr(i)=ihr(i)-1
      endif
c
      call GAPCALC(i)
c
      if(.not.single_turbo)then
         write(16,'(1x,i3,1x,3i2.2,1x,2i2.2,f6.2,
     &              1x,f7.4,a1,1x,f8.4,a1,1x
     &              f6.2,f5.2,i4,f6.3,1x,3f7.2)')
     &             i,iyr(i),imo(i),iday(i),ihr(i),nin,sec,
     &             xlat,cns,xlon,cew,e(4,i),
     &             emag(i),knobs(i)-nobswithw0,rms(i),
     &             xxx,(e(j,i),j=3,4)
      endif
c
      if(isingle.ne.0)then
         write(6,'(1x,3i2.2,1x,2i2.2,f6.2,
     &              1x,f7.4,a1,1x,f8.4,a1,1x,
     &              f6.2,f5.2,i4,f6.3,2f7.2,f6.2)')
     &             iyr(i),imo(i),iday(i),ihr(i),nin,sec,
     &             xlat,cns,xlon,cew,e(4,i),
     &             emag(i),knobs(i)-nobswithw0,rms(i),
     &             xxx,(e(j,i),j=3,4)
      endif
c
      if(ismpout.eq.1)then
         year19 = iyr(i) + 1900
	 sec10 = nint(sec * 10.)
	 lat1000 = nint(xlat * 1000.)
	 lon1000 = nint(xlon * 1000.)
	 xout = nint(xxx)
	 yout = nint(e(3,i))
	 dep10 = nint(e(4,i))
	 if (dep10.lt.0.0) then
	    dep10 = 0.0
         endif
	 mag10 = nint(emag(1) * 10.)
	 if (mag10.lt.0.0) then
	    mag10 = 0.0
         endif
	 rms100 = nint(rms(i) * 100.)
c    
c        Calculating dmin
c
         dmin = 9999
         do aa=1,knobs(i)
	    if (kpwt(i,aa).lt.4) then
	       xsta = -x(istm(aa,i),1)
	       ysta = x(istm(aa,i),2)
               dist = sqrt((xxx-xsta)**2 + (e(3,i) - ysta)**2)
	       if (dist.lt.dmin) then
	          dmin = dist
	       endif
	    endif
	 enddo
	 idmin = nint(dmin)
         write(smpline,'(i4,4i2.2,i3.3,i5.5,a1,i6.6,a1,
     &              i3.3,i2.2,''Ml'',i4.4,2i3,''000000000SEDL'',
     &              2i3,i4.4,i2)')
     &              year19,imo(i),iday(i),ihr(i),nin,sec10,
     &              lat1000,cns,lon1000,cew,
     &              dep10,mag10,nreg,xout,yout,
     &              rms100,igap(i),idmin,knobs(i)-nobswithw0
      endif
c
chrm  If running in single mode, the diagonal elements of the covariance matrix
chrm  will also added to the smpline. In the simultaneous mode the line will be
chrm  written to the file here.
c
cek
cek   this option turned off for single event mode since smp-file only in
cek   SED special format
cek
      if(isingle.eq.0.and.ismpout.eq.1) then
         write(11,'(a80)') smpline
      endif
ccc      if(.not.single_turbo) write(16,*)
      if(isingle.ne.0)then
         if(.not.single_turbo)then
          write(16,'(1x,''Event# '',i3,'' GAP = '',i3)')
     &	  isingle,igap(i)
         endif
         if(icoordsystem.eq.2.and.nreg.ge.1000)then
            if(.not.single_turbo)then
               write(16,'(1x,a32,''   L+T Nr.: '',i4)') 
     &      regionname,nreg
            endif
            write(6,'(1x,a32,''   L+T Nr.: '',i4)') 
     &      regionname,nreg
         else
            if(.not.single_turbo)then
               write(16,'(1x,a32,''   F-E Nr.: '',i4)') 
     &         regionname,nreg
            endif
            write(6,'(1x,a32,''   F-E Nr.: '',i4)') 
     &      regionname,nreg
         endif
      endif
      if(icnvout.eq.0) goto 12  ! do NOT write on file07
c
c---- output final hypocenters and travel times to file07
c
cek     write(7,'(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,2f7.2)')
cek     &          iyr(i),imo(i),iday(i),ihr(i),nin,sec,
cek     &          xlat,cns,xlon,cew,
cek     &          e(4,i),emag(i)
cek
cek next statement changed by E.Kissling, 21.12.90:
cek          added output of gap and rms of event to converted format
cek
      write(7,'(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,2f7.2,
     &         4x,i3,5x,f5.2)')
     &         iyr(i),imo(i),iday(i),ihr(i),nin,sec,
     &         xlat,cns,xlon,cew,
     &         e(4,i),emag(i),igap(i),rms(i)
c
      imin(i)=nin
      do 18 j=1,knobs(i)
   18 tt(j)=pt(j,i)-e(1,i)
      E(1,I)=SEC
c
      do j=1,knobs(i)
         phzz(j)='P'
         if(sphase(j,i).eq.1.) phzz(j)='S'
         if(sphase(j,i).eq.-1.) phzz(j)='M'
         if(sphase(j,i).eq.2.) phzz(j)='-'
      enddo
      write(7,19) (smn(j,i),PHZZ(j),kpwt(j,i),tt(j),j=1,knobs(i))
      write(7,*)
   19 format(6(a4,a1,i1,f6.2))
c
  12  continue
      if(.not.single_turbo) write(16,*)
c
      RETURN
      end ! of subr. finalhypocout
