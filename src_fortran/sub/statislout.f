c
c
c
c
c
c
      subroutine STATISLOUT
c
c     print location-output in file *.VEL , compatible for program STATISL
c
cek     implemented print output of studentized residuals (studres)
cek  29.3.95 implemented print output of station name for all phases
c             (incl. S)
c
c       principally studentized res = residual/(sigma * sqrt(1.-diagofhat))
c
c       here we print out studres= residual/sqrt(1.-diag of hat matrix)
cek  281093
c
      implicit none
      include '../inc/vel_com.inc'
c
      real sec,xlat,xlon,aar,ofd,tfd
      real erh,studres
      real erx,ery,erz
      real xstn,ystn,xhyp,yhyp,azi,tobs,tcorr
      integer nin,idmin,k,js,jd,no,jav,knobs1,i,in,kk,iazi
      integer iamx,iprisecondcard
      character ctime*20, card*90
      character*1 clay, char1, cns,cew
      character*1 q,qs,qd
c
      character*1 class(4)
      data class/'A','B','C','D'/
c
      if( (knobs(1)-nobswithw0) .lt. nvar .and.iabort.eq.0)then
         iabort=1
         if(.not.single_turbo)then
            write(16,*)'knobs(i)-nobswithw0 < nvar !!!'
            write(16,*)'Event cannot be located!!!'
         endif
         write(6,*)'knobs(i)-nobswithw0 < nvar !!!'
         write(6,*)'Event cannot be located!!!'
      endif
      if(iabort.eq.1)then
         write(2,'('' ERROR: insufficient data to locate the quake!'')')
         RETURN
      endif
c
      if(ifixsolution.eq.0)then
         write(2,'(''0 DATE  ORIGIN   TIME   LAT      LON     DEPTH '',
     &             '' MAG  NO  DM GAP  RMS   ALE D-SPR'')')
      endif
      if(ifixsolution.eq.1)then
         write(2,'(''0 DATE  ORIGIN TIME   LAT       LON     *DEPTH*'',
     &             '' MAG  NO  DM GAP  RMS   ALE D-SPR'')')
      endif
      if(ifixsolution.eq.9)then
         write(2,'(''0 DATE  ORIGIN TIME  *LAT*     *LON*    *DEPTH*'',
     &             '' MAG  NO  DM GAP  RMS   ALE D-SPR'')')
      endif
c
      sec=e(1,1)
      nin=imin(1)
      if(sec.lt.0.)then
         sec=sec+60.
         nin=nin-1
      endif
      if(sec.gt.60.)then
         sec=sec-60.
         nin=nin+1
      endif
      if(nin.lt.0)then     !  U.K. 3.Feb.87
         nin=nin+60
         ihr(1)=ihr(1)-1
      endif
c
c     convert hypocenter into degrees:
c
      if(icoordsystem.eq.2)then
         call GEOKO( -e(2,1), e(3,1), xlat, xlon , 1 )  ! calc. LAT/LON
         xlon=-xlon
      else
         call SDC( e(2,1), e(3,1), xlat, xlon , 1 )  ! calc. LAT/LON
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
      idmin=999   ! minimum distance (epicenter --> receiver)
      aar=0.0   ! absolute average weighted residual
      do k=1,knobs(1)
         delta=SQRT(  (e(2,1)-d(k,1,1))**2
     &              + (e(3,1)-d(k,2,1))**2  )
         idelta(k)=NINT(delta)
         if(w(k,1).gt.0.0)then
            if(idelta(k).lt.idmin) idmin=idelta(k)
            aar=aar+ABS(res(k,1))*w(k,1)
         endif
      enddo
      aar=aar/(knobs(1)-nobswithw0)
c
c write first line of summary card:
c
      write(2,'(1x,3i2.2,1x,i2,'':'',i2,'':'',f6.3,1x
     &          f7.4,a1,f8.4,a1,1x,f7.3,2x,
     &          f3.1,2x,i2,1x,i3,1x,i3,f5.2,f6.2,
     &          1x,f5.2)')
     &          iyr(1),imo(1),iday(1),ihr(1),nin,sec,
     &          xlat,cns,xlon,cew,e(4,1),
     &          emag(1),knobs(1)-nobswithw0,idmin,igap(1),rms(1),ale(1),
     &          spread
      write(2,'(''0  ERX  ERY  ERZ Q SQD  ADJ  IN NR  '',
     &          ''AVR   AAR  NM AVXM  SDXM IT'')')
c
c     standard deviations :
c     (square-root of diagonalelements of UNIT covariance matrix
c
      erh=sqrt(s(2)**2 + s(3)**2)
      erx=s(2)
      ery=s(3)
      erz=s(4)
      JS=4
      IF((RMS(1).LT.1.00).AND.(ERH.LE.15.0))JS=3
      IF((RMS(1).LT.0.60).AND.(ERH.LE.8.0).AND.(erz.LE.15.0))JS=2
      IF((RMS(1).LT.0.30).AND.(ERH.LE.2.0).AND.(erz.LE.6.0)) JS=1
      JD=4
      no=knobs(1)-nobswithw0
      OFD=e(4,1)      ! focal depth
      TFD=2.*e(4,1)
      IF(OFD .LT.10.) OFD=10.
      IF(TFD .LT. 30.) TFD=30.
      IF((iGAP(1).LE.180).OR.(NO.GE.4).AND.(iDMIN.LE.100)) JD=3
      IF((iGAP(1).LE.135).OR.(NO.GE.5).AND.(iDMIN.LE.TFD )) JD=2
      IF((iGAP(1).LE. 90).OR.(NO.GE.6).AND.(iDMIN.LE.OFD )) JD=1
      JAV=(JS+JD+1)/2
      Q=class(JAV)
      QS=class(JS)
      QD=class(JD)
      knobs1=0
      do i=1,knobs(1)
         if(kpwt(i,1).ge.5) knobs1=knobs1+1
      enddo
      knobs1=knobs(1)-knobs1 ! # of obs with weight less 5 !!!
      IN=0
c
c write second line of summary card
c
      write(2,'(1x,f5.1,f5.1,f5.1,1x,a1,1x,a1,''/'',a1,f6.2,1x,
     &          i2,i3,2f6.2,
     &          i3,2x,f3.1,2x,f3.1,1x,i3)')
     &                          erx,ery,erz,q,qs,qd,steplen,IN,
     &                          knobs1,avres(1),AAR,
     &                          nmag,xmagnitude,sdxmagnitude,nitt
      if(icoordsystem.eq.2)then
         if(nreg.ge.1000)then
            write(2,'(''0 L+T NR:'',i4,1x,a32,''CH-COORD.:'',
     &       f9.3,'' /'',f9.3, '' KM'')') nreg,regionname,-e(2,1),e(3,1)
         else
            write(2,'(''0 F-E NR:'',i4,1x,a32,''CH-COORD.:'',
     &       f9.3,'' /'',f9.3, '' KM'')') nreg,regionname,-e(2,1),e(3,1)
         endif
      else
         write(2,'(''0 F-E NR:'',i4,1x,a32)') nreg,regionname
      endif
      write(2,'(''0 STN  DIST AZM AIN PRMK HRMN  P-SEC  TPOBS  TPCAL '',
     &          '' -TSCOR  P-RES   P-WT IMP STURES'')')
      write(2,'(''        AMX PRX     SRMK XMAG  S-SEC  TSOBS  TSCAL '',
     &          '' -TSCOR  S-RES   S-WT IMP STURES'')')
c
      k=0
c
c  DO LOOP for EACH OBSERVATION
c
      do kk=1,knobs(1)
      k=k+1
c
c     station# are stored in array  ISTM(iobs,ievent)
c     stn(stn#)   is the station-name (a4)
c     smn(iobs,iev)  is the station-name
c     sphase(nobs,iev) : 0 for P, 1 for S, -1 for reflected P (M-phase)
c     and 2 for an s-p phase
c
      clay=' '
      if(sphase(k,1).eq.0) char1='P'
      if(sphase(k,1).eq.1) char1='S'
      if(sphase(k,1).eq.2) char1='-'
      if(sphase(k,1).eq.-1)then
         char1='P'
         clay='M'
      endif
c        Azimuth (stn --> hypoc) = 57.296*ATAN2(xhyp-xstn,yhyp-ystn)
      xstn=d(k,1,1)
      ystn=d(k,2,1)
      xhyp=e(2,1)
      yhyp=e(3,1)
      azi=57.296*ATAN2(xhyp-xstn,yhyp-ystn)
      if(azi.lt.0) azi=azi+360.
c        Azimuth (hypoc --> stn) = Azimuth (stn --> hypoc) + 180 deg
      azi=azi+180.
      azi=MOD(azi,360.)
      iazi=NINT(azi)
c
      tobs=pt(k,1)-e(1,1)  !  P observed travel time
cek
cek correct time for minuite change as previously done for origin time
cek (see above)
cek this correction added by ek 2.May1996
c
      if(tobs.lt.0.)then
         pt(k,1)=pt(k,1)+60.
      endif
      if(tobs.gt.60.)then
         pt(k,1)=pt(k,1)-60.
      endif
c
      tobs=pt(k,1)-e(1,1)  !  P observed travel time
c
      if(tcalc(k).lt.0.)then
         tcalc(k)=tcalc(k)+60.
      endif
      if(tcalc(k).gt.60.)then
         tcalc(k)=tcalc(k)-60.
      endif
c
c station corrections:
      tcorr=ptcor(istm(k,1))  
      if(nsp.eq.2.and.sphase(k,1).eq.1.0) tcorr=stcor(istm(k,1))
      if(nsp.eq.3.and.sphase(k,1).eq.1.0) tcorr=tcorr*vpvs
      if(sphase(k,1).eq.2) tcorr= 0.0
c
      card=' '
      studres=res(k,1)/sqrt(1.-drm(k,k))
      if(studres.gt.999.) studres=999.999
c
cek print P card (or S-card if no P obs for this station available
c
      iprisecondcard=0
      write(card,'(2X,A4,1x,2i4,I4,1X,3A1,I1,a1,2I2,3F7.3,F7.3,1x,
     &             f7.3,1x,f6.2,1x,f6.4,f7.3)')
     &        smn(k,1),idelta(k),iazi,iain(k),prmk(k,1),char1,prmk(k,2),
     &        kpwt(k,1),
     &        clay,ihr(1),nin,pt(k,1),tobs,tcalc(k),
     &        tcorr,res(k,1),w(k,1),drm(k,k),studres
cek      if(w(k,1).eq.0.0) write(card(61:61),'(''*'')')
      if(kpwt(k,1).eq.5)then
         card(21:80)=' ' ! there is no P arrival!
         if(k.lt.knobs(1))then
            if(smn(k+1,1).eq.smn(k,1))then
               write(card(16:19),'(i4)') iain(k+1)  ! angle of next phase (S) !!
               write(2,'(a)') card !next card will be an S
            endif
         else
            goto 99
         endif
      else
         write(2,'(a)') card
      endif
      iamx=NINT(amx(k))
      card=' '
      if(xmagni(k).ne.-13.)then  ! a magnitude was calculated for this observ.!
         write(card(1:15),'(2x,4x,i5,f4.1)') iamx, prx(k)
         write(card(26:28),'(f3.1)') xmagni(k)
      endif
      if(k.lt.knobs(1))then
chrm?   Must anything changed here for s-p phases?  
c
c  NOW  PRINT  S-PHASE 
c   
         if(smn(k+1,1).eq.smn(k,1))then
            if(sphase(k+1,1).eq.1)then
              iprisecondcard=1
              k=k+1
c
              char1='S'
              clay=' '
              tobs=pt(k,1)-e(1,1)
cek
cek correct time for minuite change as previously done for origin time
cek (see above)
cek this correction added by ek 2.May1996
c
              if(tobs.lt.0.)then
                  pt(k,1)=pt(k,1)+60.
              endif
              if(tobs.gt.60.)then
                  pt(k,1)=pt(k,1)-60.
              endif
c
              tobs=pt(k,1)-e(1,1)  !  S observed travel time
c
              if(tcalc(k).lt.0.)then
                  tcalc(k)=tcalc(k)+60.
              endif
              if(tcalc(k).gt.60.)then
                  tcalc(k)=tcalc(k)-60.
              endif
c
              tcorr=ptcor(istm(k,1))
              if(nsp.eq.2.and.sphase(k,1).eq.1.0) then
		 tcorr=stcor(istm(k,1))
              endif
              if(nsp.eq.3.and.sphase(k,1).eq.1.0) then
		 tcorr=tcorr*vpvs
              endif
cek       NOW write a following card (second phase) for same station
cek  29.3.95
              write(card(3:6),'(a4)') smn(k-1,1) !write station name for S
              write(card(21:25),'(3a1,i1,a1)')
     &                          prmk(k,1),char1,prmk(k,2),kpwt(k,1),clay
              studres=res(k,1)/sqrt(1.-drm(k,k))
              if(studres.gt.999.) studres=999.999
              write(card(30:),'(3F7.3,F7.3,1x,f7.3,1x,F6.2,1x,f6.4,
     &                          f7.3)')
     &                          pt(k,1),tobs,tcalc(k),
     &                          tcorr,res(k,1),w(k,1),drm(k,k),studres
cek              if(w(k,1).eq.0.0) write(card(61:61),'(''*'')')
            endif
         endif
      endif
cek next statements may95
      if(iprisecondcard.eq.1) then
         write(2,'(a)') card
      endif
cek      if(kpwt(k,1).lt.5) write(2,'(a)') card
 2022 continue
      if(k.eq.knobs(1)) goto 99
      enddo
c
  99  continue
      call DATETIME(ctime)  ! get date&time from system
      write(2,'(''  $$$ '',2x,''VELEST-Version ETH-11FEB92'',
     &          '' located at: '',a20)') ctime
      write(2,*)
      write(2,*)
c
      return
      end ! of subr. statislout
