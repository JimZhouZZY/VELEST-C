c
c
c
c
c
c
      subroutine STATISTICSOUT
c
c     output the statistics done in this VELEST run
c
      implicit none
      include '../inc/vel_com.inc'
c
      integer i,ii,irefr,isour,itotal,nhit,iavgap,lesseq1
      integer mge6,nobslesseq1,nobsmge6,mag,mge5l,mge5i
      integer irefl,idepi,idepl,ntot,i1,i2,i3,i4
      real depthi,depthl,res0,res1,res2,res3,res4
      real tot
      real  rlen,err,mini,maxi,xmag
      integer magnr(49), nobsnr(49)
      integer depthnri(105), depthnrl(105)
      character cstari*51, cstarl*51
c
      do i=1,105
         depthnri(i)=0
         depthnrl(i)=0
      enddo
      do i=1,49
         magnr(i)=0
         nobsnr(i)=0
      enddo
c
      do i=1,neqs
         do ii=1,nltot-1
            if(e(4,i).ge.h(ii).and.e(4,i).lt.h(ii+1))then
               ihypoclayer(ii)=ihypoclayer(ii)+1
            endif
         enddo
         if(e(4,i).ge.h(nltot)) ihypoclayer(nltot)=ihypoclayer(nltot)+1
      enddo
      rlen=0.0
      do i=1,nltot
         rlen=rlen+refraylen(i)       !  rlen  will be total refr.raylength [km]
      enddo
      do i=1,nltot
chrm  The next if else statement seems to be neccessary to avoid
chrm  floating point exceptions in case of rlen = 0
      if (rlen.gt.0.00001) then
         refraylen(i)=100.*refraylen(i)/rlen
      else 
	 refraylen(i) = 0.0
      endif
         if(hitlay(i,1).ge.1)then
            hitlay(i,2)=hitlay(i,2)/hitlay(i,1)  !  average horiz. km in layer i
            hitlay(i,3)=hitlay(i,3)/hitlay(i,1)  !  average verti. km in layer i
         endif
      enddo                                  ! hitlay(i,1) = nofHITS of layer i
      write(16,*)
      write(16,*)
      write(16,*)
      write(16,*)
      write(16,*)'RAY-STATISTICS FOR  L A S T  ITERATION'
      write(16,*)'--------------------------------------'
      write(16,*)
      write(16,*)'NHYP : nr of hypocenters in this layer'
      write(16,*)'NREF : nr of headwaves in this layer'
      write(16,*)'%len : % of "refracted km" in this layer '
      write(16,*)'       with respect to all refracted kilometers'
      write(16,*)'NHIT : nr of rays passed thru this layer'
      write(16,*)'xy-km: average horizontal ray length [km] in layer'
      write(16,*)'z-km : average  vertical  ray length [km] in layer'
      write(16,*)'RFLX : number of reflections at bottom of this layer'
      write(16,*)
      write(16,'(a)')' nlay   top ..... bottom     velocity   '//
     &                'NHYP NREF %len  NHIT xy-km  z-km  RFLX'
      write(16,*)
      irefr=0
      isour=0
      irefl=0
      do i=1,nltot-1
         nhit=NINT(hitlay(i,1))
         write(16,'(2x,i2,3x,f6.2,''...'',f6.2,'' km   '',
     &         f5.2,'' km/s'',2(1x,i4),1x,f5.1,1x,i5,2(1x,f5.1),2x,i4)')
     &              i,h(i),h(i)+thk(i),
     &              v(i),ihypoclayer(i),irefrlayer(i),refraylen(i),
     &              nhit,hitlay(i,2),hitlay(i,3),irefllayer(i)
         irefr=irefr+irefrlayer(i)
         isour=isour+ihypoclayer(i)
         irefl=irefl+irefllayer(i)
      enddo
      nhit=NINT(hitlay(nltot,1))
      write(16,'(2x,i2,3x,f6.2,''...'',5x,''  km   '',
     &      f5.2,'' km/s'',2(1x,i4),1x,f5.1,1x,i5,2(1x,f5.1))')
     &      nltot,h(nltot),
     &      v(nltot),ihypoclayer(nltot),irefrlayer(nltot),
     &      refraylen(nltot), nhit,hitlay(nltot,2),hitlay(nltot,3)
      irefr=irefr+irefrlayer(nltot)
      isour=isour+ihypoclayer(nltot)
      irefl=irefl+irefllayer(nltot)
c
      write(16,*)
      write(16,'('' Total nr of events was '',i4)') isour
      write(16,*)
      write(16,'('' Total nr of refracted rays = '',i5)') irefr
      write(16,'('' Total nr of reflected rays = '',i5)') irefl
      write(16,'('' Total nr of   other   rays = '',i5)') noheadwave
      itotal=noheadwave+irefr+irefl
      write(16,'(''                              ------'')')
      write(16,'('' Total nr of    all    rays = '',i5)') itotal
      write(16,*)
      write(16,*)'Average (absolute) error of the raytracers:'
      err=0.0
      if(noheadwave.gt.0) err=(sterr+direrr)/noheadwave
      write(16,'(1x,'' Straight and direct rays : '',f7.2,'' meters'')')
     &           err
      err=0.0
      if(irefr.gt.0) err=refrerr/irefr
      write(16,'(1x,'' Refracted           rays : '',f7.2,'' meters'')')
     &           err
      err=0.0
      if(irefl.gt.0) err=reflerr/irefl
      write(16,'(1x,'' Reflected           rays : '',f7.2,'' meters'')')
     &           err
      write(16,*)
      write(16,*)
      avhraylen=avhraylen/itotal
      write(16,*)'ALL RAYS TOGETHER:'
      write(16,'(1x,''Average horizontal ray length = '',
     &           f6.1,'' km   (Hypocenter --> Station) '')')
     &           avhraylen
      avvraylen=avvraylen/itotal
      write(16,'(1x,''Average  vertical  ray length = '',
     &           f6.1,'' km   (Deepest ray-point --> Station)'')')
     &           avvraylen
      write(16,*)
c
      write(16,*)
      write(16,*)
      write(16,*)'... and some more STATISTICS '
      write(16,*)'---------------------------- '
      write(16,*)' GAP of final epicenters:'
      write(16,*)
      write(16,'(5(''  Event# -> GAP''))')
      write(16,'(5(5x,i3,'' -> '',i3))') (i,igap(i),i=1,neqs)
      write(16,*)
      mini=361
      maxi=-1
      iavgap=0
      do i=1,neqs
         if(igap(i).gt.maxi) maxi=igap(i)
         if(igap(i).lt.mini) mini=igap(i)
         iavgap=iavgap+igap(i)
      enddo
      iavgap=NINT(float(iavgap)/float(neqs))
      write(16,'(1x,''GAPs were between '',i3,'' and '',i3)')
     &           nint(mini),nint(maxi)
      write(16,'(''      (average GAP was '',i3,'')'')') iavgap
      write(16,*)
c
      write(16,*)
      write(16,*)
      write(16,*)' MAGNITUDES of INPUT-DATA:'
      write(16,*)
      write(16,*)'Magnitude (# of events)   ***average number of obs***'
      write(16,*)
      lesseq1=0
      mge6=0
      nobslesseq1=0
      nobsmge6=0
      do i=1,neqs
         mag=NINT(emag(i)*10.)
         if(mag.le.10)then
            nobslesseq1=nobslesseq1+knobs(i)
            lesseq1=lesseq1+1
         endif
         if(mag.gt.10.and.mag.lt.60)then
            magnr(mag-10)=magnr(mag-10)+1
            nobsnr(mag-10)=nobsnr(mag-10)+knobs(i)
         endif
         if(mag.ge.60)then
            mge6=mge6+1
            nobsmge6=nobsmge6+knobs(i)
         endif
      enddo
      cstari='*********1*********2*********3*********4*********5>'
      if(lesseq1.gt.0) nobslesseq1=
     &                 NINT(float(nobslesseq1)/float(lesseq1))
      if(nobslesseq1.gt.50)nobslesseq1=51
      if(nobslesseq1.gt.0)then
         write(16,'(1x,''MAG<= 1.0 '',''('',i3,'') '',a)')
     &              lesseq1,cstari(1:nobslesseq1)
      else
         write(16,'(1x,''MAG<= 1.0 '',''('',i3,'') '',a)')
     &              lesseq1,' '
      endif
      do i=1,49
         xmag=(10.+i)/10.
         if(magnr(i).gt.0) nobsnr(i)=
     &                     NINT(float(nobsnr(i))/float(magnr(i)))
         if(nobsnr(i).gt.50)nobsnr(i)=51
         if(nobsnr(i).gt.0)then
            write(16,'(1x,''MAG = '',f3.1,'' ('',i3,'') '',a)')
     &                xmag,magnr(i),cstari(1:nobsnr(i))
         else
            write(16,'(1x,''MAG = '',f3.1,'' ('',i3,'') '',a)')
     &                xmag,magnr(i),' '
         endif
      enddo
      if(mge6.gt.0) nobsmge6=NINT(float(nobsmge6)/float(mge6))
      if(nobsmge6.gt.50)nobsmge6=51
      if(nobsmge6.gt.0)then
         write(16,'(1x,''MAG>= 6.0 '',''('',i3,'') '',a)')
     &              mge6,cstari(1:nobsmge6)
      else
         write(16,'(1x,''MAG>= 6.0 '',''('',i3,'') '',a)')
     &              mge6,' '
      endif
      write(16,*)
      write(16,*)
      write(16,*)
c
      write(16,*)
      write(16,*)
      write(16,*)' DEPTHs of INPUT-DATA and of LAST ITERATION:'
      write(16,*)' ==========================================='
      write(16,*)
      write(16,*)'      Depth       # of events  '
      write(16,*)
      lesseq1=0
      mge6=0
      mge5l=0
      mge5i=0
      do i=1,neqs
         idepi=NINT(depthsofinput(i))
         if(idepi.lt.100)then
            depthnri(idepi+5)=depthnri(idepi+5)+1
         endif
         if(idepi.ge.100)then
            mge5i=mge5i+1
         endif
         idepl=NINT(e(4,i))
         if(idepl.lt.100)then
            depthnrl(idepl+5)=depthnrl(idepl+5)+1
         endif
         if(idepl.ge.100)then
            mge5l=mge5l+1
         endif
      enddo
      cstari='.........1.........2.........3.........4.........5>'
      cstarl='*********1*********2*********3*********4*********5>'
      do i=1,105
         depthi=i-5
         depthl=i-5
         if(depthnri(i).gt.50) depthnri(i)=51
         if(depthnrl(i).gt.50) depthnrl(i)=51
         if(depthnri(i).gt.0)then
            write(16,'(1x,''DEPTH ( input ) = '',f4.0,'' km :  '',a)')
     &                depthi,cstari(1:depthnri(i))
         else
            write(16,'(1x,''DEPTH ( input ) = '',f4.0,'' km :  '',a)')
     &                depthi,' '
         endif
         if(depthnrl(i).gt.0)then
            write(16,'(1x,''DEPTH (last_IT) = '',f4.0,'' km :  '',a)')
     &                depthl,cstarl(1:depthnrl(i))
         else
            write(16,'(1x,''DEPTH (last_IT) = '',f4.0,'' km :  '',a)')
     &                depthl,' '
         endif
cc         write(16,*)
      enddo
      if(mge5i.gt.50)mge5i=51
      if(mge5l.gt.50)mge5l=51
      if(mge5i.gt.0)then
         write(16,'(1x,''DEPTH ( input ) > 100. '',''km :  '',a)')
     &              cstari(1:mge5i)
      else
         write(16,'(1x,''DEPTH ( input ) > 100. '',''km :  '',a)')
     &              ' '
      endif
      if(mge5l.gt.0)then
         write(16,'(1x,''DEPTH (last_IT) > 100. '',''km :  '',a)')
     &              cstarl(1:mge5l)
      else
         write(16,'(1x,''DEPTH (last_IT) > 100. '',''km :  '',a)')
     &              ' '
      endif
      write(16,*)
      write(16,*)
c
      write(16,*)'Residuals of the stations according to the azimuth:'
      write(16,*)'(RES  = total average residual at station)'
      write(16,*)'(RES1 = average residual of rays from 1st quadrant)'
      write(16,*)'(RES2 = average residual of rays from 2nd quadrant)'
      write(16,*)'(RES3 = average residual of rays from 3rd quadrant)'
      write(16,*)'(RES4 = average residual of rays from 4th quadrant)'
      write(16,*)
      write(16,'(1x,''Stn#  Stn     RES          RES1        '',
     &           '' RES2         RES3         RES4'')')
      do i=1,nsta
         res0=0.0
         res1=0.0
         res2=0.0
         res3=0.0
         res4=0.0
         if(stnazires(i,2).gt.0) res1=stnazires(i,1)/stnazires(i,2)
         if(stnazires(i,4).gt.0) res2=stnazires(i,3)/stnazires(i,4)
         if(stnazires(i,6).gt.0) res3=stnazires(i,5)/stnazires(i,6)
         if(stnazires(i,8).gt.0) res4=stnazires(i,7)/stnazires(i,8)
         tot=stnazires(i,2)+stnazires(i,4)+stnazires(i,6)+stnazires(i,8)
         if(tot.gt.0.0) res0=(
     &                        stnazires(i,2)*res1
     &                       +stnazires(i,4)*res2
     &                       +stnazires(i,6)*res3
     &                       +stnazires(i,8)*res4) / tot
         if(tot.gt.0.0)then
            ntot=NINT(tot)
            i1=NINT(stnazires(i,2))
            i2=NINT(stnazires(i,4))
            i3=NINT(stnazires(i,6))
            i4=NINT(stnazires(i,8))
            write(16,'(1x,i3,3x,a4,1x,5(f7.2,''('',i4,'')''))')
     &      i,stn(i),res0,ntot,res1,i1,
     &      res2,i2,res3,i3,res4,i4
         else
            write(16,'(1x,i3,3x,a4,1x,''-.-'')') i,stn(i)
         endif
      enddo
      write(16,*)
      write(16,*)
c
      RETURN
      end ! of subr. statisticsout
