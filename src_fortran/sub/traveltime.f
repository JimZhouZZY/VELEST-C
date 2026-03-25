c
c
c
c
c
c
      subroutine TRAVELTIME(i,nobs,iresflag)
c
c     computes the forward problem: all the raytracing is done, traveltimes
c     AND traveltime-derivatives are calculated. Moreover, a whole bunch
c     of statistics is done here.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Input parameter:                                            c
c                                                                 c
c     i          = ith event                                      c
c     nobs       = nr of obs for the ith event                    c
c     iresflag   =
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer i,nobs,iresflag
      include '../inc/vel_com.inc'
c
      integer ii,nittdone,ifirst,k2,nl,k1,mll,nrp,nrtn
      integer nrpdeep,ir,kk1,j,jzz
      real z,ttt,tkh,r1,r2,ster,direr,refrer,refler
      real takeoff_angle
      real dtddrefl,dtdhrefl,pobs,extrat1,extrat2
      integer do_s
      real rp(3,inrpmax),
     &     ss(inrpmax),x2(inrpmax),y2(inrpmax),z2(inrpmax)
      real v1(inltot),vsq1(inltot)
      real dtdr_s(3),dtdv_s(inltot),ttt_s
      save nITTdone, ifirst
c
c     reset statistics-parameter if iteration is new
c     or if NITT is the same but all rays are shot again (e.g. after
c     model/hypocenter-backup)
c
      if(isingle.gt.0)then
         if(kpwt(nobs,i).eq.5)then  ! this is not a real phase!!!
            if(w(nobs,i).ne.0.0) stop'TRAVELTIME>>> error in obs-wt!'
            res(nobs,i)=0.0
            do ii=1,3
               dtdr(ii)=0.0
            enddo
            RETURN
         endif
      endif
      if(isingle.eq.0)then
         if(i.eq.1.and.nobs.eq.1.and.nitt.gt.nITTdone)ifirst=1 ! 1st obs 1x done
         if(i.eq.1.and.nobs.eq.1.and.nitt.eq.nITTdone)ifirst=2
c ifirst=1   --->  first of all the obs for the first time done in this iterat.
c ifirst=2   --->  first of all the obs for the 2nd time coming in this iterat.
         if(nitt.gt.nITTdone.or.ifirst.eq.2)then
            nITTdone=nitt
            ifirst=1
            call RESETSTATIS
         endif
      endif
c
c     Reset s-p switch
c
      do_s = 1
c      write(*,*) 'DEBUG traveltime: i=',i,', nobs=',nobs,', k2=',k2
c      write(*,*) 'DEBUG e: e(1,i)=',e(1,i),', e(2,i)=',e(2,i),
c     &      ', e(3,i)=',e(3,i),', e(4,i)=',e(4,i)
c
c
c    z -- event depth
      z=e(4,i)
c    k2 -- velocity model number
      k2=iphase(nobs,i) 
   10 continue      
c     nl -- number of layers for model k2
      nl=nplay(k2)
c    k1 -- station number
      k1=istm(nobs,i)
c    v,vsq are the velocity and velocity squared
c     h(ii) is the depth to the top of layer ii
c     thk(ii) -- thickness of layer ii
c
      do ii=1,nl
         v(ii)=vp(k2,ii)
         v1(ii)=vp(k2,ii)
         if(nsp.eq.3.and.sphase(nobs,i).eq.1.0) v(ii)=vp(k2,ii)/vpvs
         if(nsp.eq.3.and.sphase(nobs,i).eq.2.0) then
	    if (do_s.eq.1) then
	       v(ii)=vp(k2,ii)/vpvs 
	    endif
	 endif   
         vsq(ii)=v(ii)**2
         vsq1(ii)=v1(ii)**2
         h(ii)=hp(k2,ii)
         thk(ii)=thkp(k2,ii)
      enddo
c
c     calculate tkh, r1,r2,delta :
c
      ttt=0.
      tkh=x(k1,3)-h(1)   ! distance from station up to top of model
      r1=(e(2,i)-x(k1,1))
      r2=(e(3,i)-x(k1,2))
c check if event too shallow
      if(itopo.gt.0)then
         if(e(4,i).lt.0.0)then  ! depth above sea-level...so near surface !
            call CHTOP(-e(2,i),e(3,i),zmin,
     &                 topo1file,topo2file) ! zmin:==surface at this point
         else
            zmin=zmininput ! depth below zero ... so zmin=zmininput
         endif
      endif
      if(i.le.neqs)then
         if(e(4,i).le.zmin)then
            if(isingle.eq.0) e(4,i)=zmin+0.1     ! earthquakes
            if(isingle.ne.0) e(4,i)=zmin+0.011   ! earthquake, single event
         endif
      else
         if(e(4,i).le.zmin) e(4,i)=zmin+0.011   ! shots
      endif
      if(ifixsolution.gt.0.and.e(4,i).le.0.0) e(4,i)=zmin+0.001 ! fix depth to
c                                                        !   min_depth allowed!!
cc      r3=(e(4,i)-x(k1,3))
      delta=sqrt(r1*r1+r2*r2)
c
c     set reflector-layer-boundary :
c
      MLL=0
      if(ireflector.gt.0)then
         if(sphase(nobs,i).eq.-1.0)then
            MLL=ireflector    ! P is reflected at bottom of layer ireflector
         endif
      endif
c
c     do the ray-tracing :
c
      call RAYPATH(1,1,1,1.,1.,1.,1.,nl,thk,h,v,vsq,
     &            e(2,i),e(3,i),e(4,i),x(k1,1),x(k1,2),x(k1,3),
     &            rp,nrp,nrtn,jl,tkj,1,ttt,MLL,ster,direr,refrer,refler,
     &            DTDDrefl,DTDHrefl)
c
      call CHECKRAYPATH(rp,nrp)
c
c
c     raypoints stored in RP [1st index: x,y,z; 2nd: raypt# (max. 2*nlay)]
c     nrp = nr of raypoints stored at moment in RP
c
      if(itopo.eq.2)then !check each raypoint whether it is below surface or not
        call RAYPOINTCHECK(rp,nrp,stn(k1))
      endif
c
c     if ray travels through the air, bend it below surface !!!
c
      if(itopo.eq.3) call BENDRAY(rp,nrp,stn(k1),v(1),ttt)
c
c     sum raytracer-errors
c
c     
c     For s-p phases the next statements should be done after the
c     p-wave calculation
c
      if (sphase(nobs,i).ne.2.0.or.do_s.eq.0) then     
         if(isingle.eq.0)then
            sterr=sterr+ster
            direrr=direrr+direr
            refrerr=refrerr+refrer
            reflerr=reflerr+refler
         endif
c      
         if(isingle.eq.0) call LAYERHIT(rp,nrpdeep,nl,nrp,mll)
c      
         if(irayout.eq.1)then
            write(13,1301) i,stn(k1),nrp
            if(icoordsystem.eq.2)then
               write(13,1302) (-rp(1,ir),rp(2,ir),rp(3,ir),ir=1,nrp)
            else
               write(13,1302) (rp(1,ir),rp(2,ir),rp(3,ir),ir=1,nrp)
            endif
 1301    format(/,' event =',i5,' station ',a4,' num. raypoints= ',i5)
 1302    format(3(2x,3f7.2))
         endif
      endif ! if .NOT. ...	 
c
c     calculate traveltime-derivatives according to the raytype :
c
      goto (5,5,6,6,5,6,5,8),nrtn   !   nrtn=8  if ray is reflected !!
c
      stop'TRAVELTIME>>> illegal nrtn from raytracer!!!'
c
  5   call TRAVDERIV('direct',
     &                nl,mll,v1,vsq1,rp,nrp,x2,y2,z2,ss,r1,r2,i,nobs)
      goto 40
  6   call TRAVDERIV('refracted',
     &                nl,mll,v1,vsq1,rp,nrp,x2,y2,z2,ss,r1,r2,i,nobs)
      goto 40
  8   call TRAVDERIV('reflected',
     &                nl,mll,v1,vsq1,rp,nrp,x2,y2,z2,ss,r1,r2,i,nobs)
 40   continue
c
c     For s-p phases, the traveltime of the s-wave was calculated. Now
c     the p-wave calculation will be performed.
c
      if (sphase(nobs,i).eq.2.0.and.do_s.eq.1) then
         do_s = 0
	 ttt_s = ttt
	 do ii=1,3
	   dtdr_s(ii) = dtdr(ii)
	 enddo  
	 do ii=1,nl
	    dtdv_s(ii) = dtdv(ii)
	 enddo   
	 if(nsp.eq.2) then
	    k2 = k2 - 1
	 endif   
	 goto 10
      endif
      if(sphase(nobs,i).eq.2.0) then
c
c        calculate the s-p traveltime derivative 
c
         do ii=1,3
            dtdr(ii)=dtdr_s(ii) - dtdr(ii)
         enddo
c
c        calculate the s-p time difference
c
         ttt = ttt_s - ttt
      endif	 
c
c     calculate residual by using appropriate station-corrections :
c
      if (sphase(nobs,i).ne.2.0) then
         pobs=pt(nobs,i)-e(1,i)
      else
         pobs=pt(nobs,i)
      endif	 
      extrat1=ptcor(k1)
      extrat2=0.0
      if(nsp.eq.2.and.sphase(nobs,i).eq.1.0)
     &extrat1=stcor(k1)
      if(nsp.eq.2.and.sphase(nobs,i).eq.2.0) then  ! s-p phases
         extrat1=stcor(k1) - ptcor(k1)
      endif	 
      if(nsp.eq.3.and.sphase(nobs,i).eq.1.0)
     &extrat1=ptcor(k1)*vpvs
      if(nsp.eq.3.and.sphase(nobs,i).eq.2.0) then ! s-p phases
         extrat1=ptcor(k1)*vpvs - ptcor(k1)
      endif
      if(nshcor.eq.0) goto 600
      if(i.le.neqs) goto 600
      kk1=map2(i-neqs)
      if(kk1.eq.0) goto 600
      do 601 j=1,nsta
      if(kk1.eq.map1(j)) goto 602
601   continue
      goto 600
c
602   extrat2=ptcor(j)
      if(nsp.eq.2.and.sphase(nobs,i).eq.1.0)
     &extrat2=stcor(j)
      if(nsp.eq.3.and.sphase(nobs,i).eq.1.0)
     &extrat2=ptcor(j)*vpvs
600   continue
      res(nobs,i)=pobs-(ttt+extrat1+extrat2)
c
      if(isingle.ne.0)then
         tcalc(nobs)=ttt
         takeoff_angle=
     &   SQRT( (rp(1,2)-rp(1,1))**2
     &        +(rp(2,2)-rp(2,1))**2 )
     &   /
     &   SQRT( (rp(1,2)-rp(1,1))**2
     &        +(rp(2,2)-rp(2,1))**2
     &        +(rp(3,2)-rp(3,1))**2 )
         takeoff_angle=57.296*ASIN(takeoff_angle)
         if( (rp(3,2)-rp(3,1)) .lt. 0.0 )then
            takeoff_angle=180.-takeoff_angle     ! ray is going upwards
         endif
         iain(nobs)=NINT(takeoff_angle)  ! with respect to positive z; downwards
      endif
c
c---- save residual according to the ray-type:
c
      if(isingle.eq.0) call RESISAVE(nrp,nrpdeep,rp,nobs,i,k1,mll)
c
c
c     if VELEST is used in single-event-mode (isingle <> 0 ),
c     set weight to 0 if abs(residual) after 2nd iteration is still > 2. sec
c     ... and normalize the weights once more !
c     --> if a reading with a reseted weight gets a residual which has
c         become small again, 'revive' it !
c
cek      if(isingle.ne.0)then
cek         if(nitt.gt.2.and.abs(res(nobs,i)).gt.2.0
cek     &               .and.igap(i).lt.250)then
cek            call REJECTOBS(i,nobs,iresflag)
cek            call GAPCALC(1)
cek         endif
cek         if(nitt.gt.2.and.abs(res(nobs,i)).lt.1.0
cek     &               .and.igap(i).lt.250)then
cek            if(w(nobs,i).eq.0.0)then
cek               call REVIVEOBS(i,nobs,iresflag)
cek               call GAPCALC(1)
cek            endif
cek         endif
cek      endif
c
      tctime(nobs,i)=ttt
      h(1)=x(k1,3)-tkh ! reset h(1) after it has been altered in subr. raypath!!
      thk(1)=thk(1)+tkh
      if(idrvout.eq.1)then
         write(21,*)'i=',i,'  nobs=',nobs,'  nitt=',nitt,'  nrtn=',nrtn
         write(21,500) ttt,pobs,res(nobs,i)
         write(21,501) (dtdr(jzz),jzz=1,3)
         write(21,502) (dtdv(jzz),jzz=1,nl)
 500     format(' ttt=',f10.5,' pobs= ',f10.5,' res= ',f10.5)
 501     format('dtdr=',3f10.5)
 502     format('dtdv=',7f10.5)
      endif
510   continue
      return
      end ! of subr. traveltime
