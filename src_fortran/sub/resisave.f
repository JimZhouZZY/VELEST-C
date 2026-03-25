c
c
c
c
c
c
      subroutine RESISAVE(nrp,nrpdeep,rp,nobs,i,k1,mll)
c
c     saves the residual according to the ray-type
c
      implicit none
      integer nrp,nrpdeep,nobs,i,k1,mll
      include '../inc/vel_com.inc'
c
      real xxx,yyy,xhyp,yhyp,xstn,ystn,azi,dist
      integer iazi
      real rp(3,inrpmax)
c
      if(rp(3,nrpdeep).eq.rp(3,nrpdeep+1))then ! ray is horizontal : refracted !
         avrefrres=avrefrres+    res(nobs,i)       ! refracted
         abrefrres=abrefrres+ABS(res(nobs,i))      ! refracted
         nrrefrres=nrrefrres+1
c
         if(lmax.eq.10)then ! refractor is Moho...
            if(irfrout.eq.1)then
c
c              write (x,y) of point A, where ray enters MOHO,
c              of point B, where ray leaves MOHO and of point between A and B
c
               xxx=rp(1,nrpdeep)
               if(icoordsystem.eq.2) xxx=-xxx
               yyy=rp(2,nrpdeep)
               write(78,'(2x,''RFR'',3f10.3,''        .1'')')
     &                    xxx,yyy,res(nobs,i)
               xxx=(rp(1,nrpdeep)+rp(1,nrpdeep+1))/2.
               if(icoordsystem.eq.2) xxx=-xxx
               yyy=(rp(2,nrpdeep)+rp(2,nrpdeep+1))/2.
               write(78,'(2x,''RFR'',3f10.3,''        .1'')')
     &                    xxx,yyy,res(nobs,i)
               xxx=rp(1,nrpdeep+1)
               if(icoordsystem.eq.2) xxx=-xxx
               yyy=rp(2,nrpdeep+1)
               write(78,'(2x,''RFR'',3f10.3,''        .1'')')
     &                    xxx,yyy,res(nobs,i)
             endif
         endif
c
      else
         if(mll.eq.0)then
            avotheres=avotheres+    res(nobs,i)    ! straight or direct
            abotheres=abotheres+ABS(res(nobs,i))   ! straight or direct
            nrotheres=nrotheres+1
         else
            avreflres=avreflres+    res(nobs,i)    ! reflected
            abreflres=abreflres+ABS(res(nobs,i))   ! reflected
            nrreflres=nrreflres+1
c
            if(irflout.eq.1)then
              if(icoordsystem.eq.2)then
                 xxx=-rp(1,nrpdeep)
              else
                 xxx=rp(1,nrpdeep)
              endif
              write(77,'(2x,''RFL'',3f10.4,''        .1'')')
     &                   -rp(1,nrpdeep),rp(2,nrpdeep),res(nobs,i)
            endif
         endif
      endif
c
c     save residual for station# k1 according to the azimuth:
c
      xhyp=rp(1,1)
      yhyp=rp(2,1)
      xstn=rp(1,nrp)
      ystn=rp(2,nrp)
c        Azimuth (stn --> hypoc) = 57.296*ATAN2(xhyp-xstn,yhyp-ystn)
      azi=57.296*ATAN2(xhyp-xstn,yhyp-ystn)
      if(azi.lt.0) azi=azi+360.
      iazi=0
      if(azi.ge.  0.0 .and. azi.lt. 90.0) iazi=1
      if(azi.ge. 90.0 .and. azi.lt.180.0) iazi=2
      if(azi.ge.180.0 .and. azi.lt.270.0) iazi=3
      if(azi.ge.270.0 .and. azi.le.360.0) iazi=4
      stnazires(k1,2*iazi-1)=stnazires(k1,2*iazi-1)+res(nobs,i)
      stnazires(k1,2*iazi)=stnazires(k1,2*iazi)+1.
c
c     type residuals as a function of focus-receiver distance:
      if(iresout.eq.1)then
         dist=(  (rp(1,1)-rp(1,nrp))**2
     &          +(rp(2,1)-rp(2,nrp))**2
     &          +(rp(3,1)-rp(3,nrp))**2  )
         dist=sqrt(dist)
         write(79,'(1x,f6.2,2x,f7.3)') dist,res(nobs,i)
      endif
c
      RETURN
      end ! of subr. resisave
