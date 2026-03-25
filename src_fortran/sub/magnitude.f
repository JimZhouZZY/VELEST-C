c
c
c
c
c
c
      subroutine MAGNITUDE(ievent,itime)   ! Urs Kradolfer, July 1987
c
c     calculate magnitude for an event with depth Z and NRP observations
c
      implicit none
      integer ievent,itime
      include '../inc/vel_com.inc'
c
      integer iseis,iscon,isdmp,isamp,ier,i
      real cormag,sconst,sdampf,voltgain,xmag
      character*2 ifilt
      character staname*4
c
      nmag=0
      xmagnitude=0.
      sdxmagnitude=0.
      do i=1,knobs(ievent)
         if(kpwt(i,ievent).gt.4) goto 40 ! weights > 4 are no readings!!!
c                              station# are stored in array  ISTM(iobs,ievent)
c                              smn(iobs,iev)  is the station-name
         write(staname,'(a4)') smn(i,ievent)
         call STINP(itime,staname,ifilt,iseis,iscon,
     &                            isdmp,isamp,cormag,ier)
         if(ier.lt.0) goto 40  ! seismo-parameters of station not found
         if(iyr(ievent).ge.1984.and.ifilt.eq.'DE') ifilt='AD'
         sconst=float(iscon)/10.
         sdampf= float(isdmp)/10.
         voltgain=float(isamp)
         delta=SQRT(  (e(2,ievent)-d(i,1,ievent))**2
     &              + (e(3,ievent)-d(i,2,ievent))**2  )
         xmag=0.
         call MUK(delta,e(4,ievent),ifilt,iseis,sconst,sdampf,
     &            voltgain,cormag,amx(i),prx(i),xmagni(i))
         if(xmagni(i).eq.-13.) goto 40
         nmag=nmag+1
         xmag=xmagni(i)
         xmagnitude=xmagnitude+xmag
         sdxmagnitude=sdxmagnitude+xmag**2
  40     continue
      enddo
      if(nmag.ne.0)then
         if(nmag.ge.2)then
            sdxmagnitude=
     &      SQRT( (nmag*sdxmagnitude-xmagnitude*xmagnitude)
     &           / (nmag*(nmag-1)) )
         else
             sdxmagnitude=0.
         endif
         xmagnitude=xmagnitude/nmag
      else
         xmagnitude=0.0
      endif
c
      return ! calculated magnitude is: AVXM +/- SDXM   and all xmagni(1...nobs)
c
      END ! of subr. magnitude
