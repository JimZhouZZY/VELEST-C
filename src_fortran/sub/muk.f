c
c
c
c
c
c
      subroutine MUK(epdist,depth,ifilter,isecpendel,seismkonst,
     &               seismdamp,voltgain,stacor,devampl,period,xmag)
c
c------------------------------- Urs Kradolfer 1984 -------------------
C
C     Berechnet die Magnitude Ml fuer Epizentraldistanz <=  700 km
C     Berechnet die Magnitude Mb fuer Epizentraldistanz > 2200 km
c
c     Input:
c             Epizentraldistanz         [km]       epdist  >=0.  {real}
c                    "                  [deg]      epdist  < 0.  {real}
c             Herdtiefe                 [km]       depth         {real}
c             filter-typ (z.b. DE, AD, PD)         filter-typ    {char}
c             Seismometer-Grenzperiode  [sec]      isecpendel    {integer}
c                  "     -Konstante     [V/(cm/s)] seismkonst    {real}
c                  "     -Daempfung     [1]        seismdamp     {real}
c             Verstaerkung (Elektronik) [1]        voltgain      {real}
c             Stationskorrektur         [Mag]      stacor        {real}
c             Max. Ablese-Amplitude     [mm]       devampl       {real}
c                 (Devco oder Plot)
c             dominierende Periode      [sec]      period        {real}
c
c     Output: Magnitude (Ml oder Mb)    [Mag]      xmag          {real}
c
c        >>   Falls die Magnitude nicht berechnet werden kann, erfolgt
c        >>   ein RETURN und der Output-Parameter xmag wird -13. gesetzt!
c
c     Zur Subr. muk gehoert ebenfalls die Subr. UFWABO und die
c     zur letzteren gehoerenden vier Complex-Functions.
c---------------------------------------------------------------------------
c
      implicit none
      real epdist,depth,seismkonst,seismdamp,voltgain
      real stacor,devampl,period,xmag
      integer isecpendel
c
      real epdistkm,ampl,waampl,boampl,delta,sigma,a,b
      integer iampltype,i
      character*2 ifilter
      real RDELT(12), SIGAR(12)
      DATA RDELT/4.0,6.,13.,16.,28.,87.,
     & 114.,120.,134.,141.,146.,180./,
     &     SIGAR/6.2,7.,7.0,5.8,6.6,7.0,7.5,
     & 7.5, 6.9, 7.1, 6.9, 6.9 /
      if(isecpendel.le.0.or.seismkonst.le.0..or.seismdamp.le.0..or.
     &   voltgain.le.0..or.devampl.le.0..or.period.le.0.)then
         xmag=-13.    ! Magnitude kann nicht berechnet werden
         return
      endif
      xmag=-13.  ! falls xmag nicht anders berechnet wird !!
      if(epdist.lt.0.)then
         epdistkm= -epdist/360.*40030     ! epdist < 0. --> [epdist] = Grad
      else
         epdistkm=sqrt(epdist**2+depth**2)
      endif
      if(epdistkm.le.2200.)iampltype=1    ! Wood-Anderson --> Ml
      if(epdistkm.gt.2200.)iampltype=2    ! Bodenbewegung --> Mb
      if(ifilter.eq.'DE'.or.ifilter.eq.'AD')
     &call ufwabo(isecpendel,seismkonst,seismdamp,voltgain,
     &            devampl,period,iampltype,ampl,ifilter)
      if(ifilter.eq.'PD')
     &call mpdr2(isecpendel,seismkonst,seismdamp,voltgain,
     &            devampl,period,iampltype,ampl)
      if(ampl.le.1.0e-10)then
         xmag=-13.  ! avoid taking LOG of non-positive value !!!
         RETURN
      endif
      if(iampltype.eq.1)waampl=ampl/2. ! p-p-Ampl. --> 0-p-Ampl.
      if(iampltype.eq.2)boampl=ampl
      if(epdistkm.gt.2200.)then
         delta=epdistkm/40030.*360.     ! delta:=epdistkm in Grad
         if(delta.gt.180.)then
            delta=180.
            sigma=sigar(12)
            goto 2
         endif
         i=0
1        i=i+1
         if(delta.gt.rdelt(i))goto 1
         ! jetzt ist delta < rdelt(i)     --> Linearisieren
         ! y=a*x+b
         a=(sigar(i)-sigar(i-1))/(rdelt(i)-rdelt(i-1))
         b=sigar(i)-a*rdelt(i)
         sigma=a*delta+b
2        xmag=log10(boampl/period) + sigma + stacor
         RETURN
      endif
      if(epdistkm.gt.0..and.epdistkm.le.60.)then
         xmag=log10(waampl) + 0.018 *epdistkm+1.77 + 0.40
      endif
      if(epdistkm.gt.60..and.epdistkm.le.700.)then
         xmag=log10(waampl) + 0.0038*epdistkm+2.62 + 0.40
      endif
      if(epdistkm.gt.1100..and.epdistkm.le.1700.)then
         xmag=log10(waampl) + 0.0029*epdistkm+3.40 + 0.40 - 2.  ! EMPIRISCH
      endif
      if(xmag.eq.-13.) RETURN     ! xmag konnte nicht berechnet werden
      xmag=xmag+stacor   ! Stationskorrektur Cs
c
      return
      end
