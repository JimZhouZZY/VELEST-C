c
c
c
c
c
c
      subroutine AVRESISTATIST
c
c     compute average residual statistics for the iteration just finished
c
      implicit none
      integer ifirstrun
      integer nrtotres
      real abtotres,avtotres,proz,oldres
c
      include '../inc/vel_com.inc'
c
      save oldres, ifirstrun
c
      write(16,*)
      if(nrotheres.gt.0) avotheres=avotheres/nrotheres
      if(nrotheres.gt.0) abotheres=abotheres/nrotheres
c
      if(nrrefrres.gt.0) avrefrres=avrefrres/nrrefrres
      if(nrrefrres.gt.0) abrefrres=abrefrres/nrrefrres
c
      if(nrreflres.gt.0) avreflres=avreflres/nrreflres
      if(nrreflres.gt.0) abreflres=abreflres/nrreflres
c
      nrtotres=nrotheres+nrrefrres+nrreflres
c
      if(nrtotres.gt.0)abtotres=(abotheres*nrotheres+
     &   abrefrres*nrrefrres+abreflres*nrreflres)/nrtotres
c
      if(nrtotres.gt.0)avtotres=(avotheres*nrotheres+
     &   avrefrres*nrrefrres+avreflres*nrreflres)/nrtotres
c
      write(16,'(1x,''After'',i3,'' iterations we got:'')') nitt
      write(16,*)'Average absolute & unweighted [and mean] residual of'
      write(16,'(1x,i5,'' straight and direct rays ='',f9.5,'' ['',f9.5,
     &           '']'')')nrotheres,abotheres,avotheres
      write(16,'(1x,i5,'' refracted           rays ='',f9.5,'' ['',f9.5,
     &           '']'')')nrrefrres,abrefrres,avrefrres
      write(16,'(1x,i5,'' reflected           rays ='',f9.5,'' ['',f9.5,
     &           '']'')')nrreflres,abreflres,avreflres
      write(16,*)
      if(ifirstrun.ne.10000001)then ! first run; no 'oldres' available...
         ifirstrun=10000001
         proz=0.0
      else ! NOT the first time in this routine; proz may be calculated
         if(ABS(oldres).gt.1.0e-10)then
            proz=100.*(abtotres-oldres)/oldres
         endif
      endif
      write(16,'(1x,i5,'' ALL                 RAYS ='',f9.5,'' ['',f9.5,
     &           '']'',5x,f7.2,'' %'')')nrtotres,abtotres,avtotres,proz
      oldres=abtotres
      write(16,*)
c
      return
      end ! of subr. avresistatist
