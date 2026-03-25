c
c
c
c
c
c
      subroutine SETUNT(nitt,invertratio,nsinv,icount,
     &                       xythet,stathet,othet,vthet,zthet,scale)
c
c     set scale factors
c
      implicit none
      integer nitt,invertratio,nsinv,icount
      real xythet,stathet,othet,vthet,zthet
      real scale(7)
c
c***  set damping so that velocities are adjusted every
c***  "invertratio" iterations
      icount=mod(nitt,invertratio)
c
c     scale(1) : origin-time
c     scale(2) : x
c     scale(3) : y
c     scale(4) : z
c     scale(5) : stacorr
c     scale(6) : veloc.model
c     scale(7) : (not used)
c
      scale(1)=1.0
      scale(7)=1.0   ! not used in this version!!!
      scale(2)=sqrt(othet/xythet)
      scale(3)=scale(2)
      scale(4)=sqrt(othet/zthet)
      scale(5)=0.0
      if(nsinv.ne.0.and.icount.eq.0) scale(5)=sqrt(othet/stathet)
      scale(6)=0.0
      if(icount.eq.0) scale(6)=sqrt(othet/vthet)
      return
      end ! of subr. setunt
