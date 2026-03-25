c
c
c
c
c
c
      subroutine CHECKSOLUTION(istopflag,better)
c
c     checks, whether the current solution is 'better' than the previous one.
c     'better' means here, that data variance has decreased. A warning is
c     output in case that although the data variance has decreased the RMS
c     didn't.
c     If datvar has decreased, the output-variable DECREASING is set to 1  .
c
      implicit none
      include '../inc/vel_com.inc'
c
      integer istopflag
      real datvar,xmsqrs2,varat1,varat2
      integer decreasing
      logical better
c
      save datvar,xmsqrs2
c
      decreasing=0
      better=.true.
c
      if(nitt.eq.0)then
         datvar=davar1
         xmsqrs2=xmsqrs1
      endif
      if(.not.single_turbo)then
         write(16,*)
         if(ibackups.eq.0)then
            write(16,'('' Iteration nr '',i2,'' obtained:'')') nitt
         else
            write(16,'(''(Iteration nr '',i2,'')   BACKUP nr '',i1,
     &                 '' obtained:'')') nitt,ibackups
         endif
         write(16,10) davar1,xmsqrs1,sqrt(xmsqrs1)
 10      format(' DATVAR=',f12.6,' mean sqrd residual= ',f12.6,
     &          '  RMS RESIDUAL= ',f12.6)
      endif
c
      if(isingle.ne.0) call GAPCALC(1)
c
      if(isingle.eq.0) call AVRESISTATIST
c
c***  form ratio of old to new variance
      varat1=datvar/davar1
c
      if(isingle.ne.0)then ! stop calc. for single event if changes get small
         if(nitt.gt.2.and.abs(datvar-davar1).lt.1e-6)then
            istopflag=1
            if(.not.single_turbo)then
               write(16,*)'Changes in datvar < 1e-6   : STOPPING...'
            endif
         endif
      endif
c*** form ratio of old to new mean sqrd residual
      varat2=xmsqrs2/xmsqrs1
c*** test to see if variance is increasing
      if(varat1.ge..99)then
         decreasing=1
         goto 1
      endif
c*** test to see if mean sqrd residual is increasing
      if(varat1.lt..99.and.varat2.ge..99)then
c
c       msqrd residual decreasing, but data variance increased:
c
         if(.not.single_turbo)then
        write(16,'('' *** WARNING: the data variance has increased.'')')
         endif
        decreasing=1
        goto 1
      endif
      if(varat2.ge..99)then
         decreasing=1
         goto 1
      endif
c
   1  if(decreasing.eq.1)then
c
c        variance is decreasing:
c
         datvar=davar1
         xmsqrs2=xmsqrs1
      else
c
c        variance increased so backup
c
         better=.false.
      endif
c
      RETURN
c
      end ! of subr. checksolution
