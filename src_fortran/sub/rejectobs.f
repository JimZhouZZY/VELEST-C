c
c
c
c
c
c
      subroutine REJECTOBS(i,nobs,iresflag)
c
c     reject the current observation for further calculation in case the
c     reading (or the seismologist who read it...) seems to be nuts:
c
      implicit none
      integer i,nobs,iresflag
      include '../inc/vel_com.inc'
c
      real wsum,xkndw
      integer knobst,iii,iobswt
      if(w(nobs,i).eq.0.0) RETURN
      if( (knobs(i)-nobswithw0) .eq. nvar ) RETURN
c
cc         w(nobs,i)=0.0039  !  = 1./(2**(2*reading_weight)) = 1./256
c
c     re-normalize weights of this event :
c
      nobswithw0=nobswithw0+1
      wsum=0.0
      knobst=knobs(i)
      do iii=1,knobst
         iobswt=kpwt(iii,i)
         if(iobswt.lt.4.and.w(iii,i).ne.0.0)then
            if(sphase(iii,i).eq.1.0.or.sphase(iii,i).eq.2.0)then
               w(iii,i)=swtfac*1.0/(2**(iobswt*2))   ! S-phase or s-p phase
            else
               w(iii,i)=       1.0/(2**(iobswt*2))   ! P- or M-phase
            endif
         else
            w(iii,i)=0.0   ! observation-weight 4 ==> don't use this arrival
         endif
         if(iii.eq.nobs) w(iii,i)=0.0
         wsum=wsum+w(iii,i)
      enddo
      xkndw=float(knobst-nobswithw0)/wsum
      do iii=1,knobst
         w(iii,i)=w(iii,i)*xkndw
      enddo
      if(.not.single_turbo)then
      write(16,*)'WARNING:'
      write(16,'('' Iteration# '',i3,'' Event#'',i4,'' Res#'',i3,
     &   '' ='',f7.2,''; ABS > 2.0 ---> weight set to zero !'')')
     &   nitt,isingle,nobs,res(nobs,i)
      endif
      write(6,'('' Iteration# '',i3,'' Event#'',i4,'' Res#'',i3,
     &   '' ='',f7.2,''; ABS > 2.0 ---> weight set to zero !'')')
     &   nitt,isingle,nobs,res(nobs,i)
      if(.not.single_turbo)then
         write(16,*)'knobs(i)   = ',knobst
         write(16,*)'nobswithw0 = ',nobswithw0
      endif
cc      write(6,*)'knobs(i)   = ',knobst
cc      write(6,*)'nobswithw0 = ',nobswithw0
cc      write(6,*)
      iresflag=1  ! inhibit a pre-stopping in subr. OUTPUT due to
                  ! increased datvar; give location a second chance.
c
      return
      end ! of subr. rejectobs
