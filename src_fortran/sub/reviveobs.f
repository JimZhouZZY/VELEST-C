c
c
c
c
c
c
      subroutine REVIVEOBS(i,nobs,iresflag)
c
c     If observationweight < 4 and weight=0.0 then set 'original' weight
c     if residual has become smaller since it has been reseted:
c
      implicit none
      integer i,nobs,iresflag
      include '../inc/vel_com.inc'
c
      real wsum,xkndw
      integer knobst,iii,iobswt
      if(w(nobs,i).ne.0.0) RETURN
      if(kpwt(nobs,i).eq.4) RETURN
c
cc         w(nobs,i)=0.0039  !  = 1./(2**(2*reading_weight)) = 1./256
c
c     re-normalize weights of this event :
c
      nobswithw0=nobswithw0-1
      wsum=0.0
      knobst=knobs(i)
      do iii=1,knobst
         iobswt=kpwt(iii,i)
         if(iii.ne.nobs.and.w(iii,i).eq.0.0)then
            w(iii,i)=0.0   ! keep zero-weights as they are !!!
         else
            if(iobswt.lt.4)then
               if(sphase(iii,i).eq.1.0.or.sphase(iii,i).eq.2.0)then
                  w(iii,i)=swtfac*1.0/(2**(iobswt*2))   ! S-phase
               else
                  w(iii,i)=       1.0/(2**(iobswt*2))   ! P- or M-phase
               endif
            else
               w(iii,i)=0.0   ! observation-weight 4 ==> don't use this arrival
            endif
         endif
         wsum=wsum+w(iii,i)
      enddo
      xkndw=float(knobst-nobswithw0)/wsum
      do iii=1,knobst
         w(iii,i)=w(iii,i)*xkndw
      enddo
      if(.not.single_turbo)then
      write(16,*)'WARNING:'
      write(16,'('' Iteration# '',i3,'' Event#'',i4,'' Res#'',i3,
     &           '' ='',f7.2,''; ABS < 1.0 ---> weight revived !'')')
     &           nitt,isingle,nobs,res(nobs,i)
      endif
      write(6,'('' Iteration# '',i3,'' Event#'',i4,'' Res#'',i3,
     &           '' ='',f7.2,''; ABS < 1.0 ---> weight revived !'')')
     &           nitt,isingle,nobs,res(nobs,i)
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
      end ! of subr. reviveobs
