c
c
c
c
c
c
      subroutine RMSDATVAR
c
c     compute RMS and DATVAR for all the events (EQs & shots)
c
      implicit none
      include '../inc/vel_com.inc'
      integer jj,i,j2,j,knobst
      real tres
c
c
c   now calculate rms for all events
c
      jj=0
      davar1=0.0
      xmsqrs1=0.0
      tres=0.0
      knobst=0
      do 11 i=1,legs
         j2=knobs(i)
         knobst=knobst+j2
         do 12 j=1,j2
            if(w(j,i).le.0.0) goto 12
            avres(i)=avres(i)+res(j,i)*w(j,i)
Cek  changed by ek to res*w*res*w::::>
            rms(i)=rms(i)+ (res(j,i)*res(j,i))*w(j,i)*w(j,i)
 12      continue
         davar1=davar1+rms(i)
         tres=tres+avres(i)
cek next statement
         if( (j2-nobswithw0) .le.1) goto 11
         rms(i)=sqrt( (rms(i)-avres(i)**2/(j2-nobswithw0) )
     &               /( (j2-nobswithw0) -1))
         avres(i)=avres(i)/(j2-nobswithw0)
 11   continue
      if(nitt.eq.0)then
         if(.not.single_turbo)then
            write(16,310) knobst
         endif
      endif
310   format(/,' Total number of observations is: ',i5,/)
      if(nitt.eq.0.and.isingle.ne.0)then
         if(.not.single_turbo)then
            write(16,*)'Number of observations with '//
     &                 'normalized weight 0.0 :', nobswithw0
            write(16,*)'knobs(i)   = ',knobst
            write(16,*)'nobswithw0 = ',nobswithw0
         endif
      endif
c
c     calculate data variance:
c
c     NVAReff is the actual number of unknowns to solve for !!!
c
      if( (knobst-nobswithw0).gt.nvareff)then
        davar1=(davar1-tres*tres/(knobst-nobswithw0))
     &          / ((knobst-nobswithw0)-nvareff)
        xmsqrs1=davar1*((knobst-nobswithw0)-nvareff)/(knobst-nobswithw0)
      else
         davar1=999.99
         xmsqrs1=999.99
      endif
c
      i f (.not.single_turbo) t h e n
      if(isingle.eq.0)then
         write(16,*)
         write(16,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         write(16,*)'Events with  | AVRES | > 1.0 SEC are suspicious !'
         write(16,*)
         j=0
         do i=1,legs
            if(abs(avres(i)).gt.1.0)then
               write(16,'('' Event# '',i3,'' >>> '',1x,3i2.2,1x,2i2.2,
     &                    ''  AVRES ='',f6.2,'' NOBS ='',i3)')
     &                    i,iyr(i),imo(i),iday(i),ihr(i),imin(i),
     &                    avres(i),knobs(i)
               j=j+1
            endif
         enddo
         if(j.gt.0)then
         write(16,*)
         write(16,*)'^^^^^^^^^ C H E C K   these events above ^^^^^^^^'
         write(16,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         else
            write(16,*)'ZERO events of this kind found! (lucky guy!!!)'
         endif
         write(16,*)
      endif
      e n d i f
c
      return
      end ! of subr. rmsdatvar
