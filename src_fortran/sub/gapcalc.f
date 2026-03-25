c
c
c
c
c
c
      subroutine GAPCALC(i)
c
c     determine GAP for one event
c
      implicit none
      integer i
      include '../inc/vel_com.inc'
      integer nofgaps,j,ig
      real xstn,ystn,xhyp,yhyp,dxstnhyp,dystnhyp
c
      integer iga(200)     !  max. 200 obs. for a single event
c
c     event: i  knobs(i)  +/-e(2,i)=x    e(3,i)=y    x=+/-x(nsta,1)  y=...,2)
c     station# are stored in array  ISTM(iobs,ievent)
c
c---- compute GAP for this event
      if(knobs(i).gt.200)then
         i f (.not.single_turbo)then
         write(16,*)' WARNING:'
         write(16,*)' Event# ',i,'   Nobs = ',knobs(i),' > 200 !!!'
         write(16,*)' Array IGA(100) is too small; redimension it in'
         write(16,*)' subr. GAPCALC'
         e n d i f
         write(6,*)' Event# ',i,'   Nobs = ',knobs(i),' > 200 !!!'
         write(6,*)' Array IGA(200) is too small; redimension it in'
         write(6,*)' subr. GAPCALC'
         stop'subr. GAPCALC >>> array IGA is too small !'
      endif
      nofgaps=0
      do j=1,knobs(i)
         if(w(j,i).gt.0.0)then
            nofgaps=nofgaps+1
            xstn=x(istm(j,i),1)
            ystn=x(istm(j,i),2)
            xhyp=e(2,i)
            yhyp=e(3,i)
cek    avoiding call to atan2(zero,zero)
cek              write(6,*) ' Gapcalc: i,nofgaps,xstn,ystn,xhyp,yhyp'
cek              write(6,'(1x,i4,2x,i3,2x,4f10.3)') i,nofgaps,
cek     +                  xstn,ystn,xhyp,yhyp
            dxstnhyp=abs(xstn-xhyp)
            dystnhyp=abs(ystn-yhyp)
            if(dxstnhyp.gt.0.0001.or.dystnhyp.gt.0.0001) then
              iga(nofgaps)=57.296*ATAN2(xstn-xhyp,ystn-yhyp)
            else
              iga(nofgaps)=359
            endif
cek
            if(iga(nofgaps).lt.0) iga(nofgaps)=iga(nofgaps)+360
         endif
      enddo
      if(nofgaps.gt.0)then
         call SORTI(iga,nofgaps)
      else
         write(6,*)'WARNING: Event-# :',i,'has zero observations!'
         if(.not.single_turbo)then
            write(16,*)'WARNING: Event-# :',i,'has zero observations!'
         endif
      endif
      igap(i)=0
      ig=iga(1)-iga(nofgaps)
      if(ig.lt.0) ig=ig+360
      if(ig.gt.igap(i)) igap(i)=ig
      do j=2,nofgaps
         ig=iga(j)-iga(j-1)
         if(ig.lt.0) ig=ig+360
         if(ig.gt.igap(i)) igap(i)=ig
      enddo
c---- IGAP(i) is the gap of this event (nr. i)
c
      return
      end ! of subr. gapcalc
