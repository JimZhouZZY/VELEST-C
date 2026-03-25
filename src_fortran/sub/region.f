c
c
c
c
c
c
      subroutine REGION (ityp,x,y,cname,nreg,
     &                   regnamfile,regkoordfile) ! Urs Kradolfer, 6.7.87
c
c     originally from Manfred Baer
c
c     modified 6.7.87, 28.6.91 Urs Kradolfer 
c
      implicit none
c
      character*(*) regnamfile, regkoordfile
      character*32 cname, place
      integer ityp,nreg, iii, iregread
      real x,y, xlat,xlon
c
      save iregread
c
      if(iregread.eq.0)then
         call REGREAD(regnamfile,regkoordfile)
         iregread=1
      endif
c
      nreg=0
c
      if(ityp.eq.1) goto 1
      if(ityp.eq.2.or.ityp.eq.3) goto 2
      stop'REGION>>> illegal coordinate-type!'
    1 continue
      place=' '
      call REGCH(x,y,place,nreg)
      if(nreg.eq.0)then
         call GEOKO(x,y,xlat,xlon,1)
         goto 22
      else
         cname=place
      endif
      goto 999
c
    2 continue
      xlat=x
      xlon=y
  22  place=' '
c     call REGNU(0,xlat,xlon,GELAT,GELON,PLACE,NREG)
      if(ityp.eq.1) iii=0
      if(ityp.eq.2) iii=0
      if(ityp.eq.3) iii=1
      call REGWORLD(iii,xlat,xlon,place,nreg)
      if(place.eq.' ')then
         cname='*****'
      else
         cname=place
      endif
      goto 999
c
  999 continue
      return
      end ! of subr. region
