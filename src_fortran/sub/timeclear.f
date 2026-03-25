c
c
c
c
c
c
      subroutine TIMECLEAR(iyr,imo,iday,ihr,imin,sec,itime)
c
c     Urs Kradolfer, Winter 1987/88
c
      implicit none
c
      integer iyr,iyr1,imo,iday,ihr,imin, juliam,itime
      real sec,sec1
c
      if(iyr.lt.1900)then
         iyr1=iyr+1900
      else
         iyr1=iyr
      endif
      itime=JULIAM(iyr1,imo,iday,ihr,imin)
   1  sec1=sec
      if(sec.lt.0)then
         sec=sec1+60.
         itime=itime-1
      endif
      if(sec.gt.60.)then
         sec=sec1-60.
         itime=itime+1
      endif
      if(sec.lt.0.0.or.sec.gt.60.0) goto 1
      call DATUM(itime,iyr1,imo,iday,ihr,imin)
      if(iyr.lt.1900)then
         iyr=iyr1-1900    ! if input-year was of form  87  give it out so !
      else
         iyr=iyr1
      endif
c
      return
c
      end ! of subr. timeclear
c
      INTEGER FUNCTION JULIAM(IYR,IMO,IDY,IHR,IMN)
C UMRECHNEN VON JAHR-MONAT-TAG-STUNDEN-MINUTEN IN MINUTEN:
C   (WENN IMN 4-BYTE INTEGER, DANN JAHR < 4000)
      implicit none
      integer iyr,imo,idy,ihr,imn
      integer KMO(12)
      integer leap,ky,km,kd,ky4,ky1,ky0,kl,l
      DATA KMO/0,31,59,90,120,151,181,212,243,273,304,334/
      DATA LEAP/1/
      KY= IYR
      KM= IMO
      KD= IDY
      IF(KM.LE.0) KM= 1
10    JULIAM= 365*KY
      KD= KMO(KM)+KD
      KY4= KY/4
      KY1= KY/100
      KY0= KY/1000
      KL= LEAP*(KY4-KY1+KY0)
      L=0
      IF(KY4*4.EQ.KY.AND.(KY1*100.NE.KY.OR.KY0*1000.EQ.KY))L= LEAP
      IF(L.NE.0.AND.KM.LT.3) KL= KL-LEAP
      JULIAM= JULIAM+KD+KL
      JULIAM= JULIAM*24+IHR
      JULIAM= JULIAM*60+IMN
      return
      END ! of integer function juliam
