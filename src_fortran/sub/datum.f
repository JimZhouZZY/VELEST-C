c
c
c
c
c
c
      subroutine DATUM(ITF,IYR,IMO,IDY,IHR,IMN)
C UMRECHNEN DES DATUMS IN MINUTEN (CF. JULIAM) IN YR-MO-DY-HR-MI
C   (MIT IMN<2**31, JAHR < 4000
c
      implicit none
      integer iyr,imo,idy,ihr,imn
c      
      integer id,l,iyr4,iyrh,iyrt,ld,i
      integer KMO(12)
      INTEGER ITF,K,KH
      DATA KMO/31,28,31,30,31,30,31,31,30,31,30,31/
      K= ITF/60
      IMN= ITF-K*60
      KH= K/24
      IHR= K-KH*24
      IYR= KH/365
5     ID= KH-IYR*365
      L= 0
      IYR4= IYR/4
      IYRH= IYR/100
      IYRT= IYR/1000
      LD= IYR4-IYRH+IYRT
      IF(IYR4*4.EQ.IYR.AND.(IYRH*100.NE.IYR.OR.IYRT*1000.EQ.IYR)) L= 1
      ID= ID-LD+L
      IF(ID.GT.0) GOTO 10
      if(id.eq.0.and.ihr.eq.0.and.imn.eq.0) then
          idy= 0
          imo= 0
          return
      endif
      IYR= IYR-1
      GOTO 5
10    KMO(2)= 28+L
      DO 20 I=1,12
      ID= ID- KMO(I)
      IF(ID.LE.0) GOTO 30
20    CONTINUE
      I=12
30    IDY= ID+KMO(I)
      IMO= I
      RETURN
      end ! of subr. datum
c
      integer function TRIMLEN(t)   ! Urs Kradolfer, June 1986
c
c     Call:    nc=TRIMLEN(char)
c
c          --> nc says, how many characters the input-string has
c              (ignoring trailing blanks!).
c
      implicit none
      character t*(*)
      do 1 trimlen=LEN(t),1,-1
    1    if(t(trimlen:trimlen).ne.' ')RETURN
      trimlen=1
      end ! of integer function trimlen
