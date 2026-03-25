c
c
c
c
c
c
      subroutine SORTI(iX,NO)   ! Urs Kradolfer, Okt. 1986 (aus HYPOSUBR.FOR )
c
c  SORTA sortiert einen Integer-Array iX so, dass iX(1) < iX(2) < ... < iX(N)
c
c  Aufruf: call SORTI(iX,N)
c
c  mit INTEGER iX(N)
c      INTEGER N       (Anzahl Elemente in X)
c
c  Array iX wird in SORTI veraendert !!
c
      implicit none
      integer i,no,mo,ko,jo
      integer iX(no),itemp
c
      if(no.eq.0) return
      MO=NO
 2    IF (MO-15) 21,21,23
 21   IF (MO-1) 29,29,22
 22   MO=2*(MO/4)+1
      GO TO 24
 23   MO=2*(MO/8)+1
 24   KO=NO-MO
      JO=1
 25   I=JO
 26   IF (iX(I)-iX(I+MO)) 28,28,27
 27   iTEMP=iX(I)
      iX(I)=iX(I+MO)
      iX(I+MO)=iTEMP
      I=I-MO
      IF (I-1) 28,26,26
 28   JO=JO+1
      IF (JO-KO) 25,25,2
 29   continue
      return
      END ! of subr. sorti
