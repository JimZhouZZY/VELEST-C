c
c
c
c
c
c
      subroutine SORTR(X,NO)   ! Urs Kradolfer, Okt. 1986 (aus HYPOSUBR.FOR )
c
c  SORTR sortiert einen Real-Array X so, dass X(1) < X(2) < ... < X(N)
c
c  Aufruf: call SORTR(X,N)
c
c  mit REAL    X(N)
c      INTEGER N       (Anzahl Elemente in X)
c
c  Array X wird in SORTR veraendert !!
c
      implicit none
      integer i,no,mo,ko,jo
      real X(no),temp
c
      MO=NO
 2    IF (MO-15) 21,21,23
 21   IF (MO-1) 29,29,22
 22   MO=2*(MO/4)+1
      GO TO 24
 23   MO=2*(MO/8)+1
 24   KO=NO-MO
      JO=1
 25   I=JO
 26   IF (X(I)-X(I+MO)) 28,28,27
 27   TEMP=X(I)
      X(I)=X(I+MO)
      X(I+MO)=TEMP
      I=I-MO
      IF (I-1) 28,26,26
 28   JO=JO+1
      IF (JO-KO) 25,25,2
 29   continue
      return
      END ! of subr. sortr
