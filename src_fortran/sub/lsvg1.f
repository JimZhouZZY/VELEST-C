c
c
c
c
c
c
      subroutine LSVG1  (A,B,DCOS,DSIN,SIG)
c
C   IMSL ROUTINE NAME   - LSVG1
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - HP1000/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE LSVDB
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      implicit none
C                                  SPECIFICATIONS FOR ARGUMENTS
      real                         A,B,DCOS,DSIN,SIG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      real                         AA,BB
C                                  FIRST EXECUTABLE STATEMENT
      double precision doubleb, doubleaa, doubler
      IF (ABS(A).LE.ABS(B)) GO TO 5
      AA=ABS(A+A)
cuk   next statement caused an underflow! u.kradolfer, 16.7.81
cuk      SIG=AA*SQRT(0.25e0+(B/AA)**2)
cuk   next statement produced a compiler warning with range checking option
cuk      SIG=AA*sngl( dSQRT(0.25d0+(dble(B)/dble(AA))**2) ) ! u.k. 7.2.92
      doubleb=DBLE(b)
      doubleaa=DBLE(aa)
      doubler=dSQRT(0.25d0+(doubleb/doubleaa)**2)
      sig=AA*sngl(doubler)
      DCOS=A/SIG
      DSIN=B/SIG
      RETURN
    5 IF (B.EQ.0.0e0) GO TO 10
      BB=ABS(B+B)
      SIG=BB*SQRT(0.25e0+(A/BB)**2)
      DCOS=A/SIG
      DSIN=B/SIG
      RETURN
   10 SIG=0.0e0
      DCOS=0.0e0
      DSIN=1.0e0
      RETURN
      END
