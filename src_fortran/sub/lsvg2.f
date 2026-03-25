c
c
c
c
c
c
      subroutine LSVG2  (DCOS,DSIN,X,Y)
c
C   IMSL ROUTINE NAME   - LSVG2
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
C
      implicit none
C                                  SPECIFICATIONS FOR ARGUMENTS
      real                         DCOS,DSIN,X,Y
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      real                         XR
C                                  FIRST EXECUTABLE STATEMENT
      XR=DCOS*X+DSIN*Y
      Y=-DSIN*X+DCOS*Y
      X=XR
      RETURN
      END
