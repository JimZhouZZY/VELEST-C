c
c
c
c
c
c
      subroutine REFLECT(nl,v,vsq,thk,jl,tkj,delta,z,mll,trefl,ain,ierr,
     &             DTDDrefl,DTDHrefl)
c
c     Urs Kradolfer, Nov. 1986
c
      implicit none
      integer nl,ierr
      real tkj,trefl,ain,dtddrefl,dtdhrefl
c
      integer nit,m
      real zbig,fac,a,sina,dellit,dddas,sqr,ddda,da,da2
      real tt,dtda,ti,dti,dtdd,sina2,dtdh
      integer jl   ! hypocenter-layer
      integer mll  ! reflection at bottom of layer MLL
      real thk(nl),v(nl),vsq(nl)
      real vsqu(100),vi(100),div(100)  ! dimension = max. number of layers allowed
      real anin    ! auftauchwinkel in rad
      real test4   ! iterate until hypoc-adjustments are less than (TEST4)**-2
      real z       ! hypocenter-depth
      real delta  ! distance epicenter-receiver
      integer i
      real sum
c
      ierr=0
      DTDDrefl=0.0
      DTDHrefl=0.0
c
      test4=0.01                  !  see HYMD01.DAT
      do i=1,nl
         vi(i)=1./v(i)
      enddo
c
c---- determine, which layer contains hypocenter
   10 sum=0.0
      do jl=1,nl
         sum=sum+thk(jl)
         if(z.le.sum)goto 19
      enddo
c---- hypocenter is in the halfspace:
      jl=nl
      tkj=z-sum
      if(jl.eq.0) goto 10
      goto 20
   19 tkj=z-sum+thk(jl)
      goto 20
   20 continue
c     hypocenter is TKJ [km] below top of hypoc-layer
c     JL is the hypocenter layer
c
      if(jl.gt.mll)then   !  hypocenter below reflector !!!
         ierr=-1
         write(6,*)'WARNING:'
         write(6,*)'subr. REFLECT >>> hypocenter is below reflector!'
         write(16,*)'WARNING:'
         write(16,*)'subr. REFLECT >>> hypocenter is below reflector!'
ccc         stop'subr. REFLECT >>> hypocenter below reflector!'
      endif
      if(mll.gt.nl)then
         stop'subr. REFLECT >>> reflector-nr > number of layers !'
      endif
c
C
C  CALCULATE REFLECTED WAVES
C
      NIT= 0
      ZBIG= TKJ
      SUM= 0.0
      DO M=1,MLL
C LAYER M IS PASSED ONLY ONCE:
         FAC= 1.
         IF(JL.GT.M) GOTO 505
C LAYER M IS PASSED TWICE
         FAC= 2.
         IF(JL.LT.M) GOTO 505
C LAYER M CONTAINS HYPOCENTER:
         FAC= (2.*THK(M)-ZBIG)/THK(M)
  505    DIV(M)= FAC*THK(M)
         VSQU(M)= VSQ(MLL)/VSQ(M)
         SUM= SUM+DIV(M)
      enddo
C ITERATE TO FIND EMERGING ANGLE FOR THE GIVEN DISTANCE:
      A= ATAN2(delta,SUM)
  515 SINA= SIN(A)
      if(sina.eq.1.00) sina=0.9999999
      SINA2= SINA*SINA
      DELLIT= 0.
      DDDAS= 0.
      DO M=1,MLL
         SQR= VSQU(M)-SINA2
c
       if(sqr.le.0)then
         write(16,*)'WARNING:                  subr. REFLECT'
         write(16,*)'probably a low velocity layer detected!'
         write(16,*)'sina=',sina
         write(16,*)'sina2=',sina2
         write(16,*)'VSQU(M):  vsqu(',m,')=',vsqu(m)
         write(16,*)'MLL=',mll
         write(16,*)'JL=',jl
         write(6,*)'sina=',sina
         write(6,*)'sina2=',sina2
         write(6,*)'VSQU(M):  vsqu(',m,')=',vsqu(m)
         write(6,*)'MLL=',mll
         write(6,*)'JL=',jl
         write(6,*)
       endif
c
c!!!!!         if(sqr.le.0.) goto 530   !auftauchwinkel gefunden
c   sina2 is always < 1 ;
c   vsqu is always > 1 if no velocity-layers are above the event-layer !!
         DDDA= DIV(M)/SQRT(SQR)
         DELLIT= DELLIT+DDDA
         DDDAS= DDDAS+VSQU(M)*DDDA/SQR
      enddo
      DELLIT= DELLIT*SINA
      DDDAS= DDDAS*COS(A)
      DA= DELLIT-delta
      DA2= DA*DA
      IF(DA2.LE.TEST4) GOTO 530
      NIT= NIT+1
cuk      IF(NIT.LT.15) GOTO 525
      IF(NIT.LT.50) GOTO 525
c      WRITE(6,524) delta
c  524 FORMAT(' REFLECTED WAVE DID NOT CONVERGE WITHIN',
c     1' 50 ITERATIONS AT DISTANCE:',F8.1,
c     2' SET WEIGHT TO 0')
      GOTO 27
  525 A= A-DA/DDDAS
      if(a.ge.1.570796) a=1.5  ! 1.5 ~ 86 deg ;   a = never .gt.  pi/2.  !!!
      GOTO 515
c---- successful iteration:
  530 TT=0.
      DTDA= 0.
      DO M=1,MLL
         TI= DIV(M)*VI(M)/SQRT(1.-SINA2/VSQU(M))
         DTI= TI/(VSQU(M)-SINA2)
         TT= TT+TI
         DTDA= DTDA+DTI
      enddo
      DTDA= DTDA*SINA*COS(A)
      DTDD=  DTDA/DDDAS
      DTDH= -(1.-V(JL)*SINA*DTDD)/(V(JL)*SQRT(1.-SINA2))
      ANIN= SINA
      GOTO 260
C------- IF NOT POSSIBLE, SET P + S WEIGHTS TO ZERO
   27 continue
      write(16,*)'WARNING:'
      write(16,*)'Reflected wave did not converge within 50 iterations!'
      write(16,*)'---> trying another ray-type ...'
      write(6,*)'Reflected wave did not converge within 50 iterations!'
      write(6,*)'sina=',sina,'  delta=',delta
      write(6,*)'da=',da,'  da2=',da2,'  test4=',test4
      write(6,*)'---> trying another ray-type ...'
      write(6,*)
      ierr=50
      return
cc      stop'STOP in subr. REFLECT >>> nit > 50'
c
  260 continue
c---- anin is emerging angle at reflector; find now takeoff-angle at source!
      ain=anin
      if(mll.gt.jl)then
         do i=mll-1,jl,-1
            ain=ain * v(i)/v(i+1)
         enddo
      endif
      trefl=tt
      DTDDrefl=dtdd
      dtdhrefl=dtdh
cc     wain=57.29578*asin(ain)
cc      write(6,*)'takeoff angle = ',wain
      return
c
      end ! of subr. reflect
