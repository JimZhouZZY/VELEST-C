c
c
c
c
c
c
      subroutine SETUPMATRIXG(neq,i)
c
c     puts one row S(i) of the data kernel into the symmetric matrix G in
c     order to accumulate the normal equations.
c
      implicit none
      integer neq,i
      include '../inc/vel_com.inc'
      integer j,nn,nni,mm,is1,kk1,ng,mg,ksta1,ksta2,k,k2
      integer nl,nf,mf,n,k1,is
c
ccc      real WK(8,2),bwk(1)
c
      do 2 j=1,nvar
   2  s(j)=0.
c
c     No trvdrv's for reseted obs:
c
      if(isingle.ne.0.and.w(i,neq).eq.0.0) goto 34
c
c---- calculate hypocenter indices
c***   neq-- event number, i--observation number
      nn=4*neq-2
      nni=nn-1
c***    nni- space for time term; nn thru mm -- 3 places for dtdr's
      if(neq.gt.neqs) nni=3*neqs+neq
      mm=nn+2
      if(neq.gt.neqs) mm=nni
      s(nni)=1.0*scale(1)      
      if (sphase(i,neq).eq.2.0) then ! new: for s-p phases (hrm 30.3.92)
         s(nni) = 0.0      ! set derivative for orig. time to zero
      endif	 
      is1=0
      if(neq.le.neqs) goto 40
      s(nni)=1.0*scale(7)
      kk1=map2(neq-neqs)
      if(nshfix.eq.1.and.kk1.ne.0) s(nni)=0.0
      if(nshcor.eq.0) goto 40
      if(kk1.le.0.or.kk1.gt.ksta) goto 40
      if(scale(5).eq.0) goto 40
c
      if(nsp.eq.3.and.sphase(i,neq).eq.1.) goto 40
      if(nsp.eq.2.and.sphase(i,neq).eq.1.) goto 43
      ng=4*neqs+nshot+nltot+1
      mg=ng-1+ksta
      ksta1=ksta
      if(kk1.gt.ksta1) goto 40
      is1=ng-1+kk1
      goto 44
43    continue
      ng=4*neqs+nshot+nltot+(ksta/2)+1
      mg=4*neqs+nshot+nltot+ksta
      ksta1=(ksta/2)
      ksta2=ksta-(ksta/2)
      if(kk1.gt.ksta2) goto 40
      is1=ng-1+kk1
44    continue
      s(is1)=1.0*scale(5)
40    continue
      k=0
      if(zadj.eq.0.0) dtdr(3)=0.0
cc      if(isingle.ne.0.and.nitt.lt.3) dtdr(3)=0.0
cc      if(isingle.ne.0.and.jgap.gt.250) dtdr(3)=0.0
      if(neq.gt.neqs) goto 3
c
cEK Dez. 1994  next statement again put in effect:
      if(ifx(neq).eq.1) dtdr(2)=0.0   ! no longer in use  U.K. 2.Oct.87
      do 4 j=nn,mm
      k=k+1
   4  s(j)=dtdr(k)*scale(k+1)
   3  continue
c---- calculate velocity indices
      if(scale(6).eq.0.) goto 34
15    k2=iphase(i,neq)
      nl=nplay(k2)
      nf=4*neqs+nshot+laysum(k2)
      mf=nf+nl-1
      k=0
      do 5 n=nf,mf
      k=k+1
      if(veladj.eq.0.0) dtdv(k)=0.0
chrm    5 s(n)=dtdv(k)*scale(6)
    5 s(n)=dtdv(k)*scale(6)/vdamp(k2,k)
c---- calculate indices for station terms
      if(scale(5).eq.0.0) goto 34
      k1=istm(i,neq)
      ksta1=ksta
      if(nsp.eq.2) ksta1=ksta/2
      if(nsp.eq.2.and.sphase(i,neq).eq.1.0) goto 16
      if(nsp.eq.3.and.sphase(i,neq).eq.1.0) goto 34
18    continue
      ng=4*neqs+nshot+nltot+1
      mg=ng-1+ksta
      if(map1(k1).eq.0) goto 34
      if(map1(k1).gt.ksta1) goto 34
      is=ng-1+map1(k1)
      goto 17
16    ng=4*neqs+nshot+nltot+(ksta/2)+1
      ksta2=ksta-(ksta/2)
      mg=ng-1+ksta2
      if(map1(k1).eq.0) goto 34
      if(map1(k1).gt.ksta2) goto 34
      is=ng-1+map1(k1)
17    continue
      s(is)=1.0*scale(5)
34    continue
cc
cc     store one line of G in GG:
cc
c
c      if(isingle.ne.0)then
c        do i1=1,4
c            gg(i,i1)=s(i1)
c        enddo
cc         if(i.eq.knobs(1))then ! matrix GG is fully calculated
cc            call LSVDF(GG,100,knobs(1),4,BWK,1,0,SV,WK,ier)
cc            call ALESUBR(SV,4,ale(1))
ccccccccc  test          ale(1)=-log10( sv(4) )   ! test test
ccccc            write(6,*) sv(1),sv(2),sv(3),sv(4), ale(1)
cc         endif
cc      endif
c
c     now accumulate normal equations (vector s contains the traveltime-derivs.)
c
      call ACCUNORMEQS(s,nvar,res(i,neq),w(i,neq),G,RHT)
c
      return
      end ! of subr. setupmatrixg
