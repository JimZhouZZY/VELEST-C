c
c
c
c
c
c
      subroutine REFLECT1(nl,v,vsq,thk,jl,tkj,delta,z,mll,trefl,ain,
     &                    ierr,DTDDrefl,DTDHrefl)
c
c     Urs Kradolfer, Nov. 1986
c
      implicit none
c      
      integer nl,ierr
      real tkj,trefl,dtddrefl,dtdhrefl
c     
      integer m,lmax,j1,l,kount
      real del, fac,depth,tklmax,vlmax,ua,ub,uasq,ubsq
      real xa,xb,dela,delb,ubdiv,x,u,usq
      real xtest,tdir
      integer jl   ! hypocenter-layer
      integer mll  ! reflection at bottom of layer MLL
      real thk(nl),v(nl),vsq(nl)
      real div(100)  ! dimension = max. number of layers allowed
      real ain     ! auftauchwinkel in rad
      real z       ! hypocenter-depth
      real delta  ! distance epicenter-receiver
      integer i
c
      ierr=0
      DTDDrefl=0.0
      DTDHrefl=0.0
c
c
      del=0.0   ! U.K. 28. Jan. 1987
c
      if(jl.gt.mll)then   !  hypocenter below reflector !!!
         write(6,*)'hypocenter-layer jl = ',jl
         write(6,*)'reflection at bottom of layer mll = ',mll
         ierr=-1
         write(6,*)'WARNING:'
         write(6,*)'subr. REFLECT1 >>> hypocenter is below reflector!'
         write(16,*)'WARNING:'
         write(16,*)'subr. REFLECT1 >>> hypocenter is below reflector!'
ccc         stop'subr. REFLECT1 >>> hypocenter below reflector!'
      endif
      if(mll.gt.nl)then
         stop'subr. REFLECT1 >>> reflector-nr > number of layers !'
      endif
c
c  determine out of real model (THK) transformed model (DIV)
c
      DO M=1,MLL
C LAYER M IS PASSED ONLY ONCE:
         FAC= 1.
         IF(JL.GT.M) GOTO 505
C LAYER M IS PASSED TWICE
         FAC= 2.
         IF(JL.LT.M) GOTO 505
C LAYER M CONTAINS HYPOCENTER:
         FAC= (2.*THK(M)-tkj)/THK(M)
  505    DIV(M)= FAC*THK(M)
      enddo
      depth=0.0
      do i=1,mll
         depth=depth+div(i)
      enddo
c
c
c  now 'new' model is established; layer-thicknesses are stored in DIV(1...mll)
c
      lmax=mll
      tklmax=div(mll)
      vlmax=v(mll)
      j1=mll-1
      do 23184l=1,j1
      if(.not.(v(l).gt.vlmax))goto 23186
      lmax=l
      tklmax=div(l)
      vlmax=v(l)
23186 continue
23184 continue
23185 continue
C CHANGE BY E.KISSLING MARCH 1984
      IF(tklmax.le.0.05) tklmax=0.05
c
c   find initial bounds on the sine of the takeoff angle
c
      ua=(v(mll)/vlmax)*delta/sqrt(delta**2+depth**2)
      ub=(v(mll)/vlmax)*delta/sqrt(delta**2+tklmax**2)
c
c   calculate horizontal travel distances
c
      uasq=ua**2
      ubsq=ub**2
C CHANGE BY E.KISSLING MARCH 1984
      IF(UBSQ.GE.1.) UBSQ=0.99999
      IF(UASQ.GE.1.) UASQ=0.99999
      xa=div(mll)*ua/sqrt(1.0-uasq)
      if(.not.(lmax.eq.mll))goto 23188
      xb=delta
      goto 23189
23188 continue
      xb=div(mll)*ub/sqrt(1.0-ubsq)
23189 continue
      dela=xa
      delb=xb
      do 23190l=1,j1
      dela=dela+div(l)*ua/sqrt(vsq(mll)/vsq(l)-uasq)
      ubdiv=sqrt(vsq(mll)/vsq(l)-ubsq)
      if(ubdiv.GT.1.e-20) GOTO 1002
      write(16,*)'WARNING:'
      write(16,1000) mll,l,lmax,vsq(mll),vsq(l),ubsq,delta,tklmax,vlmax
 1000 format(/,2x,'subr. reflect1: ',3i4,2f10.4,f15.6,3f10.5,/)
      ubdiv=1.e-20
 1002 continue
      delb=delb+div(l)*ub/sqrt(vsq(mll)/vsq(l)-ubsq)
23190 continue
23191 continue    !  NOT used... !!!
c
c   loop to find the zero of del-delta by teh method of false position
c
      do 23192kount=1,25
      if(.not.((delb-dela).lt.0.02))goto 23194
      x=0.5*(xa+xb)
      u=x/sqrt(x**2+div(mll)**2)
      usq=u**2
      goto 23193
23194 continue
      x=xa+(delta-dela)*(xb-xa)/(delb-dela)
      u=x/sqrt(x**2+div(mll)**2)
      usq=u**2
      del=x
      do 23196l=1,j1
      del=del+div(l)*u/sqrt(vsq(mll)/vsq(l)-usq)
23196 continue
23197 continue
      xtest=del-delta
      if(.not.(abs(xtest).lt.0.02))goto 23198
      goto 23193
23198 continue
      if(.not.(xtest.lt.0.0))goto 23200
      xa=x
      dela=del
      goto 23201
23200 continue
      xb=x
      delb=del
23201 continue
23192 continue
23193 continue
c
c   calculate direct ray travel time
c
c
      if(del.eq.0.0) del=x   ! U.K. 28. Jan. 1987
c
      tdir=(sqrt(x**2+div(mll)**2))/v(mll)
      do 23202l=1,j1
      tdir=tdir+div(l)*v(mll)/(vsq(l)*sqrt(vsq(mll)/vsq(l)-usq))
23202 continue
23203 continue
      tdir=tdir-(u/v(mll))*(del-delta)
c
      trefl=tdir
c---- u is sine of the 'takeoff-angle' in
c---- transformed model = emerging angle at reflector!
c---- Find now takeoff-angle at source!
      ain=u
      if(mll.gt.jl)then
         do i=mll-1,jl,-1
            ain=ain * v(i)/v(i+1)
         enddo
      endif
c
      return
      end ! of subr. reflect1
