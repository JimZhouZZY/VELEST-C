c
c
c
c
c
c
      subroutine STEPLENGTHDAMP(damp)
c
c     set step length damping to model-vector in case of LARGE VELOCITY CHANGE
c     (  | deltavelocity |   >   MAXveladjALLOWED  )
c
      implicit none
      real damp
      include '../inc/vel_com.inc'
      integer i,j11,j22,jjj
      real btemp
c
      damp=1.0
      if(icount.eq.1) goto 900
      if(scale(6).eq.0.0) goto 900
      do i=1,nmod
         j11=4*neqs+nshot+laysum(i)
         j22=j11+nplay(i)-1
         do 121 jjj=j11,j22
            if(veladj.gt.abs(b(jjj))) goto 121
            btemp=veladj/abs(b(jjj)) ! ABS[b(jjj)] > veladj --> damping enlarged
            if(btemp.lt.damp) damp=btemp ! eff. damping=smallest of all dampings
 121     continue
      enddo
c
c     apply step length damping just calculated to ALL unknowns !!! :
c
      do jjj=1,nvar
         b(jjj)=b(jjj)*damp
      enddo
900   continue
c
      return
      end ! of subr. steplengthdamp
