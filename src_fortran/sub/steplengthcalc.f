c
c
c
c
c
c
      subroutine STEPLENGTHCALC
c
c     compute the step-length STEPLEN actually applied;
c     STEPLEN is the euclidean length of the model-vector
c
      implicit none
      include '../inc/vel_com.inc'
      integer i
c
cu      real csum(4,2)
c
      steplen=0
cu      bsum=0
cu      k=0
cu      do j=1,4
cu         csum(j,1)=0.0
cu      enddo
      do i=1,nvar
cu        k=k+1
cu        if(i.le.4*neqs) csum(k,1)=csum(k,1)+b(i) ! CSUM is never used...!!??!!
cu        if(nitt.eq.1) csum(k,2)=csum(k,1)
cu        if(k.eq.4) k=0
cu        if(i.gt.4*neqs+nshot) bsum=bsum+b(i)  ! BSUM is never used... !!??!!
         steplen=steplen+b(i)*b(i)
      enddo
      steplen=sqrt(steplen)
      if(ibackups.gt.0) steplen=-steplen
      if(.not.single_turbo)then
         write(16,*)
         write(16,53) steplen
 53      format (' (Applied) Step length = ', f7.3/)
      endif
c
      return
      end ! of subr. steplengthcalc
