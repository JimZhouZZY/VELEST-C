c
c
c
c
c
c
      subroutine STOREG(k,l)
c
c---- save column k of g in gcopy
c
      implicit none
      integer k,l
      include '../inc/vel_com.inc'
      integer ii,jj,i,j
c
      ii=0
      jj=0
      do 1 i=1,nvar
      do 1 j=1,i
      jj=jj+1
      if(i.ne.k.and.j.ne.k) goto 1
      ii=ii+1
      gcopy(l,ii)=g(jj)
   1  continue
      return
      end ! of subr. storeg
