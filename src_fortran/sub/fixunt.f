c
c
c
c
c
c
      subroutine FIXUNT(b,neqs,nshot,nl,ksta,scale,
     &                  vdamp,itotmodels,inltot,nlfirst)
c
c---- restore solution to proper units
c
      implicit none
      integer neqs,nshot,nl,ksta
      real b(4*neqs+nshot+nl+ksta),scale(7)
      integer nlfirst,itotmodels,inltot
      real vdamp(itotmodels,inltot)
      integer i,j,k,l,m
      i=0
      if(neqs.le.0) goto 10
      do 1 j=1,neqs
      do 1 k=1,4
      i=i+1
    1 b(i)=b(i)*scale(k)
10    continue
      if(nshot.le.0) goto 3
      do 2 j=1,nshot
      i=i+1
    2 b(i)=b(i)*scale(7)
3     if(scale(6).eq.0.) return
      l=0
      m=1
      do 5 j=1,nl
      i=i+1
      l=l+1
      if(l.gt.nlfirst) then
         m=m+1
	 l=1
      endif
chrm      write(6,*)m,'  ',l
5     b(i)=b(i)*scale(6)/vdamp(m,l)
      if(ksta.eq.0) return
      do 4 j=1,ksta
      i=i+1
    4 b(i)=b(i)*scale(5)
      return
      end ! of subr. fixunt
