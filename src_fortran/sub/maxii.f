c
c
c
c
c
c
      subroutine MAXII(n,nx,imax,jndex)
c determine maximum-value of an integer-array
      implicit none
      integer jndex,i,n,imax, nx(n)
      jndex=1
      do 1 i=1,n
   1  if(nx(jndex).le.nx(i)) jndex=i
      imax=nx(jndex)
      return
      end ! of subr. maxni
