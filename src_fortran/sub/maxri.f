c
c
c
c
c
c
      subroutine MAXRI(n,x,xmax,jndex)
c determine maximum-value of a real-array
      implicit none
      integer i,n,jndex
      real x(n), xmax
      jndex=1
      do 1 i=1,n
   1  if(x(jndex).le.x(i)) jndex=i
      xmax=x(jndex)
      return
      end ! of subr. maxri
