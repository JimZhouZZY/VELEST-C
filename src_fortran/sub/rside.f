c
c
c
c
c
c
      subroutine RSIDE(rht,s,n,res,w)
c
c   this routine accumulates in rht the quantity (At)b
c   for the normal equations.
c
c   input:
c    s(i)-one column of the matrix (At)
c    n - the length of s(i)
c    res - one value of the vector b.
c    w  -  a weight for the weighted least sqs.
c
      implicit none
      integer n
      real rht(n),s(n),w,res
      integer i
c
      do 1 i=1,n
      if(s(i).eq.0) goto 1
      rht(i)=res*w*s(i) + rht(i)
1     continue
      return
      end ! of subr. rside
