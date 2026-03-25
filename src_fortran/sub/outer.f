c
c
c
c
c
c
      subroutine OUTER(g,s,n,w)
c
c   this routine does an outer product of a vector
c   s with itself.  This is used in least squares
c   to accumulate (At)A knowing only one
c   row of A at a time.  (At)A is a symmetric
c   matrix which is n by n.  This matrix is stored
c   in the vector g in symmetric storage mode.  Thus,
c   given the ith row and jth column of the matrix, then
c            nsym=((i*(i-1))/2) + j
c   gives the index for the element nsym of g
c   for i.ge.j.  For j less than i, the
c   element ij is identical to the element ji.
c
c   input:
c    s(i) one row of the matrix (A) .
c    n - the length of s(i)
c    w - a weight for the least squares.
c
      implicit none
      integer n
      real s(n),g((n*(n-1)/2)+n),w
      integer i,j,nsym
      real a,b
c
      do 1 i=1,n
      if(s(i).eq.0) goto 1
      a=s(i)*w
c
      do 2 j=1,i
      if(s(j).eq.0) goto 2
      b=s(j)*a
c
      nsym=((i*(i-1))/2) + j
      g(nsym)=g(nsym) + b
  2   continue
c
  1   continue
      return
      end ! of subr. outer
