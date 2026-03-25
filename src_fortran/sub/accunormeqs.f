c
c
c
c
c
c
      subroutine ACCUNORMEQS(rowofa,nvar,res,w,g,rhs)
c
c     Urs Kradolfer, 22. 4. 1987
c
c     accumulate normal equations, knowing only one row of A at a time;
c     subr. ACCUNORMEQS must be called nvar times (once for
c     each observation).
c     w is the (normalized ! ) weight of the observation.
c
c        A * x = res          INPUT: one row of A and one element of res
c
c     At*A * x = At*res
c
c       G  * x =  RHS         OUTPUT: G (nvar x nvar) and RHS (nvar)
c
      implicit none
      real g(*),rowofa(*),rhs(*),res,w
      integer nvar
      call OUTER(g,rowofa,nvar,w)
      call RSIDE(rhs,rowofa,nvar,res,w)
c
      return
      end ! of subr. accunormeqs
