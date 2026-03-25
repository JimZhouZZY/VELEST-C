c
c
c
c
c
c
      subroutine LUDECP (a,ul,n,d1,d2,ier)
c
c     function:        - cholesky decomposition of a matrix A -
c     ---------
c     A   nxn pos. def. symm. matrix, stored in symm. storage mode
c     UL  nxn lower triangular matrix,  with UL * ULt = A
c     !!! diagonal elements of UL are stored in reziprocal form !!!
c     det(A) = D1*(2.**D2)      d1, d2 is output of LUDECP
c     IER=0  --> A is pos. def.;    everything o.k.
c
      implicit none
      integer n
      integer ier
      real   d1,d2
      real   a(n*(n+1)/2),ul(n*(n+1)/2)
      real   zero,one,four,sixtn,sixth,rn,x
      integer ip,i,iq,ir,j,k,ip1
      data   zero,one,four,sixtn,sixth/0.0,1.,4.,16.,.0625/
c
      d1=one
      d2=zero
      rn=one/(n*sixtn)
      ip=1
      ier=0
      do 45 i=1,n
         iq=ip
         ir=1
         do 40 j=1,i
            x=a(ip)
            if (j .eq. 1) go to 10
            do 5  k=iq,ip1
               x=x-ul(k)*ul(ir)
               ir=ir+1
    5       continue
   10       if (i.ne.j) go to 30
            d1=d1*x
            if (a(ip)+x*rn .le. a(ip)) go to 50
   15       if (abs(d1) .le. one) go to 20
            d1=d1 * sixth
            d2=d2 + four
            go to 15
   20       if (abs(d1) .ge. sixth) go to 25
            d1=d1 * sixtn
            d2=d2 - four
            go to 20
   25       ul(ip)=one/sqrt(x)
            go to 35
   30       ul(ip)=x * ul(ir)
   35       ip1=ip
            ip=ip+1
            ir=ir+1
   40    continue
   45 continue
      go to 9005
   50 ier=129
 9000 continue
 9005 return
      end ! of subr. ludecp
