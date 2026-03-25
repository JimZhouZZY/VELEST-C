c
c
c
c
c
c
      subroutine luelmp (a,b,n,x)
c
c     function:     - elimination part of the solution of Ax=b
c     ---------
c       A  (nxn) lower triangular matrix (outout of LUDECP)
c       !!! diagonal elements of A are stored in reziprocal form !!!
c       b  n-vector
c       x  n-vector (solution of equation Ax=b)
c
      implicit none
      integer n
      real a(n*(n+1)/2),b(n),x(n)
      real zero,t
      integer ip,iw,i,im1,k,n1,ii,is,iq,kk
      data zero/0./
c                              solution of ly = b
      ip=1
      iw=0
      do 15 i=1,n
         t=b(i)
         im1=i-1
         if (iw .eq. 0) go to 9
         ip=ip+iw-1
         do 5 k=iw,im1
            t=t-a(ip)*x(k)
            ip=ip+1
    5    continue
         go to 10
    9    if (t .ne. zero) iw=i
         ip=ip+im1
   10    x(i)=t*a(ip)
         ip=ip+1
   15 continue
c                                  solution of ux = y
      n1=n+1
      do 30 i=1,n
         ii=n1-i
         ip=ip-1
         is=ip
         iq=ii+1
         t=x(ii)
         if (n.lt.iq) go to 25
         kk=n
         do 20 k=iq,n
            t=t-a(is)*x(kk)
            kk=kk-1
            is=is-kk
   20    continue
   25    x(ii)=t*a(is)
   30 continue
      return
      end ! of subr. luelmp
