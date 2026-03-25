c
c
c
c
c
c
      subroutine ALEsubr(SV,m,ale)
c
c     Urs Kradolfer, June 1st, 1987
c
c     Input : Array SV, containing the M singular values
c             m , dimension of SV
c
c     Output: ALE-value
c
c             ALE = - 1./(m-izero) * sum(i=1...m) log10( SV(i)/SV(1) )
c                   + izero*10.
c
c             where izero = # zero-singular-values
c
      implicit none
c
      integer m
      real SV(m)
      real ale,alesum
      integer i,izero
c
      if(m.eq.0)stop'ALEsubr>>> m was zero !!!'
c
      ale=0.0
      alesum=0.0
      izero=0
      do i=1,m
         ale=sv(i)/sv(1)
         if(ale.le.0.0)then
            izero=izero+1
         else
            alesum=alesum+log10(ale)
         endif
      enddo
      ale=-(alesum/float(m-izero)) + 10*izero
c
      return
      end ! of subr. alesubr
