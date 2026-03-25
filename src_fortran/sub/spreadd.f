c
c
c
c
c
c
      subroutine SPREADd(R,m,spread) ! Urs Kradolfer
c     DIRICHLET SPREAD FUNCTION
c     measures the goodness of the resolution spread, based on the
c     L2 norm of the difference between the resolution matrix and the
c     identity matrix [ Dirichlet spread function ]
c     R = I   <==>  spread(R) = 0
c
      implicit none
      integer m
      real r(m,m)
      real spread
      integer i,j
c
      spread=0.0
      do i=1,m
         do j=1,m
            if(i.eq.j) spread=spread+(r(i,j)-1)*(r(i,j)-1)
            if(i.ne.j) spread=spread+r(i,j)*r(i,j)
         enddo
      enddo
      return
      end ! of subr. SPREADd
