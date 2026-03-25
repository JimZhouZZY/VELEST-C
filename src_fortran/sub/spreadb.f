c
c
c
c
c
c
      subroutine SPREADb(R,m,spread) ! Urs Kradolfer
c     BACKUS-GILBERT SPREAD FUNCTION
c     measures the goodness of the resolution spread, based on the
c     L2 norm of the WEIGHTED difference between the resolution matrix and the
c     identity matrix [ Backus-Gilbert spread function ]
c     R = I   <==>  spread(R) = 0
c
      implicit none
      integer m
      real r(m,m)
      real spread
      integer j,i
c
      spread=0.0
      do i=1,m
         do j=1,m
            spread=spread+(i-j)*(i-j)*r(i,j)*r(i,j)
         enddo
      enddo
      return
      end ! of subr. spreadb
