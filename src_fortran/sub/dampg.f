c
c
c
c
c
c
      subroutine DAMPG
c
c     apply damping to the diagonal elements of the symmetric matrix G
c
      implicit none
      include '../inc/vel_com.inc'
      integer j,k
c
c     apply damping to diagonal elements of G-matrix :
c
      j=0
      do k=1,nvar          !   k = 1 2 3  4  5 ...
         j=j+k             !   j = 1 3 6 10 15 ...   <== diagonal elements of G
         g(j)=g(j)+othet
      enddo
c
      return
      end ! of subr. dampg
