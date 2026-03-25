c
c
c
c
c
c
      subroutine MATRTRAN(A,n,m,AT)
c     Urs Kradolfer, 28.3.1987
      implicit none
c     input:  (nxm)-matrix A
c     output: (mxn)-matrix AT ( = A transpose )
      integer n,m ,
     &        i,j
      real A(n,m),AT(m,n)
      do i=1,n
         do j=1,m
            AT(j,i)=A(i,j)
         enddo
      enddo
      return
      end ! of subr. matrtran
