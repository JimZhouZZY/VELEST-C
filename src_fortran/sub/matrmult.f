c
c
c
c
c
c
      subroutine MATRMULT(A,m,p,B,p1,n,C,m1,n1)
c     Urs Kradolfer, 28.3.1987
      implicit none
c     input:  (mxp)-matrix A
c             (p1xn)-matrix B
c     output: (m1xn1)-matrix C = A*B
      integer m,p, p1,n, m1,n1 ,
     &        i,j,k
      real A(m,p), B(p1,n), C(m1,n1) ,
     &     s
      if(m.ne.m1)stop'Matrices cannot be multiplied !'
      if(p.ne.p1)stop'Matrices cannot be multiplied !'
      if(n.ne.n1)stop'Matrices cannot be multiplied !'
      do i=1,m
         do j=1,n
            s=0.
            do k=1,p
               s=s+A(i,k)*B(k,j)
            enddo
            C(i,j)=s
         enddo
      enddo
      return
      end ! of subr. matrmult
