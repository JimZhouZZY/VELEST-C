c
c
c
c
c
c
      subroutine SVDSOLUK(Ain,Bin,m,eigmin,X,S,ale,COV,R)
c
c     Urs Kradolfer, 24. April 1987
c
c     solve NORMAL EQUATIONS  Gt*G * X = Gt*R
c                               A  * X  =  B    , so A is symmetric !!!
c
c     Input : left-hand-side A(mxm), right-hand-side  B(m),
c             eigmin : >=0  cutoff-value for solution
c                      = -1 only S (singular values) are returned
c
c     Output: solution-vector X(m), singular-(eigen-)values S(m),
c             diagonal elements of covariance-matrix COV(m),
c             resolution-matrix R(mxm)
c             Average Logarithmic Eigenvalue ALE :
c
c             ALE = 1./(m-izero) * sum(i=1...m) log10(S(i))
c
c             where izero = # zero-singular-values
c
c             before, LGCN was used !!!
c             GCN = LOG10 [ average quotient of 1stEigenvalue to ithEigenvalue
c                           / average eigenvalue ]
c
      implicit none
      integer m
      real eigmin,ale,sum
      integer mm,mm2,i,j,ier,izero,nfre
      parameter (mm=4,mm2=8)
      real Ain(mm,mm),Bin(mm),A(4,4),B(mm),
     &     X(mm),S(mm),COV(mm),R(mm,mm)
      real At(mm,mm)
      real WK(mm2,2)  ! must be dimensioned to : 2m,2
c
      if(m.gt.4)stop'SVDSOLUK>>> m.gt.4 ; redimension arrays !!!'
c
c     copy input-matrix/-vector into local ones (so the input doesn't change) :
c
      do i=1,m
         do j=1,m
            a(i,j)=ain(i,j)
         enddo
         b(i)=bin(i)
         s(i)=0.0
         cov(i)=0.0
      enddo
      do i=1,mm2
         do j=1,2
            wk(i,j)=0.0
         enddo
      enddo
c
c     do SVD (singular value decomposition) of matrix A = V * S * Vt   :
c
      call LSVDF(A,m,m,m,B,m,1,S,WK,ier)   ! IMSL-routine
c
c     determine ALE  ( Average Logarithmic Eigenvalue )
c
c      izero=0
c      if(m.eq.1)then
c         GCN=1.0
c      else
c         GCN=0.0
c         sum=s(1)
c         do i=2,m
c            if(s(i).gt.0.0)then
c               GCN=GCN+s(1)/s(i)
c               sum=sum+s(i)
c            else
c               izero=izero+1
c            endif
c         enddo
c      endif
c      sum=sum/m   ! average eigenvalue.
c      GCN=izero*100.+log10((GCN/float((m-1)-izero))/sum)
cc      ale=0.0
      izero=0
cc      do i=1,m
cc            ale=ale+log10((s(i)/s(1)))
ccc            ale=ale+log10((s(i)/s(1)) + 1e-10)
cc      enddo
cc      ale=-ale/float(m)
      call ALEsubr(s,m,ale)
c
      if(eigmin.lt.0) return
c
c     determine degree of freedom (neglect sing. vals. < cutoff-value eigmin) :
c
      nfre=0
      do i=1,m
         x(i)=0.0
         if(s(i).gt.eigmin)then
            nfre=nfre+1
         else
ccc            write(6,*)'SVDSOLUK>>> skipping singular value = ',s(i)
            s(i)=0.0        ! set eigenvalue to zero
            do j=1,m
               a(j,i)=0.0   ! set appr. COLUMN of V (=eigenvector!!!) to zero
            enddo
         endif
      enddo
c
c     calculate solution-vector X = V * inv(S) * Vt*B    :
c                            [  X = A * inv(S) *   B   ]
c
      do i=1,m
         do j=1,nfre
            x(i)=x(i)+a(i,j)*b(j)/s(j)
         enddo
      enddo
c
c     calculate diagonal elements of UNSCALED covariance matrix :
c
c     C = datvar * V * inv(inv(S)) * Vt
c
c     C(i,j)= SUM(k=1...m)  V(i,k)*V(j,k)/(S(k,k)*S(k,k))  <-- for all elements
c     C(i,i)= SUM(k=1...m)  V(i,k)*V(i,k)/(S(k,k)*S(k,k))  <-- diag. elem. only
c
      do i=1,m
         sum=0.0
         do j=1,nfre
ccc            write(6,*)'i=',i,'j=',j,'V(i,j)=',a(i,j)
            sum=sum+a(i,j)*a(i,j)/(s(j)*s(j))
         enddo
         cov(i)=sum
      enddo
c
c     compute resolution-matrix:  R = V * Vt    :
c                               [ R = A * At ]
c
      call MATRTRAN(A,m,m,At)
      call MATRMULT(A,m,m,At,m,m,R,m,m)
c
      return
      end ! of subr. SVDSOLUK
