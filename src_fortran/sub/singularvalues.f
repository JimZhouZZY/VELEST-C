c
c
c
c
c
c
      subroutine SINGULARVALUES(i)
c
c     compute the singular values (in this case: eigenvalues) of the
c     symmetric matrix G.
c
      implicit none
      integer i
      include '../inc/vel_com.inc'
      integer k,ii,jj,i1,j
c
      real A(4,4), Xsol(4)
c
      do k=1,nvar
         ii=0
         jj=0
         do i1=1,nvar
            do j=1,i1
               jj=jj+1
               if(.NOT.(i1.ne.k.and.j.ne.k))then
                  ii=ii+1
                  A(k,ii)=g(jj)
               endif
            enddo
         enddo
      enddo
c
      call SVDSOLUK(A,RHt,nvar,-1.,Xsol,SV,ale(i),COVs,Rs)
c
      if(ifixsolution.eq.1) call ALESUBR(SV,3,ale(i))
      if(ifixsolution.eq.9) call ALESUBR(SV,1,ale(i))
c
      if(.not.single_turbo)then
         write(16,*)'Singular values; iteration #',nitt
         write(16,'(1x,4(2x,f10.6))') (sv(jj),jj=1,nvar)
         write(16,*)'ALE = ',ale(i)
      endif
c
ccc      write(16,*)'SVD-solution vector:'
ccc      write(16,'(1x,4(2x,f16.6))') (Xsol(jj),jj=1,nvar)
c
      return
      end ! of subr. singularvalues
