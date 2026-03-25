c
c
c
c
c
c
      subroutine RESOLCOVAR(davari)
c
c     calculate resolution and covariance matrices.
c
      implicit none
      real davari
      include '../inc/vel_com.inc'
c
      integer mef,nef,n1,n2,n3,k,i,l,j
      real scale1,spread1,spread2,size,avresol
      real Rdiag(inva)
c
      character*80 pcard
c
c---- resolution and covariance calculations
      if(.not.single_turbo)then
         write(16,'(/////)')
         write(16,40)
   40    format(' Resolution and covariance calculations:')
         if(isingle.ne.0)then
cc       write(16,*)'CHOLESKY-decomposition:'
         write(16,*)'    RESOLUTION-matrix                   ',
     &              '    COVARIANCE-matrix'
         write(16,*)
         endif
      endif
   25 mef=4*neqs+nshot-4
      nef=4+nltot
      if(isingle.ne.0) NEF=4     ! only first 4 rows and first 4 standard-devs.
      n1=4*neqs
      n2=n1+nshot
      n3=n2+nltot
      if(iresolcalc.eq.2.and.isingle.eq.1) 
     &write(80,*)'Covariance matrice'
      do 15 k=1,nef
      do 16 i=1,nvar
   16 rht(i)=gcopy(k,i)
      l=k
      if(k.gt.4) l=k+mef
      call LUELMP(g,rht,nvar,rht)
      if(isingle.eq.0)then
         if(.not.single_turbo)then
            write(16,17) l
 17         format(1x,'resolution row',i3)
            write(16,18)  (rht(j),j=1,nvar)
         endif
         Rdiag(k)=rht(l)
      else
         write(pcard,18) (rht(j),j=1,nvar)
         do j=1,4
            Rc(l,j)=rht(j)     !  resolution-matrix
         enddo
      endif
 18   format(8f10.4)
      call LUELMP(g,rht,nvar,rht)
c---- put covariance into proper unit
      j=mod(l,4)
      if(j.eq.0) j=4
      scale1=scale(j)
      if(l.gt.n1) scale1=scale(1)
      if(l.gt.n2) scale1=scale(6)
      if(l.gt.n3) scale1=scale(5)
      do 19 j=1,nvar
ccc   19 rht(j)=rht(j)*davari*scale1
   19 rht(j)=rht(j)*scale1  ! <== UNIT covariance matrix !!!
       call FIXUNT(rht,neqs,nshot,nltot,ksta,scale,
     &             vdamp,itotmodels,inltot,nplay(1))
      if(isingle.eq.0)then
         if(.not.single_turbo)then
            write(16,20) l
 20         format(1x,'covariance row',i3)
            write(16,18)  (rht(j),j=1,nvar)
         endif
      else
         write(pcard(41:80),18) (rht(j),j=1,nvar)
         if(.not.single_turbo) write(16,*) pcard
         do j=1,4
            COVc(l,j)=rht(j)
         enddo
      endif
c      s(k)=sqrt(abs(rht(l)))  ! compiler produces warning with range checking
      s(k)=abs(rht(l))
      s(k)=sqrt(s(k))
 15   continue
c
      if(isingle.ne.0)then
         if(.not.single_turbo) write(16,*)
         call SPREADd(Rc,4,spread1)
         spread=spread1
         if(ifixsolution.eq.1)then
            do j=1,3
               do k=1,3
                  R3(j,k)=Rc(j,k)
               enddo
            enddo
            call SPREADd(R3,3,spread)
         endif
         if(ifixsolution.eq.9)then
            call SPREADd(Rc(1,1),1,spread)
         endif
         call SPREADb(Rc,4,spread2)
         size=0.0
         do j=1,4
            size=size+COVc(j,j)
         enddo
         if(.not.single_turbo)then
            write(16,'(1x,''D-Spread(R) = '',f6.3,
     &                    ''   B-G Spread(R) = '',f6.3,
     &                 ''       Size (C) = '',f6.3
     &             )') spread,spread2, size
         endif
      endif
c
      if(.not.single_turbo)then
         write(16,*)
         write(16,21)
 21      format(' standard deviation of selected model parameters:')
      endif
      if(isingle.ne.0)then
        avresol=0.0
        do j=1,nef
           avresol=avresol+Rc(j,j)
        enddo
        avresol=avresol/float(nef)
        i f (.not.single_turbo) t h e n
        write(16,*)
        write(16,*)'   OT (sec)   X (km)    Y (km)    Z (km) '
        e n d i f
        write(6,'(23x,''  OT (sec)   X (km)    Y (km)    Z (km) '')')
c        write(6,'(1x,''CHOLESKY  (othet='',f5.3,'') :'')') othet
        write(6,'('' Sigma (CHD):         '',4f10.4)') (s(j),j=1,nef)
        write(6,'('' Resolution (CHD):    '',4f10.4,x,
     &            ''D-spread ='',f6.3)') (Rc(j,j),j=1,nef) , spread
        write(6,'('' Data Variance      = '',f10.4)') davari
ccc        call SPREADd(Rs,4,spread3)
ccc        do j=1,4
ccc           COVs(j)=SQRT(COVs(j))  !  <-- standard deviation
ccc        enddo
        write(6,'('' Singular values:     '',4f10.4,5x,''ALE ='',f7.3)')
     &               (SV(j),j=1,nef) , ale(1)
ccc        write(6,'('' Sigma (SVD):         '',4f10.4)') (COVs(j),j=1,nef)
ccc        write(6,'('' Resolution (SVD):    '',4f10.4,3x,
ccc     &            ''D-spread ='',f6.3)') (Rs(j,j),j=1,nef) , spread3
      endif
      if(.not.single_turbo) write(16,22) (s(j),j=1,nef)
 22   format(8f10.4)
      if(.not.single_turbo)then
         write(16,*)
         write(16,*)'Rdiag of selected model parameters:'
         write(16,22) (Rdiag(j),j=1,nef)
         write(16,*)
         write(16,*)
      endif
c
      return
      end ! of subr. resolcovar
