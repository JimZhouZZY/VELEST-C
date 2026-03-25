c
c
c
c
c
c
      subroutine ACTUALSTATIONS
c
cek added separate report of P and S phases 20.8.98
c
      implicit none
      include '../inc/vel_com.inc'
      integer i,nofreadings,k,nobsp(ist),nobss(ist),nofreadp,nofreads
c
      do i=1,nsta
         nactualsta(i)=0
         nobsp(i)=0
         nobss(i)=0
      enddo
      nstaeff=0
      nofreadings=0
      do i=1,legs  ! = neqs+nshot
         do k=1,knobs(i)
            nactualsta( istm(k,i) ) = nactualsta( istm(k,i) ) + 1
            if(nsp.eq.2) then
             if(sphase(k,i).eq.0.0) nobsp(istm(k,i))=nobsp(istm(k,i))+1
             if(sphase(k,i).eq.1.0) nobss(istm(k,i))=nobss(istm(k,i))+1
            endif
            nofreadings=nofreadings+1
         enddo
      enddo
      if(.not.single_turbo) write(16,*)
      do i=1,nsta
         if(nactualsta(i).gt.0)then
            if(.not.single_turbo)then
               if(nsp.eq.2) then
               write(16,'('' readings for station '',a4,'' : tot='',
     &                  i4,''  P:'',i4,''  S:'',i4)')
     &                  stn(i),nactualsta(i),nobsp(i),nobss(i)
               else
               write(16,'('' readings for station '',a4,'' :'',i4)')
     &                  stn(i),nactualsta(i)
               endif
            endif
            nstaeff=nstaeff+1
         endif
         nofreadp=nofreadp+nobsp(i)
         nofreads=nofreads+nobss(i)
      enddo
      if(.not.single_turbo)then
         write(16,*)
         write(16,'('' Total number of stations with readings:'',i4)')
     &           nstaeff
         write(16,*)
         write(16,'('' Total number of readings: '',i7)') nofreadings
         write(16,'('' Total number of P readings: '',i7)') nofreadp
         write(16,'('' Total number of S readings: '',i7)') nofreads
         write(16,*)
      endif
c
      RETURN
      end ! of subr. actualstations
