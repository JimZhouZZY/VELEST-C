c
c
c
c
c
c
      subroutine FINALSTARESI
c
c     computes final station residuals for all stations used in this run.
c
      implicit none
      include '../inc/vel_com.inc'
c
      integer i,j,m,k,iz,iwarn
      real stcor1
      real del(ist)
chrm      real phz(3)
      character*1 phz(3)
      real aa(ist),bb(ist),dd(ist),ee(ist)
      integer icc(ist), iccs(ist)
      real aas(ist),bbs(ist),dds(ist),ees(ist)
      character*1 cns,cew
      character cline*80, sta*4
      character*1 phzz(ist)
      data phz(1),phz(2),phz(3)/'S','P','m'/
c
      if(isingle.eq.0)then
         if(iturbo.eq.0)then
            write(16,*)
            write(16,*)
            write(16,*)
            write(16,*)
         else
            write(16,*)'TURBO-option is set; residuals are NOT printed '
     &                 //'for each event!'
         endif
      endif
      write(16,*)
      if(isingle.eq.0)then
         if(iturbo.eq.0)then
            write(16,1)
 1          format(1h1)
         else
            goto 8001    ! do NOT output residuals for each event if TURBO set.
         endif
      endif
c
c     output residuals for each event seperately:
c
      do 2 i=1,legs
         write(16,*)
         if(imin(i).lt.0)then  !  U.K.  3.Feb.87
            imin(i)=imin(i)+60
            ihr(i)=ihr(i)-1
         endif
         if(isingle.eq.0)then
          write(16,'(1x,''Station residuals for event='',i4,
     &               3x,3i2.2,1x,2i2.2,1x,f5.2)')
     &          i,iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i)
         else
          write(16,'(1x,''Station residuals for event='',i4,
     &               3x,3i2.2,1x,2i2.2,1x,f5.2)')
     &          isingle,iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i)
         endif
c        no more active
c3       format(1h ,'Station residuals for event=',i4,
c    &               3x,3i2.2,1x,2i2.2,1x,f5.2)
         write(16,51)
51       format(1x,2('sta ph wt  res   ttime delta',5x))
         do 52 j=1,knobs(i)
            del(j)=sqrt((e(2,i)-d(j,1,i))**2 + (e(3,i)-d(j,2,i))**2)
            phzz(j)='p'
            if(sphase(j,i).eq.1.0) phzz(j)='s'
            if(sphase(j,i).eq.-1.0) phzz(j)='m'
52       continue
         write(16,53) (smn(j,i),phzz(j),kpwt(j,i),res(j,i),tctime(j,i),
     &                 del(j),j=1,knobs(i))
53       format(2(1x,a4,1x,a1,1x,i2,f7.3,2f6.2,4x))
2     continue
      write(16,*)
      if(isingle.ne.0) goto 8001
      write(16,*)
      write(16,*)
      write(16,*)
      write(16,*)
      write(16,*)
      write(16,*)
 8001 continue
c
c     output station statistics (nobs, avres, ...) :
c
      if(isingle.eq.0) then
        write(16,7777) nsp
7     format(/,1x,'sta phase nobs avres  avwres    std    wsum    ',
     &         'delay',/)
        write(16,7)
7777  format(//,2x,' station statistics, remember nsp was set to:',i2,/)
      endif
48    do m=1,nsta
         aa(m)=0.0
         bb(m)=0.0
         icc(m)=0
         dd(m)=0.0
         ee(m)=0.0
         aas(m)=0.0
         bbs(m)=0.0
         iccs(m)=0
         dds(m)=0.0
         ees(m)=0.0
      enddo
c loop 9: collect residual data from all events
      do 9 m=1,nsta
         sta=stn(m)
         do 10 i=1,legs
            k=knobs(i)
            do 11 j=1,k
               if(sta.eq.smn(j,i).and.sphase(j,i).eq.0.) goto 12    ! P
               if(sta.eq.smn(j,i).and.sphase(j,i).eq.-1.0) goto 12  ! refl. P
11          continue
            goto 110
12          continue
            aa(m)=aa(m)+res(j,i)*w(j,i)
cek         res: *w(j,i)
            bb(m)=bb(m)+res(j,i)*res(j,i)*w(j,i)*w(j,i)
            ee(m)=ee(m)+res(j,i)
            dd(m)=dd(m)+w(j,i)
            icc(m)=icc(m)+1
110         continue
            if(nsp.eq.1) goto 10
c
chrm        Doing S and S-P phases
c
            do 101 j=1,k
               if(sta.eq.smn(j,i).and.sphase(j,i).eq.1.) goto 102 ! s-phase
               if(sta.eq.smn(j,i).and.sphase(j,i).eq.2.) goto 102 ! s-p phase
101         continue
            goto 10
102         continue
            aas(m)=aas(m)+res(j,i)*w(j,i)
c           res: *w(j,i)
            bbs(m)=bbs(m)+res(j,i)*res(j,i)*w(j,i)*w(j,i)
            ees(m)=ees(m)+res(j,i)
            dds(m)=dds(m)+w(j,i)
            iccs(m)=iccs(m)+1
10       continue
9     continue
c end loop 9
c
c     Station file output:
c
      if(istaout.gt.0)then
         open(12,file=stafile,status='unknown')
         write(12,1201) fm
 1201    format(a80)
      endif
c
c begin loop 13: print output
      do 13 m=1,nsta
c
c                              write on file 12
         if(istaout.gt.0)then
            iz=ielev(m)
            if(xla(m).lt.0.0)then
               cns='S'
               xla(m)=-xla(m)
            else
               cns='N'
            endif
            if(xlo(m).lt.0.0)then
               cew='E'
               xlo(m)=-xlo(m)
            else
               cew='W'
            endif
            write(cline,fm) stn(m),xla(m),cns,xlo(m),cew,
     &                   iz,model(m),map1(m),
     &                   ptcor(m),stcor(m)
            if(cns.eq.'S') xla(m)=-xla(m)
            if(cew.eq.'E') xlo(m)=-xlo(m)
            if(m.eq.1)then
               cline(47:)='       lon,z,model,icc,ptcor,stcor'
            endif
            write(12,'(a)') cline
         endif
c                            end writing on file 12
c
         if(dd(m).eq.0.or.icc(m).lt.2) goto 1013
         bb(m)=(bb(m)-aa(m)**2/dd(m))*icc(m)/(dd(m)*(icc(m)-1))
         iwarn=0
         if(bb(m).lt.-0.1)then
           iwarn=1
           write(16,'(5x,''WARNING: Station = '',a4,
     &          5x,''!!! Variance  bb('',i3,'') = '',f7.3,'' < 0 !!'')')
     &          stn(m),m,bb(m)
cek           write(6,'(5x,''WARNING: Station = '',a4,
cek     &         5x,''!!! Variance  bb('',i3,'') = '',f7.3,'' < 0 !!'')') 
cek     &         stn(m),m,bb(m)
cek           write(6,*)
            bb(m)=0.0
         else
            bb(m)=abs(bb(m))
         endif
         if(iwarn.eq.0) then
            bb(m)=sqrt(bb(m))
         else
         endif
         aa(m)=aa(m)/dd(m)
         ee(m)=ee(m)/icc(m)
c
cek  phz(2)=  P-phases
c
cek         if(iwarn.eq.0)then
         write(16,'(1x,a4,3x,a1,i4,3f8.4,1x,2f8.4)') 
     &             stn(m),phz(2),icc(m),ee(m),aa(m),bb(m),dd(m),ptcor(m)
cek         else
cek            write(16,'(1x,a4,3x,a1,i4,2f8.4,''   ?.?.?''1x,2f8.4,
cek     &              '' <--- !!'')')
cek     &             stn(m),phz(2),icc(m),ee(m),aa(m),dd(m),ptcor(m)
cek            iwarn=0
cek         endif
1013     if(nsp.eq.1) goto 39
cek         write(16,'(1x,a4,3x,a1,i4,25x,f8.4)') 
cek     &                         stn(m),phz(1),iccs(m),dds(m)
c  print output for S-wave data  (nsp=2 or 3):
         if(dds(m).eq.0.or.iccs(m).lt.2) goto 13
         bbs(m)=(bbs(m)-aas(m)**2/dds(m))*iccs(m)/(dds(m)*(iccs(m)-1))  
         iwarn=0
         if(bbs(m).lt.-0.1)then
           iwarn=1
           write(16,'(5x,''WARNING: Station = '',a4,
     &        5x,''!!! Variance  bbs('',i3,'') = '',f7.3,'' < 0 !!'')')
     &        stn(m),m,bbs(m)
cek           write(6,'(5x,''WARNING: Station = '',a4,
cek     &         5x,''!!! Variance  bb('',i3,'') = '',f7.3,'' < 0 !!'')') 
cek     &         stn(m),m,bb(m)
cek           write(6,*)
         else
            bbs(m)=abs(bbs(m))
         endif
         if(iwarn.eq.0) then
            bbs(m)=sqrt(bbs(m))
         else
            bbs(m)=0.0
         endif
         aas(m)=aas(m)/dds(m)
         ees(m)=ees(m)/iccs(m)
c
c   (nsp=3):
         if(nsp.eq.3) then
            stcor1=ptcor(m)*vpvs
c   (nsp=2):
         else
            stcor1=stcor(m)
         endif
c
cek  phz(1)=  S-phases
c
cek         if(iwarn.eq.0)then
         write(16,'(1x,a4,3x,a1,i4,3f8.4,1x,2f8.4)') 
     &             stn(m),phz(1),iccs(m),ees(m),aas(m),bbs(m),dds(m),
     &             stcor1
cek         else
cek            write(16,'(1x,a4,3x,a1,i4,2f8.4,''   ?.?.?''1x,2f8.4,
cek     &              '' <--- !!'')')
cek     &             stn(m),phz(1),iccs(m),ees(m),aas(m),dds(m),stcor1
cek            iwarn=0
cek         endif
c
 14      format(1x,a4,3x,a1,i4,3f8.4,1x,2f8.4)
 39      continue
 13   continue
c
c end loop print one station
c
      write(16,*)
      write(16,*)
      if(istaout.gt.0)then
         write(12,*)
         close(12)
      endif
c
chrm  Output of the final model
c
      if(istaout.eq.2)then
         open(12,file='velout.mod',status='unknown')
         write(12,*)'Output model:'
c	 write(12,*)(nplay(m),m=1,nmod)
	 do m=1,nmod
            write(12,*)nplay(m)
	    do i=1,nplay(m)
	       write(12,'(f5.2,5x,f7.2,2x,f7.3,4x)')
     &              vp(m,i),hp(m,i),vdamp(m,i)
	    enddo
	 enddo
	 close(12)
      endif 
      return
      end ! of subr. finalstaresi
