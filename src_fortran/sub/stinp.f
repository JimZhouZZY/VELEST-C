c
c
c
c
c
c
      subroutine STINP(itime,stn,sfreq,isper,iscon,isdmp,isamp,scor,ier)
c tabelle fuer stationen im input file
c itime= eventtime in minuten
c
c      ifilt(l)= sfreq
c      iseis(l)= isper
c      sconst(l)=float(iscon)/10.
c      sdampf(l)= float(isdmp)/10.
c      voltgain(l)= isamp
c      cormag(l)= float(iscor)/1000.
c
      implicit none
c      
      integer isper,iscon,isdmp,isamp,ier
      real scor
      integer ifirstcall,ilin,ierr,isyr1,ismo1,isdy1,ishr1,ismin1
      integer isyr2,ismo2,isdy2,ishr2,ismin2,iscor
      character*2 sfreq
chrm      character*4 snam,stn
      character*4 stn
cNEW!!:
      character*6 snam6,stn6
      integer itime,istime1,istime2,juliam
c
      character*80 stlin(600)    ! file STLIST.DAT contains currently 600 lines
      integer maxstlin
      save stlin
      save maxstlin, ifirstcall
c
c     for the first time in this subr. --> read seismo-file into STLIN :
c
      if(ifirstcall.ne.10000001)then
         ifirstcall=10000001
         maxstlin=1
1000     read(10,'(a)',end=2000) stlin(maxstlin)
         maxstlin=maxstlin+1
         if(maxstlin.eq.600) stop'STINP>>> array stlin too small!!'
         goto 1000
2000     maxstlin=maxstlin-1
      endif
c
      stn6=stn
      ilin=0
      ierr=0
1     ilin=ilin+1
      if(ilin.gt.maxstlin) goto 40
      read(stlin(ilin),'(a6)') snam6   !!!!NEW
      if(snam6(1:4).ne.stn6(1:4)) goto 1
20    ilin=ilin+1
      if(ilin.gt.maxstlin) goto 40
      read(stlin(ilin),11)isyr1,ismo1,isdy1,ishr1,ismin1,
     &                    isyr2,ismo2,isdy2,ishr2,ismin2,
     &                    sfreq,isper,iscon,isdmp,isamp,scor
11    format(4x,2(i4,4i2,1x),a2,i2,2i3,i6,f6.2)
      if(isyr1.ne.0) then
         istime1=juliam(isyr1,ismo1,isdy1,ishr1,ismin1)
         istime2=juliam(isyr2,ismo2,isdy2,ishr2,ismin2)
         if(itime.lt.istime1.or.itime.gt.istime2) then
            goto 20
         else
30          ilin=ilin+1
            if(ilin.gt.maxstlin) goto 40
            read(stlin(ilin),11) isyr1
            if(isyr1.ne.0) goto 30
         endif
      endif
      iscor= scor*1000
      return
c
c     if station not found in stationfile with seismometer-parameters:
c
40    continue
      ierr= -6
      return
      end ! of subr. stinp
