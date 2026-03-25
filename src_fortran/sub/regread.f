c
c
c
c
c
c
      subroutine REGREAD(regnamfile,regkoordfile)
      implicit none
      integer indfe,lt50,lt25
      real xnrlon
      character*(*) regnamfile,regkoordfile
      character irname*24300
      common /regioncom/ indfe(730),lt50(111),lt25(421),
     &                   irname,xnrlon(14400)
c
      integer trimlen,nstart,nr
      integer nc
      integer i,j,k
      integer iu1, iu2
      character*40 cregion
cu      character*64 inputfile
c
cu      inputfile='/users/kradi/velest/regionsnamen.dat'
      call FREEUNIT(iu1)
cVMS      open(iu1,file=regnamfile,status='OLD',readonly)
      open(iu1,file=regnamfile,status='OLD')
      nstart=0
      do i=1,729
          indfe(i)=nstart
          read(iu1,'(i5,2x,a)') nr,cregion
          if(nr.ne.i) then
              write(6,'(''sequence error at:'',i5,a,i5)') nr,cregion,i
              stop'REGREAD>>> error reading file REGIONSNAMEN.DAT'
          endif
          nc=TRIMLEN(cregion)
          irname(nstart+1:nstart+nc)=cregion(1:nc)
          nstart=nstart+nc
      enddo
      indfe(730)=nstart
      do i=1,420
          lt25(i)=nstart
          read(iu1,'(i5,2x,a)') nr,cregion
          if(i.le.400) then
             k=i+999
          else
             k=1019+(i-401)*20
          endif
          if(k.ne.nr) then
              write(6,'(''sequence error at:'',i5,a,i5)') nr,cregion,i
              stop'REGREAD>>> error reading file REGIONSNAMEN.DAT'
          endif
          nc=TRIMLEN(cregion)
          irname(nstart+1:nstart+nc)=cregion(1:nc)
          nstart=nstart+nc
      enddo
      lt25(421)=nstart
      do i=1,110
          lt50(i)=nstart
          read(iu1,'(i5,2x,a)') nr,cregion
          if(i.le.100) then
              k=i+199
          else
              k=209+(i-101)*10
          endif
          if(k.ne.nr) then
              write(6,'(''sequence error at:'',i5,a,i5)') nr,cregion,i
              stop'REGREAD>>> error reading file REGIONSNAMEN.DAT'
          endif
          nc=TRIMLEN(cregion)
          irname(nstart+1:nstart+nc)=cregion(1:nc)
          nstart=nstart+nc
      enddo
      lt50(111)=nstart
      close(iu1)
c      write(6,'(''anzahl worte fuer regions namen:'',i6)') nstart
c  jetzt: flinn-engdahl-regionen
cu      inputfile='/users/kradi/velest/regionskoord.dat'
      call FREEUNIT(iu2)
cVMS      open(iu2,file=regkoordfile,status='OLD',readonly)
      open(iu2,file=regkoordfile,status='OLD')
      do j=1,14400,8
          read(iu2,'(8f9.3)') (xnrlon(i),i=j,j+7)
      enddo
      close(iu2)
c
      return
      end ! of subr. regread
