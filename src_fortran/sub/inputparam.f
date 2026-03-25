c
c
c
c
c
c
      subroutine INPUTPARAM ! old name was: INPUT1 ;
c                           ! this subr. needs another input-format !!
c
c     reads in the control-file (type VELEST.CMN) which contains all the
c     control-parameters for a VELEST-run
c
      implicit none
      include '../inc/vel_com.inc'
      real olat,olon,avelo,z,dx,dy,dz,rotate
      integer ifil,j,ifl,mode,icc,ml,jndex,k
c
      integer trimlen, i,n
      character*1 reflch, cns,cew
      character*40 titl
      character*80 card,line(32), titleline
      logical lexist
c
c     Open control-inputfile:
c
      inquire(file='velest.cmn',exist=lexist)
      if(.not.lexist)then
         stop'INPUTPARAM>>> control file `velest.cmn` not found!'
      endif
cVMS      open(10,file='velest.cmn',status='old',err=9910,readonly)
      open(10,file='velest.cmn',status='unknown',err=9910)
c
c     input center of coordinate system
c
      i=0
 111  read(10,'(a)',end=222) card
      if(card(1:1).eq.'*') goto 111
      i=i+1
      if(i.gt.32) stop'INPUTPARAM>>> control-file not correct!'
      line(i)=card
      goto 111
 222  continue
      if(.NOT.(i.eq.24.or.i.eq.32.))
     &    stop'INPUTPARAM>>> control-file not correct!'
      read(line(1),'(a)') titleline
      read(line(2),*) olat,olon,icoordsystem,zshift,itrial,ztrial,ised
      read(line(3),*) neqs,nshot,rotate
      read(line(4),*) isingle, iresolcalc
      read(line(5),*) dmax,itopo,zmininput,veladj,zadj,lowveloclay
      read(line(6),*) nsp, swtfac,vpvs, nmod
      read(line(7),*) othet,xythet,zthet,vthet,stathet
      read(line(8),*) nsinv,nshcor,nshfix, iuseelev,iusestacorr
      read(line(9),*) iturbo, icnvout,istaout,ismpout
      read(line(10),*) irayout,idrvout,ialeout,idspout,
     &                irflout,irfrout,iresout
      read(line(11),*) delmin,ittmax,invertratio
      read(line(12),'(a)') modelfilename
      read(line(13),'(a)') stationfilename
      read(line(14),'(a)') seismofilename
      read(line(15),'(a)') regnamfile
      read(line(16),'(a)') regkoordfile
      read(line(17),'(a)') topo1file
      read(line(18),'(a)') topo2file
      read(line(19),'(a)') phasefile
      read(line(20),'(a)') shotfile
c     Output files:
      read(line(21),'(a)') outfile
      if(outfile.eq.' ') outfile='vel.out'
      read(line(22),'(a)') velfile
      if(velfile.eq.' ') velfile='velout.vel'
      read(line(23),'(a)') cnvfile
      if(cnvfile.eq.' ') cnvfile='velout.cnv'
      read(line(24),'(a)') stafile
      if(stafile.eq.' ') stafile='velout.sta'
c
c     the next few ouputfiles are not very often used, therefore
c     either all or none of them have to be specified in the controlfile:
c
      if(i.eq.32)then
         read(line(25),'(a)') smpfile
         if(smpfile.eq.' ') smpfile='velout.smp'
         read(line(26),'(a)') rayfile
         if(rayfile.eq.' ') rayfile='velout.ray'
         read(line(27),'(a)') drvfile
         if(drvfile.eq.' ') drvfile='velout.drv'
         read(line(28),'(a)') alefile
         if(alefile.eq.' ') alefile='velout.ale'
         read(line(29),'(a)') dsprfile
         if(dsprfile.eq.' ') dsprfile='velout.dspr'
         read(line(30),'(a)') rflfile
         if(rflfile.eq.' ') rflfile='velout.rfl'
         read(line(31),'(a)') rfrfile
         if(rfrfile.eq.' ') rfrfile='velout.rfr'
         read(line(32),'(a)') resfile
         if(resfile.eq.' ') resfile='velout.res'
      endif
c
      single_turbo=.false.
      if(isingle.eq.1.and.iturbo.eq.1) single_turbo=.true.
c
c     open the main-output-file:
c
      if(.not.single_turbo)then
cVMS         open(16,file=outfile,status='new',carriagecontrol='list')
         open(16,file=outfile,status='unknown')
         write(16,'(a)') headerline(1)
         write(16,*)
         write(16,'(a)') headerline(2)
         write(16,'(a)') headerline(3)
         write(headerline(3),*)' (Authors: see source code)'
         write(16,'(a)') headerline(3)
         write(16,'(a)') headerline(2)
         write(16,*)
         write(16,*)
         write(16,*)'Title of this VELEST run:'
         write(16,*)
         n=trimlen(titleline)
         write(16,*) titleline(1:n)
         do i=1,n
            titleline(i:i)='-'
         enddo
         write(16,*) titleline(1:n)
         write(16,*)
         write(16,*)
         write(16,*)'Current array-dimensions of program VELEST:'
         write(16,*)
         write(16,'('' Max. number of '',
     &              ''EARTHQUAKES for simult. inversion IEQ = '',i3)')
     &              ieq
         write(16,'('' Max. number of '',
     &              ''SHOTS for simult. inversion    INSHOT = '',i3)')
     &              inshot
         write(16,'('' Max. number of '',
     &              ''(different) 1D-MODELS      ITOTMODELS = '',i3)')
     &              itotmodels
         write(16,'('' Max. number of '',
     &              ''LAYERS per one 1D-model        INLTOT = '',i3)')
     &              inltot
         write(16,'('' Max. number of '',
     &              ''STATIONS in stationlist           IST = '',i3)')
     &              ist
         write(16,'('' Max. number of '',
     &              ''OBSERVATIONS per event MAXOBSPEREVENT = '',i3)')
     &              maxobsperevent
         write(16,*)
         write(16,*)
      endif
c
cek 210495: check for reasonable switch combinations nmod,nsp
c
      if(nsp.lt.1) nsp=1
      if(nsp.gt.3) nsp=3
c nsp=1 P-model and p-data only
c nsp=2 P-model (+data) and S-model (+data) independently used
c nsp=3 P-model (+data) and S-model dependent on P-model (vpvs-factor used)
      if(nsp.eq.2) then
        nmod=2
      else
        nmod=1
      endif
c
cek 210495
c
c     print the input-parameters
c
      if(.not.single_turbo)then
         write(16,'(///)')
         write(16,*)'   INPUT - P A R A M E T E R S :'
         write(16,'(//)')
      write(16,'(a)')'***'
      write(16,'(a)')'***  olat       olon   icoordsystem      '//
     &              'zshift   itrial ztrial    ised'
      write(16,'(1x,f10.6,f11.6,6x,i1,10x,f7.3,7x,i1,3x,f6.2,6x,i1)')
     &          olat,olon,icoordsystem,zshift,itrial,ztrial,ised
      write(16,'(a)')'***'
      write(16,'(a)')'*** neqs   nshot    rotate'
      write(16,'(1x,i6,1x,i6,3x,f6.1)') neqs,nshot,rotate
      write(16,'(a)')'***'
      write(16,'(a)')'*** isingle   iresolcalc'
      write(16,'(7x,i1,10x,i1)') isingle, iresolcalc
      write(16,'(a)')'***'
      write(16,'(a)')'*** dmax    itopo    zmin'//
     &              '     veladj    zadj   lowveloclay'
      write(16,'(2x,f7.2,5x,i1,4x,f6.2,6x,f5.2,3x,f5.2,7x,i1)')
     &          dmax,itopo,zmininput,veladj,zadj,lowveloclay
      write(16,'(a)')'***'
      write(16,'(a)')'*** nsp    swtfac   vpvs       nmod'
      write(16,'(5x,i1,5x,f5.2,3x,f6.3,7x,i2)') nsp, swtfac, vpvs, nmod
      write(16,'(a)')'***'
      write(16,'(a)')'***   othet   xythet    zthet    vthet   stathet'
      write(16,'(4x,5(f7.2,2x))') othet,xythet,zthet,vthet,stathet
      write(16,'(a)')'***'
      write(16,'(a)')'*** nsinv   nshcor   nshfix'//
     &              '     iuseelev    iusestacorr'
      write(16,'(7x,3(i1,7x),4x,i1,12x,i1)') nsinv, nshcor, nshfix,
     &                         iuseelev,iusestacorr
      write(16,'(a)')'***'
      write(16,'(a)')'*** iturbo    icnvout   istaout   ismpout'
      write(16,'(7x,i1,9x,i1,9x,i1,9x,i1)') iturbo,
     &                                      icnvout,istaout,ismpout
      write(16,'(a)')'***'
      write(16,'(a)')'*** irayout   idrvout   ialeout   idspout'//
     &               '   irflout   irfrout   iresout'
      write(16,'(7x,i1,6(9x,i1))') irayout,idrvout,ialeout,idspout,
     &                            irflout,irfrout,iresout
      write(16,'(a)')'***'
      write(16,'(a)')'*** delmin   ittmax   invertratio'
      write(16,'(3x,f6.3,6x,i2,9x,i2)') delmin, ittmax, invertratio
      write(16,'(a)')'***'
      write(16,'(a)')'*** Modelfile:'
      write(16,'(a)') modelfilename(1:trimlen(modelfilename))
      write(16,'(a)')'***'
      write(16,'(a)')'*** Stationfile:'
      write(16,'(a)') stationfilename(1:trimlen(stationfilename))
      write(16,'(a)')'***'
      write(16,'(a)')'*** Seismofile:'
      write(16,'(a)') seismofilename(1:trimlen(seismofilename))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with region names:'
      write(16,'(a)') regnamfile(1:trimlen(regnamfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with region coordinates:'
      write(16,'(a)') regkoordfile(1:trimlen(regkoordfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File #1 with topo data:'
      write(16,'(a)') topo1file(1:trimlen(topo1file))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File #2 with topo data:'
      write(16,'(a)') topo2file(1:trimlen(topo2file))
      write(16,'(a)')'***'
      write(16,'(a)')'*** DATA INPUT files:'
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with Earthquake data:'
      write(16,'(a)') phasefile(1:trimlen(phasefile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with Shot data:'
      write(16,'(a)') shotfile(1:trimlen(shotfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** OUTPUT files:'
      write(16,'(a)')'***'
      write(16,'(a)')'*** Main print output file:'
      write(16,'(a)') outfile(1:trimlen(outfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with single event locations:'
      write(16,'(a)') velfile(1:trimlen(velfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with final hypocenters in *.cnv format:'
      write(16,'(a)') cnvfile(1:trimlen(cnvfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with new station corrections:'
      write(16,'(a)') stafile(1:trimlen(stafile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with summary cards (e.g. for plotting):'
      write(16,'(a)') smpfile(1:trimlen(smpfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with raypoints:'
      write(16,'(a)') rayfile(1:trimlen(rayfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with derivatives:'
      write(16,'(a)') drvfile(1:trimlen(drvfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with ALEs:'
      write(16,'(a)') alefile(1:trimlen(alefile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with Dirichlet spreads:'
      write(16,'(a)') dsprfile(1:trimlen(dsprfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with reflection points:'
      write(16,'(a)') rflfile(1:trimlen(rflfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with refraction points:'
      write(16,'(a)') rfrfile(1:trimlen(rfrfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with residuals:'
      write(16,'(a)') resfile(1:trimlen(resfile))
      write(16,'(a)')'***'
      write(16,'(///)')
      endif
c
c     remarks concerning OLAT & OLON:
c     Both olat and olon are input in decimal-degrees.
c     For LAT N(orth) and LON W(est) their values are positive.
c     For LAT S(outh) and LON E(ast) their values must be negative!
c
c     For more details see comments on top of this program-source!
c
      zmin=zmininput
c
      if(.not.single_turbo)then
         write(16,*)'Origin of cartesian coordinates:'
         if(olon.le.0.0)then
           write(16,'(6x,f10.5,2hN ,f11.5,1hE)')olat,-olon
         else
           write(16,'(6x,f10.5,2hN ,f11.5,1hW)')olat,olon
         endif
         if(olon.eq.0.0.and.olat.eq.0.0)then
            write(16,*)'Origin BERNE is taken!'
            write(16,*)'OLAT=46.95240 n OLON=7.439583 e'
         endif
         write(16,*)
         write(16,*)' X,Y-AXES rotated clockwise from North'
         write(16,*)'          Rotation angle (rotate):'
         write(16,'(12x,7hrotate=,1x,f6.1,1x,5hdegr.)')rotate
         write(16,*)
c
         if(itrial.gt.0)then
            write(16,*)
            write(16,*)' ============================================'
            write(16,*)' trial epicenter ~ earliest station'
            write(16,*)' trial depth = ztrial = ',ztrial
            write(16,*)' ============================================'
            write(16,*)
         endif
c
         if(icoordsystem.eq.2)then
            write(16,*)'SWISS COORDINATES will be used instead of the',
     &                 ' short distance conversion!'
            write(16,*)'origin BERNE: x=600.km, y=200.km'
            write(16,*)'x is positive towards east,',
     &                 ' y is positive towards north'
ccc            write(16,*)'---> no rotation can be performed !!!'
         else
            write(16,*)'icoordsystem = ',icoordsystem
            write(16,*)'normal SHORT DISTANCE CONVERSION will be made'
            write(16,*)'x is positive towards west,',
     &                 ' y is positive towards north'
         endif
      endif
c
c     use setorg to set up coordinate system
c
      ifil=0
      if(.not.single_turbo) ifil=16
c
      if(icoordsystem.ne.2) call SETORG(olat,olon,rotate,ifil)
c
c  now origin is setup.
c
4     format(f5.2,i5,4f5.2,f6.2,3i1,5f6.2)
c
      if(dmax.eq.0.)then
         if(.not.single_turbo)then
            write(16,*)'WARNING: dmax was zero ! ... set to 150 km !'
         endif
         dmax=150.
      endif
c
      rmsmin=0.0
      if(invertratio.le.0) invertratio=999  ! do not invert for sta-corr & model
      if(.not.single_turbo)then
         if(iuseelev.eq.0)then
            write(16,*)
            write(16,*)'Station-elevations internally set to ZERO !'
            write(16,*)'(but correctly printed in file12)'
            write(16,*)
         endif
         if(iusestacorr.eq.0)then
            write(16,*)
            write(16,*)'Station-corrections set to ZERO !'
            write(16,*)'(if you do NOT invert for station-corrections,'
            write(16,*)' then these 0.0-values are printed in file12 !)'
            write(16,*)
         endif
      endif
c***  nsinv=  0 no inversion for station corrections
c***          1 inversion for station corrections
c
c***  legs is the total number of shots and quakes
cek
cek   legs and neqs MUST be 1 in single_event mode
c
      if(isingle.ne.1) then
         legs=neqs+nshot
      else
         legs=1
         neqs=1
         nshot=0
      endif
c
      if(.not.single_turbo)then
         write(16,'(///)')
         write(16,*)'   INPUT - M O D E L :'
         write(16,'(//)')
         write(16,*)'Model(s) read from file :'
         write(16,*) modelfilename
         write(16,*)
      endif
      close(10)
c
c   ***********
c NOW READ model file
c   ***********
c
cVMS      open(10,file=modelfilename,status='old',err=9911,readonly)
      open(10,file=modelfilename,status='old',err=9911)
c
      ireflector=0
cek
      if(nsp.eq.3) then
         write(16,7788)
7788  format(//,2x,'attention nsp=3! >s-data used for p-model too!',//)
      endif
c
cek  modifications for new model file
      read(10,'(a40)') titl
      if(.not.single_turbo) write(16,'(a)') 
     +                 ' model file title: ',titl
cek      read(10,*) (nplay(j),j=1,nmod)
8     format(10i5)
      do 14 i=1,nmod
      if(.not.single_turbo)then
         write(16,17) i
17       format(1x,'Velocity structure for model',i4,' :')
         write(16,*)
         write(16,11)
11       format(1h ,'layer    vel   depth   vdamp  reflector')
      endif
      titl=' '
      read(10,'(i3)') nplay(i)
      do 9 j=1,nplay(i)
       if(j.eq.1)then
         read(10,1212) vp(i,j),hp(i,j),vdamp(i,j),reflch,titl
1212     format(f5.2,5x,f7.2,2x,f7.3,3x,a1,1x,a40)
       else
         read(10,12) vp(i,j),hp(i,j),vdamp(i,j),reflch
12       format(f5.2,5x,f7.2,2x,f7.3,3x,a1)
       endif
      if(reflch.ne.' ')then
         if(reflch.eq.'m'.or.reflch.eq.'M')then
            reflchar=reflch
            ireflector=j
            if(.not.single_turbo)then
               if(vp(i,ireflector).gt.8.0)then
                  write(16,*)'WARNING ::: velocity ABOVE reflector is',
     &                       ' greater than 8.0 km/s  .....'
                  write(6,*)'WARNING ::: velocity ABOVE reflector is',
     &                       ' greater than 8.0 km/s  .....'
                  write(6,*)
               endif
            endif
         else
            if(.not.single_turbo)then
              write(16,*)'WARNING:'
              write(16,*)'Reflector indicated in velocity-model is',
     &                   ' marked with a : ',reflch
              write(16,*)'Only   m   or   M   are allowed!'
              write(16,*)'Reflected phases will be ignored (wrong mark)'
              write(6,*)'Reflected phases will be ignored (wrong mark)'
              write(6,*)
            endif
         endif
      endif
      if(j.eq.1)then
         if(.not.single_turbo)then
            write(16,'(2x,i2,3x,2f7.2,1x,f7.2,1x,a1,1x,a40)') j,
     &                vp(i,j),hp(i,j),vdamp(i,j),reflch,titl
         endif
         titl=' '
      else
         if(.not.single_turbo)then
            write(16,13) j,vp(i,j),hp(i,j),vdamp(i,j),reflch
13          format(2x,i2,3x,2f7.2,1x,f7.2,1x,a1)
         endif
      endif
c
      if(j.gt.1.and.hp(i,j).lt.0.0)then
         write(6,*)'WARNING:'
         write(6,*)'Only top of the first layer can be negative !!'
         write(6,*)
         if(.not.single_turbo)then
            write(16,*)'WARNING:'
            write(16,*)'Only top of the first layer can be negative !!'
         endif
      endif
c
    9 continue
c
c    calculate and print average velocities of the model i :
c
      if(.not.single_turbo)then
        ifl=1
        write(16,*)
        write(16,*)'Calculation of average velocity starts at layer # ',
     &             ifl
        avelo=0
        do j=ifl+1,nplay(i)
           avelo=avelo + ( hp(i,j)-hp(i,j-1) ) * vp(i,j-1)
           write(16,'('' z ='',f6.2,'' ...'',f6.2,''    v = '',f4.2,
     &               ''    AVEL = '',f4.2,'' downto z = '',f6.2)')
     &           hp(i,j-1),hp(i,j),vp(i,j-1),avelo/(hp(i,j)-hp(i,ifl)),
     &           hp(i,j)
        enddo
        write(16,*)
c
        ifl=2
        write(16,*)
        write(16,*)'Calculation of average velocity starts at layer # ',
     &             ifl
        avelo=0
        do j=ifl+1,nplay(i)
           avelo=avelo + ( hp(i,j)-hp(i,j-1) ) * vp(i,j-1)
           write(16,'('' z ='',f6.2,'' ...'',f6.2,''    v = '',f4.2,
     &               ''    AVEL = '',f4.2,'' downto z = '',f6.2)')
     &           hp(i,j-1),hp(i,j),vp(i,j-1),avelo/(hp(i,j)-hp(i,ifl)),
     &           hp(i,j)
        enddo
        write(16,*)
        write(16,*)
      endif
c
14    continue
      i f (ireflector.ne.0) t h e n
         if(.not.single_turbo)then
         write(16,*)
        write(16,'(1x,''Phases in the input-phaselist marked with a  '',
     &            a1,/,'' are treated as reflections from the bottom '',
     &            ''of layer nr. '',i2)') reflchar,ireflector
            if(lowveloclay.eq.1)then
               write(16,*)'Switch LOWVELOCLAY is set to 1, but'//
     &                    ' reflected phases are allowed to occur'
               write(16,*)'This is improper (no low velocity-layers '//
     &                    'allowed for reflected waves!!)'
               lowveloclay=0
               write(16,*)'LOWVELOCLAY now set to 0 '
               write(6,*)'Switch LOWVELOCLAY is set to 1, but'//
     &                    ' reflected phases are allowed to occur'
               write(6,*)'This is improper (no low velocity-layers '//
     &                    'allowed for reflected waves!!)'
               write(6,*)'LOWVELOCLAY now set to 0 '
            endif
         endif
      e n d i f
c
c
c   ***********
c NOW READ station file
c   ***********
c
c    read in station data:
c
      if(.not.single_turbo)then
         write(16,'(///)')
         write(16,*)'   INPUT - S T A T I O N S :'
         write(16,'(//)')
         write(16,*)'Station-parameters read from file :'
         write(16,*) stationfilename
         write(16,*)
      endif
      close(10)
cVMS      open(10,file=stationfilename,status='old',err=9912,readonly)
      open(10,file=stationfilename,status='old',err=9912)
      read(10,1) fm
1     format(a80)
      if(.not.single_turbo)then
         write(16,*)
         write(16,2)
2        format(1h ,3x,' stn latitude longitude elev',
     &   3x,'x      y      z',5x,'ptcor stcor model icc')
      endif
      nsta=0
c
10    nsta=nsta+1
      read(10,fm) stn(nsta),xla(nsta),cns,xlo(nsta),cew,
     &            ielev(nsta),mode,icc,
     &            ptcor(nsta),stcor(nsta)
      call CASEFOLD(cns)
      call CASEFOLD(cew)
      if(cns.eq.'S') xla(nsta)=-xla(nsta)
      if(cew.eq.'E') xlo(nsta)=-xlo(nsta)
c
      if(stn(nsta).eq.' ') goto 41
cek      if(mode.eq.0) mode=1
      mode=1
c
      if(nsp.eq.2) then
         model(2*nsta-1)=mode
         model(2*nsta)=mode+1
      else
         model(nsta)=mode
      endif	 
      z=-ielev(nsta)/1000.
      if(icoordsystem.eq.2)then
         call GEOKO(dx,dy,xla(nsta),-xlo(nsta),-1) ! calc. cart. coord.
         dx=-dx
      else
         call SDC(dx,dy,xla(nsta),xlo(nsta),-1) ! calc. cart. coord.
      endif
      dz=z
      x(nsta,1)=dx
      x(nsta,2)=dy
      if(iuseelev.eq.1)then
         x(nsta,3)=dz
      else
         x(nsta,3)=0.     ! station-elevation not used !
      endif
      if(iusestacorr.eq.0)then
         ptcor(nsta)=0.0  ! initial (from input) station-corrections not used !
         stcor(nsta)=0.0  !   "          "               "            "    "
      endif
      if(.not.single_turbo)then
         if(icoordsystem.eq.2) dx=-dx
         if(cns.eq.'S') xla(nsta)=-xla(nsta)
         if(cew.eq.'E') xlo(nsta)=-xlo(nsta)
         write(16,27) nsta,stn(nsta),xla(nsta),cns,xlo(nsta),cew,
     &                ielev(nsta),dx,dy,dz,
     &                ptcor(nsta),stcor(nsta),mode,icc
27       format(1x,i3,1x,a4,f7.4,a1,f8.4,a1,1x,i5,3f7.2,1x,2f6.2,2i4)
         if(cns.eq.'S') xla(nsta)=-xla(nsta)
         if(cew.eq.'E') xlo(nsta)=-xlo(nsta)
      endif
c
c   an icc of 0 holds the station delay fixed.
c   the highest nonzero icc station has its p delay
c   held fixed but its s delay is allowed to float.
c
      map1(nsta)=icc
      goto 10
41    continue
      close(10)  ! stationfile
c
c     read in seismo data (seismometer-specifications and so on):
c
      if(.not.single_turbo)then
         if(seismofilename.ne.' ')then
            write(16,*)
            write(16,*)'Seismo-parameters read from file :'
            write(16,*) seismofilename
            write(16,*)
            write(16,*)'This file contains station-informations which'
            write(16,*)'are used for magnitude-determination'
            write(16,*)'(seismometer-constant,filterparameters, etc.)'
            write(16,*)
            write(16,*)
         else
            write(16,*)
            write(16,*)'NO Seismo-parameter-file specified'
            write(16,*)'--> NO magnitudes calculated!'
            write(16,*)
            write(16,*)
            write(16,*)
         endif
      endif
      if(seismofilename.ne.' ')
cVMS     &open(10,file=seismofilename,status='old',err=9913,readonly)
     &open(10,file=seismofilename,status='old',err=9913)
c     the seismo-parameters will be read in later; we must know the eventtime
c     before we can pick the time-dependent (!) seismo-parameters !!!
c
c***  nsta is the number of stations read in
      nsta=nsta-1
c
      call MAXII(nsta,map1,ml,jndex) ! determine MAX (=ml) of map1 (icc-values)
c***  ksta is the number of total station corrections
c***    to invert for
      ksta=0
      if(nsinv.ne.0)then  !  invert for station corrections
c                            the last station for p  will be held fixed.
         ksta=ml-1
         if(nsp.eq.2) ksta=2*ksta + 1
      endif
c
c    find the total number of layers
c
c**  laysum(i) is an index for the location of the first
c***  layer of model i
      laysum(1)=1
      laysum(2)=nplay(1)+1
      if(nmod.lt.3) goto 50
      do 42 i=3,nmod
42    laysum(i)=laysum(i-1)+nplay(i-1)
50    continue
c***  nltot is the total number of velocity layers
c***  to invert for
      nltot=0
      do 18 i=1,nmod
      nltot=nltot+nplay(i)
18    continue
c
      if(.not.single_turbo)then
         write(16,*)
         write(16,19) neqs,nshot,nltot
      endif
 1619 format(10x,3i20)
19    format(3x,'neqs=',i5,'  nshot=',i5,'  nltot=',i5,/)
c
      do 25 i=1,itotmodels     ! max. number of 1D-models
      do 25 j=1,inltot         ! max. number of layers per model
25    thkp(i,j)=0.0
      do 22 i=1,nmod
      if(nplay(i).eq.1) goto 22
      k=nplay(i)-1
      do 23 j=1,k
23    thkp(i,j)=hp(i,j+1)-hp(i,j)
22    continue
      RETURN
c
 9910 call OPENERROR('inputparam','control-input-file FOR010')
      return
 9911 call OPENERROR('inputparam','model-input-file (FOR010)')
      return
 9912 call OPENERROR('inputparam','station-input-file (FOR010)')
      return
 9913 call OPENERROR('inputparam','seismo-input-file (FOR010)')
      return
c
      end ! of subr. inputparam
