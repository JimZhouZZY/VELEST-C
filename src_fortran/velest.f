      program velest
c
c     version 3.3  E.Kissling, Institute of Geophysics, ETH Zurich
c
c                  21. August 1998
c
c  this version is running as: 
c                           
c                               -  velest.f for SUN-Unix 
c                               -  velestUNIXhp.f   for HP-Unix
c
c
c    the program is based on the version 3.0 HRM,UK,EK,and W.Ellsworth
c    of 28.oct93, and on version 3.1 EK 21April95.
c    (see PhD thesis by U.Kradolfer 1989, and by H.R. Maurer 1993
c     at Institute of Gephysics, ETH Zurich, Switzerland)
c
c *********************************************************************
c
c  This version of VELEST can be used to simultaneously invert a large
c  number of earthquake travel time data for hypocentral and velocity
c  model (P and S velocities) parameters. Alternatively (isingle=1) it
c  may also be used as a single event location program.
c
c  This version may read from input data file in the converted (*.cnv)
c  format (ised=0), in a velest-archive type (*arcvel) format (ised=1),
c  and in the new SED (Swiss Seismological Service) format (ised=2).
c  VELEST optionally writes output data files in converted and summary
c  card format and in singleevent mode writes a location file for
c  archive purposes (print output file).
c
c
cEK addition to 3.0:
c  9.nov93  version 3.01 added subr. inputarcvel (ised=1)
c  10.nov93              modified subr. inputparam for model file
c  10.nov93              modified  "        "      for print ouput in
c                                          single event location mode
c  15.nov93              modified in subr. traveltime no calls to 
c                                 routines rejectobs and reviveobs
c  16. nov93             modified in subr. gapcalc to avoid division
c                        zero by zero for call to atan2
c  23. nov93             modified subr inputcnv to recognize end of file
c
c  28. nov93             testing singleevent and multievent mode and
c                        corrections for print output in subr. inputparam
c
c  30. nov93             corrected error in inputparam for isingle and
c                        no of eqs .ne.0:  neqs and legs MUST be one
c                        for isingle=1
c
c implemented dez 1994 EK:   read rotate angle of coordinate system
c                                        as in earlier versions
c                              subr. inputparam read neqs,nshots,rotate
c
c  also dez94 EK:          if ifx(i)=0 normally
c                                   =1 to inhibit adjustement of (rotated)
c                                      y coordinate (see subr. MATRIX)
c                                   (i.e. dtdr(2)=0.0 if ifx(i)=1)
c
c
c EK version 3.3 as of 21August1998:  implemented and improved display of
c                          P and S phase usage and station delays
c                                          
cEK  still to be implemented:  flat Earth approximation
c
cEK
c
c   *vel.out file for STATIS written in subr. STATISLOUT
c
c------------------------------------------------------------------------------
c-----START OF SOURCE-CODE OF PROGRAM VELEST (MAIN PROGRAM) -------------------
c------------------------------------------------------------------------------
c
c
      implicit none
      include 'inc/vel_com.inc'
c
      integer i4enabletest
      integer juliam
      integer nittc,istopflag,k,j,i,l,iresflag,ii
      integer ier,iminold,itime
      real damp,xlat,xlon,cpmintot,d1,d2
      real ggti_v(16),gtg_v(16)
      logical better
c
      character cline*80
c
      character*20 ctime
      real cpusec
c      integer errx,erry,errz
c
c--------------------------------------------
c
c     reset internal cpu-timer:
c
      call CPUTIMER(cpusec)
c
c test to AVOID runs of velest with compilation-option   NOi4
      i4enabletest=1234567890
c end of this 'mini-test'
c
c     save report of this VELEST-run (write it ev. to file16 in subr. INPUTPARAM
c
      call DATETIME(ctime)  ! get date&time from system
      write(headerline(1),'(1x,''>>> Start of program VELEST at '',a20,
     &              '' <<<'')') ctime
      write(6,'(1x,''>>> Start of program VELEST at '',a20,
     &          '' <<<'')') ctime
      write(6,*)
c
      write(headerline(2),*)'::::::::::::::::::::::::::::::::::::::'
      write(headerline(3),*)' V E L E S T  - Version : January 3, 1995'
c
c      reset...
c
      nitt=0   ! number of iterations made
      nittc=0  !   "          "       for wich CPU is calculated
      nvar=0   ! number of unknowns to solve for
      istopflag=0
c
c     input control-parameters, model and stations
c
      call INPUTPARAM
      if(.not.single_turbo) write(16,*)'~~~ input parameters read.'
c
c     if a file containing all the raypaths (ray-points) of all events
c     is desired, open the file (unit=13):
c
      if(irayout.eq.1) open(13,file=rayfile,
     &                      status='unknown')
      if(idrvout.eq.1) open(21,file=drvfile,
     &                      status='unknown')
c
c     file07 - final hypocenters and travel times ( *.CNV format, for next it.)
c     file11 - summary cards of final hypocenters (for plotting)
c
      if(icnvout.eq.1) open(7,file=cnvfile,
     &                      status='unknown')
      if(ismpout.eq.1) open(11,file=smpfile,
     &                      status='unknown')
      if(isingle.ne.0) open(2,file=velfile,
     &                      status='unknown')  !  statisL-compatible output 
      if(ialeout.eq.1) open(75,file=alefile,
     &                      status='unknown')  !  x,y, ALE of all events
      if(idspout.eq.1) open(76,file=dsprfile,
     &                      status='unknown')  !  x,y, DSPR of all events
      if(irflout.eq.1) open(77,file=rflfile,
     &                      status='unknown')  !  x,y, resid. of reflected ray
      if(irfrout.eq.1) open(78,file=rfrfile,
     &                      status='unknown')  !  x,y, resid. of refracted ray
      if(iresout.eq.1) open(79,file=resfile,
     &                      status='unknown')  !  dist, residual of ray
c
      k=0
ccc      call getsymbol('zirvax')  ! not used in this version
c
 1080 continue     ! START OF THE LOOP OVER  O N E   E V E N T FOR SINGLE-TASK !
c
      if(k.eq.-1)then ! means: end of input-file detected !!!
         close(2)
cek
c    no smp-file for single_event_mode since smp only in SED special format
c    use *.cnv file for plotting etc.
cek
	 if(ismpout.eq.1.and.isingle.eq.1) close(11)
         call DATETIME(ctime)  ! get date&time from system
         call CPUTIMER(cpusec)
         if(.not.single_turbo)then
         write(16,'(x,''>>>   End of program VELEST at '',a20,
     &           '' <<<'','' [ CPU-sec : '',f7.1,'']'')') ctime,cpusec
         endif
         write(6,'(x,''>>>   End of program VELEST at '',a20,
     &           '' <<<'','' [ CPU-sec : '',f7.1,'']'')') ctime,cpusec
         stop'...end...(VELEST was running in SINGLE-EVENT-OPTION)'
      endif
c
c----------------------------------------------------------------------------
c
c     do forward problem first:
c
      nitt=0
c
c     do scaling of damping-factors (and determine ICOUNT):
c
      if(.not.single_turbo)
     &              write(16,*)'~~~ set units for first iteration ...'
c                 1 --> so icount will not be 0 !!!
      call SETUNT(1,invertratio,nsinv,icount,
     &              xythet,stathet,othet,vthet,zthet,scale)
c
c     determine nr of unknowns NVAR & nr of equations LIP to be solved :
c
      if(.not.single_turbo)
     &           write(16,*)'~~~ determine number of unknowns ...'
      call DETNOFUNKNOWNS
c
c     reset G-matrix and RHT-vector :
c
      if(.not.single_turbo) write(16,*)'~~~ reset G and RHT ...'
      do j=1,kvar
         g(j)=0.0
      enddo
      do k=1,nvar
         rht(k)=0.0
      enddo
c
c     for each event do: reset avres&rms, input of the data,
c                        raytracing and accumulate normal equations
c                        ( = setup of matrix G and vector RHT ) :
c
      if(irayout.eq.1) rewind(13) ! rewind file (if desired) with raypaths
      if(idrvout.eq.1) rewind(21) ! rewind file (if desired) with raypaths
      if(.not.single_turbo)
     &        write(16,*)'~~~ input data, raytracing, setup G ...'
      if(.not.single_turbo)then
         write(16,'(/////)')
         write(16,'(''    I N P U T - D A T A  '',/)')
         if(nsp.eq.2) then
         write(16,'(1x,'' eq    origin-time    latitude longitude '',
     &        ''depth    x      y      z    mag ifxOBStot obsP obsS'')')
         else
         write(16,'(1x,'' eq    origin-time    latitude longitude '',
     &                 ''depth    x      y      z    mag ifxOBS'')')
         endif
         endif
      do i=1,legs
         avres(i)=0.
         rms(i)=0.
c        it's the starting: input trial hypocenter and phase-data
         call INPUTDATA(i)
         k=knobs(i)
         if(k.eq.-1)then
            goto 1080   ! means: end of input-file (data) detected
         endif
c
c        for each arrival do raytracing, residual-calc. and setup G & RHT :
c
         do l=1,k
            call TRAVELTIME(i,l,iresflag)
            call SETUPMATRIXG(i,l)
         enddo
         if(isingle.ne.0)then
            if( (knobs(i)-nobswithw0) .lt. nvar )then
               iabort=1
               if(.not.single_turbo)then
                  write(16,*)'knobs(i)-nobswithw0 < nvar !!!'
                  write(16,'('' knobs(i)='',i4,'' nobswithw0='',i4,
     +            '' nvar='',i2)') knobs(i),nobswithw0,nvar
                  write(16,*)'Event cannot be located!!!'
               endif
               write(6,*)'knobs(i)-nobswithw0 < nvar !!!'
               write(6,'('' knobs(i)='',i4,'' nobswithw0='',i4,
     +            '' nvar='',i2)') knobs(i),nobswithw0,nvar
               write(6,*)'Event cannot be located!!!'
               goto 98989
            endif
            if(.not.single_turbo)then
               write(16,*)'~~~ compute singular values of G ...'
            endif
            call SINGULARVALUES(i)
         endif
      enddo
c
c     now all the events are read-in; now check, how many stations
c     appear in the input-data (= nofstainput)
c
      call ACTUALSTATIONS
c
c     save selected columns of g for resolution calculations
c
      if(.not.single_turbo)then
         write(16,*)
         write(16,*)'~~~ store parts of G ...'
      endif
      do k=1,4
         call STOREG(k,k)
      enddo
      if(scale(6).ne.0.0)then
         i=4*neqs+nshot+1
         j=i+nltot-1
         ii=4
         do k=i,j
            ii=ii+1
            call STOREG(k,ii)
         enddo
      endif
c
c     all events: raytracing done, residuals calculated
c                 and normal equations accumulated
c
c     calculate rms & data variance for all events:
c
      if(.not.single_turbo) write(16,*)'~~~ compute RMS and DATVAR ...'
      call RMSDATVAR
c
      if(ittmax.eq.0)then
         if(.not.single_turbo)then
            write(16,*)'WARNING:  ITTmax=0  -->  no iteration is made'
         endif
         write(6,*)'WARNING:  ITTmax=0  -->  no iteration is made'
         goto 9999
      endif
c
c     calculate GAP, average-residuals of ray-types and check,
c     whether variance/msqrd res have decreased (this is
c     always the case if NITT=0 !) or not (in the
c     latter case do hypocenter/model backup):
c
      if(.not.single_turbo) write(16,*)'~~~ save initial datvar ...'
      call CHECKSOLUTION(istopflag,better)
c
c
 10   continue     !  START OF THE LOOP FOR ONE   I T E R A T I O N
c
c     reset CPU-timer for this iteration
c
      call DATETIME(ctime)  ! get date&time from system
      call CPUTIMER(cpusec)
      if(isingle.eq.0)then
         write(6,'(1x,'' finished Iteration #'',i3,''  at'',
     &             3x,a20,5x,''CPU-sec = '',f7.1)')
     &             nittc,ctime,cpusec
      endif
      nittc=nittc+1
c
c
      nitt=nitt+1
c
c     print header of this iteration :
c
      cline='----------------------------------------'
     &    //'----------------------------------------'
      if(.not.single_turbo)then
         write(16,'(///,a,///)') cline
         write(16,'(11x,''ITERATION no '',i2,/)') nitt
      endif
c
c     do scaling of damping-factors (and determine ICOUNT) for this iteration:
c
      if(.not.single_turbo)then
         write(16,*)'~~~ set units for this iteration ...'
      endif
      call SETUNT(nitt,invertratio,nsinv,icount,
     &                 xythet,stathet,othet,vthet,zthet,scale)
c
c     set nr of unknowns back to nr for this iteration
c
      if(.not.single_turbo)then
        write(16,*)'~~~ determine number of unknowns for this iter. ...'
      endif
      call DETNOFUNKNOWNS
c
c     print out the damping-factors to be used in this iteration :
c
      if(isingle.eq.0)then
         if(.not.single_turbo)then
            write(16,*)'Damping factors to be used in this iteration:'
            write(16,'(1x,''  Othet='',f8.3,''   XYthet  ='',f8.3,
     &                    ''   Zthet='',f7.3)') othet,xythet,Zthet
            write(16,'(1x,''STAthet='',f8.3,''   Vthet   ='',f8.3,/)')
     &                 stathet,vthet
         endif
      endif
c
c     apply damping to diagonal elements of G-matrix :
c
      if(.not.single_turbo) write(16,*)'~~~ damp G ...'
      call DAMPG
c
c     copy right-hand-side of normal equations (RHT=At*RES) to vector B :
c
      do k=1,nvar
         b(k)=rht(k)
      enddo
c
c
c     Cholesky-decomposition of matrix G :
c
      if(.not.single_turbo) write(16,*)'~~~ solve G * b = RHT ...'
         write(6,'(''DEBUG LUDECP(pre): nitt='',i4,'' nvar='',i6,
     &          '' kvar='',i10,'' g1='',1pe12.4,'' gdiag='',1pe12.4,
     &          '' rht1='',1pe12.4)')
     &          nitt,nvar,kvar,g(1),g((nvar*(nvar+1))/2),rht(1)
      call LUDECP(g,g,nvar,d1,d2,ier)
      if(ier.ne.0)then   ! matrix not positive definit !!!
         if(.not.single_turbo)then
            write(16,'('' WARNING: error in ludecp     ier='',i5)')ier
         endif
         write(6,'('' WARNING: error in ludecp     ier='',i5)')ier
         STOP'error in (damped) matrix G !!!'
      endif
c
c     solve normal equations (simple elimination; G is now lower triangular) :
c
      call LUELMP(g,b,nvar,b)
c
c     put solution into proper units :
c
      if(.not.single_turbo) write(16,*)'~~~ fix units ...'
      call FIXUNT(b,neqs,nshot,nltot,ksta,scale,
     &            vdamp,itotmodels,inltot,nplay(1))
c
c     check velocity-change;
c     if large, apply step-length-damping to adjustment-vector !
c
      if(.not.single_turbo)then
         write(16,*)'~~~ apply eventually step-length damping ...'
      endif
      call STEPLENGTHDAMP(damp)
c
c     adjust hypocenters, velocities & stationcorrections :
c
      if(.not.single_turbo)then
       write(16,*)'~~~ adjust model (hypocenters, stacorr & veloc.) ...'
      endif
      call ADJUSTMODEL(damp)
c
c     calculate total step-length of this iteration and print it out:
c
      if(.not.single_turbo)then
         write(16,*)'~~~ compute effective step-length ...'
      endif
      call STEPLENGTHCALC
c
c
 1010  continue   ! come here, if backup has been performed   !
c
c     do forward problem for solution of iteration# NITT:
c
c
c     do scaling of damping-factors (and determine ICOUNT) for NEXT iteration:
c
      if(ibackups.lt.4)then
         if(.not.single_turbo)then
            write(16,*)'~~~ set units for next iteration ...'
         endif
         call SETUNT(nitt+1,invertratio,nsinv,icount,
     &                      xythet,stathet,othet,vthet,zthet,scale)
      else
         if(.not.single_turbo)then
            write(16,*)'~~~ set units for this iteration ...'
         endif
         call SETUNT(nitt,invertratio,nsinv,icount,
     &                    xythet,stathet,othet,vthet,zthet,scale)
      endif
c
c     determine nr of unknowns NVAR & nr of equations LIP to be solved for NEXT:
c
      if(.not.single_turbo)then
       if(ibackups.lt.4)then
        write(16,*)'~~~ determine number of unknowns for next iter. ...'
       else
        write(16,*)'~~~ determine number of unknowns for this iter. ...'
       endif
      endif
      call DETNOFUNKNOWNS
c
c     reset G-matrix and RHT-vector :
c
      if(.not.single_turbo) write(16,*)'~~~ reset G and RHT ...'
      do j=1,kvar
         g(j)=0.0
      enddo
      do k=1,nvar
         rht(k)=0.0
      enddo
c
c     for each event do: reset avres&rms, input (if 1st iteration),
c                        raytracing and accumulate normal equations
c                        ( = setup of matrix G and vector RHT ) :
c
      if(irayout.eq.1) rewind(13) ! rewind file (if desired) with raypaths
      if(idrvout.eq.1) rewind(21) ! rewind file (if desired) with raypaths
      if(.not.single_turbo) write(16,*)'~~~ raytracing, setup G ...'
      do i=1,legs
         avres(i)=0.
         rms(i)=0.
         k=knobs(i)
c
c        for each arrival do raytracing, residual-calc. and setup G & RHT :
c
         do l=1,k
            call TRAVELTIME(i,l,iresflag)
            call SETUPMATRIXG(i,l)
c
chrm    Calculate the matrix G for the data resolution matrix
	    if(isingle.ne.0)then
	       gg2(l,1) = 1 * w(l,1)
	       gg2(l,2) = dtdr(1) * w(l,1)
	       gg2(l,3) = dtdr(2) * w(l,1)
	       gg2(l,4) = dtdr(3) * w(l,1)
	    endif
         enddo
c
         if(isingle.ne.0)then
            if( (knobs(i)-nobswithw0) .lt. nvar )then
               iabort=1
               if(.not.single_turbo)then
                  write(16,*)'knobs(i)-nobswithw0 < nvar !!!'
                  write(16,*)'Event cannot be located!!!'
               endif
               write(6,*)'knobs(i)-nobswithw0 < nvar !!!'
               write(6,*)'Event cannot be located!!!'
               goto 98989
            endif
            if(.not.single_turbo)then
               write(16,*)'~~~ compute singular values of G ...'
            endif
            call SINGULARVALUES(i)
         endif
      enddo
c
c     save selected columns of g for resolution calculations
c
c     Store parts of matrix G here (valid for next iteration), before SETUNT
c     for this iteration!
c
      if(.not.single_turbo) write(16,*)'~~~ store parts of G ...'
      do k=1,4
         call STOREG(k,k)
      enddo
      if(scale(6).ne.0.0)then
         i=4*neqs+nshot+1
         j=i+nltot-1
         ii=4
         do k=i,j
            ii=ii+1
            call STOREG(k,ii)
         enddo
      endif
c
c     do scaling of damping-factors (and determine ICOUNT) for CURRENT iter.:
c
      if(.not.single_turbo)then
         write(16,*)'~~~ set units back for current iteration ...'
      endif
      call SETUNT(nitt,invertratio,nsinv,icount,
     &                 xythet,stathet,othet,vthet,zthet,scale)
c
c     set nr of unknowns back to nr for this iteration
c
      if(.not.single_turbo)then
        write(16,*)'~~~ determine number of unknowns for this iter. ...'
      endif
      call DETNOFUNKNOWNS
c
c     all events: raytracing done, residuals calculated
c                 and normal equations accumulated
c
c     calculate rms & data variance for all events:
c
      if(.not.single_turbo)then
         write(16,*)'~~~ compute RMS and DATVAR ...'
      endif
      call RMSDATVAR
c
c     calculate GAP, average-residuals of ray-types and check,
c     whether variance/msqrd res have decreased (this is
c     always the case if NITT=1 !) or not (in the
c     latter case do hypocenter/model backup):
c
      if(.not.single_turbo)then
         write(16,*)'~~~ check solution (better? worse?) ...'
      endif
      call CHECKSOLUTION(istopflag,better)
c
c better=.false. if variance inc'd and residual inc'd, so do readjustments
c if backup was done, but less than 4 times, do forward problem again !!
c
      if(ibackups.lt.4)then
         if(better)then
            if(ibackups.lt.4) ibackups=0
         else
            if(ibackups.eq.0)then
               if(.not.single_turbo)then
                  call NITTOUTPUT(damp)
               endif
            endif
            if(ibackups.lt.4)then
               if(.not.single_turbo)then
                  write(16,'('' msqrd res increased'',
     &                       '' --> hypocenter and model backup.'')')
                  write(16,*)'~~~ solution is worse --> backup ...'
               endif
               call BACKUP
               ibackups=ibackups+1
               call STEPLENGTHCALC
               goto 1010
            endif
         endif
      else
         if(.not.single_turbo)then
            write(16,*)'4 times backup made!'
         endif
      endif
c
c     come here, if no backup made or 4 backups made
c
      if(.not.single_turbo)then
         write(16,*)'~~~ iteration done; output results ...'
      endif
      if(.not.single_turbo) call NITTOUTPUT(damp)
c
      if(.not.single_turbo) write(16,*)'~~~ another iteration? ...'
c
c  IF 4 TIMES 'HALF ADJUSTMENTS MADE', SOLUTION IS FINAL --> STOP ITERATIONS
c
      if(ibackups.eq.4)then
         if(.not.single_turbo)then
            write(16,*)'4 TIMES BACKUP MADE --> sTOP ITERATIONS.'
         endif
         goto 9999
      endif
c
c     IF  STEP-LENGTH < DELMIN ( MINIMUM STEP-LENGTH ) --> STOP ITERATIONS :
c
      if(steplen.gt.0)then  ! otherwise we are during backups...
         if(steplen.lt.delmin)then
            if(.not.single_turbo)then
               write(16,*)'STEP-LENGTH < DELMIN  --> stop iterations.'
            endif
            goto 9999
         endif
      endif
c
c     IF  #ITERATIONS = ITTMAX  --> STOP ITERATIONS :
c
      if(nitt.eq.ittmax)then
         if(.not.single_turbo)then
            write(16,*)'NITT = ITTMAX  --> stop iterations.'
         endif
         goto 9999
      endif
c
c     SINGLE-EVENT-MODE: IF CHANGES IN DATVAR < 1.E-6 --> STOP ITERATIONS :
c
      if(isingle.ne.0.and.istopflag.eq.1)then
         if(.not.single_turbo)then
            write(16,*)'Changes in datvar < 1.e-6  --> stop iterations.'
         endif
         goto 9999
      endif
c
c     DO NEXT ITERATION :
c
      ibackups=0
      GOTO 10
c
c----------------------------------------------------------------------------
c
c     FINAL solution reached!   output final solution of all events,
c                               all residuals and station-corrections:
c
 9999 continue
c
      call CPUTIMER(cpusec)
      if(isingle.eq.0)then
         write(6,'(1x,'' finished Iteration #'',i3,''  at'',
     &             3x,a20,5x,''CPU-sec = '',f7.1)')
     &             nittc,ctime,cpusec
      endif
      nittc=nittc+1
c
      if(.not.single_turbo)then
         write(16,*)'~~~ final solution reached ...'
         write(16,*)
         if(iresolcalc.gt.0)then
         else
            write(16,*)'~~~ damp G & Chol.-Decomp.: NOT made '//
     &                 '(switch iresolcalc is NOT set !)'
         endif
      endif
      if(iresolcalc.gt.0)then
c
c        G was calculated for next iteration, so compute nvar, kvar etc. for
c        next iteration as well:
         if(.not.single_turbo) 
     &      write(16,*)'~~~ set units for next iteration ...'
         call SETUNT(nitt+1,invertratio,nsinv,icount,
     &                      xythet,stathet,othet,vthet,zthet,scale)
         if(.not.single_turbo) 
     &      write(16,'(a)')' ~~~ determine number of unknowns '//
     &                          'for next iter. ...'
         call DETNOFUNKNOWNS
c
         if(.not.single_turbo) 
     &      write(16,*)'~~~ damp G ...'
         call DAMPG   ! necessary for subr. RESOLCOVAR !!!
         if(.not.single_turbo)then
            write(16,*)'~~~ Cholesky decomposition of G ...'
         endif
         call LUDECP(g,g,nvar,d1,d2,ier)   ! necessary for subr. RESOLCOVAR !
         if(ier.ne.0)then   ! matrix not positive definit !!!
            if(.not.single_turbo)then
              write(16,'('' WARNING: error in ludecp     ier='',i5)')ier
            endif
            write(6,'('' WARNING: error in ludecp     ier='',i5)')ier
            STOP'error in (damped) matrix G !!!'
         endif
      endif
c
      if(isingle.ne.0.and.icoordsystem.eq.2)then
       if(.not.single_turbo)then
          write(16,*)'~~~ compute magnitude ...'
       endif
       do i=1,legs
        if(iyr(i).lt.100)then
           iminold=JULIAM(iyr(i)+1900,imo(i),iday(i),ihr(i),imin(i))
        else
           iminold=JULIAM(iyr(i),imo(i),iday(i),ihr(i),imin(i))
        endif
        call TIMECLEAR(iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i),
     &                                                      itime)
        do j=1,knobs(i)
           pt(j,i)=pt(j,i)+(iminold-itime)*60.
        enddo
        call MAGNITUDE(i,itime)
        emag(i)=xmagnitude
c       standard-deviation of emag(i)=sdxmagnitude
c       the magnitudes for each station are now stored on XMAGNI(1...nobs)
       enddo
      endif
c
      if(.not.single_turbo)then
         write(16,*)
         if(isingle.ne.0)then
            write(16,*)'Damping factors to be used in this iteration:'
            write(16,'(1x,''  Othet='',f8.3,''   XYthet  ='',f8.3,
     &                    ''   Zthet='',f8.3)') othet,xythet,Zthet
            write(16,'(1x,''STAthet='',f8.3,''   Vthet   ='',f8.3)')
     &                 stathet,vthet
            write(16,*)
            write(16,*)'Final solution of event#',isingle
            write(16,*)'------------------------------------'
         else !         print final rms-values of all the events
            write(16,59) (rms(i),i=1,legs)
 59         format(' Final rms values: ',10f5.2,/,(19x,10f5.2))
            write(16,'(//)')
         endif
      endif
c
      if(isingle.ne.0)then
         if(icoordsystem.eq.2)then
            call REGION(1,-e(2,1),e(3,1),regionname,nreg,
     &                                   regnamfile,regkoordfile)
         else
            call SDC(e(2,1),e(3,1),xlat,xlon,1) ! calc. LAT/LON
            call REGION(2,xlat,-xlon,regionname,nreg,
     &                       regnamfile,regkoordfile) ! Subr. REGION needs LON E
         endif                                        ! as input (E = positive)
      endif
      if(.not.single_turbo)then
         write(16,*)'~~~ output final hypocenters ...'
      endif
      call FINALHYPOCOUT
      if(.not.single_turbo)then
       if(isingle.ne.0)then
       write(16,'('' Singular values:     '',4f10.4,5x,''ALE ='',f7.3)')
     &                (SV(j),j=1,4) , ale(1)
          write(16,'(1x,''data variance = '',f10.6)') davar1
       endif
      endif
c
      if(.not.single_turbo)then
         write(16,*)'~~~ output final station residuals ...'
         call FINALSTARESI
      endif
      if(iresolcalc.gt.0)then
         if(.not.single_turbo)then
           write(16,*)'~~~ output resolution- & covariance matrices ...'
         endif
         call RESOLCOVAR(davar1)
      endif
c      if (s(2).gt.99.9) then
c         errx = 999
c      else
c      endif
c     if (s(3).gt.99.9) then
c        erry = 999
c      else
c        erry = nint(s(3)*10.)
c      endif
c      if (s(4).gt.99.9) then
c         errz = 999
c      else
c         errz = nint(s(4)*10.)
c      endif
c      write(smpline(46:54),'(3i3)'),errx,erry,errz
cek     if(isingle.eq.0) then
cek         write(11,'(a80)') smpline
cek      endif
      if(ialeout.eq.1)then
         write(75,'(2x,''ALE'',3f10.4,''        .1'')')
     &              -e(2,1),e(3,1), ale(1)
      endif
      if(idspout.eq.1)then
         write(76,'(2x,''DSP'',3f10.4,''        .1'')')
     &              -e(2,1),e(3,1), spread
      endif
c
c     output all the statistics-stuff...
c
      if(isingle.eq.0)then
         write(16,*)'~~~ output statistics ...'
         call STATISTICSOUT
      endif
c
c     close output-files:
c
      if(iturbo.eq.0)then
         close(7)
cek         if(ismpout.eq.1.and.isingle.eq.0) close(11)
      endif
      if(irayout.eq.1)then
         close(13)
         close(21)
      endif
c
c     output cpu-statistics for whole VELEST-run:
c
      if(.not.single_turbo) write(16,*)
      call CPUTIMER(cpusec)
      if(isingle.eq.0)then
         write(16,*)'CPU - statistics:'
         write(16,*)
      endif
      if(.not.single_turbo)then
         cpmintot=cpusec/60.
         write(16,'(1x,''TOTAL CPU-sec      '',3x,'' ='',f7.1,5x,
     &              ''(CPU-minutes:'',f7.3,'')'')') cpusec,cpmintot
         write(16,*)
      endif
c
      if(isingle.eq.0)then
         call DATETIME(ctime)  ! get date&time from system
         write(16,'(x,''>>>   End of program VELEST at '',a20,
     &           '' <<<'','' [ CPU-sec : '',f7.1,'']'')') ctime,cpusec
         write(6,'(1x,''>>>   End of program VELEST at '',a20,
     &           '' <<<'','' [ CPU-sec : '',f7.1,'']'')') ctime,cpusec
      endif
c
c
c     in case of 'single-event-mode' output NITT, GAP,
c                      SPREAD of resolution and final hypocenter in
c                      statisL-compatible
c
98989 continue   ! return address for 'insufficient data' !!!
      if(isingle.gt.0)then
         call DATETIME(ctime)  ! get date&time from system
         if(iabort.ne.1)
     &      write(6,'(1x,''GAP = '',i3,''   NITT = '',i2,
     &                ''   D-Spread ='',f5.2)')igap(1),nitt,spread
         write(6,*)
         if(.not.single_turbo)then
            write(16,*)'end of event# ',isingle
            write(16,*)
            write(16,*)
         endif
c
c        output statisL/SED-compatible summary in file02   ( *.VEL )
c
         if(.not.single_turbo)then
         write(16,*)'~~~ output solution in HYPO71-compatible format...'
         endif
c
chrm     Calculation of the data resolution matrix
c        gg2 = Matrix G
c        ggt = G transposed
c	 gtg = G * G transposed
c	 ggti = Inverse of gtg
c        ggg = Generalized inverse
c        drm = data resolution matrix
c
chrm     Setup ggt
c
	 call matrtran(gg2,maxobsperevent,4,ggt)
c
chrm     Setup gtg
c
	 call matrmult(ggt,4,maxobsperevent,gg2,maxobsperevent,4,gtg,4,4)
c
chrm     Damping of gtg ...
c
	 gtg(1,1) = gtg(1,1) + othet
	 gtg(2,2) = gtg(2,2) + xythet
	 gtg(3,3) = gtg(3,3) + xythet
	 gtg(4,4) = gtg(4,4) + zthet
c
chrm     Convert gtg to a vectorized form
c
	 do i=1,4
	    do k=1,4
	       gtg_v((i-1)*4 + k) = gtg(i,k)
	    enddo
	 enddo
	 call matrinv(4,gtg_v,ggti_v)
c
chrm     Reconvert the vectorized matrix ggti_v to a 4 x 4 matrix
c
	 do i=1,4
	    do k=1,4
	       ggti(i,k) = ggti_v((i-1)*4 + k)
	    enddo
	 enddo
c
chrm     Calculate the generalized inverse
c
	 call matrmult(ggti,4,4,ggt,4,maxobsperevent,ggg,4,
     &              maxobsperevent)
c
chrm     Calculate the data resolution matrix
c
	 call matrmult(gg2,maxobsperevent,4,ggg,4,maxobsperevent,drm,
     &             maxobsperevent,maxobsperevent)
c
chrm     Data resolution matrix calculation finished!
c
cek   print output *out.VEL for singleevent locations written by
c     subr. STATISLOUT:
c
         call STATISLOUT
c
c        reset parameters...
c
         iabort=0
         isingle=isingle+1
         nitt=0
         nvar=0
         icount=0
         istopflag=0
         iresflag=0
         nittc=0
         ibackups=0
         ale(1)=0.0
         do i=1,3
            isconstrain(i)=0
         enddo
         iconstrain(1)=0
         emag(1)=0.0
         nmag=0
         do i=1,knobs(1)
            xmagni(i)=0.0
            amx(i)=0.0 ! not necessary in this version;
            prx(i)=0.0 ! just in case, AMX&PRX will be used for *.CNV files also
            istm(i,1)=0
         enddo
         ifixsolution=0
         zmin=zmininput
         goto 1080 ! goto NEXT EVENT (SINGLE-EVENT-MODE !)
      endif
c
c     come here if NOT 'single-event-mode' to terminate program...
c
      close(8) ! close data-input-file ( *.CNV )
c
      goto 99999
c
  998 call OPENERROR('main program','data-input-file FOR008')
c
99999 continue
      END
