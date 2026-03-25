c
c
c
c
c
c
      subroutine RAYPATH(nx,ny,nz,x,y,z,vel,nl,thk,d,v,vsq,xe,ye,ze,
     & xr,yr,zr,rp,nrp,nrtn,jl,tkj,itype,ttt,MLL,sterr,direrr,refrerr,
     & reflerr,DTDDrefl,DTDHrefl)
c
c       raypath attempts to approximate the least time path for a
c  seismic ray traveling in a medium with velocities specified at
c  points on a 3-dimensional grid.  for very small delta a straight
c  line path is chosen.  for larger delta a layered velocity model
c  approximating the region between the event and the receiver is
c  constructed and the least time path in that layered model is chosen
c  both direct and refracted rays are considered and the possibility
c  of low velocity layers is allowed.
c       the 3 - dimensional grid should be specified by the
c  intersections of 3 sets of planes normal to the x,y and z coordinate
c  axes.  the intersection of the planes with the axes are specified
c  by x(i), y(j), and z(k), and the spacing of these may be
c  completely arbitrary.
c       subr. layers produces layer velocities from the
c  velocities of the 3 - dimensional model between the event and
c  the receiver, and because these will differ for every event -
c  receiver pair, layers must be called every time raypath is
c  called (unless delta is very small).  the number, thicknesses,
c  and depths of the layers should depend only on the array z(k),
c  and need only to be calculated once for a given velocity model.
c  this is why nl, thk(l), and d(l) are included as input parameters
c  for raypath.  when calling layers, the event and receiver must
c  be located within the horizontal confines of the grid, that is,
c  x(1) <= xe, xr <= x(nx) and y(1) <= ye, yr <= y(ny).
c       subrs. refract and direct determine refracted and direct
c  ray travel times and other parameters in a relatively standard
c  fashion, but are sufficiently general to account for low velocity
c  layers.
c       once the nature of the fastest ray  has been determined, the
c  calculation of a set of raypoints along the path is straightforward
c  and is performed by subr. strpath, refpath, or dirpath.  these
c  points consist of the event and receiver coordinates and, in the cases
c  of refpath and dirpath, every intermediate point at which the ray
c  crosses a boundary between layers.
c       because of the extensive branching in raypath, a parameter
c  to index the the 7 return statements in the main subr.,
c  nrtn is included and takes on values from 1 to 7 (to 8 with extension
c  of reflected waves ! U.K. 11/1986).
c
c  input:
c           nx,ny,nz - number of grid planes in x,y,z directions
c     x(i),y(j),z(k) - gridpoint coordinates
c         vel(i,j,k) - velocity at grid point x(i),y(j),z(k)
c
c                 nl - number of layers (should equal nz)
c             thk(l) - thickness of layer l
c               d(l) - depth to top of layer l
c
c           xe,ye,ze - event coordinates
c           xr,yr,zr - receiver coordinates
c
c  output:
c              v(l) - velocity of layer l
c            vsq(l) = v(l) ** 2
c               nrp - number of raypoints
c    rp(1,2,or 3,i) - coordinates of raypoint i
c
c  important internal arrays and variables:
c
c            delta - horizontal distance between source and receiver
c            depth - depth of event below receiver
c               jl - number of event layer
c              tkj - depth of event from top of event layer
c
c               kk - refracting layer for fastest refracted ray
c           didjkk - critical distance for refraction in layer kk
c             tref - refracted ray travel time
c           xovmax - an upper bound on delta for which the direct ray
c                       can be the first arrival
c
c           salpha - sine of the takeoff angle of the direct ray
c            deljl - horizontal travel distance of the
c                        direct ray in the event layer
c             tdir - direct ray travel time
c             nrtn - index to point of return from raypath
c
c  subrs. called:
c
c          layers      - constructs  averaged layered model
c
c          stpath      - straight raypath
c          strpath     - NOT USED !! straight raypath from event to receiver
c
c          direct      - for the direct ray:
c                          salpha, deljl, tdir
c          dirpath     - direct ray path
c
c          reflect1     - for reflected ray
c          reflect     - NOT USED !! for reflected ray
c          reflectpath - reflected ray path
c
c          refract     - for the fastest refracted ray:
c                         kk, didjkk, tref, xovmax
c          tiddid      - travel-time intercept & crit. distance for a ray
c          refpath     - refracted ray path
c
      implicit none
      integer nx,ny,nz,nl,nrp,nrtn,jl,itype,mll
      real xe,ye,ze,xr,yr,zr,tkj,ttt,sterr
      real direrr,refrerr,reflerr,dtddrefl,dtdhrefl
      integer l,ierr,kk
      real depth,delta,tdir,salpha,deljl,trefl
      real ain,tref,didjkk,xovmax
      real x(nx),y(ny),z(nz),vel(nx,ny,nz)
      real d(nl),thk(nl),v(nl),vsq(nl)
cek max. nr of raypoints =50
      real rp(3,200)
c
1313  continue   ! restart here in case of an error in REFLECT !!
c
      sterr=0.0
      direrr=0.0
      refrerr=0.0
      reflerr=0.0
c
      if(ze.lt.d(1)) stop'RAYPATH>>> earthquake ABOVE model!'
      if(zr.lt.d(1)) stop'RAYPATH>>> station ABOVE model!'
c
      depth=ze-zr
      delta=sqrt((xr-xe)**2+(yr-ye)**2)
      l=nl
      continue
23107 if(.not.(d(l).ge.ze))goto 23108
c
c   determine event layer
c
      l=l-1
      goto 23107
23108 continue
      jl=l
      d(1)=zr
c   adjust for nonzero depth of reciever
c
      thk(1)=d(2)-zr
      tkj=ze-d(jl)
      if(.not.(delta.lt..05))goto 23105
      call stpath(xe,ye,ze,xr,yr,zr,rp,nrp,ttt,nl,jl,tkj,v,d,thk,sterr)
c
c   if delta is small, set a straight path
c
      nrtn=1
      return
23105 continue
c
c   assign averaged layer velocities
c
      if(itype.eq.1) go to 1
      call layers(nx,ny,nz,x,y,z,vel,xe,ye,xr,yr,nl,v,vsq)
1     continue
      if(.not.(jl.eq.nl))goto 23109
c
c   consider only the direct ray if the event is in the half space
c
      call direct1(nl,v,vsq,thk,jl,tkj,delta,depth,tdir,salpha,deljl)
      call dirpath(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,jl,tkj,
     2 salpha,deljl,rp,nrp,direrr)
      nrtn=2
      ttt=tdir
      return
23109 continue
c
c
      if(mll.le.0)goto 23110
c     MLL > 0  : a reflection layer is specified; ray MUST be reflected!
ctest      call reflect(nl,v,vsq,thk,jl,tkj,delta,depth,mll,trefl,ain,ierr,
ctest     &             DTDDrefl,DTDHrefl)
      call reflect1(nl,v,vsq,thk,jl,tkj,delta,depth,mll,trefl,ain,ierr,
     &             DTDDrefl,DTDHrefl)
      if(ierr.ne.0)then
         write(6,*)'trying another ray-type...'
         write(16,*)'trying another ray-type...'
         MLL=0      ! maybe it works with another ray-type...
         goto 1313  ! retry
      endif
      call reflectpath(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,jl,tkj,
     &                 ain,mll,rp,nrp,reflerr)
c
      ttt=trefl
      nrtn=8   ! reflected wave
      RETURN
23110 continue
c
c   find the first refracted ray to arrive
c
      call refract(nl,v,vsq,thk,jl,tkj,delta,kk,tref,didjkk,xovmax)
      if(.not.(delta.gt.xovmax))goto 23111
c   for delta greater than xovmax
c   1st arrival is the refracted ray
c
      call refpath(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,tkj,jl,
     2 kk,didjkk,rp,nrp,refrerr)
      nrtn=3
      ttt=tref
      return
23111 continue
      if(.not.(jl.eq.1))goto 23113
      tdir=(sqrt(tkj**2+delta**2))/v(1)
      if(.not.(tref.lt.tdir))goto 23115
c
c   if refracted ray is 1st arrival
c
      call refpath(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,tkj,jl,
     2 kk,didjkk,rp,nrp,refrerr)
      nrtn=4
      ttt=tref
      return
23115 continue
c
c   if direct ray is 1st arrival
c
      call stpath(xe,ye,ze,xr,yr,zr,rp,nrp,ttt,nl,jl,tkj,v,d,thk,sterr)
      nrtn=5
      ttt=tdir
      return
23113 continue
c
c   for event below layer 1
c
      call direct1(nl,v,vsq,thk,jl,tkj,delta,depth,tdir,salpha,deljl)
      if(.not.(tref.lt.tdir))goto 23117
c
c   if refracted ray is 1st arrival
c
      call refpath(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,tkj,jl,
     2 kk,didjkk,rp,nrp,refrerr)
      ttt=tref
      nrtn=6
      return
23117 continue
c
c   if direct ray is 1st arrival
c
      call dirpath(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,jl,tkj,
     2 salpha,deljl,rp,nrp,direrr)
      nrtn=7
      ttt=tdir
c      write(6,'(A,I3,A,F10.6,A,I3,A,I3,A,F10.6)')
c     & 'DEBUG raypath: nrtn=',nrtn,', ttt=',ttt,', nrp=',nrp,
c     & ', jl=',jl,', tkj=',tkj
c      write(6,'(A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F10.6)')
c     & 'DEBUG raypath: delta=',delta,', depth=',depth,
c     & ', xe=',xe,', ye=',ye,', ze=',ze
c      write(6,'(A,F10.6,A,F10.6,A,F10.6)')
c     & 'DEBUG raypath: xr=',xr,', yr=',yr,', zr=',zr
      return
      end ! of subr. raypath
