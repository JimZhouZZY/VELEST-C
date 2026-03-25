c
c
c
c
c
c
      subroutine DIRECT1(nl,v,vsq,thk,jl,tkj,delta,depth,tdir,u,x)
c
c       for the direct seismic ray from an event to a receiver in
c  a layered velocity structure, direct predicts the travel time, the
c  sine of the takeoff angle, and the horizontal distance of travel in
c  the event layer.  the receiver must be located at the top of layer
c  1 and the event must be located below layer 1.  low velocity
c  layers are permitted.
c       to find the takeoff angle of the ray, a numerical approach
c  is required.  the basic scheme adopted here is the method of false
c  position.  (see acton, 1970, 'numerical methods that work,' for
c  example.)  first, the characteristics of the fastest layer
c  between the event and the surface are determined.  these permit
c  placing definite lower and upper bounds, ua and ub, on the
c  sine of the takeoff angle.  in turn, lower and upper bounds, xa
c  and xb, on the horizontal travel distance in the event layer are
c  determined.  the total horizontal travel distance for a ray with
c  with horizontal travel distance x in the event layer is denoted
c  by del, and the zero of del - delta is obtained by using xa and
c  xb as initial guesses for x in the method of false position
c  from x and tkj, the depth of the event below the top of the event
c  layer, the sine of the takeoff angle, u , is calculated.
c       from u and x, tdir is found by summing the travel time in
c  each layer.  finally, a slight correction to tdir is made, based
c  on the misfit between the final del and delta.
c
c  input:     nl - number of layers
c           v(l) - velocity of layer l
c         vsq(l) = v(l) ** 2
c         thk(l) - thickness of layer l
c             jl - event layer
c            tkj - depth of event from top of event layer
c          delta - horizontal distance between event and receiver
c          depth - depth of event
c
c  output:  tdir - direct ray travel time
c              u - sine of the takeoff angle
c              x - horizontal travel distance in the event layer
c
c
c  find the fastest layer, lmax, above and including jl
c
      implicit none
      integer nl,jl
      real tkj,delta,depth,tdir,u,x
      real v(nl),vsq(nl),thk(nl)
c
      real del,tklmax,vlmax,ua,ub,uasq,ubsq,usq
      real xa,xb,dela,delb,ubdiv,xtest
      integer lmax,j1,l,kount
      del=0.0   ! U.K. 28. Jan. 1987
c
      lmax=jl
      tklmax=tkj
      vlmax=v(jl)
      j1=jl-1
      do 23184l=1,j1
      if(.not.(v(l).gt.vlmax))goto 23186
      lmax=l
      tklmax=thk(l)
      vlmax=v(l)
23186 continue
23184 continue
23185 continue
C CHANGE BY E.KISSLING MARCH 1984
      IF(tklmax.le.0.05) tklmax=0.05
c
c   find initial bounds on the sine of the takeoff angle
c
      ua=(v(jl)/vlmax)*delta/sqrt(delta**2+depth**2)
      ub=(v(jl)/vlmax)*delta/sqrt(delta**2+tklmax**2)
c
c   calculate horizontal travel distances
c
      uasq=ua**2
      ubsq=ub**2
C CHANGE BY E.KISSLING MARCH 1984
      IF(UBSQ.GE.1.) UBSQ=0.99999
      IF(UASQ.GE.1.) UASQ=0.99999
      xa=tkj*ua/sqrt(1.0-uasq)
      if(.not.(lmax.eq.jl))goto 23188
      xb=delta
      goto 23189
23188 continue
      xb=tkj*ub/sqrt(1.0-ubsq)
23189 continue
      dela=xa
      delb=xb
      do 23190l=1,j1
      dela=dela+thk(l)*ua/sqrt(vsq(jl)/vsq(l)-uasq)
      ubdiv=sqrt(vsq(jl)/vsq(l)-ubsq)
      if(ubdiv.GT.1.e-20) GOTO 1002
      write(16,*)'WARNING:'
      write(16,1000) jl,l,lmax,vsq(jl),vsq(l),ubsq,delta,tklmax,vlmax
 1000 format(/,2x,'subr. direct1: ',3i4,2f10.4,f15.6,3f10.5,/)
      ubdiv=1.e-20
 1002 continue
      delb=delb+thk(l)*ub/sqrt(vsq(jl)/vsq(l)-ubsq)
23190 continue
23191 continue    !  NOT used... !!!
c
c   loop to find the zero of del-delta by teh method of false position
c
      do 23192kount=1,25
      if(.not.((delb-dela).lt.0.02))goto 23194
      x=0.5*(xa+xb)
      u=x/sqrt(x**2+tkj**2)
      usq=u**2
      goto 23193
23194 continue
      x=xa+(delta-dela)*(xb-xa)/(delb-dela)
      u=x/sqrt(x**2+tkj**2)
      usq=u**2
      del=x
      do 23196l=1,j1
      del=del+thk(l)*u/sqrt(vsq(jl)/vsq(l)-usq)
23196 continue
23197 continue
      xtest=del-delta
      if(.not.(abs(xtest).lt.0.02))goto 23198
      goto 23193
23198 continue
      if(.not.(xtest.lt.0.0))goto 23200
      xa=x
      dela=del
      goto 23201
23200 continue
      xb=x
      delb=del
23201 continue
23192 continue
23193 continue
c
c   calculate direct ray travel time
c
c
      if(del.eq.0.0) del=x   ! U.K. 28. Jan. 1987
c
      tdir=(sqrt(x**2+tkj**2))/v(jl)
      do 23202l=1,j1
      tdir=tdir+thk(l)*v(jl)/(vsq(l)*sqrt(vsq(jl)/vsq(l)-usq))
23202 continue
23203 continue
      tdir=tdir-(u/v(jl))*(del-delta)
      return
      end ! of subr. direct1
