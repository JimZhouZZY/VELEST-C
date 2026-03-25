c
c
c
c
c
c
      subroutine REFRACT (nl,v,vsq,thk,jl,tkj,delta,
     &kk,tref,didjkk,xovmax)
c
c       for refracted rays in a layered earth model, refract
c  determines the fastest travel time, tref, the layer
c  in which the fastest ray is refracted, kk, the
c  critical distance for refraction in that layer,
c  didjkk, and an upper bound on delta for which a
c  direct ray can be a first arrival, xovmax.  refract
c  allows for the possibility of low velocity layers.
c       note that there may not be a refracted ray, either because
c  all layers below the event layer are low velocity layers or
c  because for all layers below the event layer which are not low
c  velocity layers the critical distance exceeds delta.  in such
c  cases tref, didjkk, and xovmax are set very large, kk is set to
c  zero, and refract returns to the calling program.
c
c  input:  nl - number of layers
c        v(l) - velocity of layer l
c      vsq(l) - v(l) ** 2
c      thk(l) - thickness of layer l
c          jl - event layer
c         tkj - depth of event in event layer
c       delta - horizontal distance between event and receiver
c
c  output:   kk - refracting layer for fastest refracted ray
c          tref - travel time of fastest refracted ray
c        didjkk - critical distance for refraction in layer kk
c        xovmax - an upper bound on delta for which the direct ray can
c                       be the first arrival
c  internal arrays:
c
c       tr(m) - travel time for refraction in layer m
c     tinj(m) - traveltime intercept
c      tid(m) - terms in travel time intercept which are
c                     independent of tkj
c     didj(m) - critical distance
c      did(m) - terms in critical distance which are
c                     independent of tkj
c
c
c  call subr. tiddid to evaluate tid(m) and
c  did(m), the terms in the travel time intercept and
c  critical distance for a ray refracted in layer m
c  which are independent of tkj.
c
c
c  determine tref, kk, didjkk
c
      implicit none
c
      integer nl,jl,kk,j1,m,lx,m1,l,jx
      real tkj,delta,tref,didjkk,xovmax
      real sqt,tim
      real v(nl),vsq(nl),thk(nl)
      real tid(100),did(100),tinj(100),didj(100),tr(100)
      call tiddid(jl,nl,v,vsq,thk,tid,did)
      tref=100000.
      j1=jl+1
      do 23151m=j1,nl
      if(.not.(tid(m).eq.100000.))goto 23153
      tr(m)=100000.
      goto 23154
23153 continue
      sqt=sqrt(vsq(m)-vsq(jl))
      tinj(m)=tid(m)-tkj*sqt/(v(m)*v(jl))
      didj(m)=did(m)-tkj*v(jl)/sqt
      tr(m)=tinj(m)+delta/v(m)
      if(.not.(didj(m).gt.delta))goto 23155
      tr(m)=100000.
23155 continue
23154 continue
      if(.not.(tr(m).lt.tref))goto 23157
      tref=tr(m)
      kk=m
      didjkk=didj(kk)
23157 continue
23151 continue
23152 continue
c
c   if there is no refracted ray:
c
      if(.not.(tref.eq.100000.))goto 23159
      didjkk=100000.
      xovmax=100000.
      kk=0
      return
23159 continue
c
c   if there is a refracted ray, determine xovmax:
c   find lx, the 1st layer below the event layer which
c   is not a low velocity layer
c
      m=jl+1
      continue
23161 if(.not.(tid(m).eq.100000.))goto 23162
      m=m+1
      goto 23161
23162 continue
      lx=m
c
c   check whether the event is in the 1st layer
c
      if(.not.(jl.eq.1))goto 23163
      xovmax=tinj(lx)*v(lx)*v(1)/(v(lx)-v(1))
      return
23163 continue
      m=jl
c
c   find jx, the 1st layer above and including the event
c   layer which is not a low velocity layer
c
      continue
23165 continue
      tid(m)=0.
      m1=m-1
      do 23168l=1,m1
      if(.not.(vsq(m).le.vsq(l)))goto 23170
      tid(m)=100000.
      goto 23171
23170 continue
      sqt=sqrt(vsq(m)-vsq(l))
      tim=thk(l)*sqt/(v(l)*v(m))
      tid(m)=tid(m)+tim
23171 continue
23168 continue
23169 continue
      m=m-1
c
c  decide whether or not jx=1 and calculate xovmax
c
23166 if(.not.(tid(m+1).lt.100000..or.m.eq.1))goto 23165
23167 continue
      if(.not.(tid(m+1).lt.100000.))goto 23172
      jx=m+1
      xovmax=(tinj(lx)-tid(jx))*v(lx)*v(jx)/(v(lx)-v(jx))
      goto 23173
23172 continue
c
c   jx=1
c
      xovmax=tinj(lx)*v(lx)*v(1)/(v(lx)-v(1))
23173 continue
      return
      end ! of subr. refract
