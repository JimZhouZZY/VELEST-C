c
c
c
c
c
c
      subroutine LAYERS (nx,ny,nz,x,y,z,vel,xe,ye,xr,yr,nl,v,vsq)
c
c  Subr. layers converts a 3 - dimensional velocity model to a
c  1 - dimensional layered model by averaging slownesses between
c  specified event and receiver positions.  the input
c  model must consist of velocities specified at grid
c  points at the intersections of 3 orthogonal grid planes.
c  slownesses are averaged along a segment joining the
c  horizontal event and receiver positions in each grid plane
c  specified by z(k).
c
c  input:
c           nx,ny,nz - numbers of grid planes
c     x(i),y(j),z(k) - grid point coordinates
c         vel(i,j,k) - velocity at point x(i),y(j),z(k)
c              xe,ye - horizontal event coordinates
c              xr,yr - horizontal receiver coordinates
c                 nl - number of layers (should equal nz)
c
c  output:
c           v(l) - velocity of layer l
c         vsq(l) = v(l) ** 2
c
      implicit none
      integer nx,ny,nz,nl
      real xe,ye,xr,yr
      real x(nx),y(ny),z(nz),vel(nx,ny,nz),v(nl),vsq(nl)
      real w(4),u(100)
      integer i,ie,j,je,ir,jr,imin,jmin,numseg,numpts,l
      integer ipts,ipoint,jpoint
      real a,b,c,d,denom,vpt1,vpt2,vpt3,vpt4,vpoint
      real xseg,yseg,xpoint,ypoint
c
c   locate the event and the reciever within horizontal
c   grid rectangles
c
      i=1
      continue
23119 if(.not.(x(i).lt.xe))goto 23120
      i=i+1
      goto 23119
23120 continue
      ie=i-1
      if(.not.(ie.eq.0))goto 23121
      ie=1
23121 continue
      j=1
      continue
23123 if(.not.(y(j).lt.ye))goto 23124
      j=j+1
      goto 23123
23124 continue
      je=j-1
      if(.not.(je.eq.0))goto 23125
      je=1
23125 continue
      i=1
      continue
23127 if(.not.(x(i).lt.xr))goto 23128
      i=i+1
      goto 23127
23128 continue
      ir=i-1
      if(.not.(ir.eq.0))goto 23129
      ir=1
23129 continue
      j=1
      continue
23131 if(.not.(y(j).lt.yr))goto 23132
      j=j+1
      goto 23131
23132 continue
      jr=j-1
      if(.not.(jr.eq.0))goto 23133
      jr=1
23133 continue
      imin=min0(ie,ir)
      jmin=min0(je,jr)
c
c   choose the number of points to use for averaging slownesses
c   and calculates x and y increments between points
c
      numseg=iabs(ir-ie)+iabs(jr-je)+1
      numpts=numseg+1
      xseg=(xr-xe)/(numseg*1.0)
      yseg=(yr-ye)/(numseg*1.0)
      do 23135l=1,nl
      u(l)=0
23135 continue
23136 continue
c
c   loop to sum slownesses in each layer, u(l), at points
c   between the event and the reciever
c
      do 23137ipts=1,numpts
      xpoint=xe+(ipts-1)*xseg
      ypoint=ye+(ipts-1)*yseg
c
c   locate (xpoint,ypoint) within a grid rectangle
c
      i=imin
      continue
23139 if(.not.(xpoint.gt.x(i)))goto 23140
      i=i+1
      goto 23139
23140 continue
      ipoint=i-1
      if(.not.(ipoint.eq.0))goto 23141
      ipoint=1
23141 continue
      j=jmin
      continue
23143 if(.not.(ypoint.gt.y(j)))goto 23144
      j=j+1
      goto 23143
23144 continue
      jpoint=j-1
      if(.not.(jpoint.eq.0))goto 23145
      jpoint=1
23145 continue
c
c   choosing weighting factors for interpolation within
c   a horizontal rectangle on the velocity grid
c
      a=xpoint-x(ipoint)
      b=x(ipoint+1)-xpoint
      c=ypoint-y(jpoint)
      d=y(jpoint+1)-ypoint
      denom=(a+b)*(c+d)
      w(1)=b*d/denom
      w(2)=a*d/denom
      w(3)=b*c/denom
      w(4)=a*c/denom
c
c   loop to interpolate for the velocity, vpoint, at
c   (xpoint,ypoint) in each layer and add its
c   reciprocal to the sum of the slownesses in that layer
c
      do 23147l=1,nl
      vpt1=w(1)*vel(ipoint,jpoint,l)
      vpt2=w(2)*vel(ipoint+1,jpoint,l)
      vpt3=w(3)*vel(ipoint,jpoint+1,l)
      vpt4=w(4)*vel(ipoint+1,jpoint+1,l)
      vpoint=vpt1+vpt2+vpt3+vpt4
      u(l)=u(l)+1/vpoint
23147 continue
23148 continue
23137 continue
23138 continue
c
c   loop to calculate the velocity from the sum of the
c   slownesses for each layer
c
      do 23149l=1,nl
      v(l)=(numpts*1.0)/u(l)
      vsq(l)=v(l)**2
23149 continue
23150 continue
      return
      end ! of subr. layers
