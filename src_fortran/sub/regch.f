c
c
c
c
c
c
      subroutine REGCH(x,y,place,nreg)
      implicit none
      real x,y
      integer nreg
      integer indfe,lt50,lt25,nx,ny,nrpu,i1,i2
      real xnrlon,dx,dy,xx,yy
      character irname*24300
      COMMON /REGIONcom/ INDFE(730),LT50(111),LT25(421),
     &                   irname,XNRLON(14400)
      character place*32
c
      nreg=0
      IF(X.LT.480.) RETURN
      IF(X.GE.865.) RETURN
      IF(Y.LE. 62.) RETURN
      IF(Y.GT.302.) RETURN
      dx=17.5
      dy=12.0
      xx=X-480.
      yy=302.-Y
      nx=xx/dx
      NY=yy/dy
      nrpu=ny*20+nx
      nreg=nrpu
      if(nx.ge.20) then
        nrpu=400+ny
        nreg=ny*20 + 19
      endif
      i1=lt25(nrpu+1)+1
      i2=lt25(nrpu+2)
      nreg=nreg+1000
      place=irname(i1:i2)
      return
      END ! of subr. regch
