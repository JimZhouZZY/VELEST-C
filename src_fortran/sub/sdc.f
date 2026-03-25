c
c
c
c
c
c
      subroutine SDC(x,y,xlat,xlon,i)     ! 7. Jan. 1988
c
c     SDC needs subroutines SETORG, DIST and REDIST
c
c
      implicit none
      real x,y,xlat,xlon
      integer i
c
      real olat,olon,aa,bb,bc,sint,cost,rotate
      integer icoordsystem
      common/GEO_COORSYSTEM/olat,olon,
     &                      rearth,ellip,rlatc,rad,
     &                      aa,bb,bc,sint,cost,rotate, icoordsystem
      double precision rearth,ellip,rlatc,rad
c
c
cc      call setorg()  ! is called separately !!
c
      if(.not.(i.eq.-1.or.i.eq.1)) stop'SDC>>> specify conversion !!!'
c
      if(i.eq.1)  call redist(x,y,xlat,xlon)
      if(i.eq.-1) call dist(xlat,xlon,x,y)
      write(*,*) 'DEBUG: sdc xlat = ', xlat
      write(*,*) 'DEBUG: sdc xlon = ', xlon
      write(*,*) 'DEBUG: sdc i = ', i
      write(*,*) 'DEBUG: sdc x = ', x
      write(*,*) 'DEBUG: sdc y = ', y
c
      end ! of subr. SDC
