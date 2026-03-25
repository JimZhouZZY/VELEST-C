c
c
c
c
c
c
      subroutine GEOKO(x,y,xlat,xlon,i)    ! Urs Kradolfer, 24.3.87
c                         ->         1
c                         <-        -1
c
c Conversion of Swiss-coordinates X,Y [km] to LAT,LON [deg]    i = +1
c Conversion of LAT,LON [deg] to Swiss-coordinates X,Y [km]    i = -1
c
      implicit none
      real x,y,xlat,xlon
      integer i
      real seichmy
      if(.not.(i.eq.1.or.i.eq.-1)) stop'GEOKO>>> specify conversion !!!'
c
      if(i.eq.1)  call EBELL(x,y,xlon,xlat,seichMY)
      if(i.eq.-1) call ELLEB(xlon,xlat,x,y)
c
      end
