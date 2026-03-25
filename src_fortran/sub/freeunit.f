c
c
c
c
c
c
      subroutine FREEUNIT(iunit)
      implicit none
c
c     Purpose: Get a free FORTRAN-unit
c
c     Output: iunit  a free (file-)unit number
c
c     Urs Kradolfer, 16.7.90
c
      integer iunit
      logical lopen
c
      do iunit=10,999
         if(iunit.eq.999) stop'FREEUNIT>>> no free unit found!'
         inquire(unit=iunit,opened=lopen)
         if(.not.lopen) RETURN
      enddo
c
      RETURN
c
      end ! of subr. freeunit
