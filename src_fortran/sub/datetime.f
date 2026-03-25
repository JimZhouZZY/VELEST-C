c
c
c
c
c
c
      subroutine DATETIME(dattim)       ! Urs Kradolfer, 16.9.91
c
c     return UNIX-date and time
c
c     runs on: SunOS , BSD UNIX
c
c     Call:    call DATETIME(dattim) 
c
c     returned value: dattim  is a character*20
      implicit none
      character datum*26, dattim*20
      character ctime
      integer(KIND=2) it
      
c external changed to intrinsic problem was solved on goo.gl/e1PcnK 

      intrinsic ctime !$pragma C(ctime)
      
      it=time8()
      datum=ctime(it)
      dattim=datum(5:24)
      end ! of subroutine datetime
