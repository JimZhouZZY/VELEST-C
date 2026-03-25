c
c
c
c
c
c
      subroutine CPUTIMER(cpusec)   ! Urs Kradolfer, 16.9.91
      implicit none
      real cpusec, icpu
      
      call cpu_time(icpu)
      cpusec=icpu
      
      end ! of subroutine cputimer
