c
c
c
c
c
c
      subroutine OPENERROR(subr,char)
c
      implicit none
      include '../inc/vel_com.inc'
c
      character*(*) subr,char
      write(6,*)'WARNING:'
      write(6,*)'SUBROUTINE :',subr,'    ERROR OPENING FILE: ',char
      if(.not.single_turbo)then
         write(16,*)'WARNING:'
         write(16,*)'SUBROUTINE :',subr,'    ERROR OPENING FILE: ',char
      endif
      stop'Openerror; program VELEST stopped !'
      end ! of subr. openerror
