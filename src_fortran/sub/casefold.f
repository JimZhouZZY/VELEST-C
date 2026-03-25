c
c
c
c
c
c
      subroutine CASEFOLD(cn) ! Urs Kradolfer, 20. Feb. 1990
c
c Input : character-string of any length
c Output: same character-string, but all letters are CAPITAL now
c         (other characters are not changed)
c
c   --> this subr. needs the function TRIMLEN !
c
c Call  : call CASEFOLD(charstring)
c
      implicit none
      character cn*(*)
      integer i, trimlen,ilen, ival
      ilen=trimlen(cn)
      do i=1,ilen
         ival=ICHAR(cn(i:i))
         if(ival.ge.97.and.ival.le.122)then
            cn(i:i)=CHAR(ival-32)
         endif
      enddo
      RETURN
      end ! of subr. casefold
