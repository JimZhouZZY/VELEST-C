c
c
c
c
c
c
      subroutine DETNOFUNKNOWNS
c
      implicit none
      include '../inc/vel_com.inc'
      integer lipeff
c
c     determine the number of unknowns NVAR and
c     calculate the number of equations LIP to be solved :
c
      if(icount.ne.0)then
         lip=4*neqs + nshot +1 ! do NOT invert for velocity parameters
         lipeff=lip
      else
         lip=4*neqs + nshot + nltot + ksta +1 ! invert for velocity parameters
         lipeff=4*neqs + nshot + nltot + nstaeff +1 ! invert for vel. parameters
         if(nsinv.eq.0) lip=4*neqs + nshot + nltot + 1 !do NOT inv. for sta-corr
      endif
c     NVAR is the number of unknowns to solve for :
      nvar=lip-1
c     NVAREFF is the REAL number of unknowns to solve for!!!
      nvareff=lipeff-1
c     KVAR is the number of elements on or above the main diagonal of G=(At(A))
      kvar=nvar*lip/2              !  cuk !BELOW!! cuk
c
c
      if(.not.single_turbo)then
         write(16,*)
         write(16,'('' Number of unknowns (for array-indexing): '',
     &              '' nvar = '',i4)') nvar
         write(16,'('' Number of effective unknowns        : '',
     &              '' nvareff = '',i4)') nvareff
         write(16,'('' Number of elements on/below main diagonal '',
     &              ''of matrix G = At*A : kvar = '',i7)') kvar
         write(16,*)
      endif
c
      return
      end ! of subr. detnofunknowns
