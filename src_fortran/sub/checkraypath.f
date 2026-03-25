c
c
c
c
c
c
      subroutine CHECKRAYPATH(rp,nrp)
c
c
      implicit none
      include '../inc/vel_com.inc'
c
c     do not allow two identical raypoints !
c     (first two raypoints may be identical; e.g. if hypocenter is on
c      layer-boundary)
c
      real rp(3,inrpmax)
      integer nrp
c
      if(rp(1,1).eq.rp(1,2).and.
     &   rp(2,1).eq.rp(2,2).and.
     &   rp(3,1).eq.rp(3,2))then
chrm         do j=1,nrp-1
chrm            rp(1,j)=rp(1,j+1)          ! delete first ray-point
chrm            rp(2,j)=rp(2,j+1)          ! and move elements in array RP 'down'
chrm            rp(3,j)=rp(3,j+1)
chrm         enddo
chrm         nrp=nrp-1                     ! nr_of_raypoints is now smaller by 1 !!
chrm     move hypocenter 10 cm away from the layer boundary (upwards)
         rp(3,1) = rp(3,1) - 0.0001        	                                
         RETURN
      endif
c
      RETURN
      end ! of subr. checkraypath
