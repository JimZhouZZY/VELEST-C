c
c
c
c
c
c
      subroutine RESETSTATIS
c
c     called by subr. TRAVELTIME. The statistics-variables are reset here.
c
      implicit none
      include '../inc/vel_com.inc'
      integer ihitl,jhitl
c
      if(irflout.eq.1) rewind(77)
      if(irfrout.eq.1) rewind(78)
      if(iresout.eq.1) rewind(79)
c
      do ihitl=1,nltot
         irefllayer(ihitl)=0
         irefrlayer(ihitl)=0
         refraylen(ihitl)=0.0
         hitlay(ihitl,1)=0.0
         hitlay(ihitl,2)=0.0
         hitlay(ihitl,3)=0.0
      enddo
      noheadwave=0     ! nr of straight & direct waves
      avhraylen=0.0    ! average horiz. raylength
      avvraylen=0.0    !    "    vertic.    "
      sterr=0.0        ! raytracer-errors
      direrr=0.0       !    "
      refrerr=0.0      !    "
      reflerr=0.0      !    "
      avrefrres=0.0    ! average residual for this ray-type
      avotheres=0.0    !    "
      avreflres=0.0    !    "
      abrefrres=0.0    ! average absolute residual for this ray type
      abotheres=0.0    !    "
      abreflres=0.0    !    "
      nrrefrres=0.0    ! nr of residuals for this ray-type
      nrotheres=0.0    !    "
      nrreflres=0.0    !    "
      do ihitl=1,nsta
         do jhitl=1,8
            stnazires(ihitl,jhitl)=0.0
         enddo
      enddo
c
      RETURN
      end ! of subr. resetstatis
