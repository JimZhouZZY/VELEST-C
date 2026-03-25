c
c
c
c
c
c
      subroutine LAYERHIT(rp,nrpdeep,nl,nrp,mll)
c
c     counts the number of hits in each layer.
c
      implicit none
      include '../inc/vel_com.inc'
c
      integer nrpdeep,nl,nrp,mll
      real rp(3,inrpmax)
      real rpmax,avraydepth
      integer j,jlay
c
c
c     find deepest layer hit by this ray (= refraction layer or eventlayer
c                                           for direct wave)
      rpmax=-999.
      do j=1,nrp
         if(rp(3,j).gt.rpmax)then
            rpmax=rp(3,j)  ! = max. depth of ray
            nrpdeep=j   ! index of deepest raypoint
         endif
      enddo
      rpmax=rpmax+0.000001   ! avoid raypoint on layer-boundary
      lmax=1
      do j=2,nl
         if(h(j).gt.rpmax)then   !  h(j)= top of layer j [km]
            lmax=j-1  ! lmax = layer number for depth rpmax
            goto 4000
         endif
      enddo
      lmax=nl  ! lowest layer
 4000 continue
      if(rp(3,nrpdeep).eq.rp(3,nrpdeep+1))then   ! ray is horizontal
         irefrlayer(lmax)=irefrlayer(lmax)+1       ! headwave in layer LMAX
         refraylen(lmax)=refraylen(lmax)+
     &   SQRT( (rp(1,nrpdeep)-rp(1,nrpdeep+1))**2 +
     &         (rp(2,nrpdeep)-rp(2,nrpdeep+1))**2 )
      else
         if(mll.eq.0)then
            noheadwave=noheadwave+1
         else
            irefllayer(lmax-1)=irefllayer(lmax-1)+1
c                                           ! -1 because of rpmax=rpmax+0.000001
         endif
      endif
      avhraylen=avhraylen+
     & SQRT( (rp(1,1)-rp(1,nrp))**2 +
     &       (rp(2,1)-rp(2,nrp))**2 )
      avvraylen=avvraylen+
     & ABS( rp(3,nrpdeep)-rp(3,nrp) )
c
      do j=2,nrp
         avraydepth=( rp(3,j)+rp(3,j-1) ) / 2.
         avraydepth=avraydepth+0.000001    ! avoid raydepth on layer-boundary
         do jlay=1,nl
            if(avraydepth.ge.h(jlay).and.
     &         avraydepth.lt.(h(jlay)+thk(jlay)) )then
               hitlay(jlay,1)=hitlay(jlay,1)+1
               hitlay(jlay,2)=hitlay(jlay,2)+
     &                        SQRT( (rp(1,j)-rp(1,j-1))**2 +
     &                              (rp(2,j)-rp(2,j-1))**2 )
               hitlay(jlay,3)=hitlay(jlay,3)+
     &                        ABS( rp(3,j)-rp(3,j-1) )
            endif
         enddo
      enddo
      RETURN
      end ! of subr. layerhit
