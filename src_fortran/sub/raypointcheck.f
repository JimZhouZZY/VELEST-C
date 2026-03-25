c
c
c
c
c
c
      subroutine RAYPOINTCHECK(rp,nrp,staname)
c
c     checks, whether all the raypoints with z<0 are below the actual
c     surface of the topography (of Switzerland).
c
      implicit none
      include '../inc/vel_com.inc'
c
      integer nrp
      real rp(3,inrpmax)
      character*4 staname
c
      integer j
      real zzz,dzzz
      do j=1,nrp
         if(rp(3,j).lt.0.0)then
            call CHTOP(-rp(1,j),rp(2,j),zzz,topo1file,topo2file)
            dzzz=rp(3,j)-zzz
            if(dzzz.lt.0.0)then
               write(6,'(1x,''ray in the air... ! rp3='',f6.3,
     &         '' ZZ='',f6.3,'' dz='',f6.3,'' rp# ='',i2,'' nrp='',i2,
     &         '' STN='',a4,''i '',i4)')
     &         rp(3,j),zzz,dzzz,j,nrp,staname,isingle
            endif
         endif
      enddo
c
      RETURN
      end ! of subr. raypointcheck
