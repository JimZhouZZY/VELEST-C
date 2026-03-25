c
c
c
c
c
c
      subroutine CHTOP(xx,yy,zk,
     &                 topo1file,topo2file)   !  Dummy routine
cEK
c     this routine for a given pair of coordinates (xx,yy) provides
c     the elevation (zk in km) by reading it off from topo-arrays
cEK
c     If necessary, replace it with appropriate routine and 
c     topo information
c
      implicit none
      real xx,yy,zk
      character*(*) topo1file,topo2file
c
c
      zk=0.0
      r e t u r n
      end
