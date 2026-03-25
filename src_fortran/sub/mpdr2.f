c
c
c
c
c
c
      subroutine MPDR2(isecpendel,seismkonst,seismdamp,
     &               voltgain,pdrampl,period,iampltyp,ampl)
c
c------------------------------------ Kaspar G. Renggli 1984 ---------------
C
c     bodenamliptude fuer pdr2-stationen des nw-ch-netzes
c     parameter wie bei subr. muk, ausser pdrampl in [1/10 um/sec]
c---------------------------------------------------------------------------
c  
      implicit none
      integer isecpendel,iampltyp
      real seismkonst,seismdamp,voltgain,pdrampl,period,ampl
      complex G,gseis,ghp,glp2,glp2a,glp3,j,jom
      real zpi,ts,s,hs,om,om2,omg1,omg2,omg3,omg4,omg5,twa
      real omw,hw,vwa
      data zpi/6.28319/ j/(0.,1.)/
      TS=isecpendel
      S=seismkonst
      hs=seismdamp
      om=zpi/period
      jom=j*om
      om2=om**2
c-----Uebertragungsfunktionen fuer Seismometer und Filter des PDR-2-Systems---
c
c------------ Seismometer als Wegaufnehmer (nur Mechanik) -------------
      omg1=zpi/TS
      gseis=jom*om2/(omg1**2+2*jom*omg1*hs-om2)
c------------- 2 Hochpaesse 1. Ordnung  =>  2. Ordnung, fc=0.3 Hz -----
      omg2=zpi*0.3
      ghp=1/(1-2*j*omg2/om-omg2**2/om2)
c------------- Butterworth Tiefpass 2. Ordnung, fc=30.Hz --------------
      omg3=zpi*30.
      glp2=1/(1+1.4142*j*om/omg3-om2/omg3**2)
c----------- Butterworth Tiefpass 3.Ordnung, fc=25 Hz (Discriminator)----
      omg4=zpi*25.
      glp3=1/(1+2*j/omg4-2*om2/omg4**2-jom*om2/omg4**3)
c----------- Butterworth Tiefpass 2.Ordnung, fc=24 Hz (Antialiasing) ----
      omg5=zpi*24.
      glp2a=1/(1+1.4142*jom/omg5-om2/omg5**2)
c----------- Uebertragungsfunktion des gesamten Systems -----------------
      G=gseis*ghp*glp2*glp3*glp2a
c---- Deconvolution mit der Uebertragungsfunktion eines Wood-Anderson-Seism.
      Twa=0.8
      omw=zpi/Twa
      hw=0.78
      Vwa=2800
      G=(om2 *Vwa/(omw**2+2*jom*omw*hw-om2))/G
      ampl=cabs(G) ! [sec] Amplitudengang = Betrag der kompl. Transferfkt.
      ampl=ampl*pdrampl/10000.  ! [sec*1/10um/sec*1/1000=mm]
      return
      end
