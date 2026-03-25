c
c
c
c
c
c
      subroutine ufwabo(isecpendel,seismkonst,seismdamp,voltgain,
     &                  devampl,period,iampltype,ampl,ifilter)
c---- subr. UFWABO ---------------------------------------------------
c
c     Berechnet die komplexe Uebertragungsfunktion im Frequenzbereich
c     fuer das System 'Bodenbewegung --> Devco-Auswertetisch' bzw.
c     fuer das System 'Bodenbewegung --> ADC (Mini)'
c
c     Options: iampltype = 1  --> [ampl] = mm(Wood-Anderson)/mm(Devco)
c              iampltype = 2  --> [ampl] = Mikron(Bodenbewegung)/mm(Devco)
c
c              ifilter   = 'DE' --> devampl von Develocorder
c              ifilter   = 'AD' --> devampl von ADC
c
c     Zur Subr. UFWABO gehoeren ebenfalls die Complex-Functions:
c     HP1, STS373, DEV und B6.
c--------------------------------------------------------------------------
      implicit none
      integer isecpendel,iampltype
      real seismkonst,seismdamp,voltgain,devampl
      real period,ampl
c
      real zpi,t02,t01,ts,s,hs,om,om2,omg2,tau2
      real om0,gain,twa,o0w,hw,vwa,amplresp,omg1,tau1
      integer iline,i
      character*(*) ifilter
      complex g1,g2,G,j,B6,HP1,STS373,DEV,jom
      complex gseismb
      data zpi/6.28319/j/(0.,1.)/
c
c---- Spezifikation der Hochpaesse auf der Gesamtuebertragung (PREampl. -->
c---- ADC) wie mit Programm Calap berechnet (bzw. verifiziert)
c
      iline=isecpendel
      if(iline.eq.1)then
         T02=11.9  ! 1-sec-Linie
         T01= 3.7  !     ''
      endif
      if(iline.eq.2)then
         T02= 12.7  ! 2-sec-Linie
         T01=112.   !     ''
      endif
      TS=isecpendel
      S=seismkonst
      hs=seismdamp
c
c---- Berechnung der Uebertragungsfunktion {Bodenbewegung --> ADC} ---
c
      om=zpi/period
      jom=j*om
      om2=om**2
c
c--------------------------- HP 2. Ordnung ----------------------------
c--------------------------- zus.gesetzt aus 2 HP 1. Ordnung ----------
      omg2=zpi*1./T02
      tau2=1./omg2
c
      g2=jom*tau2
      g2=(g2/(1.+g2))**2
c--------------------------- HP 1. Ordnung ----------------------------
      omg1=zpi*1./T01
      tau1=1./omg1
c
      g1=jom*tau1
      g1=g1/(1.+g1)
c----- Seismometer ----------------------------------------------------
100   om0=zpi/TS
      ! fuer Bodenbewegung ------------------------------------------
      gseismb=om2*jom*S/(om0**2+2*jom*om0*hs-om2)    ! [V/cm]
c----------------------------------------------------------------------
c
c---- Die Uebertragungsfunktionen wurden urspruenglich via ADC berechnet
c     (WK Kradi ~1983) und demzufolge ist der anti-alias-Filter schon
c     in der Uebertragungsfunktion mit inbegriffen. Um die Uebertragungs-
c     funktionen nur bis zum Diskriminator zu berechnen, muss dieser
c     Filter rueckgaengig gemacht werden:
c---- Inverse Uebertragung des DC-Unterdrueckungs-HP vor ADC
c
      G=g1*g2*gseismb/HP1(om)
c---- (damit sind wir beim Diskriminator)
c
      if(ifilter.eq.'DE')then
c
c        Uebertragung {Diskriminator --> Vergroesserungsschirm}
c
         G=G*STS373(om)*DEV(om)
                     ! Bis jetzt [mm Devco/mm Bodenbew.]
c
      endif
c
      if(ifilter.eq.'AD')then
c
c        Uebertragung {Diskriminator --> ADC }
c
         G=G*HP1(om)*B6(om)
                       ! Bis jetzt [volt/cm Bodenbewegung]
         gain=16.      ! [plot mm / volt]   (Ablesung ist so skaliert)
         G=G/10.       ! [volt/mm Bodenbewegung]
         G=G*gain      ! [mm plot/mm Bodenbewegung]
c
      endif
c
c---- Umrechnung in Mikron/Devco-Millimeter  bzw.
c---- Deconvolution mit der Uebertragungsfunktion eines Wood-Anderson-Seism.
c
      i=iampltype        !!! Ablesung = devco oder plot
      if(i.eq.2)then
         G=G/1000.     ! jetzt: mm Ablesung/um Bodenbewegung
         goto 999
      endif
      if(i.eq.1)then
         Twa=0.8
         o0w=zpi/Twa
         hw=0.78
         Vwa=2800
         G=G/(om2 *Vwa/(o0w**2+2*jom*o0w*hw-om2))
C                                             ! [1]     nur Wegvergr. !!
                     ! jetzt: mm Ablesung/mm Wood-Anderson !
         goto 999
      endif
c
999   amplresp=cabs(G)   ! Amplitudengang := Betrag der kompl. Transferfkt.
c
      ampl=devampl/(voltgain*amplresp)     ! erst hier reziproke Einheiten !!
c
c fuer IAMPLTYPE =          1          |           2
c dann [ampl]    = [mm Wood-Anderson]  |  [um Bodenbewegung]
c
      RETURN
      end
c
c
c------------------- B6: Bessel-Filter (TP) ------------------------------
      complex function B6(om)
      implicit none
      complex j,gbn1,gbn2,gbn3
      real om,zpi,fg1,omg1
c
      data j/(0.,1.)/,zpi/6.28319/
      fg1=12.     ! Grenzfrequenz 12 Hz
      omg1=om/(zpi*fg1)
c
      gbn1= 1. + 1.2217*j*omg1 - 0.3887*omg1*omg1
      gbn2= 1. + 0.9686*j*omg1 - 0.3505*omg1*omg1
      gbn3= 1. + 0.5131*j*omg1 - 0.2756*omg1*omg1
c
      b6=1./(gbn1*gbn2*gbn3)
c
      end
c
c------------------- HP 0.1 Hz 1. Ordnung (vor ADC) [DC-Unterdrueckung]-
      complex function HP1(om)
      implicit none
      real om
      complex j
      real zpi,fg,tau
      data j/(0.,1.)/zpi/6.28319/
      Fg=0.1 ! [Hz]     Grenzfrequenz
      tau=1./(zpi*Fg)
c
      HP1=j*om*tau
      HP1=HP1/(1.+HP1)
c
      return
      end
c
c------------------- TP STS-373 (Bessel 4. Ordnung, 5 Hz ) ----------
      complex function STS373(om)
      implicit none
      real om
c      implicit InTeGeR*4 (i-n)
      complex    j,gbn1,gbn2, o
      real zpi,fg1,omg1
      data j/(0.,1.)/zpi/6.28319/
      fg1=5.
      omg1=zpi*fg1
      o=om/omg1
c
      gbn1=1.+1.3397*j*o-0.4889*o**2
      gbn2=1.+0.7743*j*o-0.3890*o**2
c
      STS373=1/(gbn1*gbn2)
c
      return
      end
c
c------------------- Develocorder & 12-fache Schirmvergroesserung -----
      complex function DEV(om)
      implicit none
      real om
      complex j
      real zpi,fg,tau
      data j/(0.,1.)/zpi/6.28319/
                                 ! 1.6 cm/V MIT 12-f. Vergr. am Tisch !!
      Fg=15. ! [Hz]              ! Galvanometer      krit. Daempfung
      tau=1./(zpi*Fg)            ! (TP)              2. Ordnung
c
      DEV=1.6*( 1./(1.+j*om*tau) )**2           ! [cm/V]
c
c!!!!!!!! bei neuen functions diese in Haupt-Subroutine deklarieren  !!!
      return
      end
