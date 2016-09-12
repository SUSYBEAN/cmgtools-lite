CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC Figure out where bin boundaries are for one histogram

      SUBROUTINE getBins(histo)

      IMPLICIT NONE

      CHARACTER*128 binfile
      INTEGER histo, i, binStat

      INCLUDE 'histos.f'

CCCCC file 17 should be open in init_hist
      read(17,*) dummystring2, binfile, LOWEDGE(histo),
     &    HIGHEDGE(histo), ISHISTON(histo)
      IF (ichar(binfile(1:1)).ge.48.and.ichar(binfile(1:1)).le.57) THEN ! # of bins is a number
         read(binfile, fmt='(I3)') numBins(histo) ! convert to number
c         eqlSpace(histo) = .true.
         IF (numBins(histo).gt.BINPHIST) THEN ! too many, we'd segfault
            write(*,*) "Warning: too many bins, setting to maximum.  "
     &               //"To change maximum, see histos.f and recompile."
            numBins(histo) = BINPHIST
         ENDIF
         DO i=0,numBins(histo)
            binBounds(i,histo) = lowEdge(histo) + 
     &         (highEdge(histo)-lowEdge(histo))*i/numBins(histo)
         ENDDO
      ELSE 
c         eqlSpace(histo) = .false.
         numBins(histo) = 0
         OPEN(unit=18, file=binfile, status='old')
         i = 0
         binStat = 1
         DO WHILE (binStat.ge.0.and.numBins(histo).le.binphist)
            READ(18, *, iostat=binStat) binBounds(i,histo)
            i = i + 1
         ENDDO
         numBins(histo) = i - 2
         lowEdge(histo) = binBounds(0,histo)
         highEdge(histo) = binBounds(numBins(histo),histo)
      ENDIF

      RETURN
      END      


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC Initialize histograms

      subroutine init_hist(histofile)

      CHARACTER*128 histofile
      INTEGER COUNTER1

      include 'histos.f'

      open(unit=17,file=histofile,status='old')
      read(17,*)
      read(17,*)

CCC   Get data for all histograms from histogram input file

      HISTNAME(H_ZPT) = 'Z/W pT'
      HISTNAME(H_ZRAP) = 'Z/W rapidity'
      HISTNAME(H_INV_MASS) = 'Q_ll Invariant Mass'
      HISTNAME(H_LPT) = 'l-/lep. pT' 
      HISTNAME(H_LETA) = 'l-/lep. eta'
      HISTNAME(H_ALPT) = 'l+/neu. pT'
      HISTNAME(H_ALETA) = 'l+/neu. eta'
      HISTNAME(H_J1PT) = 'jet 1 pT'
      HISTNAME(H_J1ETA) = 'jet 1 eta'
      HISTNAME(H_J2PT) = 'jet 2 pT'
      HISTNAME(H_J2ETA) = 'jet 2 eta'
      HISTNAME(H_PHTPT) = 'phot. pT'
      HISTNAME(H_PHTETA) = 'phot. eta'
      HISTNAME(H_BEAM_THRUST) = 'beam thrust'
      HISTNAME(H_DRLL) = 'dR_lep,lep'
      HISTNAME(H_DRJ1L) = 'dR_jet1,l-/lep.'
      HISTNAME(H_DRJ1AL) = 'dR_jet1,l+/neu.'
      HISTNAME(H_DRJ2L) = 'dR_jet2,l-/lep.'
      HISTNAME(H_DRJ2AL) = 'dR_jet2,l+/neu.'
      HISTNAME(H_DRJJ) = 'dR_jet,jet'
      HISTNAME(H_DRPHTL) = 'dR_phot,lep'
      HISTNAME(H_H_T) = 'H_T'
      HISTNAME(H_M_T) = 'M_T'
      do COUNTER1=1,H_M_T ! all "normal" histograms
         CALL getBins(COUNTER1)
      enddo

CCC   Forward-backward Asymmetry

      HISTNAME(H_AFB) = 'A_FB numer'
      CALL getBins(H_AFB)
      HISTNAME(H_AFB_DEN) = 'A_FB denom'
      NUMBINS(H_AFB_DEN) = NUMBINS(H_AFB)
      LOWEDGE(H_AFB_DEN) = LOWEDGE(H_AFB)
      HIGHEDGE(H_AFB_DEN) = HIGHEDGE(H_AFB)
      ISHISTON(H_AFB_DEN) = ISHISTON(H_AFB)

      DO i=0,numBins(H_AFB)
         binBounds(i,H_AFB_DEN) = binBounds(i,H_AFB)
      ENDDO

CCC   COLLINS-SOPER ANGLES AND MOMENTS
      read(17,*)
      read(17,*)

      HISTNAME(H_A_0) = 'A_0 numer'
      CALL getBins(H_A_0)

      HISTNAME(H_A_1) = 'A_1 numer'
      NUMBINS(H_A_1) = NUMBINS(H_A_0)
      LOWEDGE(H_A_1) = LOWEDGE(H_A_0)
      HIGHEDGE(H_A_1) = HIGHEDGE(H_A_0)
      ISHISTON(H_A_1) = ISHISTON(H_A_0)

      HISTNAME(H_A_2) = 'A_2 numer'
      NUMBINS(H_A_2) = NUMBINS(H_A_0)
      LOWEDGE(H_A_2) = LOWEDGE(H_A_0)
      HIGHEDGE(H_A_2) = HIGHEDGE(H_A_0)
      ISHISTON(H_A_2) = ISHISTON(H_A_0)

      HISTNAME(H_A_3) = 'A_3 numer'
      NUMBINS(H_A_3) = NUMBINS(H_A_0)
      LOWEDGE(H_A_3) = LOWEDGE(H_A_0)
      HIGHEDGE(H_A_3) = HIGHEDGE(H_A_0)
      ISHISTON(H_A_3) = ISHISTON(H_A_0)

      HISTNAME(H_A_4) = 'A_4 numer'
      NUMBINS(H_A_4) = NUMBINS(H_A_0)
      LOWEDGE(H_A_4) = LOWEDGE(H_A_0)
      HIGHEDGE(H_A_4) = HIGHEDGE(H_A_0)
      ISHISTON(H_A_4) = ISHISTON(H_A_0)

      HISTNAME(H_A_DEN) = 'A_i denom'
      NUMBINS(H_A_DEN) = NUMBINS(H_A_0)
      LOWEDGE(H_A_DEN) = LOWEDGE(H_A_0)
      HIGHEDGE(H_A_DEN) = HIGHEDGE(H_A_0)
      ISHISTON(H_A_DEN) = ISHISTON(H_A_0)

      DO i=0,numBins(H_A_0)
         binBounds(i,H_A_1) = binBounds(i,H_A_0)
         binBounds(i,H_A_2) = binBounds(i,H_A_0)
         binBounds(i,H_A_3) = binBounds(i,H_A_0)
         binBounds(i,H_A_4) = binBounds(i,H_A_0)
         binBounds(i,H_A_DEN) = binBounds(i,H_A_0)
      ENDDO

      HISTNAME(H_PHI_CS) = 'phi CS'
      CALL getBins(H_PHI_CS)

      HISTNAME(H_PHI_CS_NORM) = 'phi CS norm'
      NUMBINS(H_PHI_CS_NORM) = 1
      LOWEDGE(H_PHI_CS_NORM) = LOWEDGE(H_PHI_CS)
      HIGHEDGE(H_PHI_CS_NORM) = HIGHEDGE(H_PHI_CS)
      ISHISTON(H_PHI_CS_NORM) = ISHISTON(H_PHI_CS)

      HISTNAME(H_COSTH_CS) = 'costh CS'
      CALL getBins(H_COSTH_CS)

      HISTNAME(H_COSTH_CS_NORM) = 'cth CS norm'
      NUMBINS(H_COSTH_CS_NORM) = 1
      LOWEDGE(H_COSTH_CS_NORM) = LOWEDGE(H_COSTH_CS)
      HIGHEDGE(H_COSTH_CS_NORM) = HIGHEDGE(H_COSTH_CS)
      ISHISTON(H_COSTH_CS_NORM) = ISHISTON(H_COSTH_CS)

      HISTNAME(H_DPHI) = 'Delta_Phi'
      CALL getBins(H_DPHI)

CCC   If Histogram is added, need to change number of histograms
CCC   in ../include/histos.f

c      HISTNAME() = ''
c      CALL getBins(#)

CC    Smoothing parameters

      read(17,*)
      read(17,*)
      read(17,*) dummystring2, smoothing
      read(17,*) dummystring2, smcut
      read(17,*)
      read(17,*) dummystring2, histmethod

      close(17)

CCC   TOTBINS is not used
C      TOTBINS=0
C      COUNTER1=1
C      DO WHILE(COUNTER1.le.NUMHIST)
C         TOTBINS=TOTBINS+NUMBINS(COUNTER1)
C         COUNTER1=COUNTER1+1
C      ENDDO
C      COUNTER1=1

      call clearhist()

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC Fill histograms

      subroutine fill_hist(part,weight,pdfvec)
      implicit none
      include 'histos.f'
      include 'pdfcontrol.f'
      include 'cutcontrol.f'
      include 'particles.f'
    
      double precision weight
      double precision etadiff, phidiff
      double precision pi,pi2

      double precision z_pt,z_Y,invm_Z,QtauB

      logical isjet1ok,isjet2ok,isphtok

ccc   Collins-Soper related Variables   cccccccc
      double precision costh_CS,phi_CS
c      double precision plot_phi_CS
c      double precision tanphi_CS
      double precision tanphiCS_N,tanphiCS_D
      double precision newrand ! built from tmprand

c      integer ib,ibb
cccccccccccccccccccccccccccccccccccccccccccccccc

      double precision m2,scm
      integer ptype
      COMMON/params/ m2,scm
      COMMON/progmode/ ptype

CCCCC Array holding the result for each PDF set

      double precision pdfvec
      dimension pdfvec(NPDFerror)

      double precision fillValue(NHIST),wweight(NHIST)
      double precision remain,smearwgt,reweight
      INTEGER thisBin,cc,i,j

      PARAMETER(pi = 3.1415926536d0)
      PARAMETER(pi2 = 3.1415927d0)
CCC --------------------------------------------------------------
CCC ---- Compute variables for the histograms defined in init_hist      
CCC
CCC ---------------------- Definition of particle 4-momenta -----------------
CCC   Particle number:
CCC     1. Z-boson  (or W-boson)
CCC     2. Lepton  (or Charged Lepton)
CCC     3. Anti-lepton  (or Neutrino/Missing Energy)
CCC     4. Parton 1
CCC     5. Parton 2
CCC     6. Photon
CCC -------------------------------------------------------------------------

      isjet1ok = .false.
      isjet2ok = .false.
      isphtok = .false.

      invm_Z = dsqrt(part(1,1)**2-part(1,2)**2-
     &     part(1,3)**2-part(1,4)**2)
      z_pt = dsqrt(part(1,2)**2+part(1,3)**2)
      z_Y = 0.5d0*dlog(dabs((part(1,1)+part(1,4))
     & /(part(1,1)-part(1,4))))

CCC   Collins-Soper Variables

CCC   Initialize the lepton + anit-lepton 4-vector,
CCC   and Collins & Soper angles

      IF(isHistOn(H_A_0).gt.0.or.isHistOn(H_PHI_CS).gt.0.or.isHistOn(H_COSTH_CS).gt.0) THEN

c         do j=1,4
c            lepsum(j)=part(2,j)+part(3,j)
c         enddo

         if ((ptype.eq.2.or.ptype.eq.4)
     &      .or.((ptype.eq.1.or.ptype.eq.3)
     &       .and.z_Y.gt.zY_CScutoff)) then ! theta defined relative to beam1
            costh_CS = (
     &           2d0*(part(2,4)*part(3,1)-part(2,1)*part(3,4))/
     &           dsqrt(part(1,1)**2-part(1,2)**2-
     &           part(1,3)**2-part(1,4)**2)/
     &           dsqrt(part(1,1)**2-part(1,4)**2)
     &           )
         elseif ((ptype.eq.1.or.ptype.eq.3)
     &          .and.(z_Y.lt.(-zY_CScutoff))) then
CCC      theta goes with beam 2
            costh_CS = (
     &           2d0*(-part(2,4)*part(3,1)-part(2,1)*(-part(3,4)))/
     &           dsqrt(part(1,1)**2-part(1,2)**2-
     &           part(1,3)**2-part(1,4)**2)/
     &           dsqrt(part(1,1)**2-part(1,4)**2)
     &           )
         else
            costh_CS = 999d0 ! undefined/cut, doesn't get binned
         endif

CCC      Need both numerator and denominator signs for fortran datan2
         tanphiCS_N = dsqrt(part(1,1)**2-part(1,4)**2)*
     &        (-part(1,3)*(part(2,2)-part(3,2))+
     &        part(1,2)*(part(2,3)-part(3,3)))
         tanphiCS_D = dsqrt(part(1,1)**2-part(1,2)**2-
     &        part(1,3)**2-part(1,4)**2)*
     &        (part(1,2)*(part(2,2)-part(3,2))+
     &        part(1,3)*(part(2,3)-part(3,3)))

         if ((ptype.eq.2.or.ptype.eq.4).or.
     &       ((ptype.eq.1.or.ptype.eq.3).and.
     &        (dabs(z_Y).gt.zY_CScutoff))) then ! .and.
c     &        (z_pt.gt.1d-8)) then
CC      First digit of tmprand used in hard/soft coin flip, rescale
            newrand = dabs(2*tmprand-1)
            if (z_pt.lt.1d-1) then
               phi_CS = pi*(2*newrand-1)
            else
               phi_CS = sign(datan2(tanphiCS_N,
     &              tanphiCS_D), newrand-0.5d0)
            endif
         else
            phi_CS = 999d0
         endif

      ENDIF

CCC   Default reweight for all histogram bins

      reweight = 1d0 ! it'll be recycled after wweight(cc) assigned

CCC   Z pT

      fillValue(H_ZPT) = z_pt
      wweight(H_ZPT) = reweight

CCC   Z rapidity

      fillValue(H_ZRAP) = z_Y
      wweight(H_ZRAP) = reweight

CCC   lepsum invariant mass

      fillValue(H_INV_MASS) = invm_Z
      wweight(H_INV_MASS) = reweight

CCC   l- pT

      fillValue(H_LPT) = dsqrt(part(2,2)**2+part(2,3)**2)
      wweight(H_LPT) = reweight

CCC   l- rapidity

      fillValue(H_LETA) = 0.5d0*dlog(dabs((part(2,1)+part(2,4))
     & /(part(2,1)-part(2,4))))
      wweight(H_LETA) = reweight

CCC   l+ pT

      fillValue(H_ALPT) = dsqrt(part(3,2)**2+part(3,3)**2)
      wweight(H_ALPT) = reweight

CCC   l+ rapidity

      fillValue(H_ALETA) = 0.5d0*dlog(dabs((part(3,1)+part(3,4))
     & /(part(3,1)-part(3,4))))
      wweight(H_ALETA) = reweight

CCC   parton 1 pT

      fillValue(H_J1PT) = dsqrt(part(4,2)**2+
     &     part(4,3)**2)
      if (fillValue(H_J1PT).ge.minptj4obs) isjet1ok = .true.
      if (isjet1ok) then
         wweight(H_J1PT) = reweight
      else
         wweight(H_J1PT) = 0d0
      endif

CCC   parton 1 rapidity

      fillValue(H_J1ETA) = 0.5d0*dlog(dabs((part(4,1)+part(4,4))
     & /(part(4,1)-part(4,4))))
      if (isjet1ok) then
         wweight(H_J1ETA) = reweight
      else
         wweight(H_J1ETA) = 0d0
      endif

CCC   parton 2 pT

      fillValue(H_J2PT) = dsqrt(part(5,2)**2+
     &     part(5,3)**2)
      if (fillValue(H_J2PT).ge.minptj4obs) isjet2ok = .true.
      if (isjet2ok) then
         wweight(H_J2PT) = reweight
      else
         wweight(H_J2PT) = 0d0
      endif

CCC   parton 2 rapidity

      fillValue(H_J2ETA) = 0.5d0*dlog(dabs((part(5,1)+part(5,4))
     & /(part(5,1)-part(5,4))))
      if (isjet2ok) then
         wweight(H_J2ETA) = reweight
      else
         wweight(H_J2ETA) = 0d0
      endif

CCC   photon pT

      fillValue(H_PHTPT) = dsqrt(part(6,2)**2+
     &     part(6,3)**2)
      if (fillValue(H_PHTPT).ge.minptr4obs) isphtok = .true.
      if (isphtok) then
         wweight(H_PHTPT) = reweight
      else
         wweight(H_PHTPT) = 0d0
      endif

CCC   photon rapidity

      fillValue(H_PHTETA) = 0.5d0*dlog(dabs((part(6,1)+part(6,4))
     & /(part(6,1)-part(6,4))))
      if (isphtok) then
         wweight(H_PHTETA) = reweight
      else
         wweight(H_PHTETA) = 0d0
      endif

CCC   Beam thrust, for the case when the final-state 
CCC   rapidity is reconstructed.  Follows 1005.4060
CCC   and computes Q*tau_B

      IF(isHistOn(H_BEAM_THRUST).gt.0) THEN

cc    zrap already computed above
c      zrap = 0.5d0*dlog(dabs((part(1,1)+part(1,4))
c     & /(part(1,1)-part(1,4))))
      
      QtauB = 0d0
      if (dabs(part(4,1)).gt.dabs(part(4,4))) then
cc    parton pT/rapidities already computed above
c         etapart1 = 0.5d0*dlog(dabs((part(4,1)+part(4,4))
c     &     /(part(4,1)-part(4,4))))
c         ptpart1 = dsqrt(part(4,2)**2+
c     &     part(4,3)**2) 
         QtauB = QtauB+fillValue(H_J1PT)*dexp(-dabs(fillValue(H_J1ETA)-fillValue(H_ZRAP)))
      endif
      
      if (dabs(part(5,1)).gt.dabs(part(5,4))) then
c         etapart2 = 0.5d0*dlog(dabs((part(5,1)+part(5,4))
c     &     /(part(5,1)-part(5,4))))
c         ptpart2 = dsqrt(part(5,2)**2+
c     &     part(5,3)**2)
         QtauB = QtauB+fillValue(H_J2PT)*dexp(-dabs(fillValue(H_J2ETA)-fillValue(H_ZRAP)))
      endif

      fillValue(H_BEAM_THRUST) = QtauB
      wweight(H_BEAM_THRUST) = reweight

      ENDIF

CCC   lepton & anti-lepton deltaR separation

      IF(isHistOn(H_DRLL).gt.0) THEN

      etadiff = 0.0d0
      phidiff = 0.0d0

c      etadiff = 0.5d0*dlog(dabs((part(3,1)+part(3,4))
c     & /(part(3,1)-part(3,4))))-
c     & 0.5d0*dlog(dabs((part(2,1)+part(2,4))
c     & /(part(2,1)-part(2,4))))

      etadiff = fillValue(H_ALETA) - fillValue(H_LETA)

      phidiff = dabs(datan2(part(3,3),part(3,2))-
     & datan2(part(2,3),part(2,2)))
      if (phidiff.gt.pi) then
         phidiff = 2d0*pi-phidiff
      endif

      fillValue(H_DRLL) = dsqrt(etadiff**2+phidiff**2)
      wweight(H_DRLL) = reweight

      ENDIF

CCC   lepton & jet1 deltaR separation

      IF(isHistOn(H_DRJ1L).gt.0) THEN

      etadiff = 0.0d0
      phidiff = 0.0d0

c      etadiff = 0.5d0*dlog(dabs((part(2,1)+part(2,4))
c     & /(part(2,1)-part(2,4))))-
c     & 0.5d0*dlog(dabs((part(4,1)+part(4,4))
c     & /(part(4,1)-part(4,4))))

      etadiff = fillValue(H_J1ETA)-fillValue(H_LETA)

      phidiff = dabs(datan2(part(2,3),part(2,2))-
     & datan2(part(4,3),part(4,2)))
      if (phidiff.gt.pi) then
         phidiff = 2d0*pi-phidiff
      endif

      fillValue(H_DRJ1L) = dsqrt(etadiff**2+phidiff**2)
      if (isjet1ok) then
         wweight(H_DRJ1L) = reweight
      else
         wweight(H_DRJ1L) = 0d0
      endif

      ENDIF

CCC   anti-lepton & jet1 deltaR separation

      IF(isHistOn(H_DRJ1AL).gt.0) THEN

      etadiff = 0.0d0
      phidiff = 0.0d0

c      etadiff = 0.5d0*dlog(dabs((part(3,1)+part(3,4))
c     & /(part(3,1)-part(3,4))))-
c     & 0.5d0*dlog(dabs((part(4,1)+part(4,4))
c     & /(part(4,1)-part(4,4))))

      etadiff = fillValue(H_J1ETA) - fillValue(H_ALETA)

      phidiff = dabs(datan2(part(3,3),part(3,2))-
     & datan2(part(4,3),part(4,2)))
      if (phidiff.gt.pi) then
         phidiff = 2d0*pi-phidiff
      endif

      fillValue(H_DRJ1AL) = dsqrt(etadiff**2+phidiff**2)
      if (isjet1ok) then
         wweight(H_DRJ1AL) = reweight
      else
         wweight(H_DRJ1AL) = 0d0
      endif

      ENDIF

CCC   lepton & jet2 deltaR separation

      IF(isHistOn(H_DRJ2L).gt.0) THEN

      etadiff = 0.0d0
      phidiff = 0.0d0

c      etadiff = 0.5d0*dlog(dabs((part(2,1)+part(2,4))
c     & /(part(2,1)-part(2,4))))-
c     & 0.5d0*dlog(dabs((part(5,1)+part(5,4))
c     & /(part(5,1)-part(5,4))))

      etadiff = fillValue(H_J2ETA) - fillValue(H_LETA)

      phidiff = dabs(datan2(part(2,3),part(2,2))-
     & datan2(part(5,3),part(5,2)))
      if (phidiff.gt.pi) then
         phidiff = 2d0*pi-phidiff
      endif

      fillValue(H_DRJ2L) = dsqrt(etadiff**2+phidiff**2)
      if (isjet2ok) then
         wweight(H_DRJ2L) = reweight
      else
         wweight(H_DRJ2L) = 0d0
      endif

      ENDIF

CCC   anti-lepton & jet2 deltaR separation

      IF(isHistOn(H_DRJ2AL).gt.0) THEN

      etadiff = 0.0d0
      phidiff = 0.0d0

c      etadiff = 0.5d0*dlog(dabs((part(3,1)+part(3,4))
c     & /(part(3,1)-part(3,4))))-
c     & 0.5d0*dlog(dabs((part(5,1)+part(5,4))
c     & /(part(5,1)-part(5,4))))

      etadiff = fillValue(H_J2ETA) - fillValue(H_ALETA)

      phidiff = dabs(datan2(part(3,3),part(3,2))-
     & datan2(part(5,3),part(5,2)))
      if (phidiff.gt.pi) then
         phidiff = 2d0*pi-phidiff
      endif

      fillValue(H_DRJ2AL) = dsqrt(etadiff**2+phidiff**2)
      if (isjet2ok) then
         wweight(H_DRJ2AL) = reweight
      else
         wweight(H_DRJ2AL) = 0d0
      endif

      ENDIF

CCC   jet-jet deltaR separation

      IF(isHistOn(H_DRJJ).gt.0) THEN

      etadiff = 0.0d0
      phidiff = 0.0d0

c      etadiff = 0.5d0*dlog(dabs((part(4,1)+part(4,4))
c     & /(part(4,1)-part(4,4))))-
c     & 0.5d0*dlog(dabs((part(5,1)+part(5,4))
c     & /(part(5,1)-part(5,4))))

      etadiff = fillValue(H_J2ETA) - fillValue(H_J1ETA)

      phidiff = dabs(datan2(part(4,3),part(4,2))-
     & datan2(part(5,3),part(5,2)))
      if (phidiff.gt.pi) then
         phidiff = 2d0*pi-phidiff
      endif

      fillValue(H_DRJJ) = dsqrt(etadiff**2+phidiff**2)
      if ((isjet1ok).and.(isjet2ok)) then
         wweight(H_DRJJ) = reweight
      else
         wweight(H_DRJJ) = 0d0
      endif

      ENDIF

CCC   lepton & photon deltaR separation

      IF(isHistOn(H_DRPHTL).gt.0) THEN

      etadiff = fillValue(H_PHTETA)-fillValue(H_LETA)
      phidiff = dabs(datan2(part(2,3),part(2,2))-
     & datan2(part(6,3),part(6,2)))
      if (phidiff.gt.pi) then
         phidiff = 2d0*pi-phidiff
      endif
c     tmporarily fill it with photon-lepton separation for now
      fillValue(H_DRPHTL) = dsqrt(etadiff**2+phidiff**2)

      etadiff = fillValue(H_PHTETA) - fillValue(H_ALETA)
      phidiff = dabs(datan2(part(3,3),part(3,2))-
     & datan2(part(6,3),part(6,2)))
      if (phidiff.gt.pi) then
         phidiff = 2d0*pi-phidiff
      endif
c     fill it the mini of photon-lepton and photon-antilepton separation
      fillValue(H_DRPHTL) = min(dsqrt(etadiff**2+phidiff**2),fillValue(H_DRPHTL))

      if (isphtok) then
         wweight(H_DRPHTL) = reweight
      else
         wweight(H_DRPHTL) = 0d0
      endif

      ENDIF

CCC   H_H_T

      IF(isHistOn(H_H_T).gt.0) THEN

      fillValue(H_H_T) = dsqrt(part(2,2)**2+part(2,3)**2)+
     &     dsqrt(part(3,2)**2+part(3,3)**2)+
     &     dsqrt(part(4,2)**2+part(4,3)**2)+
     &     dsqrt(part(5,2)**2+part(5,3)**2)+
     &     dsqrt(part(6,2)**2+part(6,3)**2)
      wweight(H_H_T) = reweight

      ENDIF

CCC   Transverse Mass

      IF(isHistOn(H_M_T).gt.0) THEN

      fillValue(H_M_T) = dsqrt(
     &     2d0 * dsqrt( (part(2,2)**2+part(2,3)**2) *
     &     (part(3,2)**2+part(3,3)**2) ) -
     &     2d0 * (part(2,2)*part(3,2) + part(2,3)*part(3,3))
     &     )
      wweight(H_M_T) = reweight

      ENDIF

CCC   Forward backward asymmetry vs invmass

      IF(isHistOn(H_AFB).gt.0) THEN

      fillValue(H_AFB) = fillValue(H_INV_MASS)
      fillValue(H_AFB_DEN) = fillValue(H_INV_MASS)

      if ((ptype.eq.1).or.(ptype.eq.3)) then
         if ((z_Y*(part(2,4)*part(3,1)-part(2,1)*part(3,4))).ge.0d0) then
            wweight(H_AFB) = reweight
         else
            wweight(H_AFB) = -reweight
         endif
      else
         if ((part(2,4)*part(3,1)-part(2,1)*part(3,4)).ge.0d0) then
            wweight(H_AFB) = reweight
         else
            wweight(H_AFB) = -reweight
         endif
      endif
      wweight(H_AFB_DEN) = reweight

      ENDIF

CCC   COLLINS-SOPER MOMENTS
CCC   These get filled unnormalized, can only normalize after sectors combined

      IF(isHistOn(H_A_0).gt.0) THEN

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   A_i vs Z pT @ invmass_CS1 +/- res_CS1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

ccc   H_A_0
      fillValue(H_A_0) = fillValue(H_ZPT)
ccc   H_A_1
      fillValue(H_A_1) = fillValue(H_A_0)
ccc   H_A_2
      fillValue(H_A_2) = fillValue(H_A_0)
ccc   H_A_3
      fillValue(H_A_3) = fillValue(H_A_0)
ccc   H_A_4
      fillValue(H_A_4) = fillValue(H_A_0)
ccc   A_i Denominator
      fillValue(H_A_DEN) = fillValue(H_A_0)

         if ((ptype.eq.2.or.ptype.eq.4).or.
     &       ((ptype.eq.1.or.ptype.eq.3).and.
     &        (dabs(z_Y).gt.zY_CScutoff))) then
ccc   H_A_0
         wweight(H_A_0) = reweight*(
     &        20d0/3d0*0.5d0*(1d0-3d0*costh_CS**2)+
     &        2d0/3d0)

         if (z_pt.gt.1d0) then
ccc   H_A_1
            wweight(H_A_1) = reweight*5d0*
     &           2d0*costh_CS*dsqrt(1d0-costh_CS**2)*
     &           dcos(phi_CS)
ccc   H_A_2
            wweight(H_A_2) = reweight*10d0*
     &           (1d0-costh_CS**2)*dcos(2d0*phi_CS)
ccc   H_A_3
            wweight(H_A_3) = reweight*4d0*dsqrt(1d0-costh_CS**2)*dcos(phi_CS)
         else
            wweight(H_A_1) = 0d0
            wweight(H_A_2) = 0d0
            wweight(H_A_3) = 0d0
         endif

ccc   H_A_4
         wweight(H_A_4) = reweight*4d0*costh_CS
ccc   A_i Denominator
         wweight(H_A_DEN) = reweight
      else
         wweight(H_A_0) = 0d0
         wweight(H_A_1) = 0d0
         wweight(H_A_2) = 0d0
         wweight(H_A_3) = 0d0
         wweight(H_A_4) = 0d0
         wweight(H_A_DEN) = 0d0
      endif

      ENDIF

CCC   Phi as measured in Collins-Soper Frame

      if ((ptype.eq.2.or.ptype.eq.4).or.
     &       ((ptype.eq.1.or.ptype.eq.3).and.
     &        (dabs(z_Y).gt.zY_CScutoff))) then
         fillValue(H_PHI_CS) = phi_CS
         wweight(H_PHI_CS) = reweight
      else
         fillValue(H_PHI_CS_NORM) = 999d0
         wweight(H_PHI_CS_NORM) = 0d0
      endif

      fillValue(H_PHI_CS_NORM) = fillValue(H_PHI_CS)
      wweight(H_PHI_CS_NORM) = wweight(H_PHI_CS)

CCC   Cos(theta)_CS as measured in Collins-Soper Frame

      if ((ptype.eq.2.or.ptype.eq.4).or.
     &       ((ptype.eq.1.or.ptype.eq.3).and.
     &        (dabs(z_Y).gt.zY_CScutoff))) then
         fillValue(H_COSTH_CS) = costh_CS
         wweight(H_COSTH_CS) = reweight
      else
         fillValue(H_COSTH_CS) = 999d0
         wweight(H_COSTH_CS) = 0d0
      endif

      fillValue(H_COSTH_CS_NORM) = fillValue(H_COSTH_CS)
      wweight(H_COSTH_CS_NORM) = wweight(H_COSTH_CS)

CCC   lepton & anti-lepton Delta Phi separation

      IF(isHistOn(H_DPHI).gt.0) THEN

      phidiff = 0.0d0

ccc   3.1415926536d0
      phidiff = dabs(datan2(part(3,3),part(3,2))-
     & datan2(part(2,3),part(2,2)))
      if (phidiff.gt.pi2) then
         phidiff = 2d0*pi2-phidiff
      endif

      fillValue(H_DPHI) = phidiff
      wweight(H_DPHI) = reweight
      ENDIF


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCC --------------------------------------------------------------
CCC   Find correct bin, fill "ENTRIES" for later use

      do cc=1,NHIST
      IF(ISHISTON(cc).gt.0) THEN
      IF((fillValue(cc).GE.LOWEDGE(cc)).AND.(fillValue(cc).
     &  LT.HIGHEDGE(cc))) THEN
c         IF(eqlSpace(cc)) THEN
c            whichBin=(fillValue(cc)-LOWEDGE(cc))/
c     &      (HIGHEDGE(cc)-LOWEDGE(cc))*NUMBINS(cc)
         
c            thisBin=int(whichBin)
c            remain = whichBin - thisBin ! should be between 0 and 1
c            thisBin = thisBin + 1
c         ELSE ! make this default method if no slower than eqlSpace method
            i = 1
            DO WHILE (i.le.numBins(cc) .and.
     &                fillValue(cc).ge.binBounds(i,cc))
               i = i + 1
            ENDDO
            thisBin = i
            remain = (fillValue(cc)-binBounds(thisBin-1,cc))/
     &               (binBounds(thisBin,cc)-binBounds(thisBin-1,cc))
c         ENDIF

         IF((remain.LT.smcut).AND.(smoothing.ne.0)
     &                       .AND.(thisBin.GT.1)) THEN
C        smoothing is on and we need to smooth to the left
           smearwgt = 5d-1*(1-(remain/smcut)**smoothing)**smoothing
           IF(mod(ISHISTON(cc),2).eq.1) THEN ! regular histogram on
           reweight = wweight(cc)*smearwgt
           ENTRIES(0,thisBin-1,cc) = ENTRIES(0,thisBin-1,cc)
     &      + reweight*weight
             do i=1,NPDFerror
               ENTRIES(i,thisBin-1,cc) = ENTRIES(i,thisBin-1,cc)
     &          + reweight*pdfvec(i)
             enddo
           reweight = wweight(cc)*(1-smearwgt)
           ENTRIES(0,thisBin,cc) = ENTRIES(0,thisBin,cc)
     &      + reweight*weight
             do i=1,NPDFerror
               ENTRIES(i,thisBin,cc) = ENTRIES(i,thisBin,cc)
     &          + reweight*pdfvec(i)
             enddo
           ENDIF
           IF((ISHISTON(cc).ge.2).and.(ISHISTON(cc).le.3)) THEN ! cumulants on (each cumulate all left bins)
           reweight = wweight(cc)*smearwgt
           CUMULANTS(0,thisBin-1,cc) = CUMULANTS(0,thisBin-1,cc)
     &      + reweight*weight
             do i=1,NPDFerror
               CUMULANTS(i,thisBin-1,cc) = CUMULANTS(i,thisBin-1,cc)
     &          + reweight*pdfvec(i)
             enddo
           reweight = wweight(cc) ! every bin to right gets full weight in cumulant
           do j=thisBin,NUMBINS(cc)
           CUMULANTS(0,j,cc) = CUMULANTS(0,j,cc)
     &      + reweight*weight
             do i=1,NPDFerror
               CUMULANTS(i,j,cc) = CUMULANTS(i,j,cc)
     &          + reweight*pdfvec(i)
             enddo
           enddo
           ELSEIF((ISHISTON(cc).ge.4).and.(ISHISTON(cc).le.5)) THEN ! reverse-cumulants on (each cumulate all right bins)
           reweight = wweight(cc)*(1-smearwgt)
           CUMULANTS(0,thisBin,cc) = CUMULANTS(0,thisBin,cc)
     &      + reweight*weight
             do i=1,NPDFerror
               CUMULANTS(i,thisBin,cc) = CUMULANTS(i,thisBin,cc)
     &          + reweight*pdfvec(i)
             enddo
           reweight = wweight(cc)
           do j=1,thisBin-1
           CUMULANTS(0,j,cc) = CUMULANTS(0,j,cc)
     &      + reweight*weight
             do i=1,NPDFerror
               CUMULANTS(i,j,cc) = CUMULANTS(i,j,cc)
     &          + reweight*pdfvec(i)
             enddo
           enddo
           ENDIF

         ELSEIF((remain.GT.(1-smcut)).AND.(smoothing.ne.0)
     &                              .AND.(thisBin.LT.NUMBINS(cc))) THEN
C        smoothing is on and we need to smooth to the right
           smearwgt = 5d-1*(1-((1-remain)/smcut)**smoothing)**smoothing
           IF(mod(ISHISTON(cc),2).eq.1) THEN ! regular histogram on
           reweight = wweight(cc)*(1-smearwgt)
           ENTRIES(0,thisBin,cc) = ENTRIES(0,thisBin,cc)
     &      + reweight*weight
             do i=1,NPDFerror
               ENTRIES(i,thisBin,cc) = ENTRIES(i,thisBin,cc)
     &          + reweight*pdfvec(i)
             enddo
           reweight = wweight(cc)*smearwgt
           ENTRIES(0,thisBin+1,cc) = ENTRIES(0,thisBin+1,cc)
     &      + reweight*weight
             do i=1,NPDFerror
               ENTRIES(i,thisBin+1,cc) = ENTRIES(i,thisBin+1,cc)
     &          + reweight*pdfvec(i)
             enddo
           ENDIF
           IF((ISHISTON(cc).ge.2).and.(ISHISTON(cc).le.3)) THEN ! cumulants on (each cumulate all left bins)
           reweight = wweight(cc)*(1-smearwgt)
           CUMULANTS(0,thisBin,cc) = CUMULANTS(0,thisBin,cc)
     &      + reweight*weight
             do i=1,NPDFerror
               CUMULANTS(i,thisBin,cc) = CUMULANTS(i,thisBin,cc)
     &          + reweight*pdfvec(i)
             enddo
           reweight = wweight(cc)
           do j=thisBin+1,NUMBINS(cc)
           CUMULANTS(0,j,cc) = CUMULANTS(0,j,cc)
     &      + reweight*weight
             do i=1,NPDFerror
               CUMULANTS(i,j,cc) = CUMULANTS(i,j,cc)
     &          + reweight*pdfvec(i)
             enddo
           enddo
           ELSEIF((ISHISTON(cc).ge.4).and.(ISHISTON(cc).le.5)) THEN ! reverse-cumulants on (each cumulate all right bins)
           reweight = wweight(cc)*smearwgt
           CUMULANTS(0,thisBin+1,cc) = CUMULANTS(0,thisBin+1,cc)
     &      + reweight*weight
             do i=1,NPDFerror
               CUMULANTS(i,thisBin+1,cc) = CUMULANTS(i,thisBin+1,cc)
     &          + reweight*pdfvec(i)
             enddo
           reweight = wweight(cc) ! every bin to left gets full weight in cumulant
           do j=1,thisBin
           CUMULANTS(0,j,cc) = CUMULANTS(0,j,cc)
     &      + reweight*weight
             do i=1,NPDFerror
               CUMULANTS(i,j,cc) = CUMULANTS(i,j,cc)
     &          + reweight*pdfvec(i)
             enddo
           enddo
           ENDIF

         ELSE ! smoothing is off or we are in unsmoothed region
           reweight = wweight(cc)
           IF(mod(ISHISTON(cc),2).eq.1) THEN ! regular histogram on
           ENTRIES(0,thisBin,cc)=ENTRIES(0,thisBin,cc)+reweight*weight
             do i=1,NPDFerror
               ENTRIES(i,thisBin,cc)=ENTRIES(i,thisBin,cc)
     &          + reweight*pdfvec(i)
             enddo
           ENDIF
           IF((ISHISTON(cc).ge.2).and.(ISHISTON(cc).le.3)) THEN
           do j=thisBin,NUMBINS(cc)
           CUMULANTS(0,j,cc)=CUMULANTS(0,j,cc)+reweight*weight
             do i=1,NPDFerror
               CUMULANTS(i,j,cc)=CUMULANTS(i,j,cc)
     &          + reweight*pdfvec(i)
             enddo
           enddo
           ELSEIF((ISHISTON(cc).ge.4).and.(ISHISTON(cc).le.5)) THEN
           do j=1,thisBin
           CUMULANTS(0,j,cc)=CUMULANTS(0,j,cc)+reweight*weight
             do i=1,NPDFerror
               CUMULANTS(i,j,cc)=CUMULANTS(i,j,cc)
     &          + reweight*pdfvec(i)
             enddo
           enddo
           ENDIF
         ENDIF
      ENDIF
      ENDIF
      ENDDO

CCC   Check code to see if we still use pdfeigens anywhere, remove later

c      if (NPDFerror.gt.0) then
         do i=1,NPDFerror
            pdfeigens(i)=pdfeigens(i)+pdfvec(i)
         enddo
c      endif
 
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC Clears histograms


      SUBROUTINE clearhist
      
      include 'histos.f'
      include 'pdfcontrol.f'

      INTEGER COUNTER1,COUNTER2,COUNTER3

      do COUNTER1=1,NHIST
      if (mod(ISHISTON(COUNTER1),2).eq.1) THEN ! no need to clear dead histograms
      do COUNTER2=1,NUMBINS(COUNTER1)
      do COUNTER3=0,NPDFerror
         ENTRIES(COUNTER3,COUNTER2,COUNTER1)=0d0
      enddo
      enddo
      endif
      enddo

ccc cumulant clearing separate loop for speed
      do COUNTER1=1,NHIST
      if (ISHISTON(COUNTER1).ge.2) THEN ! no need to clear dead cumulants
      do COUNTER2=1,NUMBINS(COUNTER1)
      do COUNTER3=0,NPDFerror
         CUMULANTS(COUNTER3,COUNTER2,COUNTER1)=0d0
      enddo
      enddo
      endif
      enddo

      do COUNTER1=1,NPDFerror
         pdfeigens(COUNTER1)=0d0
      enddo
      
      RETURN
      END
