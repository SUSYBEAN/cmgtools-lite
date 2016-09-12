
CCC ------- HISTOGRAM INFORMATION -----------------------

      INTEGER NHIST,BINPHIST,NUMBINSTOT,NPDFS
      parameter (NHIST=36) ! Number of histograms 
      parameter (BINPHIST=30) ! Max bins per histogram
      parameter (NUMBINSTOT=2160) ! NHIST*BINPHIST*2
      parameter (NPDFS=204) ! Max number of PDF eigenvectors (the same as NSMAX in pdfcontrol.f)

CCC   Parameters to give meaning to histogram numbers, to prevent coding errors

      INTEGER H_ZPT, H_ZRAP, H_INV_MASS, H_LPT, H_LETA, H_ALPT,H_ALETA,
     &        H_J1PT, H_J1ETA, H_J2PT, H_J2ETA, H_PHTPT, H_PHTETA,
     &        H_BEAM_THRUST, H_DRLL, H_DRJ1L, H_DRJ1AL, H_DRJ2L, 
     &        H_DRJ2AL, H_DRJJ, H_DRPHTL, H_H_T, H_M_T,
     &        H_AFB, H_AFB_DEN,
     &        H_A_0, H_A_1, H_A_2, H_A_3, H_A_4, H_A_DEN,
     &        H_PHI_CS, H_PHI_CS_NORM, H_COSTH_CS, H_COSTH_CS_NORM,
     &        H_DPHI

      PARAMETER (H_ZPT=1,H_ZRAP=2,H_INV_MASS=3,
     &           H_LPT=4,H_LETA=5,H_ALPT=6,H_ALETA=7,
     &           H_J1PT=8,H_J1ETA=9,H_J2PT=10,H_J2ETA=11,
     &           H_PHTPT=12,H_PHTETA=13,
     &           H_BEAM_THRUST=14,
     &           H_DRLL=15,H_DRJ1L=16,H_DRJ1AL=17,
     &           H_DRJ2L=18,H_DRJ2AL=19,H_DRJJ=20,
     &           H_DRPHTL=21,
     &           H_H_T=22,H_M_T=23,
     &           H_AFB=24,H_AFB_DEN=25,
     &           H_A_0=26,H_A_1=27,H_A_2=28,H_A_3=29,H_A_4=30,
     &           H_A_DEN=31,
     &           H_PHI_CS=32,H_PHI_CS_NORM=33,
     &           H_COSTH_CS=34,H_COSTH_CS_NORM=35,
     &           H_DPHI=36)
      
      character*128 dummystring2

      INTEGER NUMBINS(NHIST)
      DOUBLE PRECISION LOWEDGE(NHIST), HIGHEDGE(NHIST)
      character*60 HISTNAME(NHIST)
      INTEGER ISHISTON(NHIST) ! 0 for none, 1 for hist, 2 for cumulant, 3 for both hist and cumulant
                              ! 4 for reverse-cumulant, 5 for both hist and reverse-cumulant
c      LOGICAL eqlSpace(NHIST) ! true for equally spaced bins

      DOUBLE PRECISION binBounds(0:BINPHIST,NHIST) ! edges of bins
      
      DOUBLE PRECISION ENTRIES(0:NPDFS,BINPHIST,NHIST)
      DOUBLE PRECISION CUMULANTS(0:NPDFS,BINPHIST,NHIST)
      DOUBLE PRECISION ENTcum(0:NPDFS,0:NUMBINSTOT) ! cur hist avg
      DOUBLE PRECISION ENTsig(0:NPDFS,0:NUMBINSTOT) ! cur hist sigma^2
      DOUBLE PRECISION oldcum(0:NPDFS,0:NUMBINSTOT) ! tot hist avg
      DOUBLE PRECISION oldsig(0:NPDFS,0:NUMBINSTOT) ! tot hist sigma^2
      DOUBLE PRECISION wgtsum ! sum of Vegas weights so far
      INTEGER nevit ! number of evaluations this iteration

ccc method of combine iterations for histogram bins
      INTEGER histmethod

ccc   DOUBLE PRECISION zY_CScutoff
      
      COMMON/HISTO/ ENTRIES,CUMULANTS,LOWEDGE,HIGHEDGE,ENTcum,ENTsig,
     & oldcum,oldsig,binBounds,wgtsum,nevit,ISHISTON,dummystring2
     & ,histmethod

      COMMON/HISTINFO/ NUMBINS,HISTNAME!,eqlSpace

c      logical bFillHisto
c      common/FillHisto/bFillHisto

      double precision tmprand
      COMMON/randstore/tmprand

      integer smoothing
      double precision smcut
      COMMON/smoothing/ smcut, smoothing
