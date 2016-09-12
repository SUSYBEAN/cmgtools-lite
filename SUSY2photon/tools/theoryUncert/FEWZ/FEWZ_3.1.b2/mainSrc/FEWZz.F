      double precision function goalerr(isect)
      implicit none
      integer isect
      include 'sectortimingz.f'
      goalerr=reltime(isect)
      return
      END

CCCC  Initializes Z-specific parameters
      subroutine inithelp(inputfile)
      implicit none
      include 'cutcontrol.f'
      include 'declar.f'
      include 'parameters.f'
      include 'pdfcontrol.f'

      double precision pi
      parameter (pi=3.141592653589793D0)

      double precision Etot
      double complex cw2
      double precision tmpin ! for reading a real number for a complex variable

C      double precision GammaZll

CCC ====   INPUT PARAMETERS =================

      open(unit=14,file=inputfile,status='old')
      read(14,*)
      read(14,*) dummystring, Etot
      read(14,*)
      read(14,*) dummystring, muf
      read(14,*) dummystring, mur
      read(14,*)
      read(14,*) dummystring, ptype
      read(14,*)
      read(14,*)
      read(14,*) dummystring, alphaqed0
      read(14,*) dummystring, alphaqed
      read(14,*) dummystring, Gmu
      read(14,*)
      read(14,*) dummystring, Ml
      read(14,*) dummystring, Mw
      read(14,*) dummystring, GammaW
      read(14,*) dummystring, Mz
      read(14,*) dummystring, GammaZ
C      read(14,*) dummystring, GammaZll
      read(14,*) dummystring, Mt
      read(14,*) dummystring, Mh
      read(14,*)
      read(14,*)
      read(14,*)
      read(14,*) dummystring, whichscheme
      read(14,*) dummystring, tmpin
                              sw2 = tmpin
      read(14,*) dummystring, Qu
      read(14,*) dummystring, Qd
      read(14,*) dummystring, Ql
      read(14,*) dummystring, tmpin
                              gvu = tmpin
      read(14,*) dummystring, tmpin
                              gvd = tmpin
      read(14,*) dummystring, tmpin
                              gvl = tmpin
      read(14,*) dummystring, gau
      read(14,*) dummystring, gad
      read(14,*) dummystring, gal
      read(14,*)
      read(14,*)
      read(14,*) dummystring, epsrel
      read(14,*) dummystring, epsabs
      read(14,*) dummystring, nstart
      read(14,*) dummystring, nincrease
      read(14,*) dummystring, maxeval
      read(14,*) dummystring, seed
      read(14,*)
      read(14,*) dummystring, porder 
      read(14,*) dummystring, porder_ew
      read(14,*) dummystring, zpole
C      read(14,*) dummystring, mlower
C      read(14,*) dummystring, mupper
      read(14,*) dummystring, EWflag
      read(14,*) dummystring, photflag
      read(14,*)
      read(14,*) dummystring, invmassmin
      read(14,*) dummystring, invmassmax
      read(14,*) dummystring, tmassmin
      read(14,*) dummystring, tmassmax
      invmassmin = max(invmassmin,tmassmin)
      invmassmax = min(invmassmax,Etot)
      mlower = invmassmin
      mupper = invmassmax
      read(14,*) dummystring, ZpTmin
      read(14,*) dummystring, ZpTmax
      read(14,*) dummystring, ZYmin
      read(14,*) dummystring, ZYmax
      read(14,*) dummystring, leppTmin
      read(14,*) dummystring, leppTmax
      read(14,*) dummystring, aleppTmin
      read(14,*) dummystring, aleppTmax
      read(14,*) dummystring, softleppTmin
      read(14,*) dummystring, softleppTmax
      read(14,*) dummystring, hardleppTmin
      read(14,*) dummystring, hardleppTmax
      read(14,*)
      read(14,*) dummystring, lepabseta
      read(14,*) dummystring, lepetamin
      read(14,*) dummystring, lepetamax
      read(14,*)
      read(14,*) dummystring, alepabseta
      read(14,*) dummystring, alepetamin
      read(14,*) dummystring, alepetamax
      read(14,*)
      read(14,*) dummystring, softlepabseta
      read(14,*) dummystring, softlepetamin
      read(14,*) dummystring, softlepetamax
      read(14,*)
      read(14,*) dummystring, hardlepabseta
      read(14,*) dummystring, hardlepetamin
      read(14,*) dummystring, hardlepetamax

CCC PHOTON CUTS FOLLOW
      read(14,*)
      read(14,*) dummystring, lrdRrecb
      read(14,*) dummystring, minptr4obs
      read(14,*) dummystring, maxetar4obs
      read(14,*)
      read(14,*) dummystring, minnumpht
      read(14,*) dummystring, maxnumpht

CCC JET CUTS FOLLOW
      read(14,*)
      read(14,*)
      read(14,*) dummystring, jetalgo
      read(14,*) dummystring, jjdRmin
      read(14,*) dummystring, Rsep
      read(14,*) dummystring, minptj4obs
      read(14,*) dummystring, maxetaj4obs

      read(14,*)
      read(14,*) dummystring, minnumjet
      read(14,*) dummystring, maxnumjet
      read(14,*) dummystring, jet1ptmin

CCC ISOLATION CUTS FOLLOW
      read(14,*)
      read(14,*) dummystring, lldRmin
      read(14,*) dummystring, lldPhimin
      read(14,*) dummystring, lldPhimax
      read(14,*) dummystring, ljdRmin
      read(14,*) dummystring, lrdRmin
CCC COLLINS-SOPER
      read(14,*)
      read(14,*)
      read(14,*) dummystring, zY_CScutoff
CCC
      read(14,*)
      read(14,*)
      read(14,*) dummystring, pdfset
      read(14,*) dummystring, NPDFerror
      read(14,*)
      read(14,*)
      read(14,*) dummystring, whichas
      read(14,*)
      read(14,*) dummystring, pdfCL
      close(14)

CCC   Internal notations===========================
      scm = Etot**2
      tau = (0.5d0*Mz)**2/scm 
      Gamma = 1.0d0
      Gammagg = 1.0d0

CCC   Internal options
CCC   Choose the scale to fit PDFs for efficiency and whether to use changing factorization scale
      if (muf.gt.0d0) then
         mufflag = 0d0
         fitscale = muf
      else
         if (muf.eq.0d0) then
            mufflag = 1d0 ! 0 means no rescale factor
         else
            mufflag = dabs(muf) ! store the rescaling factor of the dynamical scale
         endif
         if ((mlower.le.Mz).and.(mupper.ge.Mz)) then
            fitscale = Mz*mufflag
         elseif (mupper.le.Mz) then
            if (mupper.ge.(mlower+1d0/mufflag)) then
               fitscale = mupper*mufflag-1d0
            else
               fitscale = mupper*mufflag
            endif
         elseif (mlower.ge.Mz) then
            fitscale = mlower*mufflag+1d0
         endif
         muf = fitscale
      endif
CCC   Whether to use changing renormalzation scale
      if (mur.gt.0d0) then
         murflag = 0d0
      else
         if (mur.eq.0d0) then
            murflag = 1d0 ! 0 means no rescale factor
         else
            murflag = dabs(mur) ! store the rescaling factor of the dynamical scale
         endif
         if (mufflag.gt.0d0) then ! just temporary assignment of mur for correction ratio
            mur = muf/mufflag*murflag
         else
            mur = muf
         endif
      endif
CCC   Whether to use DIS scheme for NLO QED evolution of quark PDF
      DISonQED = 1

CCC   Check Inputs

CCC   Determine which subset of EW correction is turned off
CCC   EWflag = (XXXX)_2 is a integer number in binary form (with X=0,1)
CCC   The 4th, 3rd, 2nd and 1st binary digit corresponds to Weak, ISR*FSR, ISR, FSR respectively
CCC   Examples:
CCC       EWflag = 0 means all on
CCC       EWflag = 1 = (0001) means FSR is off
CCC       EWflag = 11 = (1011) means Weak, ISR, FSR are off
      if (EWflag.ge.16) EWflag = mod(EWflag,16) ! make sure the binary number doesn't have 5th or higher digit

CCC   Determine input scheme
CCC   Note that alphaqed must be real because it's a physical parameter
CCC   otherwise IR pole won't cancel b/c photon cannot be on-shell if alphaqed is complex
      if (whichscheme.ne.0) then ! 1. Gmu scheme 2. AlphaMz scheme
         Mz2 = Mz**2-dcmplx(0d0,1d0)*GammaZ*Mz
         Mw2 = Mw**2-dcmplx(0d0,1d0)*GammaW*Mw
         sw2 = 1d0-Mw2/Mz2
         if (whichscheme.eq.1) then
c            alphaqed = dsqrt(2d0)/pi*Gmu*dble(zabs(Mw2*sw2))
            alphaqed = dsqrt(2d0)/pi*Gmu*dble(Mw**2*(1d0-Mw**2/Mz**2))
         endif
         !ec = 2d0*dsqrt(pi*alphaqed) ! list here for reference
         !gw = ec/dsqrt(sw2) ! list here for reference
         !wc = gw/dsqrt(1d0-sw2)/4d0 ! list here for reference
         Ql = -1d0 !*ec
         Qu = 2d0/3d0 !*ec
         Qd = -1d0/3d0 !*ec
         gvl = (-1d0/4d0+sw2)*4d0 !*wc
         gal = 1d0 !*wc
         gvu = (1d0/4d0-2d0/3d0*sw2)*4d0 !*wc
         gau = -1d0 !*wc
         gvd = (-1d0/4d0+sw2/3)*4d0 !*wc
         gad = 1d0 !*wc
      else ! manual scheme most parameters rely on inputs
         if (EWflag.lt.8) EWflag = EWflag+8 ! force weak correction off
         Mz2 = Mz**2
         Mw2 = Mw**2
      endif

CCC   Always have photon-mediated channel if weak correction is on 
      if ((porder_ew.gt.0).and.(EWflag.lt.8)) photflag = 0

CCC   Lepton mass too small, set to zero
      if (Ml.le.1d-8) Ml = 0d0

CCC   Define also the relative normalization of the Z,W cross-sections

      cw2 = 1d0-sw2
C     Assume all propagator are Z propagator squared
C     Ngg and Ngz adjust for Z+photon interference propagator and photon propagator squared
C     For the interference terms M_z*xM_g+M_zxM_g*=2xRe(M_zxM_g*), i.e. only take conjugate of the photon part
      Ngg = 1d0 ! the constant part of Ngg <--- not necessary actually
      Ngz = (1d0-4d0*sw2)/8d0/sw2/cw2 ! the constant part of Z and pht propagator, the final factor is Ngz2
C     turn off contribution with photon as DY mediator for LO and higher order QCD corr.
      if (photflag.eq.1) then
         Ngg = 0d0 ! <--- not necessary
         Ngz = 0d0
      endif
CCC   Disable using experimental value of Z leptonic decay partial width: GammaZll, since we removed it from input
C      if ((porder_ew.le.0).and.(GammaZll.ne.0d0)) then
C         Nzz = 3.0d0/abs(sw2*cw2)/alphaqed*GammaZll/Mz/16.0d0
C      else ! turn off using experimental Z->ll partial width if EW correction is on OR GammaZll is set to 0d0
C         GammaZll = Mz*(abs(gvl)**2+abs(gal)**2)*alphaqed/abs(sw2*cw2)/48d0 ! overriding GammaZll just in case it is used somewhere else
C         Nzz = (abs(gvl)**2+abs(gal)**2)/(abs(sw2*cw2)**2)/256d0
C      endif
      Nzz = dble(gvl*conjg(gvl)+gal*gal) ! gal is real
     $      /dble(sw2*conjg(sw2)*cw2*conjg(cw2))/256d0 ! use conjg instead of abs just to be safe

CCC=====================================================
CCC   Initialize EW corrections by
CCC   translating variables read in here
CCC   from parameter.f to constcommon_ew.f
      call InitEWSeq()

      END

CCCC  Writes Z-specific output
      subroutine endinghelp(i, tdiff)
      implicit none
      include 'pdfcontrol.f'
      include 'parameters.f'
      include 'histos.f'
      include 'declar.f'
      include 'cutcontrol.f'
      include 'sectorcontrol.f'

      double precision pi
      parameter (pi=3.141592653589793D0)

      integer i,counter1
      double precision tdiff

      character*8 pdfclstr

C      double precision GammaZll

      write(21,*) 'CMS collision Energy = ',dsqrt(scm)
      if (ptype.eq.1) then
         write(21,*) 'Collider type = pp'
      else
         write(21,*) 'Collider type = ppbar'
      endif

      write(21,*) '======================================'
      if (mufflag.gt.0d0) then
         if (mufflag.eq.1d0) then
            write(21,*) 'Factorization scale                = dynamical'
         else
            write(21,*) 'Factorization scale                = dynamical'
     &                  //' with factor ',mufflag
         endif
      else
         write(21,*) 'Factorization scale                = ',muf
      endif
      if (murflag.gt.0d0) then
         if (murflag.eq.1d0) then
            write(21,*) 'Renormalization scale              = dynamical'
         else
            write(21,*) 'Factorization scale                = dynamical'
     &                  //' with factor ',murflag
         endif
      else
         write(21,*) 'Renormalization scale              = ',mur
      endif
      write(21,*) '=================================='
      if ((porder_ew.eq.0).and.(porder.eq.0)) then 
         write(21,*) 'Leading Order Cross-Section'
      else if ((porder_ew.eq.0).and.(porder.eq.1)) then 
         write(21,*) 'NLO QCD Cross-Section'
      else if ((porder_ew.eq.0).and.(porder.ge.2)) then 
         write(21,*) 'NNLO QCD Cross-Section'
      else if ((porder_ew.ge.1).and.(porder.eq.0)) then
         write(21,*) 'NLO EW Cross-Section'
      else if ((porder_ew.ge.1).and.(porder.eq.1)) then
         write(21,*) 'NLO QCD + NLO EW Cross-Section'
      else if ((porder_ew.ge.1).and.(porder.ge.2)) then
         write(21,*) 'NNLO QCD + NLO EW Cross-Section'
      endif

      write(21,*) '======================================'
      if (whichscheme.eq.0) then
         write(21,*) 'EW input parameters: manual input'
      elseif (whichscheme.eq.1) then
         write(21,*) 'EW input parameters: Gmu scheme'
      elseif (whichscheme.eq.2) then
         write(21,*) 'EW input parameters: alpha(Mz) scheme'
      endif
      if (whichscheme.eq.1) then ! use Fermi constant in Gmu scheme only
         write(21,*) 'Fermi constant (1/GeV^2)           = ',Gmu
      else
         write(21,*) 'Alpha QED                          = ',alphaqed
      endif
      write(21,*) '======================================'
      write(21,*) 'Lepton mass (GeV)                  = ',Ml
      if (whichscheme.ne.0) then ! use W mass all schemes except for manual input scheme
         write(21,*) 'W mass (GeV)                       = ',Mw
         write(21,*) 'W width (GeV)                      = ',GammaW
      endif
      write(21,*) 'Z mass (GeV)                       = ',Mz
      write(21,*) 'Z width (GeV)                      = ',GammaZ
      if (whichscheme.ne.0) then ! use W mass all schemes except for manual input scheme
         write(21,*) 'Top mass (GeV)                     = ',Mt
         write(21,*) 'Higgs mass (GeV)                   = ',Mh
      endif
C      write(21,*) 'Z->ll partial width                = ',GammaZll
      if (whichscheme.eq.0) then ! use the following in manual input scheme only
         write(21,*) 'sin^2(theta)                       = ',dble(sw2)
         write(21,*) 'up quark charge                    = ',Qu
         write(21,*) 'down quark charge                  = ',Qd
         write(21,*) 'lepton chage                       = ',Ql
         write(21,*) 'up quark vector coupling           = ',dble(gvu)
         write(21,*) 'down quark vector coupling         = ',dble(gvd)
         write(21,*) 'lepton vector coupling             = ',dble(gvl)
         write(21,*) 'up quark axial coupling            = ',gau
         write(21,*) 'down quark axial coupling          = ',gad
         write(21,*) 'lepton axial coupling              = ',gal
      endif

      write(21,*) '======================================'
      if (murflag.gt.0d0) then
         write(21,*) 'Strong coupling at Mz              = ',asZ
      else
         write(21,*) 'Strong coupling                    = ',as*pi
      endif
      if ((NPDFerror.gt.0).and.
     &    (pdfset.eq.'MSTW2008NLO'.or.pdfset.eq.'MSTW2008NNLO')) then ! handle additional MSTW parameters
         if (pdfCL.eq.1) then
            pdfclstr = '(68 CL)'
         else
            pdfclstr = '(90 CL)'
         endif
         if (whichas.eq.1) then
            write(21,*) 'PDF set = ',pdfset(1:12), ' + alphaS ', pdfclstr
         else
            write(21,*) 'PDF set = ',pdfset(1:12), pdfclstr
         endif
      else
         write(21,*) 'PDF set = ',pdfset
      endif

      write(21,*) '======================================'     
      if (zpole.eq.1) then
         write(21,*) 'Z pole focus ON'
      else
         write(21,*) 'Z pole focus OFF'
         write(21,*) 'Lower M limit for non-Zpole     = ',mlower
         write(21,*) 'Upper M limit for non-Zpole     = ',mupper
      endif

      if (photflag.eq.1) then
         write(21,*) '======================================'
         write(21,*) 'Photon interference turned off:  Z only'
      endif

      write(21,*) '======================================'
      write(21,*) 'Lepton-pair invariant mass minimum = ',invmassmin
      write(21,*) 'Lepton-pair invariant mass maximum = ',invmassmax
      write(21,*) 'Transverse mass minimum            = ',tmassmin
      write(21,*) 'Transverse mass maximum            = ',tmassmax
      write(21,*) 'Z pT minimum                       = ',ZpTmin
      write(21,*) 'Z pT maximum                       = ',ZpTmax
      write(21,*) 'Z rapidity minimum                 = ',ZYmin
      write(21,*) 'Z rapidity maximum                 = ',ZYmax
      write(21,*) 'Lepton pT minimum                  = ',leppTmin
      write(21,*) 'Lepton pT maximum                  = ',leppTmax
      write(21,*) 'Anti-lepton pT minimum             = ',aleppTmin
      write(21,*) 'Anti-lepton pT maximum             = ',aleppTmax
      if (lepabseta.eq.1) then
         write(21,*) 'Lepton eta absolute value?         =  ', .true.
      else
         write(21,*) 'Lepton eta absolute value?         =  ', .false.
      endif
ccc      write(21,*) 'Lepton eta absolute value?         =  ',lepabseta
      write(21,*) 'Lepton pseudorapidity minimum      = ',lepetamin
      write(21,*) 'Lepton pseudorapidity maximum      = ',lepetamax
      if (alepabseta.eq.1) then
         write(21,*) 'Anti-lepton eta absolute value?    =  ', .true.
      else
         write(21,*) 'Anti-lepton eta absolute value?    =  ', .false.
      endif
ccc      write(21,*) 'Anti-lepton eta absolute value?    =  ',alepabseta
      write(21,*) 'Anti-lepton pseudorapidity minimum = ',alepetamin
      write(21,*) 'Anti-lepton pseudorapidity maximum = ',alepetamax

      write(21,*) '======================================'
      write(21,*) 'Photon Recombination Cone Size     =   ',lrdRrecb
      write(21,*) 'Minimum pT for Observable Photon   = ',minptr4obs
      write(21,*) 'Maximum eta for Observable Photon  = ',maxetar4obs
      if (minnumpht.gt.porder_ew) then
         write(21,*) 'NOTICE: Minimum Number of Phton is ',minnumpht
         write(21,*) 'while Pert. Order in EW is N^',porder_ew,' LO'
         write(21,*) '======================================'
      endif
      write(21,*) 'Minimum Number of Photon           =  ',minnumpht
      write(21,*) 'Maximum Number of Photon           =  ',maxnumpht

      write(21,*) '======================================'
      write(21,*) 'Jet merging algorithm              =   ',jetalgo
      write(21,*) 'Jet algorithm cone size (deltaR)   = ',jjdRmin
      write(21,*) 'parton-parton Rsep (cone algo only)= ',Rsep
      write(21,*) 'Minimum pT for Observable Jets     = ',minptj4obs
      write(21,*) 'Maximum eta for Observable Jets    = ',maxetaj4obs
      if (minnumjet.gt.porder) then
         write(21,*) 'NOTICE: Minimum Number of Jets is ',minnumjet
         write(21,*) 'while Pert. Order in QCD is N^',porder,' LO'
         write(21,*) '======================================'
      endif
      write(21,*) 'Minimum Number of Jets             =  ',minnumjet
      write(21,*) 'Maximum Number of Jets             =  ',maxnumjet
      write(21,*) 'Leading jet pT minimum             = ',jet1ptmin

      write(21,*) '======================================'
      write(21,*) 'Lep-Anti-lep deltaR minimum        = ',lldRmin
      write(21,*) 'Lep-Anti-lep deltaPhi minimum      = ',lldPhimin
      write(21,*) 'Lep-Anti-lep deltaPhi maximum      = ',lldPhimax
      write(21,*) 'Lep-Jet deltaR minimum             = ',ljdRmin
      write(21,*) 'Lep-Photon deltaR minimum          = ',lrdRmin
      write(21,*) '======================================'
      write(21,*) 'Rapidity of Z cutoff for CS angles = ',zY_CScutoff
      write(21,*) ''
      write(21,*) '===========   VEGAS PARAMETERS   ========'     
      write(21,*) ''
      write(21,*) 'Maximum number of evaluations      = ', maxeval
      write(21,*) 'Actual number of evaluations       = ', neval
      write(21,*) 'Requested relative precision (%)   = ',epsrel
      write(21,*) 'Requested absolute precision       = ',epsabs
      write(21,*) 'Total integration time (sec)       = ',tdiff
      write(21,*) 'Random number seed                 = ', seed
      write(21,*) ''
      write(21,*) '===========   RESULT    ================='     
      write(21,*) ''
      if (i.eq.1) then
      write(21,*) 'Sigma (pb)                  = ',integral(1)
      write(21,*) 'Error (pb)                  = ',error(1)
      write(21,*) 'chi^2/iteration             = ',prob(1)
      else ! write out for each PDF eigenvecto
      do counter1=1,NPDFerror
      write(21,*) 'Sigma (pb)                  = ',oldcum(counter1,0)
      write(21,*) 'Error (pb)                  = '
     & ,dsqrt(oldsig(counter1,0))
      write(21,*) ''
      end do
      end if

      END
