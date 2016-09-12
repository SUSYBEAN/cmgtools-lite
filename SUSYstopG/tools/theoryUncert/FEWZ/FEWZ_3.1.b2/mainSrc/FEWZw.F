      double precision function goalerr(isect)
      implicit none
      integer isect
      include 'sectortimingw.f'
      goalerr=reltime(isect)
      return
      END

CCCC  Initializes W-specific parameters
      subroutine inithelp(inputfile)
      implicit none
      include 'cutcontrol.f'
      include 'declar.f'
      include 'parameters.f'
      include 'pdfcontrol.f'

      double precision pi
      parameter (pi=3.141592653589793D0)

      double precision Etot

      double precision tmpin

      double precision GammaWlv
      common /ioparam/ GammaWlv ! parameters to be shared between inithelp() and endinghelp()

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
      read(14,*) dummystring, alphaqed
      read(14,*) dummystring, Gmu
      read(14,*)
      read(14,*) dummystring, Mw
      read(14,*) dummystring, GammaW
      read(14,*) dummystring, GammaWlv
      read(14,*) dummystring, tmpin
                              sw2 = tmpin
      read(14,*) dummystring
      read(14,*) dummystring, V11
      read(14,*) dummystring, V12
      read(14,*) dummystring, V22
      read(14,*) dummystring, V13
      read(14,*) dummystring, V23
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
      porder_ew = 0 ! need to give it a default value
      read(14,*) dummystring, zpole
c      read(14,*) dummystring, mlower
c      read(14,*) dummystring, mupper
      read(14,*)
      read(14,*) dummystring, invmassmin
      read(14,*) dummystring, invmassmax
      read(14,*) dummystring, tmassmin
      read(14,*) dummystring, tmassmax
      invmassmin = max(invmassmin,tmassmin)
      invmassmax = min(invmassmax,Etot)
      mlower = max(1d-1,invmassmin)
      mupper = min(Etot-1d-1,invmassmax)
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
      read(14,*) dummystring, lldPhimin
      read(14,*) dummystring, lldPhimax
      read(14,*) dummystring, ljdRmin
      zY_CScutoff = 0d0
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

CCC   Define Z cuts not used for W
      Mz = dsqrt(dble(Mw**2)/(1d0-dble(sw2)))
C      GammaZ = 1d0
C      GammaZll = 1d0
C      Qu = 1d0
C      Qd = 1d0
C      Ql = 1d0
C      gvu = 1d0
C      gvd = 1d0
C      gvl = 1d0
C      gau = 1d0
C      gad = 1d0
C      gal = 1d0
C      photflag = 0

      alepabseta = 1
      alepetamin = 0d0
      alepetamax = 1d10
      softlepabseta = 1
      softlepetamin = 0d0
      softlepetamax = 1d10
      hardlepabseta = 1
      hardlepetamin = 0d0
      hardlepetamax = 1d10
      lldRmin = 0d0
      zY_CScutoff = 0d0
      minnumpht = 0
      maxnumpht = 10

CCC   Internal notations===========================
      scm = Etot**2
      tau = (0.5d0*Mw)**2/scm

CCC   Internal options
CCC   Choose the scale to fit PDFs for efficiency
      if (muf.gt.0d0) then
         mufflag = 0d0
         fitscale = muf
      else
         if (muf.eq.0d0) then
            mufflag = 1d0 ! 0 means no rescale factor
         else
            mufflag = dabs(muf) ! store the rescaling factor of the dynamical scale
         endif
         if ((mlower.le.Mw).and.(mupper.ge.Mw)) then
            fitscale = Mw*mufflag
         elseif (mupper.le.Mw) then
            if (mupper.ge.(mlower+1d0/mufflag)) then
               fitscale = mupper*mufflag-1d0
            else
               fitscale = mupper*mufflag
            endif
         elseif (mlower.ge.Mw) then
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

CCC   Define also the relative normalization of the Z,W  cross-sections

      Nww = 3d0/4d0/dble(sw2)/alphaqed*GammaWlv/Mw
C      GammaWlv = Mw*alphaqed/sw2/12d0 ! if we wanna override the manual input GammaWlv
C      Nww = 1/sw2**2/16d0


c      V11=0.97410D0**2
c      V12=0.2196D0**2
c      V13=1.0d0-V11-0.2196D0**2
c      V22=0.996D0**2
c      V23=1.0d0-V12-0.996D0**2
      V11=V11**2
      V12=V12**2
      V13=V13**2
      V22=V22**2
      V23=V23**2

      END

CCCC  W-specific output

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

      double precision GammaWlv
      common /ioparam/ GammaWlv ! parameters to be shared between inithelp() and endinghelp()

      write(21,*) 'CMS collision Energy = ',dsqrt(scm)
      if (ptype.eq.1) then
         write(21,*) 'Collider type = pp, W- mode'
      elseif (ptype.eq.3) then
         write(21,*) 'Collider type = pp, W+ mode'
      elseif (ptype.eq.2) then
         write(21,*) 'Collider type = ppbar, W- mode'
      else
         write(21,*) 'Collider type = ppbar, W+ mode'
      endif

      write(21,*) '======================================'
      if (mufflag.gt.0d0) then
         if (mufflag.eq.1d0) then
            write(21,*) 'Factorization scale                = dynamical'
         else
            write(21,*) 'Factorization scale                = dynamical'
     &                  //' with factor',mufflag
         endif
      else
         write(21,*) 'Factorization scale                = ',muf
      endif
      if (murflag.gt.0d0) then
         if (murflag.eq.1d0) then
            write(21,*) 'Renormalization scale              = dynamical'
         else
            write(21,*) 'Factorization scale                = dynamical'
     &                  //' with factor',murflag
         endif
      else
         write(21,*) 'Renormalization scale              = ',mur
      endif
      write(21,*) '=================================='
      if (porder.eq.0) then 
         write(21,*) 'Leading Order Cross-Section'
      else if (porder.eq.1) then 
         write(21,*) 'NLO Cross-Section'
      else if (porder.eq.2) then 
         write(21,*) 'NNLO Cross-Section'
      endif

      write(21,*) '======================================'
      write(21,*) 'Alpha QED                          = ',alphaqed
ccc      write(21,*) 'Fermi constant (1/GeV^2)           = ',Gmu
      write(21,*) '======================================'
      write(21,*) 'W mass (GeV)                       = ',Mw
      write(21,*) 'W width (GeV)                      = ',GammaW
C      write(21,*) 'Z mass (GeV)                       = ',Mz
C      write(21,*) 'Z width (GeV)                      = ',GammaZ
      write(21,*) 'W->lv partial width                = ',GammaWlv
      write(21,*) 'sin^2(theta)                       = ',dble(sw2)
      write(21,*) 'Vud                                = ',sqrt(V11)
      write(21,*) 'Vus=Vcd                            = ',sqrt(V12)
      write(21,*) 'Vcs                                = ',sqrt(V22)
      write(21,*) 'Vub                                = ',sqrt(V13)
      write(21,*) 'Vcb                                = ',sqrt(V23)
C      write(21,*) 'up quark charge                    = ',Qu
C      write(21,*) 'down quark charge                  = ',Qd
C      write(21,*) 'lepton chage                       = ',Ql
C      write(21,*) 'up quark vector coupling           = ',gvu
C      write(21,*) 'down quark vector coupling         = ',gvd
C      write(21,*) 'lepton vector coupling             = ',gvl
C      write(21,*) 'up quark axial coupling            = ',gau
C      write(21,*) 'down quark axial coupling          = ',gad
C      write(21,*) 'lepton axial coupling              = ',gal

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
         write(21,*) 'W pole focus ON'
      else
         write(21,*) 'W pole focus OFF'
      endif

      write(21,*) '======================================'
      write(21,*) 'Lepton-pair invariant mass minimum = ',invmassmin
      write(21,*) 'Lepton-pair invariant mass maximum = ',invmassmax
      write(21,*) 'Transverse mass minimum            = ',tmassmin
      write(21,*) 'Transverse mass maximum            = ',tmassmax
      write(21,*) 'W pT minimum                       = ',ZpTmin
      write(21,*) 'W pT maximum                       = ',ZpTmax
      write(21,*) 'W rapidity minimum                 = ',ZYmin
      write(21,*) 'W rapidity maximum                 = ',ZYmax
      write(21,*) 'Charged lepton pT minimum          = ',leppTmin
      write(21,*) 'Charged lepton pT maximum          = ',leppTmax
      write(21,*) 'Missing pT minimum                 = ',aleppTmin
      write(21,*) 'Missing pT maximum                 = ',aleppTmax
      if (lepabseta.eq.1) then
         write(21,*) 'Lepton eta absolute value?         =  ', .true.
      else
         write(21,*) 'Lepton eta absolute value?         =  ', .false.
      endif
ccc      write(21,*) 'Lepton eta absolute value?         =  ',lepabseta
      write(21,*) 'Lepton pseudorapidity minimum      = ',lepetamin
      write(21,*) 'Lepton pseudorapidity maximum      = ',lepetamax

      write(21,*) '======================================'
      write(21,*) 'Jet merging algorithm              =   ',jetalgo
      write(21,*) 'Jet algorithm cone size (deltaR)   = ',jjdRmin
      write(21,*) 'parton-parton Rsep (cone algo only)= ',Rsep
      write(21,*) 'Minimum pT for Observable Jets     = ',minptj4obs
      write(21,*) 'Maximum eta for Observable Jets    = ',maxetaj4obs

      if (minnumjet.gt.porder) then
         write(21,*) 'NOTICE: Minimum Number of Jets is ',minnumjet
         write(21,*) 'while Perturbation Order is N^',porder,' LO'
         write(21,*) '======================================'
      endif
      
      write(21,*) 'Minimum Number of Jets             =  ',minnumjet
      write(21,*) 'Maximum Number of Jets             =  ',maxnumjet
      write(21,*) 'Leading jet pT minimum             = ',jet1ptmin
      write(21,*) '======================================'
C      write(21,*) 'Lep-Anti-lep deltaR minimum        = ',lldRmin
      write(21,*) 'Lep-neutrino deltaPhi minimum      = ',lldPhimin
      write(21,*) 'Lep-neutrino deltaPhi maximum      = ',lldPhimax
      write(21,*) 'Lep-Jet deltaR minimum             = ',ljdRmin
      write(21,*) '======================================'
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
