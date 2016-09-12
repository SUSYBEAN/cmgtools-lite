
CCCCC Normalizations for gamma,Z,interference
      double precision norm,Ngg,Nzz,Nww
      double complex Ngz
      COMMON/normalization/norm,Ngg,Nzz,Nww,Ngz

CCCCC Mass and coupling parameters
      double precision Ml,Mw,Mz,Mh,Mt,GammaZ,GammaW,Qu,Qd,
     & Ql,gau,gad,gal,Gmu,alphaqed,alphaqed0
      double complex Mw2,Mz2,sw2,gvu,gvd,gvl
      COMMON/couplings/ Ml,Mw,Mz,Mh,Mt,GammaZ,GammaW,Qu,Qd,
     & Ql,gau,gad,gal,Gmu,alphaqed,alphaqed0,
     & Mw2,Mz2,sw2,gvu,gvd,gvl
    
CCCCC CKM parameters
      double precision V11,V12,V13,V22,V23
      COMMON/ckm/ V11,V12,V13,V22,V23

CCCCC Resonance parameters
      double precision Gamma,Gammagg
      COMMON/decay/ Gamma,Gammagg

CCCCC Strong coupling parameters
      double precision asZ,as
      COMMON/asorder/ asZ,as

CCCCC Scale choice, perturbative order, collider type, input coupling scheme
      INTEGER porder,porder_ew,ptype
      INTEGER whichscheme,photflag,DISonQED,EWflag
      DOUBLE PRECISION mur,muf,mufflag,murflag
      DOUBLE PRECISION m2,Scm
      COMMON/params/ m2,Scm
      COMMON/perturbative/porder,porder_ew
      COMMON/scale/ mur,muf,mufflag,murflag
      COMMON/progmode/ ptype
      COMMON/schmparam/ whichscheme,photflag,DISonQED,EWflag

CCCCC Integration region specification: Z-pole or not,
CCCCC photon on/off
      double precision mlower,mupper
      INTEGER zpole
      COMMON /intregion/ mlower,mupper,zpole
