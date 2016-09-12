      integer whichas, pdfCL, NPDFerror
      double precision fitscale
      COMMON/pdfcontrol/ fitscale, whichas, pdfCL, NPDFerror
      integer NSMAX
      parameter(NSMAX=204) ! the same as NPDFS in histos.f
      integer NUMFIT
      parameter(NUMFIT=1500)

      character*20 suffix
      character*60 pdfset
      COMMON/pdfset/ pdfset,suffix

CCCCC Array holding the result for each PDF set

      double precision pdfeigens
      dimension pdfeigens(NSMAX)
      COMMON/PDFanalysis/ pdfeigens

CCCCC PDF interpolation parameters
      COMMON /phparam/ xxph, AAPH, BBPH, CCPH, hhph
      double precision xxph, AAPH, BBPH, CCPH, hhph
      dimension xxph(0:NSMAX,NUMFIT), AAPH(0:NSMAX,NUMFIT),
     & BBPH(0:NSMAX,NUMFIT), CCPH(0:NSMAX,NUMFIT)
      COMMON /gparam/ xxg, AAG, BBG, CCG, hhg
      double precision xxg, AAG, BBG, CCG, hhg
      dimension xxg(0:NSMAX,NUMFIT), AAG(0:NSMAX,NUMFIT), 
     & BBG(0:NSMAX,NUMFIT), CCG(0:NSMAX,NUMFIT)
      COMMON /uvparam/ xxuv, AAUV, BBUV, CCUV, hhuv
      double precision xxuv, AAUV, BBUV, CCUV, hhuv
      dimension xxuv(0:NSMAX,NUMFIT), AAUV(0:NSMAX,NUMFIT), 
     & BBUV(0:NSMAX,NUMFIT), CCUV(0:NSMAX,NUMFIT)
      COMMON /dvparam/ xxdv, AADV, BBDV, CCDV, hhdv
      double precision xxdv, AADV, BBDV, CCDV, hhdv
      dimension xxdv(0:NSMAX,NUMFIT), AADV(0:NSMAX,NUMFIT),
     & BBDV(0:NSMAX,NUMFIT), CCDV(0:NSMAX,NUMFIT)
      COMMON /usparam/ xxus, AAUS, BBUS, CCUS, hhus
      double precision xxus, AAUS, BBUS, CCUS, hhus
      dimension xxus(0:NSMAX,NUMFIT), AAUS(0:NSMAX,NUMFIT),
     & BBUS(0:NSMAX,NUMFIT), CCUS(0:NSMAX,NUMFIT)
      COMMON /dsparam/ xxds, AADS, BBDS, CCDS, hhds
      double precision xxds, AADS, BBDS, CCDS, hhds
      dimension xxds(0:NSMAX,NUMFIT), AADS(0:NSMAX,NUMFIT), 
     & BBDS(0:NSMAX,NUMFIT), CCDS(0:NSMAX,NUMFIT)
      COMMON /svparam/ xxsv, AASV, BBSV, CCSV, hhsv
      double precision xxsv, AASV, BBSV, CCSV, hhsv
      dimension xxsv(0:NSMAX,NUMFIT), AASV(0:NSMAX,NUMFIT), 
     & BBSV(0:NSMAX,NUMFIT), CCSV(0:NSMAX,NUMFIT)
      COMMON /svparam/ xxss, AASS, BBSS, CCSS, hhss
      double precision xxss, AASS, BBSS, CCSS, hhss
      dimension xxss(0:NSMAX,NUMFIT), AASS(0:NSMAX,NUMFIT), 
     & BBSS(0:NSMAX,NUMFIT), CCSS(0:NSMAX,NUMFIT)
      COMMON /chmparam/ xxchm, AACHM, BBCHM, CCCHM, hhchm
      double precision xxchm, AACHM, BBCHM, CCCHM, hhchm
      dimension xxchm(0:NSMAX,NUMFIT), AACHM(0:NSMAX,NUMFIT), 
     & BBCHM(0:NSMAX,NUMFIT), CCCHM(0:NSMAX,NUMFIT)
      COMMON /btmparam/ xxbtm, AABTM, BBBTM, CCBTM, hhbtm
      double precision xxbtm, AABTM, BBBTM, CCBTM, hhbtm
      dimension xxbtm(0:NSMAX,NUMFIT), AABTM(0:NSMAX,NUMFIT), 
     & BBBTM(0:NSMAX,NUMFIT), CCBTM(0:NSMAX,NUMFIT)


CCC   For ABKM/GJR/MSTW sets, need to track whether to reweight each piece   
CCC   by the alpha_s value appropriate to the given PDF eigenset
CCC   Track whether we're dealing with LO, NLO, NNLO piece

      INTEGER CoeffFuncOrder
      DOUBLE PRECISION asRwt(NSMAX,0:2)
      COMMON /asreweight/ asRwt,CoeffFuncOrder
      SAVE /asreweight/

      COMMON /fordyi/ tau
      double precision tau

CCCCC PDF calling parameters
      double precision uval1,dval1,usea1,dsea1,sval1,ssea1,chm1,btm1,
     &                 glu1,!phot1,
     &                 uval2,dval2,usea2,dsea2,sval2,ssea2,chm2,btm2,
     &                 glu2!,phot2
      COMMON /pdfs/ uval1,dval1,usea1,dsea1,sval1,ssea1,chm1,btm1,glu1
     &              ,!phot1
     &              uval2,dval2,usea2,dsea2,sval2,ssea2,chm2,btm2,glu2
     &              !,phot2
c comment out photon in common block because photons are used in EW correction module only, which does not use getPDFs()


CCCCC Holds last x called for particular set, so we know when to use old value
      INTEGER NPDFCACHE
      PARAMETER(NPDFCACHE=4)
      DOUBLE PRECISION lastX(NPDFCACHE,0:NSMAX)

CCCCC Last PDF calls for each eigenvector
      INTEGER NPDFTYPE
      PARAMETER(NPDFTYPE=9) !instead of 10 because we don't use getPDFs() to cache photon in EW correction module 
      DOUBLE PRECISION pdfStore(NPDFTYPE,NPDFCACHE,0:NSMAX)
      COMMON /pdfCache/ pdfStore, lastX
      SAVE /pdfCache/

CCCCC Luminosity cache
      INTEGER NLUMCACHE
      PARAMETER(NLUMCACHE=32)
      DOUBLE PRECISION lastLumx1(NLUMCACHE),
     &                 lastLumx2(NLUMCACHE),
     &                 lastFlag(NLUMCACHE),
     &                 lastLum(NLUMCACHE),
     &                 lumins(0:NSMAX,NLUMCACHE)
      COMMON /lumCache/ lumins,lastLumx1,lastLumx2,lastFlag,lastLum
      SAVE /lumCache/

