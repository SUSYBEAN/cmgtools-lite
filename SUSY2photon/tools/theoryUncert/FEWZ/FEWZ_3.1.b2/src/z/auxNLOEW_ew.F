CCC -------------------------------------------------------------------

CCC Translate variables from FEWZ to here

      SUBROUTINE TRANSVAR(alf_c,Ml_c,Mw_c,GAMw_c,Mz_c,GAMz_c,
     &                    Mh_c,Mt_c,
     &                    geu_c,ged_c,gel_c,
     &                    gvu_c,gvd_c,gvl_c,gau_c,gad_c,gal_c,
     &                    smax_c,smin_c,prej0_c)
      IMPLICIT NONE

      include 'cutcontrol.f'
      include 'parameters.f'

      double precision alf_c,Ml_c,Mw_c,GAMw_c,Mz_c,GAMz_c,
     &                 Mh_c,Mt_c,
     &                 geu_c,ged_c,gel_c,
     &                 smax_c,smin_c,prej0_c
      double complex gvu_c,gvd_c,gvl_c,gau_c,gad_c,gal_c

      double precision elecpl
      double complex weakcpl
ccc to include experimental Z->ll partial width to compare with original FEWZ LO results (for test only)
c      double precision weakcpl2,weakadj
c      common/weakadjst/ weakadj

      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      Ml_c = Ml
      Mw_c = Mw
      Mz_c = Mz
      GAMw_c = GammaW
      GAMz_c = GammaZ
      Mh_c = Mh
      Mt_c = Mt

      alf_c = alphaqed
      elecpl = dsqrt(4d0*PI*alphaqed)
      weakcpl = sqrt(alphaqed*PI/4d0/sw2/(1d0-sw2)) !elecpl/sqrt(sw2*(1d0-sw2))/4d0
c      weakcpl2 = dsqrt(12d0*PI*GammaZll/Mz/(abs(gvl)**2+abs(gal)**2))
c      weakadj = weakcpl2/abc(weakcpl) ! since it always appear as weakcpl*conjg(weakcpl) for where weakadj is needed
      geu_c = Qu*elecpl !real
      ged_c = Qd*elecpl !real
      gel_c = Ql*elecpl !real
      gvu_c = gvu*weakcpl !possible complex
      gvd_c = gvd*weakcpl !possible complex
      gvl_c = gvl*weakcpl !possible complex
      gau_c = gau*weakcpl !possible complex
      gad_c = gad*weakcpl !possible complex
      gal_c = gal*weakcpl !possible complex

      smax_c = min(DBLE(scm),DBLE(mupper**2))
      smin_c = max(DBLE(mlower**2),DBLE(4d0*Ml**2))
      smax_c = max(DBLE(smax_c),DBLE(smin_c))

CCC a factor of norm is multiplied in crsection.F: prej0(original) = 1/36*0.3893792922d9 = prej0(new)*norm w/ norm = 4/9*PI*alphaqed**2*0.3893792922d9
C         prej0_c = PI/elecpl**4
CCC stop using the norm adjustment, define prej0 to be the original one since we modified the way NLOEWVAL is called
      prej0_c = 1d0/36d0*0.3893792922d9

      RETURN
      END

CCC -------------------------------------------------------------------

CCC Pass particle info variables from here to FEWZ

      SUBROUTINE PassVar(PART_c,part)
      IMPLICIT NONE

      DOUBLE PRECISION PART_c(3,0:3)
      include 'particles.f'
      include 'cutcontrol.f'

      INTEGER j
      DOUBLE PRECISION ptr

      do j=1,4
         part(1,j) = PART_c(1,j-1) + PART_c(2,j-1)
         part(2,j) = PART_c(1,j-1)
         part(3,j) = PART_c(2,j-1)
         part(4,j) = 0d0
         part(5,j) = 0d0
         part(6,j) = PART_c(3,j-1)
      enddo

c      if (part(6,1).lt.ptbuf) then ! ptr < part(6,1), no need to test for ptr
c         do j=1,4
c            part(6,j) = 0d0
c         enddo
c      else
c         ptr = dsqrt(part(6,2)**2+part(6,3)**2)
c         if (ptr.lt.ptbuf) then
c            do j=2,3
c               part(6,j) = 0d0
c            enddo
c         endif
c      endif

      RETURN
      END

CCC -------------------------------------------------------------------

CCC Pass particle info variables from here to FEWZ

      SUBROUTINE ReversePart(part)
      IMPLICIT NONE

      include 'particles.f'

      INTEGER i,j

      do i=1,NPART
      do j=4,4 !j=2,4
         part(i,j) = -part(i,j)
      enddo
      enddo

      RETURN
      END

CCC -------------------------------------------------------------------

      SUBROUTINE PDFSeqAll(x1,x2,z1,z2,convtype,channel)
      IMPLICIT NONE

CCC  The argument 'convtype' of pdfarry denotes types: 1.no convoluion; 2. convolution with parton 1; 3. convolution with parton 2
CCC  The argument 'channel': 1. qi-qbi 2. phot-phot ...
      integer channel,convtype,jcvt,iset
      double precision x1,x2,z1,z2,y1,y2,q,rsf

      include 'parameters.f' ! to read ptype
      include 'pdflumcommon_ew.f' ! part of pdfcontrol.f is included in pdflumcommon_ew.f already

      double precision uval1,usea1,dval1,dsea1,sval1,ssea1,chm1,btm1,phot1,
     &                 uval2,usea2,dval2,dsea2,sval2,ssea2,chm2,btm2,phot2
      double precision pdfval(1:10)

      double precision retpdfph,retpdfg,retpdfuv,retpdfdv,retpdfus,retpdfds,
     & retpdfsv,retpdfss,retpdfchm,retpdfbtm
      external retpdfph,retpdfg,retpdfuv,retpdfdv,retpdfus,retpdfds,
     &         retpdfsv,retpdfss,retpdfchm,retpdfbtm      

CCC To save the scale low limit allowed by PDF set for determining the right scale for the convolution case
      double precision qmin,qfct
      data qmin/-1d0/
      save qmin,qfct
      if (qmin.le.0d0) then ! make sure it's only done once
         qmin = 1d-1 ! small enough to be lower than the allowed limit
         CALL CheckPDFRange(-1d0,qmin) ! obtain the allowed limit
         qfct = min(murflag,mufflag)
         if (qfct.eq.0d0) qfct = murflag+mufflag ! in case that one of the flag is 0
         if (qfct.gt.0d0) qmin = qmin/qfct
         qmin = qmin*mufflag !! worry about murflag to be consistent with QCD LO, but only care about factorization scale for EW corr.
      endif

CCC Unlike QCD correction, can only use s_qq for dynamical scale
CCC because s_ll would break degeneracy of FSR collinear radiation
CCC which means PDF with convolution should be called with scale s/z1, s/z2 or s/z1/z2
      if (mufflag.gt.0) muf = mufflag*dsqrt(m2)
      q = muf ! in case of constant scale choice
      do jcvt = 1,max(1,convtype) !convtype should be between [1,4]
         if (jcvt.eq.1) then
            if (mufflag.gt.0d0) q = max(qmin,muf)
            y1 = x1
            y2 = x2
         elseif (jcvt.eq.2) then
            if (mufflag.gt.0d0) q = max(qmin,muf/z1)
            y1 = x1/z1
            y2 = x2
         elseif (jcvt.eq.3) then
            if (mufflag.gt.0d0) q = max(qmin,muf/z2)
            y1 = x1
            y2 = x2/z2
         elseif (jcvt.eq.4) then
            if (mufflag.gt.0d0) q = max(qmin,muf/z1/z2)
            y1 = x1/z1
            y2 = x2/z2
         endif
         rsf = y1*y2 ! rescale factor rsf is because that the obtained PDF value are actually x*f(x)
         do iset = 0,NPDFerror
            if (channel.eq.1) then
               if (dabs(fitscale-q).le.0.5d0) then
                  uval1 = retpdfuv(y1,iset)
                  dval1 = retpdfdv(y1,iset)
                  usea1 = retpdfus(y1,iset)
                  dsea1 = retpdfds(y1,iset)
                  sval1 = retpdfsv(y1,iset)
                  ssea1 = retpdfss(y1,iset)
                  chm1  = retpdfchm(y1,iset)
                  btm1  = retpdfbtm(y1,iset)
C                  glu1  = retpdfg(y1,iset)
                  usea2 = retpdfus(y2,iset)
                  dsea2 = retpdfds(y2,iset)
                  ssea2 = retpdfss(y2,iset)
                  chm2  = retpdfchm(y2,iset)
                  btm2  = retpdfbtm(y2,iset)
C                  glu2  = retpdfg(y2,iset)
               else
                  CALL ObtainPDFs(y1,q,iset,pdfval,0)
                  uval1 = pdfval(1) !uv
                  dval1 = pdfval(2) !dv
                  usea1 = pdfval(3) !us
                  dsea1 = pdfval(4) !ds
                  sval1 = pdfval(5) !sv
                  ssea1 = pdfval(6) !ss
                  chm1 = pdfval(7) !chm
                  btm1 = pdfval(8) !btm
C                  glu1 = pdfval(9) !g
C                  phot1 = pdfval(10) !ph
                  CALL ObtainPDFs(y2,q,iset,pdfval,0)
                  usea2 = pdfval(3) !us
                  dsea2 = pdfval(4) !ds
                  sval2 = pdfval(5) !sv
                  ssea2 = pdfval(6) !ss
                  chm2 = pdfval(7) !chm
                  btm2 = pdfval(8) !btm
               endif
               if (ptype.eq.1) then
                  uubplum(jcvt,iset) = ((uval1+usea1)*usea2+chm1*chm2)/rsf
                  ddbplum(jcvt,iset) = ((dval1+dsea1)*dsea2+(sval1+ssea1)*ssea2+btm1*btm2)/rsf
                  uubmlum(jcvt,iset) = uubplum(jcvt,iset)
                  ddbmlum(jcvt,iset) = ddbplum(jcvt,iset)
               else !elseif (ptype.eq.2) then
                  uval2 = retpdfuv(y2,iset)
                  dval2 = retpdfdv(y2,iset)
                  sval2 = retpdfsv(y2,iset)
                  uubplum(jcvt,iset) = ((uval1+usea1)*(uval2+usea2)+chm1*chm2)/rsf
                  ddbplum(jcvt,iset) = ((dval1+dsea1)*(dval2+dsea2)+(sval1+ssea1)*(sval2+ssea2)+btm1*btm2)/rsf
                  uubmlum(jcvt,iset) = (usea1*usea2+chm1*chm2)/rsf
                  ddbmlum(jcvt,iset) = (dsea1*dsea2+ssea1*ssea2+btm1*btm2)/rsf
               endif
               if (jcvt.le.1) then ! separately store b PDF luminosity for only the case w/ no convolution b/c it's  only used in weak corrections
                  bbbplum(jcvt,iset) = btm1*btm2/rsf
                  bbbmlum(jcvt,iset) = bbbplum(jcvt,iset)
               endif
            elseif ((channel.eq.2).and.(jcvt.le.1)) then ! photon PDF luminosity is only needed for the case with no convolution
               if (dabs(fitscale-q).le.0.5d0) then
                  phot1 = retpdfph(y1,iset)
                  phot2 = retpdfph(y2,iset)
               else
                  CALL ObtainPDFs(y1,q,iset,pdfval,1)
                  phot1 = pdfval(10) !ph
                  CALL ObtainPDFs(y2,q,iset,pdfval,1)
                  phot2 = pdfval(10) !ph
               endif
               phphplum(jcvt,iset) = phot1*phot2/rsf
               phphmlum(jcvt,iset) = phphplum(jcvt,iset)
            endif
         enddo
      enddo

      if (mufflag.gt.0d0) muf = max(qmin,muf)
C      if (murflag.gt.0d0) mur = murflag*max(qmin/mufflag,dsqrt(m2)) ! don't care about mur here since it's not used for EW corrections 

      RETURN
      END

CCC -------------------------------------------------------------------

      SUBROUTINE InitEWSeq()
      IMPLICIT NONE
      include 'constcommon_ew.f'

      DOUBLE COMPLEX IM
      PARAMETER (IM=(0.D0,1.D0))
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      DATA dels/0.1D1/
      DATA delc/0.1D1/
      DATA Ieps/(0.D0,0.1D-59)/

CCC need to make dels,delc,Ieps internal common block

      CALL TransVar(alf,Ml,Mw,GAMw,Mz,GAMz,Mh,Mt,
     &              geu,ged,gel,gvu,gvd,gvl,gau,gad,gal,
     &              smax,smin,prej0)

      Ep = dsqrt(Scm)/2d0
      Mzsq = Mz**2-IM*GAMz*Mz
      Mwsq = Mw**2-IM*GAMw*Mw
      cw2 = Mwsq/Mzsq
      sw2 = 1d0-cw2

CCC take coupling conjugate and make combinations in advance
CCC gvu*geu,gvd*ged,gvl*gel is not our concern because they correspond to photon/Z interference terms
CCC in which case the Z diagrams are multiplied by conjugate of photon diagrams
CCC and the real part of the product is taken later
      gvusq = dble(gvu*conjg(gvu)) !replace gvu**2
      gvdsq = dble(gvd*conjg(gvd)) !replace gvd**2
      gvlsq = dble(gvl*conjg(gvl)) !replace gvl**2
      gausq = dble(gau*conjg(gau)) !replace gau**2
      gadsq = dble(gad*conjg(gad)) !replace gad**2
      galsq = dble(gal*conjg(gal)) !replace gal**2
      gvugau = dble(gvu*conjg(gau)) !replace gvu*gau
      gvdgad = dble(gvd*conjg(gad)) !replace gvd*gad
      gvlgal = dble(gvl*conjg(gal)) !replace gvl*gal

CCC The following is for contol of gauge subset: ISR, FSR, ISR*FSR
      GeuSq = geu**2
      GedSq = ged**2
      GelSq = gel**2
      GeuGel = geu*gel
      GedGel = ged*gel
CCC Whether to turn off ISR*FSR
      IF (mod(EWflag/(2**(3-1)),2).eq.1) THEN ! take the 3rd digit of the binary number
         GeuGel = 0d0
         GedGel = 0d0
      ENDIF
CCC Whether to turn off ISR
      IF (mod(EWflag/(2**(2-1)),2).eq.1) THEN ! take the 2nd digit of the binary number
         GeuSq = 0d0
         GedSq = 0d0
      ENDIF
CCC Whether to turn off FSR
      IF (mod(EWflag/(2**(1-1)),2).eq.1) THEN ! take the 1st digit of the binary number
         GelSq = 0d0
      ENDIF

CCC Turn photon mediated channel off
CCC photonic correction is still on since GeuSq,GedSq,GelSq,GeuGel,GedGel are used instead
      IF (photflag.eq.1) THEN
         geu = 0d0
         ged = 0d0
         gel = 0d0
      ENDIF

      sllmin = smin
      sllmax = smax
      sqqmin = smin
      sqqmax = Scm

CCC initialize renormalization constant for weak correction
      CALL InitConstw()

      RETURN
      END

CCC -------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION NLOEWVAL(X,N,W)
      IMPLICIT NONE

      INTEGER N
      DOUBLE PRECISION X(N),W

      include 'parameters.f'

      DOUBLE PRECISION NLOEWVAL_mls,NLOEWVAL_msv,NLOEWVAL_weak,LOEWVAL_rr
      EXTERNAL NLOEWVAL_mls,NLOEWVAL_msv,NLOEWVAL_weak,LOEWVAL_rr

      NLOEWVAL = 0d0

      IF (alphaqed0.ne.0d0) NLOEWVAL = NLOEWVAL+LOEWVAL_rr(X,N,W)

!      IF (mod(EWflag/(2**(4-1)),2).eq.0) NLOEWVAL = NLOEWVAL+NLOEWVAL_weak(X,N,W)
      IF (EWflag.lt.8) NLOEWVAL = NLOEWVAL+NLOEWVAL_weak(X,N,W)

      IF (Ml.eq.0d0) THEN
         NLOEWVAL = NLOEWVAL+NLOEWVAL_mls(X,N,W)
      ELSE
         NLOEWVAL = NLOEWVAL+NLOEWVAL_msv(X,N,W)
      ENDIF
      
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

