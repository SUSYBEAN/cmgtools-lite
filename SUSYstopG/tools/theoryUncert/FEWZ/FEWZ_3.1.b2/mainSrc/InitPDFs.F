
      subroutine InitPDFs(ttmpdir)
      implicit none

      character*128 ttmpdir

      include 'parameters.f'
      include 'pdfcontrol.f'

      character*128 dir
      character*150 tmpprefix
      integer i,j
      logical pdferrflag

CCC   MSTW2008 parameters
      DOUBLE PRECISION distance,tolerance,mCharm,mBottom,alphaSQ0,
     & alphaSMZ
      INTEGER alphaSorder,alphaSnfmax
      COMMON/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
      double precision xg, GetOnePDF, ALPHAS, SSC
      external GetOnePDF, ALPHAS, SSC

CCC ABKM09 parameters
      double precision pdfs(0:8),dpdfs(0:8,25)

CCC GJR functions 
      INTEGER igjr
      double precision GJR08VFNSalphas,JR09VFNNLOalphas

CCC CT12 functions
      double precision CT12Alphas

CCC PDF grid file directory to pass for PDF calling (not used for LHAPDF)
      character*150 prefix
      common/griddir/ prefix

CCC   Stores alpha_S values for each eigenvector
      DOUBLE PRECISION asEigens(0:NSMAX)

      dir = ttmpdir
      if (NPDFerror.gt.0) then ! we use NPDFerror as flag to turn on or off PDF errors
         pdferrflag = .true.
      else
         pdferrflag = .false.
      endif

CCC  === CHANGE FOR NEW PDF'S ==========================
CCC  Value of alpha strong at the mass of the Z
CCC  These values should be consistent with the fits
CCC  of the pdf's. 

      if (pdfset .EQ. 'MRST2004QED') then

         prefix = dir(1:len_trim(dir))//"/dat/MRST/mrst2004qed" ! mrst2004qedn.dat for neutron
         asZ  =  0.1189D0
         NPDFerror = 0

      elseif (pdfset .EQ. 'MRST2006NNLO') then 

         prefix = dir(1:len_trim(dir))//"/dat/MRST/mrst2006nnlo"
         asZ  =  0.1191D0
         NPDFerror = 30

      elseif (pdfset .EQ. 'MSTW2008LO') then

         prefix = dir(1:len_trim(dir))//"/dat/Grids/mstw2008lo"
         xg = GetOnePDF(prefix,0,0.1D0,10.D0,0) ! dummy call
         CALL INITALPHAS(alphaSorder,1.D0,1.D0,alphaSQ0,
     &        mCharm,mBottom,1.D10)
         asZ = alphaSMZ
         NPDFerror = 40

      elseif (pdfset .EQ. 'MSTW2008NLO') then  

         prefix = dir(1:len_trim(dir))//"/dat/Grids/mstw2008nlo"
         xg = GetOnePDF(prefix,0,0.1D0,10.D0,0) ! dummy call
         CALL INITALPHAS(alphaSorder,1.D0,1.D0,alphaSQ0,
     &        mCharm,mBottom,1.D10)
         asZ = alphaSMZ
         NPDFerror = 40
         if (whichas.gt.0) then
            NPDFerror = 204
         endif

      elseif (pdfset .EQ. 'MSTW2008NNLO') then  

         prefix = dir(1:len_trim(dir))//"/dat/Grids/mstw2008nnlo"
         xg = GetOnePDF(prefix,0,0.1D0,10.D0,0) ! dummy call
         CALL INITALPHAS(alphaSorder,1.D0,1.D0,alphaSQ0,
     &        mCharm,mBottom,1.D10)
         asZ = alphaSMZ
         NPDFerror = 40
         if (whichas.gt.0) then
            NPDFerror = 204
         endif

      elseif (pdfset .EQ. 'CTEQ65') then

         prefix = dir(1:len_trim(dir))//"/dat/CTEQGrids/"
         asZ  =  0.118D0
         NPDFerror = 40

      elseif (pdfset .EQ. 'CTEQ66') then

         prefix = dir(1:len_trim(dir))//"/dat/CTEQGrids/"
         asZ  =  0.118D0
         NPDFerror = 44

      elseif (pdfset .EQ. 'CTEQ6L1') then

         prefix = dir(1:len_trim(dir))//"/dat/CTEQGrids/"
         asZ  =  0.130D0
         NPDFerror = 0

      elseif (pdfset .EQ. 'CTEQ10') then

         prefix = dir(1:len_trim(dir))//"/dat/CT10Grids/"
         asZ  =  0.118D0
         NPDFerror = 52

      elseif (pdfset .EQ. 'CTEQ10W') then

         prefix = dir(1:len_trim(dir))//"/dat/CT10Grids/"
         asZ  =  0.118D0
         NPDFerror = 52

      elseif (pdfset .EQ. 'CTEQ12NLO') then

         prefix = dir(1:len_trim(dir))//"/dat/CT12Grids/"
         tmpprefix = prefix(1:len_trim(prefix))//"ct10wn.00.pds"
         call SetCT12(tmpprefix)
         asZ  =  CT12Alphas(Mz)
         NPDFerror = 52

      elseif (pdfset .EQ. 'CTEQ12NNLO') then

         prefix = dir(1:len_trim(dir))//"/dat/CT12Grids/"
         tmpprefix = prefix(1:len_trim(prefix))//"ct10wnn.00.pds"
         call SetCT12(tmpprefix)
         asZ  =  CT12Alphas(Mz)
         NPDFerror = 50

      elseif (pdfset .EQ. 'ABKM09NLO') then

         prefix = dir(1:len_trim(dir))//"/dat/a09/"        
         call a09(prefix,0.1d0,Mz**2,pdfs,dpdfs,5,-1,1)
         asZ = pdfs(0)
         NPDFerror = 25

      elseif (pdfset .EQ. 'ABKM09NNLO') then    

         prefix = dir(1:len_trim(dir))//"/dat/a09/"     
         call a09(prefix,0.1d0,Mz**2,pdfs,dpdfs,5,-1,2)
         asZ = pdfs(0)        
         NPDFerror = 25

      elseif (pdfset .EQ. 'GJR08LO') then

         prefix = dir(1:len_trim(dir))//"/dat/GJR08VFNS/"
         call GJR08VFNSinit(prefix,15)
         asZ = GJR08VFNSalphas(Mz**2)
         NPDFerror = 0

      elseif (pdfset .EQ. 'GJR08NLO') then

         prefix = dir(1:len_trim(dir))//"/dat/GJR08VFNS/"
         call GJR08VFNSinit(prefix,0)
         asZ = GJR08VFNSalphas(Mz**2)
         NPDFerror = 26

      elseif (pdfset .EQ. 'JR09NNLO') then

         prefix = dir(1:len_trim(dir))//"/dat/JR09VFNS/"
         call JR09VFNNLOinit(prefix)
         asZ = JR09VFNNLOalphas(Mz**2,0)
         NPDFerror = 26

      elseif (pdfset .EQ. 'NNPDF20') then

         prefix = dir(1:len_trim(dir))//
     &            "/dat/NNPDF20/NNPDF20_100.grid"
         call InitNNPDFwrap(prefix,NPDFerror)
         asZ  =  0.119D0

      elseif (pdfset .EQ. 'NNPDF21lo') then

         prefix = dir(1:len_trim(dir))//
     &            "/dat/NNPDF21/NNPDF21_lo_as_0119_100.LHgrid"
         call initNNPDFset(prefix,NPDFerror)
         call numberNNPDF(NPDFerror)
         asZ  =  0.119D0

      elseif (pdfset .EQ. 'NNPDF21nlo') then

         prefix = dir(1:len_trim(dir))//
     &            "/dat/NNPDF21/NNPDF21_100.LHgrid"
         call initNNPDFset(prefix,NPDFerror)
         call numberNNPDF(NPDFerror)
         asZ  =  0.119D0

      elseif (pdfset .EQ. 'NNPDF21nnlo') then

         prefix = dir(1:len_trim(dir))//
     &            "/dat/NNPDF21/NNPDF21_nnlo_100.LHgrid"
         call initNNPDFset(prefix,NPDFerror)
         call numberNNPDF(NPDFerror)
         asZ  =  0.119D0

      endif

      if (pdferrflag) NPDFerror = 0 ! overide the NPDFerror to turn off PDF error

CCC in case the scale is not dynamical
      CALL ObtainAlphaS() ! to fill asRwt(i,j) and return as

      END

CCC ===============================================================================================

CCC Call with negative value to disable the check of it, i.e. CheckPDFRange(-1d0,q) will only check q
      SUBROUTINE CheckPDFRange(x,q) ! freeze the x or q at their minimum allowed value
      implicit none

      include 'pdfcontrol.f'

      double precision x,q
      double precision qmin, xmin
      integer errcount
      data errcount/-1/
      save qmin, xmin

      if (errcount.lt.0) then ! initialization
        errcount = 0
        if (pdfset.eq.'MRST2004QED') then
           qmin = 1.25d0 ! minimum scale the PDF set supports
           xmin = 1d-5 ! minimum Bjorken scale the PDF set supports 
        elseif (pdfset.eq.'MRST2006NNLO') then
           qmin = 1d0 ! minimum scale the PDF set supports
           xmin = 1d-6 ! minimum Bjorken scale the PDF set supports 
        elseif ((pdfset.eq.'MSTW2008LO').or.(pdfset.eq.'MSTW2008NLO').or.(pdfset.eq.'MSTW2008NNLO')) then
           qmin = 1d0 ! minimum scale the PDF set supports
           xmin = 1d-6 ! minimum Bjorken scale the PDF set supports 
        elseif (pdfset.eq.'CTEQ6L1') then
           qmin = 1.69d0 ! minimum scale the PDF set supports
           xmin = 1d-6 ! minimum Bjorken scale the PDF set supports 
        elseif (pdfset.eq.'CTEQ65') then
           qmin = 1.69d0 ! minimum scale the PDF set supports
           xmin = 1d-7 ! minimum Bjorken scale the PDF set supports 
        elseif ((pdfset.eq.'CTEQ66')
     &         .or.(pdfset.eq.'CTEQ10').or.(pdfset.eq.'CTEQ10W')
     &         .or.(pdfset.eq.'CTEQ12NLO').or.(pdfset.eq.'CTEQ12NNLO')) then
           qmin = 1.69d0 ! minimum scale the PDF set supports
           xmin = 1d-8 ! minimum Bjorken scale the PDF set supports 
        elseif ((pdfset.eq.'ABKM09NLO').or.(pdfset.eq.'ABKM09NNLO')) then
           qmin = 0.8d0 ! minimum scale the PDF set supports
           xmin = 1d-7 ! minimum Bjorken scale the PDF set supports 
        elseif (pdfset .EQ. 'GJR08LO') then
           qmin = 0.3d0 ! minimum scale the PDF set supports
           xmin = 1d-9 ! minimum Bjorken scale the PDF set supports 
        elseif (pdfset .EQ. 'GJR08NLO') then
           qmin = 0.5d0 ! minimum scale the PDF set supports
           xmin = 1d-9 ! minimum Bjorken scale the PDF set supports 
        elseif (pdfset .EQ. 'JR09NNLO') then
           qmin = 0.55d0 ! minimum scale the PDF set supports
           xmin = 1d-9 ! minimum Bjorken scale the PDF set supports 
        elseif ((pdfset .EQ. 'NNPDF20')
     &          .or.(pdfset .EQ. 'NNPDF21lo')
     &          .or.(pdfset .EQ. 'NNPDF21nlo')
     &          .or.(pdfset .EQ. 'NNPDF21nnlo')) then
           qmin = 2d0 ! minimum scale the PDF set supports
           xmin = 1d-9 ! minimum Bjorken scale the PDF set supports 
        endif
      endif

      if ((q.gt.0).and.(q.lt.qmin)) then
C         print *,"Warning: q =",q," for PDF is too small, "//
C     &           "possible incorrect result for small mass region"
         q = qmin
      endif
      if ((x.gt.0.).and.(x.lt.xmin)) then
         if (errcount.ge.100) then
         print *,"Warning: found x for PDF too small for 100 times, "//
     &           "possible incorrect result for large rapidity region"
         errcount = 0
         else
         errcount = errcount+1
         endif
         x = xmin
      endif

      END

CCC ===============================================================================================

      SUBROUTINE ObtainAlphaS()
      implicit none

      INCLUDE 'parameters.f'
      INCLUDE 'pdfcontrol.f'

      character*150 tmpprefix
      integer i, j, asgrid
      double precision prevmur, prevas
      data prevmur/0d0/

CCC   MSTW2008 parameters
      DOUBLE PRECISION distance,tolerance,mCharm,mBottom,alphaSQ0,
     & alphaSMZ
      INTEGER alphaSorder,alphaSnfmax
      COMMON/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
      double precision xg, GetOnePDF, ALPHAS, SSC
      external GetOnePDF, ALPHAS, SSC

CCC ABKM09 parameters
      double precision pdfs(0:8),dpdfs(0:8,25)

CCC GJR functions 
      INTEGER igjr
      double precision GJR08VFNSalphas,JR09VFNNLOalphas

CCC CT12 functions
      double precision CT12Alphas

CCC PDF grid file directory to pass for PDF calling (not used for LHAPDF)
      character*150 prefix
      common/griddir/ prefix

CCC Stores alpha_S values for each eigenvector
      DOUBLE PRECISION asEigens(0:NSMAX)

CCC Constants
      double precision pi, mufr2
      PARAMETER(pi=3.141592654d0)
      data mufr2/0d0/

CCC Saved IO
      save prevmur, prevas, mufr2

CCC ================================

      if (mur.eq.prevmur) then
         as = prevas
         return ! nothing else will touch asRwt(i,j) 
      endif

CCC ================================

      if (mufr2.le.0) then ! only run it once in the beginning
         if (((murflag*mufflag).eq.0d0).and.((murflag+mufflag).ne.0d0)) then ! only one of them is dynamical scale
            if ((pdfset.eq.'MSTW2008LO').or.(pdfset.eq.'MSTW2008NLO').or.(pdfset.eq.'MSTW2008NNLO')) then
            print *,'Warning: mu_f/mu_r must be constant for MSTW PDF'//
     &              ', automatically adjusted both to be both dynamical'
            murflag = murflag+mufflag
            mufflag = murflag+mufflag
            endif
         else
            mufr2 = muf**2/mur**2
         endif
      endif

CCC ================================

      if (pdfset .EQ. 'MRST2004QED') then
         as = SSC(mur,asZ,porder)
         asEigens(0) = as

      elseif (pdfset .EQ. 'MRST2006NNLO') then 
         as = SSC(mur,asZ,porder)
         DO i=0,NPDFerror
            asEigens(i) = as
         ENDDO

      elseif (pdfset .EQ. 'MSTW2008LO') then
CCC Use constant mufr2
C         CALL INITALPHAS(alphaSorder,mufr2,1.D0,alphaSQ0,
C     &           mCharm,mBottom,1.D10)
         as = ALPHAS(mur)/pi
         DO i=0,NPDFerror
            asEigens(i) = as
         ENDDO

      elseif (pdfset .EQ. 'MSTW2008NLO') then
         if (NPDFerror.gt.40) then ! need to consider alphaS variation
            do asgrid=0,4
               if (asgrid.eq.0) then ! choose suffix of filename for each variation of alphaS
                  suffix=""
               elseif (asgrid.eq.1) then
                  if (pdfCL.eq.1) then
                     suffix="_asmz-68cl"
                  else
                     suffix="_asmz-90cl"
                  endif
               elseif (asgrid.eq.2) then
                  if (pdfCL.eq.1) then
                     suffix="_asmz-68clhalf"
                  else
                     suffix="_asmz-90clhalf"
                  endif
               elseif (asgrid.eq.3) then
                  if (pdfCL.eq.1) then
                     suffix="_asmz+68clhalf"
                  else
                     suffix="_asmz+90clhalf"
                  endif
               elseif (asgrid.eq.4) then
                  if (pdfCL.eq.1) then
                     suffix="_asmz+68cl"
                  else
                     suffix="_asmz+90cl"
                  endif
               endif
               tmpprefix = prefix(1:len_trim(prefix))
     &              //suffix(1:len_trim(suffix))
               xg = GetOnePDF(tmpprefix,0,0.1D0,10.D0,0) ! dummy call
               CALL INITALPHAS(alphaSorder,mufr2,1.D0,alphaSQ0,
     &              mCharm,mBottom,1.D10)
               as = ALPHAS(mur)/pi
               DO i=asgrid*41,asgrid*41+40
                  asEigens(i) = as
               ENDDO
            enddo
         else ! only use the best fit alphaS
CCC Use constant mufr2
C            CALL INITALPHAS(alphaSorder,mufr2,1.D0,alphaSQ0,
C     &           mCharm,mBottom,1.D10)
            as = ALPHAS(mur)/pi
            DO i=0,NPDFerror
               asEigens(i) = as
            ENDDO
         endif
CCC Fill MSTWas(asgrid) array
C         as = asEigens(0)
C         MSTWas(0) = as
C         if (NPDFerror.gt.40) then
C         do asgrid=1,4
C            MSTWas(asgrid) = asEigens(asgrid*41)
C         enddo
C         endif

      elseif (pdfset .EQ. 'MSTW2008NNLO') then
         if (NPDFerror.gt.40) then ! need to consider alphaS variation
            do asgrid=0,4
               if (asgrid.eq.0) then ! choose suffix of filename for each variation of alphaS
                  suffix=""
               elseif (asgrid.eq.1) then
                  if (pdfCL.eq.1) then
                     suffix="_asmz-68cl"
                  else
                     suffix="_asmz-90cl"
                  endif
               elseif (asgrid.eq.2) then
                  if (pdfCL.eq.1) then
                     suffix="_asmz-68clhalf"
                  else
                     suffix="_asmz-90clhalf"
                  endif
               elseif (asgrid.eq.3) then
                  if (pdfCL.eq.1) then
                     suffix="_asmz+68clhalf"
                  else
                     suffix="_asmz+90clhalf"
                  endif
               elseif (asgrid.eq.4) then
                  if (pdfCL.eq.1) then
                     suffix="_asmz+68cl"
                  else
                     suffix="_asmz+90cl"
                  endif
               endif
               tmpprefix = prefix(1:len_trim(prefix))
     &              //suffix(1:len_trim(suffix))
               xg = GetOnePDF(tmpprefix,0,0.1D0,10.D0,0) ! dummy call
               CALL INITALPHAS(alphaSorder,mufr2,1.D0,alphaSQ0,
     &              mCharm,mBottom,1.D10)
               as = ALPHAS(mur)/pi
               DO i=asgrid*41,asgrid*41+40
                  asEigens(i) = as
               ENDDO
            enddo
         else ! only use the best fit alphaS
CCC Use constant mufr2
C            CALL INITALPHAS(alphaSorder,mufr2,1.D0,alphaSQ0,
C     &           mCharm,mBottom,1.D10)
            as = ALPHAS(mur)/pi
            DO i=0,NPDFerror
               asEigens(i) = as
            ENDDO
         endif
CCC Fill MSTWas(asgrid) array
C         as = asEigens(0)
C         MSTWas(0) = as
C         if (NPDFerror.gt.40) then
C         do asgrid=1,4
C            MSTWas(asgrid) = asEigens(asgrid*41)
C         enddo
C         endif

      elseif (pdfset .EQ. 'CTEQ65') then
         as = SSC(mur,asZ,porder)
         DO i=0,NPDFerror
            asEigens(i) = as
         ENDDO

      elseif (pdfset .EQ. 'CTEQ66') then
         as = SSC(mur,asZ,porder)
         DO i=0,NPDFerror
            asEigens(i) = as
         ENDDO

      elseif (pdfset .EQ. 'CTEQ6L1') then
         as = SSC(mur,asZ,porder)
         asEigens(0) = as

      elseif (pdfset .EQ. 'CTEQ10') then
         as = SSC(mur,asZ,porder)
         DO i=0,NPDFerror
            asEigens(i) = as
         ENDDO

      elseif (pdfset .EQ. 'CTEQ10W') then
         as = SSC(mur,asZ,porder)
         DO i=0,NPDFerror
            asEigens(i) = as
         ENDDO

      elseif (pdfset .EQ. 'CTEQ12NLO') then
         as = CT12Alphas(mur)/pi
         DO i=0,NPDFerror
            asEigens(i) = as
         ENDDO

      elseif (pdfset .EQ. 'CTEQ12NNLO') then
         as = CT12Alphas(mur)/pi
         DO i=0,NPDFerror
            asEigens(i) = as
         ENDDO

      elseif (pdfset .EQ. 'ABKM09NLO') then
         call a09(prefix,0.1d0,mur**2,pdfs,dpdfs,5,-1,1)
         as = pdfs(0)/pi
         do i=0,NPDFerror
            asEigens(i) = pdfs(0)+dpdfs(0,i)
         enddo  

      elseif (pdfset .EQ. 'ABKM09NNLO') then    
         call a09(prefix,0.1d0,mur**2,pdfs,dpdfs,5,-1,2)
         as = pdfs(0)/pi         
         do i=0,NPDFerror
            asEigens(i) = pdfs(0)+dpdfs(0,i)
         enddo

      elseif (pdfset .EQ. 'GJR08LO') then
         as = GJR08VFNSalphas(mur**2)/pi
         NPDFerror = 0
         asEigens(0) = as
c         GJRas(0) = as

      elseif (pdfset .EQ. 'GJR08NLO') then
         call GJR08VFNSinit(prefix,0)
         as = GJR08VFNSalphas(mur**2)/pi
         asEigens(0) = as
c         GJRas(0) = as
         do i=1,NPDFerror
            if ((i/2*2).eq.i) then
               igjr=-i/2
               call GJR08VFNSinit(prefix,igjr)
               asEigens(i) = GJR08VFNSalphas(mur**2)/pi
            else
               igjr=(i+1)/2
               call GJR08VFNSinit(prefix,igjr)
               asEigens(i) = GJR08VFNSalphas(mur**2)/pi
            endif
         enddo

      elseif (pdfset .EQ. 'JR09NNLO') then
         as = JR09VFNNLOalphas(mur**2,0)/pi
         asEigens(0) = as
c         GJRas(0) = as
         do i=1,NPDFerror
            if ((i/2*2).eq.i) then
               igjr=-i/2
            else
               igjr=(i+1)/2
            endif
            asEigens(i) = JR09VFNNLOalphas(mur**2,igjr)/pi
         enddo

      elseif (pdfset .EQ. 'NNPDF20') then
         as = SSC(mur,asZ,porder)
         DO i=0,NPDFerror
            asEigens(i) = as
         ENDDO

      elseif (pdfset .EQ. 'NNPDF21lo') then
         as = SSC(mur,asZ,porder)
         DO i=0,NPDFerror
            asEigens(i) = as
         ENDDO

      elseif (pdfset .EQ. 'NNPDF21nlo') then
         as = SSC(mur,asZ,porder)
         DO i=0,NPDFerror
            asEigens(i) = as
         ENDDO

      elseif (pdfset .EQ. 'NNPDF21nnlo') then
         as = SSC(mur,asZ,porder)
         DO i=0,NPDFerror
            asEigens(i) = as
         ENDDO

      endif

      DO i=0,2 ! LO to NNLO reweightings
         DO j=1,NPDFerror
            asRwt(j,i) = (asEigens(j)/asEigens(0))**i
         ENDDO
      ENDDO

      prevmur = mur
      prevas = as

      END
