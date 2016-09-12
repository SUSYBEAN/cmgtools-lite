****************************************************************
*
*     initNNPDFset.f
*
*     This routine initializes the NNPDF
*     .LHgrid files 
*     Same format as LHAPDF but does not
*     require linking to the library,
*     can be  used from an external wrapper
*
*     Analogous to the LHAPDF routine initPDFset,
*     with the extra output variable for the number
*     of replicas.
*
*     Settings are the same as those
*     from the NNPDF2.0 family of  .LHgrid files
*
*     By default, the .LHgrid files are located
*     in results/projects/lhapdf/
*
************************************************************

      subroutine initNNPDFset(wrapfile)

*     Loads the array of evolved PDF grids for the
*     LHAPDF pdf interpolation

      implicit none
      
      integer MXREP,NREP
      parameter(MXREP=1e3)
      common/CNREP/NREP

      integer NPX,NPQ2,NPL,IX,IQ2
      parameter(NPX=100,NPQ2=50)
      parameter(NPL=5000)

      integer IPDF,ILINES,NXTMP,NQ2TMP
      integer NLINES,IDUM,IREP

      double precision Q2MIN,Q2MAX,XPDFMIN,XPDFMAX,XCH,Q2CH
      parameter(Q2MAX=1d8,Q2CH=4d0)
      parameter(XPDFMIN=1d-9,XPDFMAX=1d0,XCH=1D-1)
*
      double precision XG(NPX),Q2G(NPQ2)
      double precision XPDFEV(NPX,NPQ2,-6:6,0:MXREP)
      common/CPDFGR/XPDFEV,XG,Q2G,Q2MIN,IX,IQ2
      double precision mc,mb,mt,Q_ASREF
      double precision ALPHAS(0:MXREP),Q20
      common/CLHAPARS/ALPHAS,Q_ASREF,mc,mb,mt,Q20
      character*10 FLVSCHM,PTORDC
      common/CEVOLCHAR/FLVSCHM,PTORDC
     
      double precision XMINTMP,XMAXTMP,Q2MINTMP,Q2MAXTMP
      double precision ADUM, LAMBDA4, LAMBDA5

      character*150 wrapfile
      character*8 path
      character*50 char1, char2, char3

      logical VERBOSE

*     ----------------------------------------------------

*     Set verbose to false to avoid any
*     output during the initialization of the NNPDF set
      VERBOSE=.true.

*     Read the NNPDF .LHgrid file

c      write(*,*)"wrapfile=",wrapfile
      open(unit=10,status="old",file=wrapfile(1:150))

*     Read description of the .LHgrid file
      read(10,*) char1,char2
      if(char1.ne."Version") then
         write(6,*) "Problem with .LHgrid file format"
         write(6,*) "Version = ",char1, " ",char2
         call exit(-10)
      else
         if(VERBOSE) then
         write(6,*) "******************************************** "
         write(6,*) " PDFset = ",wrapfile
         write(6,*) " LHAPDF version = ",char2
         endif
      endif

      read(10,*) char1
      if(char1.ne."Description:") then
         write(6,*) "Problem with .LHgrid file format"
         write(6,*) "Description = ",char1
         call exit(-10)
      else
         NLINES=10
         do ILINES=1,NLINES
            read(10,*) char1
           if(VERBOSE) write(6,*) char1
         enddo
         write(6,*) "   --------------------------- "
      endif

      read(10,*) char1
      if(char1.ne."Alphas:") then
         write(6,*) "Problem with .LHgrid file format"
         write(6,*) "Alphas: = ",char1
         call exit(-10)
      else
         read(10,*) FLVSCHM,PTORDC,char3
         if(VERBOSE) then
            write(6,*) "Flavour scheme = ",FLVSCHM
            write(6,*) "Perturbative order = ",PTORDC
         endif
         read(10,*) idum, Q_ASREF,mc,mb,mt
         if(VERBOSE) then
            write(6,441) " Q_ASREF = ",Q_ASREF, " GeV"
            write(6,442) " mc, mb, mt [GeV] = ",mc,mb,mt
         endif
      endif

      read(10,*) char1
      if(char1.ne."MinMax:") then
         write(6,*) "Problem with .LHgrid file format"
         write(6,*) "MinMax: = ",char1
         call exit(-10)
      else
         read(10,*) idum, idum
         read(10,*) XMINTMP,XMAXTMP,Q2MINTMP,Q2MAXTMP
      endif

      read(10,*) char1
      if(char1.ne."QCDparams:") then
         write(6,*) "Problem with .LHgrid file format"
         write(6,*) "QCDparams: = ",char1
         call exit(-10)
      else
         read(10,*) NREP, idum
         read(10,*) LAMBDA4, LAMBDA5
         write(6,*) "NREP = ",NREP
      endif

      read(10,*) char1
      if(char1.ne."Parameterlist:") then
         write(6,*) "Problem with .LHgrid file format"
         write(6,*) "Parameterlist: = ",char1
         call exit(-10)
      else
         read(10,*) char1, idum, idum
         if(wrapfile(1:5).eq."toyLH") then
            IREP=0
            read(10,*) ALPHAS(IREP)
         else
            do IREP=0,NREP
               read(10,*) ALPHAS(IREP)
            enddo
         endif
      endif

      read(10,*) char1
      if(char1.ne."Evolution:") then
         write(6,*) "Problem with .LHgrid file format"
         write(6,*) "Evolution: = ",char1
         call exit(-10)
      else
         read(10,*) char1, Q20, adum
         read(10,*) char1
         read(10,*) idum, idum
         if(VERBOSE) write(6,441) " Q20 =",Q20, " GeV2"
      endif

      if(VERBOSE) then
         write(6,*) "********************************************** "
      endif

 441  format(a,f8.2,a)
 442  format(a,f12.3,1x,f12.3,1x,f12.2)
      
*     Read and check grids
*     Read the grid in x
      read(10,*) NXTMP
      if(NXTMP.ne.NPX) then
         write(6,*) "Invalid value of NX!"
         call exit(-10)
      endif
      do IX=1,NPX
         read(10,*) XG(IX)
      enddo
*     Read the grid in Q2
      read(10,*) NQ2TMP
      read(10,*) ADUM
      if(NQ2TMP.ne.NPQ2) then
         write(6,*) "Invalid value of NQ2!"
         call exit(-10)
      endif
      do IQ2=1,NPQ2
         read(10,*) Q2G(IQ2)
      enddo

*     Read the evolved xpdf grid for each replica
CCC MODIFY: REMOVE THIS CHECK TO ALLOW MY OWN FORMAT FOR THE FILE NAME
      if(wrapfile(1:5).ne."XXXXX") then
*     Read the number of replicas
      read(10,*) NREP
         do IREP=0,NREP
            do IX=1,NPX
               do IQ2=1,NPQ2
                  read(10,*) ( XPDFEV(IX,IQ2,IPDF,IREP), IPDF=-6,6,1 )
               enddo
            enddo
         enddo
      elseif(wrapfile(1:5).eq."toyLH") then
         IREP=0
         do IX=1,NPX
            do IQ2=1,NPQ2
               read(10,*) ( XPDFEV(IX,IQ2,IPDF,IREP), IPDF=-6,6,1 )
            enddo
         enddo
      else
         write(6,*) "Invalid wrapper type!"
         call exit(-10)
      endif
      close(10)

      if(VERBOSE) then
         write(6,*) "Succeful initialization of ",
     1        wrapfile
      endif
*
      return
      end


* -----------------------------------------------


      subroutine numberNNPDF(NMEM)
      implicit none

      integer MXREP,NREP,NMEM
      parameter(MXREP=1e3)
      common/CNREP/NREP

      NMEM = NREP

      if(NMEM.lt.0.or.NMEM.gt.MXREP) then
         write(6,*) "Invalid number of replicas"
         call exit(-10)
      endif

      return
      end

*     ---------------------------------------
      
