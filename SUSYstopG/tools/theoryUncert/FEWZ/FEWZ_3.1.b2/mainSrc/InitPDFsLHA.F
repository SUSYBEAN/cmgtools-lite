      SUBROUTINE InitPDFs(ttmpdir)
      IMPLICIT NONE

      character*128 ttmpdir ! dummy argument, unneeded for LHA

      INCLUDE 'parameters.f'
      INCLUDE 'pdfcontrol.f'

      INTEGER i,j

      DOUBLE PRECISION alphasPDF
      EXTERNAL alphasPDF ! LHA alphas return

      CALL InitPDFsetByName(pdfset) ! initializes LHA pdf set with this name
      if (NPDFerror.gt.0) then ! we use NPDFerror as flag to turn on or off PDF errors
         NPDFerror = 0
      else
         CALL numberPDF(NPDFerror) ! gets number of eigenvectors in this set
      endif
      CALL InitPDF(0) ! loads best fit member
      asZ = alphasPDF(Mz) ! gets central value of alpha_S at Mz

      CALL ObtainAlphaS() ! to fill asRwt(i,j) and return as

      RETURN
      END

CCC ===============================================================================================

CCC Call with negative value to disable the check of it, i.e. CheckPDFRange(-1d0,q) will only check q
      SUBROUTINE CheckPDFRange(x,q) ! freeze the x or q at their minimum allowed value
      implicit none

      include 'pdfcontrol.f'

      double precision x,q
      integer i
      double precision tmpval, qmin, xmin, qmax, xmax
      integer errcount
      data errcount/-1/
      save qmin, xmin, qmax, xmax ! not using qmax and xmax for now

      if (errcount.lt.0) then
         errcount = 0
         call GetXmax(0,xmax)
         call GetXmin(0,xmin)
         call GetQ2max(0,qmax)
         call GetQ2min(0,qmin)
         qmax = dsqrt(qmax)
         qmin = dsqrt(qmin)
         do i=1,NPDFerror
            call GetXmax(i,tmpval)
            xmax = min(xmax,tmpval)
            call GetXmin(i,tmpval)
            xmin = max(xmin,tmpval)
            call GetQ2max(i,tmpval)
            qmax = min(qmax,dsqrt(tmpval))
            call GetQ2min(i,tmpval)
            qmin = max(qmin,dsqrt(tmpval))
         enddo
      endif

      if ((q.gt.0).and.(q.lt.qmin)) then
C         print *,"Warning: q =",q,"is too small for PDF calling, "//
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

      DOUBLE PRECISION pi
      PARAMETER(pi=3.141592654d0)

      INTEGER i,j

CCC   Stores alpha_S values for each eigenvector
      DOUBLE PRECISION asEigens(NSMAX)

      DOUBLE PRECISION alphasPDF
      EXTERNAL alphasPDF ! LHA alphas return

      double precision prevmur, prevas
      data prevmur/0d0/
      save prevmur, prevas

      if (mur.eq.prevmur) then
         as = prevas
         return ! nothing else will touch asRwt(i,j) 
      endif

      CALL InitPDF(0) ! loads best fit member
      as = alphasPDF(mur)/pi ! gets central value of alpha_S

      DO i=1,NPDFerror
         CALL InitPDF(i)
         asEigens(i) = alphasPDF(muf)/pi
         DO j=0,2
            asRwt(i,j) = (asEigens(i)/as)**j
         ENDDO
      ENDDO

      prevmur = mur
      prevas = as

      RETURN
      END
