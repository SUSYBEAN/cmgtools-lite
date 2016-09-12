      subroutine ObtainPDFs(xval,qval,iset,pdfval,onlyphtflag) ! onlyphtflag is for saving time if only need photon PDF
      implicit none
      integer i, iset, onlyphtflag
      double precision x, q, xval, qval, pdfval(1:10)
      integer previset
      double precision prevx, prevq, prevpdfval(1:10)
      data previset/-1/
      data prevx/0d0/
      data prevq/0d0/
      data prevpdfval/10*0d0/

      include 'pdfcontrol.f'

CCC For LHAPDF
      logical has_photon
      logical photyes
      double precision f(-6:6),f2

CCC Variables to save for the next time
      save prevx, prevq, previset, prevpdfval

CCC     See if we only ask for the photon PDF ===========================================

      photyes = has_photon()
      if ((onlyphtflag.gt.0).and.(.not.photyes)) then
         do i=1,10
            pdfval(i) = 0d0
         enddo
         return
      endif

CCC     Return 0 for irreasonable x =======================================

      if ((xval.le.0).or.(xval.ge.1)) then
          do i=1,10
             pdfval(i) = 0d0
          enddo
          return
      endif
      x = xval
      q = qval
      call CheckPDFRange(x,q)

CCC     Save the effort of running it again if the last run already yielded the same result =======

      if ((xval.eq.prevx).and.(qval.eq.prevq).and.(iset.eq.previset)) then
         do i=1,10
            pdfval(i) = prevpdfval(i)
         enddo
         return
      endif

CCC     Do not re-initialize if it is the same PDF eigenvector as the last time ==============

      if (iset.ne.previset) call InitPDF(iset)

CCCC    Call PDF set ===============================================

      if (photyes) then
          CALL evolvePDFphoton(x,q,f,f2)
      else
          CALL evolvePDF(x,q,f)
          f2 = 0d0
      endif

      pdfval(1) = f(2) - f(-2) !uval
      pdfval(2) = f(1) - f(-1) !dval
      pdfval(3) = f(-2) !usea
      pdfval(4) = f(-1) !dsea
      pdfval(5) = f(3) - f(-3) !sval
      pdfval(6) = f(-3) !ssea
      pdfval(7) = f(4) !chm
      pdfval(8) = f(5) !bot
      pdfval(9) = f(0) !glu
      pdfval(10) = f2 !phot

      prevx = xval
      prevq = qval
      previset = iset
      do i=1,10
         prevpdfval(i) = pdfval(i)
      enddo

      return
      end
