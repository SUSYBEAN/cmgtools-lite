      double precision function Fjet(var,ndim,z,uxa,s13,s23,s14,
     & s24,s34,s1v,s2v,s3v,s4v,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,
     & s3l2,s4l2,lum,sc1,sc2,Tcoeff)
c      implicit double precision(a-h,o-z)
      IMPLICIT NONE

      include 'parameters.f'
      include 'pdfcontrol.f'
      INCLUDE 'particles.f'

CCCCC Kinematics and event weight
      double precision Tcoeff,cweight,fjetjac
      double precision z,s13,s23,s14,s24,s34,s1v,s2v,s3v,
     & s4v,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,
     & s3l2,s4l2,uxa,x1cons,x2cons,sc1,sc2

CCCCC Required VEGAS parameters
      INTEGER ndim
      double precision var(ndim)

CCCCC External functions
      double precision bookPart,lumin
      INTEGER findLumin
      external bookPart,lumin,findLumin

      double precision fpflag
      integer iset,i

CCC   Array holding the result for each PDF set

      double precision pdfvec
      dimension pdfvec(NSmax)

      double precision x1,x2,lum,ttcoeff,tau2z

CCCCC Integration region parameters
      
      double precision BWjac,tau2
      COMMON /BW/ BWjac,tau2
      SAVE /BW/

      Fjet = 0.0d0

      IF(Tcoeff.ne.0d0) THEN ! don't bother if we're already zero
      
CCC   Cut on the z value

      if (z.gt.tau2) then

         tau2z = tau2/z
         
CCC   Note that in the qg case, I call this procedure with the argument
CCC   1/ux.  This really means flip ux, not var(6), as I do in the
CCC   argument.  Fix that here.  Note that uxa=var(6)<1; if it's
CCC   greater, that means I called it with 1/ux

CCC   Note that I'm generating the variable x1

CCC   For the x1,x2 that enter constraints, multiply by sc1,sc2, the 
CCC   convolution scaleing variables

C         uxmax = 1.0d0
C         uxmin = tau2/z
         
CCC   fpflag denotes that we must make parton 1 the pbar; for Tevatron

         fpflag = -999.0d0
C         uxjac = uxmax-uxmin
         if (uxa.le.(1.0d0)) then
C            ux = uxjac*dabs(uxa)+uxmin  
C            ux = uxa
C            uxold = ux
C            x1 = ux
            x1 = (tau2z)**dabs(1d0-uxa)
            x2 = tau2z/x1
            x1cons = x1*sc1
            x2cons = x2*sc2
         else
            fpflag = 999.0d0
C            ux = uxjac/uxa+uxmin
C            uxold = ux
C            ux = 1.0d0/ux
C            x1 = uxold
            x1 = (tau2z)**dabs(1d0-1/uxa)
            x2 = tau2z/x1
            x1cons = x2*sc2
            x2cons = x1*sc1
         endif
         
CCC   There is an x1*x2 absorbed by the PDF function calls

C         fjetjac = 1.0d0/z/x1*uxjac*BWjac
         fjetjac = -dlog(tau2z)/z*BWjac

         i = findLumin(x1,x2,fpflag,lum)
         
CCC   Call the luminosity function for central PDF value

         Fjet = fjetjac*lumins(0,i)
         cweight = Fjet*Tcoeff

CCC   Call the remaining PDF eigenvectors and store their weight in an array

         ttcoeff = Tcoeff*fjetjac
         do iset = 1,NPDFerror
            pdfvec(iset) = asRwt(iset,CoeffFuncOrder)*ttcoeff
     &                     *lumins(iset,i)
         enddo

         if ((ptype.eq.1).or.(ptype.eq.2)) then
CCC	 W-, regular lepton/neutrino order

         call kinvar2partmom(x1cons,x2cons,s13,
     &    s23,s14,s24,s34,s1v,s2v,s3v,s4v,s1l1,s2l1,s1l2,s2l2,
     &    s3l2,s4l2,part)

         else ! W+, flip lepton and neutrino in constraint call

         call kinvar2partmom(x1cons,x2cons,s13,
     &    s23,s14,s24,s34,s1v,s2v,s3v,s4v,s1l2,s2l2,s1l1,s2l1,
     &    s3l1,s4l1,part)

         endif

CCC      Throw 4-momenta and weights to constraint/histograms

         Fjet = Fjet*bookPart(part,cweight,pdfvec)

      endif

      endif

CCC   Call the histogramming routine

c      IF (Fjet.ne.0d0) THEN ! only call if nonzero weight (this procedure is expensive)
c         call fill_hist(part,cweight,pdfvec)
c      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
      END
