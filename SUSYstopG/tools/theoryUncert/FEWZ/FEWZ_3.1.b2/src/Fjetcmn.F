CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Purges Fjet buffer to fill histograms and get ready for next eval

      SUBROUTINE purgeFJBuf()
      IMPLICIT NONE
      INCLUDE 'fjetbuf.f'

      DATA bufPos /0/ ! make sure it's zero at the first time it's used

      INTEGER i,j,k
      DO i=1,bufPos
         DO j=1,NPART
            DO k=1,4
               part(j,k) = partBuf(j,k,i)
            ENDDO
c            write(*,*) part(j,1), part(j,2), part(j,3), part(j,4)
         ENDDO
c         write(*,*) weightBuf(0,i)
         DO j=1,NPDFerror
            newPdfVec(j) = weightBuf(j,i)
         ENDDO
         IF (weightBuf(0,i).ne.0d0) THEN
C            CALL doJetAlgo(part)
            if (partflgBuf(i).gt.0) CALL modifyPart(part,partflgBuf(i))
            CALL fill_hist(part,weightBuf(0,i),newPdfVec)
         ENDIF
      ENDDO

      bufPos = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Perform jet algorithm on particle 4-momentum based on the flag

      SUBROUTINE modifyPart(part,flag)
      IMPLICIT NONE
      INCLUDE 'particles.f'

      INTEGER flag
      INTEGER j

CCC   flag: 1. merge two partons; 2. set parton 1 to be zero, and parton 2 to be new parton 1 (the hard jet);
CCC         3. set parton 2 to zero; 4. set both partons to be zero

      DO j=1,4
         IF (flag.eq.1) THEN
            part(4,j) = part(4,j) + part(5,j) ! merge partons
         ELSEIF (flag.eq.2) THEN
            part(4,j) = part(5,j) ! set first parton down beam
         ELSEIF (flag.eq.4) THEN
            part(4,j) = 0d0
         ENDIF
         part(5,j) = 0d0 ! parton 2 will be set down beam anyway
      ENDDO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Returns result of constraint for inputs; calls histogramming

      DOUBLE PRECISION FUNCTION bookPart(part,cweight,pdfvec)
      IMPLICIT NONE
      INCLUDE 'fjetbuf.f'

      DOUBLE PRECISION cweight,pdfvec(NPDFerror)

      DOUBLE PRECISION constraint
      EXTERNAL constraint

      DOUBLE PRECISION oldPart(NPART,4)
      INTEGER i,j,k
      LOGICAL found,jaflag

      found = .false.
      i = 0

      DO WHILE (.not.found.and.i.lt.bufPos)
         i = i+1
         j = 2 ! Z momentum is redundant, start with leptons
         k = 1
         DO WHILE (part(j,k).eq.partBuf(j,k,i).and..not.found)
            k = k+1
            IF (k.gt.4) THEN
               j = j+1
               IF (j.gt.NPART) THEN ! finished checking this set of particles, matched
                  found = .true.
               ELSE
                  k = 1 
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      IF (found) THEN ! we found a match, add this weight to the buffers
         IF (constrBuf(i)) THEN ! this guy passed constraints when we checked
            weightBuf(0,i) = weightBuf(0,i) + cweight
            DO j=1,NPDFerror
               weightBuf(j,i) = weightBuf(j,i) + pdfvec(j)
            ENDDO
            bookPart = 1d0
         ELSE ! failed constraint, return zero to Fjet
            bookPart = 0d0 
         ENDIF
      ELSE ! not found, fill buffer
         bufPos = bufPos+1
         IF (bufPos.gt.NFJETBUF) THEN ! ran out of storage, push first guy out
            IF (weightBuf(0,1).ne.0d0) THEN ! first stored parton has histograms to fill
               DO k=1,4
                  DO j=1,NPART
                     oldPart(j,k) = partBuf(j,k,1)
                  ENDDO
               ENDDO
               DO j=1,NPDFerror
                  newPdfVec(j) = weightBuf(j,1)
               ENDDO
C               CALL doJetAlgo(oldPart)
               if (partflgBuf(1).gt.0) CALL modifyPart(oldPart,partflgBuf(1))
               CALL fill_hist(oldPart,weightBuf(0,1),newPdfVec)
            ENDIF
            DO i=1,NFJETBUF-1
               constrBuf(i) = constrBuf(i+1)
               DO j=0,NPDFerror
                  weightBuf(j,i) = weightBuf(j,i+1)
                  partflgBuf(i) = partflgBuf(i+1)
               ENDDO
               DO k=1,4
                  DO j=1,NPART
                     partBuf(j,k,i) = partBuf(j,k,i+1)
                  ENDDO
               ENDDO
            ENDDO
            bufPos = NFJETBUF
         ENDIF
CCC   Record new
         DO k=1,4
            DO j=1,NPART
               partBuf(j,k,bufPos) = part(j,k)
            ENDDO
         ENDDO
         CALL doJetAlgo(part,jaflag) ! make sure jets are in order before constraint
         partflgBuf(bufPos) = jaflag ! store the particle operation flag
         IF (constraint(part).gt.0d0) THEN ! assume everything positive passed cuts
            constrBuf(bufPos) = .true.
            weightBuf(0,bufPos) = cweight
            DO j=1,NPDFerror
               weightBuf(j,bufPos) = pdfvec(j)
            ENDDO
            bookPart = 1d0 ! return passing cuts to Fjet
         ELSE
            constrBuf(bufPos) = .false.
            DO j=0,NPDFerror
               weightBuf(j,bufPos) = 0d0
            ENDDO
            bookPart = 0d0
         ENDIF
      ENDIF

      RETURN
      END           


CCC   Gets invariant mass squared, Breit-Wigner Jacobian, stores in common blocks

      SUBROUTINE getM2(mass, width, v)
      IMPLICIT NONE
      DOUBLE PRECISION mass, width, v
      INCLUDE 'parameters.f'
      double precision s0,sinf,y0,yinf,yjac,yvar,BWjac,tau2
      COMMON /BW/ BWjac,tau2
      SAVE /BW/

CCC   Define the parameters needed for integration
CCC   over the Breit-Wigner shape

CCC   Limits of m2 integration

      s0 = mlower**2/mass**2
      sinf = mupper**2/mass**2

      if (zpole.eq.1) then

         y0 = mass/width*datan(mass/width*(s0-1.0d0))
         yinf = mass/width*datan(mass/width*(sinf-1.0d0))
         yjac = yinf-y0
         yvar = yjac*dabs(v)+y0
         m2 = mass**2*(1+width/mass*dtan(width/mass*yvar))
         BWjac = yjac/mass**2

      else

         if ((s0.le.1d-6).or.(mass.eq.Mw)) then

            yjac = sinf-s0
            yvar = yjac*v+s0
            m2 = yvar*mass**2
            BWjac = yjac/mass**2/((yvar-1)**2+width**2/mass**2)

         else

C     No Z pole focus perform variable transformation to adapt to behavior of 1/m2**3
C     The first 1/m2 is from luminosity factor and the second one is from partonic cross section
            yinf = 1d0/sinf**2
            y0 = 1d0/s0**2
            yjac = yinf-y0
            yvar = yjac*v+y0
            m2 = mass**2/dsqrt(yvar)
            yjac = m2**3*(y0-yinf)/mass**4/2d0
            BWjac = yjac/((m2-mass**2)**2+width**2*mass**2)

         endif

      endif

      tau2 = m2/scm
      CALL getmuFR()

      RETURN
      END

CCC   Gets invariant mass squared, Breit-Wigner Jacobian, stores in common blocks

      SUBROUTINE getmuFR()
      IMPLICIT NONE
      DOUBLE PRECISION q
      INCLUDE 'parameters.f'

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
      endif

C     Use dynamical scale if flagged
      if (qfct.gt.0d0) then
         q = max(qmin,dsqrt(m2))
         if (mufflag.gt.0d0) then
            muf = q*mufflag
         endif
         if (murflag.gt.0d0) then
            mur = q*murflag
            CALL ObtainAlphaS()
         endif
      endif

      RETURN
      END

