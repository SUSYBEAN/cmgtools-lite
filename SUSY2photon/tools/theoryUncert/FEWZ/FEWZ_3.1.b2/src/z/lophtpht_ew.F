      DOUBLE PRECISION FUNCTION LOEWVAL_rr(X,N,W)
      IMPLICIT NONE

      include 'constcommon_ew.f'
      include 'eventcommon_ew.f'
      include 'pdfcontrol.f'
      include 'particles.f'

      INTEGER N
      DOUBLE PRECISION X(N),W

      DOUBLE PRECISION Fi,FVec(NSMAX)
      DOUBLE PRECISION CmFct,Fsum,Pcut,Mcut
      LOGICAL DoneOnce
      INTEGER iset

      DOUBLE PRECISION constraint
      DOUBLE PRECISION CS0z_rr

      Fsum = 0d0

CCC 2->2 process (LO, photon+photon->lepton+antilepton)

      CALL EventSeq2to2_rr(X) ! generate events and prepare pdf functions
      CALL PreSeq2to2_rr(CmFct) ! generate common factor (contain phase space, flux and jacobian factors)

      CALL InfoSeq2to2_rr(part) ! translate to particle 4-momentum information
      Pcut = constraint(part) ! test cuts
      IF (Pcut.ne.0d0) THEN
         CALL CrsecSeqz_rr() ! calculate component of tree amplitude
         Fi = CS0z_rr(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0z_rr(1,iset)
               FVec(iset) = FVec(iset)*CmFct*W
            ENDDO
         ENDIF
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      CALL ReversePart(part) ! convert to parity reversed event for quark and antiquark from the other hadrons
      Mcut = constraint(part) ! test cuts
      IF (Mcut.ne.0d0) THEN
         if (Pcut.eq.0d0) CALL CrsecSeqz_rr() ! calculate component of tree amplitude only if not calculated before
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! ! calculate it if pp-bar collider (for pp collider it is exact the same result) or not done before
            Fi = CS0z_rr(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0z_rr(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif ! always safe to remove the 'if' and run them every time because it is just a hack for efficiency in the case of pp collider
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

      LOEWVAL_rr = Fsum

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE EventSeq2to2_rr(X)
      IMPLICIT NONE
      DOUBLE PRECISION X(4)
      include 'constcommon_ew.f'
      include 'eventcommon_ew.f'
      LOGICAL CIsBad

      lam1 = X(1)
      lam2 = X(2)
      lam3 = X(3)
C      lam4 = X(4)

      phi = 0
      s = 1d0/dsqrt((1d0/smax**2-1d0/smin**2)*lam1+1d0/smin**2)
      limthy = (Scm-s)/(Scm+s)
      Exp2Y = 2d0/(1d0-(2d0*lam2-1d0)*limthy)-1d0
      BtQ = dsqrt(1d0-4d0*Ml**2/s)
      if (Ml.eq.0d0) then
         CsTh = 2d0*lam3-1d0
         BCsTh = BtQ*CsTh
      else
         if (dble(Ml**2/s).lt.1d-12) then
            AthClc = (1d0+BtQ)/Ml**2*s/2d0
         else
            AthClc = (1d0+BtQ)/(1d0-BtQ)
         end if
         if ((dble(AthClc).eq.0d0).or.(CIsBad(AthClc)))
     #      AthClc = (1+BtQ)/Ml**2*s/2d0
         AthClc = dlog(AthClc)/2d0
         BCsTh = dtanh(AthClc*(2d0*lam3-1))
         CsTh = BCsTh/BtQ
      endif
      BSnTh = BtQ*dsqrt(1d0-CsTh**2)
      x1 = dsqrt(s/Scm*Exp2Y)
      x2 = s/Scm/x1
      Eq = dsqrt(s)/2d0
      El = Eq
C      cz1 = (1-x1)*lam4**2
C      cz2 = (1-x2)*lam4**2

      m2 = s
      CALL PDFSeqAll(x1,x2,1d0,1d0,1,2)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE InfoSeq2to2_rr(part)
      IMPLICIT NONE
      include 'constcommon_ew.f'
      include 'eventcommon_ew.f'
      include 'particles.f'

      INTEGER i,j
      DOUBLE PRECISION t1,t2,t3,t10,t20

      DO i=4,NPART
      DO j=1,4
         part(i,j) = 0d0
      ENDDO
      ENDDO

      t1 = x1+x2
      t2 = x1-x2
      t3 = Ep/2d0 !Ep*El/Eq/2d0
      t10 = t1*BCsTh
      t20 = t2*BCsTh
      part(2,1) = (t1+t20)*t3
      part(2,2) = 0d0 !sin(DBLE(phi))*El*BSnTh
      part(2,3) = El*BSnTh !cos(DBLE(phi))*El*BSnTh
      part(2,4) = (t2+t10)*t3
      part(3,1) = (t1-t20)*t3
      part(3,2) = -part(2,2)
      part(3,3) = -part(2,3)
      part(3,4) = (t2-t10)*t3

      DO j=1,4
         part(1,j) = part(2,j) + part(3,j)
      ENDDO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE PreSeq2to2_rr(CmFct)
      IMPLICIT NONE
      include 'constcommon_ew.f'
      include 'eventcommon_ew.f'

      DOUBLE PRECISION CmFct
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      j1 = (s**3)*(1d0/smin**2-1d0/smax**2)/2d0
      j2 = 2d0*limthy/(1d0-((Exp2Y-1d0)/(Exp2Y+1d0))**2d0)
      j1n2 = 1d0/Scm
      if (Ml.eq.0d0) then
         j3 = 2d0
      else
         j3 = (1d0-BCsTh**2)*AthClc*2d0/BtQ
      endif
      flxfct = 1d0/s/2d0
      PSfct = BtQ/PI/16d0
      j0 = prej0*flxfct*PSfct*j1n2*j1*j2*j3

CCC no spin and color average for incoming photons: *4*9
CCC but polarization average and identical bosonic particle factor 1/2!: /4/2
CCC adjust the factor here: *9/2
      CmFct = j0*9d0/2d0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CrsecSeqz_rr()
      IMPLICIT NONE
      include 'constcommon_ew.f'
      include 'eventcommon_ew.f'
      DOUBLE PRECISION Arr 
      COMMON/CRSrr/ Arr ! common block shared between CS0z_rr() and CrsecSeqz_rr()

      DOUBLE PRECISION t1,t2,t3

      t1 = Ml**2/s
      t2 = 1d0-BCsTh**2
      t3 = 1d0+BCsTh**2
      Arr = -16d0*(32d0*t1**2-8d0*t1*t2-t2*t3)/t2**2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION FUNCTION CS0z_rr(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION rlum
      include 'parameters.f' ! directly using parameters from the main common block of FEWZ
      include 'pdflumcommon_ew.f' ! directly using PDF luminosity variable from the common block for EW corr. module in FEWZ
      DOUBLE PRECISION Arr 
      COMMON/CRSrr/ Arr ! common block shared between CS0z_rr() and CrsecSeqz_rr()

      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      if (whichlum.eq.1) then
         rlum = phphplum(1,iset)
      else
         rlum = phphmlum(1,iset)
      endif
      CS0z_rr = rlum*Arr * (16d0*alphaqed0**2*PI**2) !gel**4

      RETURN
      END

