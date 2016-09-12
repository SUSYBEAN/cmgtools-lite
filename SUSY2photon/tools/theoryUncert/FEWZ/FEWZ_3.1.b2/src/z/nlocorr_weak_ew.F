      DOUBLE PRECISION FUNCTION NLOEWVAL_weak(X,N,W)
      IMPLICIT NONE

      include 'constcommon_ew.f'
      include 'eventcommon_ew.f'
      include 'pdfcontrol.f'
      include 'particles.f'

      INTEGER N
      DOUBLE PRECISION X(N),W

      DOUBLE PRECISION Fi,FVec(NSMAX)
      DOUBLE PRECISION CmFct,Fsum,Pcut,Mcut
      INTEGER iset

      DOUBLE PRECISION constraint
      DOUBLE PRECISION CS0w

      NLOEWVAL_weak = 0d0
      Fsum = 0d0

CCC Perturbation order

      IF (porder_ew.le.0) RETURN

CCC 2->2 process

      CALL EventSeq2to2(X)
      CALL PreSeq2to2(CmFct)

      CALL InfoSeq2to2(1,part)
      Pcut = constraint(part)
      IF (Pcut.ne.0d0) THEN
         CALL CrsecSeqw() ! calculate component of virtual weak amplitude
         Fi = CS0w(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0w(1,iset)
               FVec(iset) = FVec(iset)*CmFct*W
            ENDDO
         ENDIF
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      CALL ReversePart(part) ! convert again to parity reversed event for filling histogram
      Mcut = constraint(part)
      IF (Mcut.ne.0d0) THEN
         if (Pcut.eq.0d0) then
            CALL CrsecSeqw()
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then
            Fi = CS0w(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0w(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

      NLOEWVAL_weak = Fsum

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE EventSeq2to2(X)

      IMPLICIT DOUBLE PRECISION (t, s)

      DOUBLE PRECISION X(4)

      include 'constcommon_ew.f'
      include 'eventcommon_ew.f'

      lam1 = X(1)
      lam2 = X(2)
      lam3 = X(3)
      lam4 = X(4)

      phi = 0
      t1 = smax**2
      t3 = smin**2
      t4 = 1/t3
      t8 = sqrt((1/t1-t4)*lam1+t4)
      s = 1/t8
      limthy = (Scm-s)/(Scm+s)
      Exp2Y = 2/(1-(2*lam2-1)*limthy)-1
      CsTh = 2*lam3-1
      t19 = Ml**2
      BtQ = sqrt(1-4*t19/s)
      BCsTh = BtQ*CsTh
      t24 = CsTh**2
      SnTh = sqrt(1-t24)
      BSnTh = BtQ*SnTh
      t28 = s/Scm
      x1 = sqrt(t28*Exp2Y)
      x2 = t28/x1
      t31 = sqrt(s)
      Eq = t31/2
      El = Eq
      t33 = lam4**2
      cz1 = (1-x1)*t33
      cz2 = (1-x2)*t33
      j4 = 2*lam4

      m2 = s
      IF (mod(EWflag/2,2).eq.1) THEN ! don't call PDF with convolution if QED ISR is off
      CALL PDFSeqAll(x1,x2,1d0-cz1,1d0-cz2,1,1)
      ELSE
      CALL PDFSeqAll(x1,x2,1d0-cz1,1d0-cz2,3,1)
      ENDIF

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE InfoSeq2to2(PARITY,part)

      IMPLICIT DOUBLE PRECISION (t, s)

      INTEGER PARITY
      DOUBLE PRECISION PARTICLE(3,0:3)

      include 'constcommon_ew.f'
      include 'eventcommon_ew.f'
      include 'particles.f'

      IF (Ml.eq.0d0) THEN ! for the projected 2->2 event from counterterms of 2->3 event
      BtQ = 1
      BCsTh = CsTh
      BSnTh = SnTh
      ENDIF

      IF (PARITY.eq.1) THEN
      t2 = (x1+x2)*Ep
      t4 = El/Eq
      t5 = t2*t4
      t7 = (x1-x2)*Ep
      t8 = t4*BCsTh
      t9 = t7*t8
      PARTICLE(1,0) = t5/2+t9/2
      t11 = sin(DBLE(phi))
      PARTICLE(1,1) = t11*El*BSnTh
      t13 = cos(DBLE(phi))
      PARTICLE(1,2) = t13*El*BSnTh
      t15 = t7*t4
      t16 = t2*t8
      PARTICLE(1,3) = t15/2+t16/2
      PARTICLE(2,0) = t5/2-t9/2
      PARTICLE(2,1) = -PARTICLE(1,1)
      PARTICLE(2,2) = -PARTICLE(1,2)
      PARTICLE(2,3) = t15/2-t16/2
      PARTICLE(3,0) = 0.D0
      PARTICLE(3,1) = 0.D0
      PARTICLE(3,2) = 0.D0
      PARTICLE(3,3) = 0.D0
      ELSE IF (PARITY.eq.-1) THEN
      t2 = (x1+x2)*Ep
      t4 = El/Eq
      t5 = t2*t4
      t7 = (x1-x2)*Ep
      t8 = t4*BCsTh
      t9 = t7*t8
      PARTICLE(1,0) = t5/2+t9/2
      t11 = sin(DBLE(phi))
      t13 = t11*El*BSnTh
      PARTICLE(1,1) = -t13
      t14 = cos(DBLE(phi))
      t16 = t14*El*BSnTh
      PARTICLE(1,2) = -t16
      t17 = t7*t4
      t18 = t2*t8
      PARTICLE(1,3) = -t17/2-t18/2
      PARTICLE(2,0) = t5/2-t9/2
      PARTICLE(2,1) = t13
      PARTICLE(2,2) = t16
      PARTICLE(2,3) = -t17/2+t18/2
      PARTICLE(3,0) = 0.D0
      PARTICLE(3,1) = 0.D0
      PARTICLE(3,2) = 0.D0
      PARTICLE(3,3) = 0.D0
      END IF

      CALL PassVar(PARTICLE,part)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE PreSeq2to2(CmFct)

      IMPLICIT DOUBLE PRECISION (t, s)

      include 'constcommon_ew.f'
      include 'eventcommon_ew.f'

      DOUBLE PRECISION CmFct

      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)
      t1 = s**2
      t3 = smin**2
      t5 = smax**2
      j1 = t1*s*(1/t3-1/t5)/2
      t10 = (Exp2Y-1)**2
      t12 = (Exp2Y+1)**2
      j2 = 2*limthy/(1-t10/t12)
      j1n2 = 1/Scm
      j3 = 2
      flxfct = 1/s/2
      PSfct = BtQ/PI/16
      j0 = 2*prej0*flxfct*PSfct*j1n2*j1*j2

      CmFct = j0

      RETURN
      END
