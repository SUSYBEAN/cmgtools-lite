      DOUBLE PRECISION FUNCTION NLOEWVAL_msv(X,N,W)
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
      INTEGER rlflag
      INTEGER iset

      DOUBLE PRECISION constraint
      DOUBLE PRECISION CS0z_msv,CS0v_msv,CS0s_msv,CS0c_msv,CS0n_msv,
     &                 CS0n1_msv,CS0n2_msv,CS0n3_msv,CS0n4_msv,CS0n5_msv

      NLOEWVAL_msv = 0d0
      Fsum = 0d0

CCC 2->2 process (LO)

      IF ((EWflag.ge.8).or.(porder_ew.le.0)) THEN ! only regenerate events if weak correction is not on, to avoid call PDFs again
      CALL EventSeq2to2(X) ! generate events and prepare pdf functions
      ENDIF
      CALL PreSeq2to2(CmFct) ! generate common factor (contain phase space, flux and jacobian factors)

      CALL InfoSeq2to2(1,part) ! translate to particle 4-momentum information
      Pcut = constraint(part) ! test cuts
      IF (Pcut.ne.0d0) THEN
         CALL CrsecSeqz_msv() ! calculate component of tree amplitude
         Fi = CS0z_msv(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0z_msv(1,iset)
               FVec(iset) = FVec(iset)*CmFct*W
            ENDDO
         ENDIF
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      CALL ReversePart(part) ! convert to parity reversed event for quark and antiquark from the other hadrons
      Mcut = constraint(part) ! test cuts
      IF (Mcut.ne.0d0) THEN
         if (Pcut.eq.0d0) CALL CrsecSeqz_msv() ! calculate component of tree amplitude only if not calculated before
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! ! calculate it if pp-bar collider (for pp collider it is exact the same result) or not done before
            Fi = CS0z_msv(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0z_msv(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif ! always safe to remove the 'if' and run them every time because it is just a hack for efficiency in the case of pp collider
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

CCC Perturbation order

      NLOEWVAL_msv = Fsum
      IF (porder_ew.le.0) RETURN
      IF ((EWflag.eq.7).or.(EWflag.eq.15)) RETURN ! save the trouble if QED ISR, FSR, ISRxFSR all off

CCC 2->2 process

C the followings are done in LO 2->2
C      CALL EventSeq2to2(X) ! skip to use the existing info from LO calculation
C      CALL PreSeq2to2(CmFct) ! skip to use the existing info from LO calculation

C      CALL InfoSeq2to2(1,part) ! skip to use the existing info from LO calculation
      CALL ReversePart(part) ! Done in LO already, reverse particle info back just for filling histogram
C      Pcut = constraint(part)
      IF (Pcut.ne.0d0) THEN
         CALL CrsecSeqv_msv() ! calculate component of virtual amplitude
         CALL CrsecSeqs_msv() ! calculate component of soft amplitude
         CALL CrsecSeqc_msv() ! calculate component of collinear amplitude
         Fi = CS0v_msv(1,0)
         Fi = Fi+CS0s_msv(1,0)
         Fi = Fi+CS0c_msv(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0v_msv(1,iset)+CS0s_msv(1,iset)+CS0c_msv(1,iset)
               FVec(iset) = FVec(iset)*CmFct*W
            ENDDO
         ENDIF
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      CALL ReversePart(part) ! convert again to parity reversed event for filling histogram
C      Mcut = constraint(part)
      IF (Mcut.ne.0d0) THEN
         if (Pcut.eq.0d0) then
            CALL CrsecSeqv_msv()
            CALL CrsecSeqs_msv()
            CALL CrsecSeqc_msv()
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then
            Fi = CS0v_msv(2,0)
            Fi = Fi+CS0s_msv(2,0)
            Fi = Fi+CS0c_msv(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0v_msv(2,iset)+CS0s_msv(2,iset)+CS0c_msv(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

CCC 2->3 process: hacked extra part

      sqqmax = smax
      sllmax = smin
      sllmin = 4d0*Ml**2

      CALL EventSeq2to3_msv(X)
      CALL PreSeq2to3_msv(CmFct)
      DoneOnce = .false.

      CALL InfoSeq2to3_msv(1,part)
      CALL doPhtRecb(part,rlflag)
      IF (rlflag.gt.0) then
         Pcut = constraint(part)
      ELSE ! no photon recombination, don't pass it
         Pcut = 0d0
      ENDIF
      IF (Pcut.ne.0d0) THEN
         CALL CrsecSeqn_msv()
         DoneOnce = .true.
         Fi = CS0n_msv(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n_msv(1,iset)
               FVec(iset) = FVec(iset)*CmFct*W
            ENDDO
         ENDIF
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      CALL ReversePart(part) ! convert to parity reversed event for quark and antiquark from the other hadrons
      IF (rlflag.gt.0) then
         Mcut = constraint(part)
      ELSE ! no photon recombination, don't pass it
         Mcut = 0d0
      ENDIF
      IF (Mcut.ne.0d0) THEN
         if (Pcut.eq.0d0) then
            CALL CrsecSeqn_msv() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n_msv(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n_msv(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

      sqqmax = Scm
      sllmax = smax
      sllmin = smin

CCC 2->3 process: standard part

      CALL EventSeq2to3_msv(X)
      CALL PreSeq2to3_msv(CmFct)
      DoneOnce = .false.

      CALL InfoSeq2to3_msv(1,part)
      CALL doPhtRecb(part,rlflag)
      Pcut = constraint(part)
      IF (Pcut.ne.0d0) THEN
         CALL CrsecSeqn_msv() ! calculate component of real radiation amplitude
         DoneOnce = .true.
         Fi = CS0n_msv(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n_msv(1,iset)
               FVec(iset) = FVec(iset)*CmFct*W
            ENDDO
         ENDIF
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      CALL ReversePart(part) ! convert to parity reversed event for quark and antiquark from the other hadrons
      Mcut = constraint(part) ! test cuts
      IF (Mcut.ne.0d0) THEN
         if (Pcut.eq.0d0) then
            CALL CrsecSeqn_msv() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n_msv(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n_msv(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

C --- all subtraction terms start here

      IF (dble(sllup).ge.dble(sqq)) THEN ! do it here for efficiency though cuts can do it
      CALL InfoSeq2to3_msv(2,part) 
      Pcut = constraint(part)
      CALL InfoSeq2to3_msv(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_msv()
            DoneOnce = .true.
         ENDIF
         Fi = CS0n1_msv(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n1_msv(1,iset)
               FVec(iset) = FVec(iset)*CmFct*W
            ENDDO
         ENDIF
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      CALL ReversePart(part) ! convert to parity reversed event for quark and antiquark from the other hadrons
      Mcut = constraint(part) ! test cuts
      IF (Mcut.ne.0d0) THEN
         if ((.not.DoneOnce).and.(Pcut.eq.0d0)) then
            CALL CrsecSeqn_msv() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n1_msv(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n1_msv(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      ENDIF

      IF (dble(sllup).ge.dble(sqq)) THEN ! do it here for efficiency though cuts can do it
      CALL InfoSeq2to3_msv(3,part) 
      Pcut = constraint(part)
      CALL InfoSeq2to3_msv(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_msv()
            DoneOnce = .true.
         ENDIF 
         Fi = CS0n2_msv(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n2_msv(1,iset)
               FVec(iset) = FVec(iset)*CmFct*W
            ENDDO
         ENDIF
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      CALL ReversePart(part) ! convert to parity reversed event for quark and antiquark from the other hadrons
      Mcut = constraint(part) ! test cuts
      IF (Mcut.ne.0d0) THEN
         if ((.not.DoneOnce).and.(Pcut.eq.0d0)) then
            CALL CrsecSeqn_msv() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n2_msv(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n2_msv(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      ENDIF

      IF (dble(sllup).ge.dble(sqq)) THEN ! do it here for efficiency though cuts can do it
      CALL InfoSeq2to3_msv(4,part) 
      Pcut = constraint(part)
      CALL InfoSeq2to3_msv(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_msv()
            DoneOnce = .true.
         ENDIF 
         Fi = CS0n3_msv(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n3_msv(1,iset)
               FVec(iset) = FVec(iset)*CmFct*W
            ENDDO
         ENDIF
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      CALL ReversePart(part) ! convert to parity reversed event for quark and antiquark from the other hadrons
      Mcut = constraint(part) ! test cuts
      IF (Mcut.ne.0d0) THEN
         if ((.not.DoneOnce).and.(Pcut.eq.0d0)) then
            CALL CrsecSeqn_msv() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n3_msv(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n3_msv(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      ENDIF

      CALL InfoSeq2to3_msv(5,part)
      Pcut = constraint(part)
      CALL InfoSeq2to3_msv(0,part) 
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_msv()
            DoneOnce = .true.
         ENDIF
         Fi = CS0n4_msv(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n4_msv(1,iset)
               FVec(iset) = FVec(iset)*CmFct*W
            ENDDO
         ENDIF
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      CALL ReversePart(part) ! convert to parity reversed event for quark and antiquark from the other hadrons
      Mcut = constraint(part) ! test cuts
      IF (Mcut.ne.0d0) THEN
         if ((.not.DoneOnce).and.(Pcut.eq.0d0)) then
            CALL CrsecSeqn_msv() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n4_msv(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n4_msv(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

      CALL InfoSeq2to3_msv(6,part) 
      Pcut = constraint(part)
      CALL InfoSeq2to3_msv(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_msv()
            DoneOnce = .true.
         ENDIF
         Fi = CS0n5_msv(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n5_msv(1,iset)
               FVec(iset) = FVec(iset)*CmFct*W
            ENDDO
         ENDIF
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      CALL ReversePart(part) ! convert to parity reversed event for quark and antiquark from the other hadrons
      Mcut = constraint(part) ! test cuts
      IF (Mcut.ne.0d0) THEN
         if ((.not.DoneOnce).and.(Pcut.eq.0d0)) then
            CALL CrsecSeqn_msv() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n5_msv(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n5_msv(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

      NLOEWVAL_msv = Fsum

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C      SUBROUTINE EventSeq2to2_msv(X) ! replaced by EventSeq2to2(X) in nlocorr_weak_ew.F

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE EventSeq2to3_msv(X)

      IMPLICIT DOUBLE PRECISION (t, s)

      DOUBLE PRECISION X(6)

      include 'constcommon_ew.f'
      include 'eventcommon_ew.f'

      LOGICAL CIsBad
      DOUBLE COMPLEX spf6
      DOUBLE PRECISION epsMl
      PARAMETER (epsMl=1d-12)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)
      lam1 = X(1)
      lam2 = X(2)
      lam3 = X(3)
      lam4 = X(4)
      lam5 = X(5)
      lam6 = X(6)

      phir = 0
      t1 = smax**2
      t3 = smin**2
      t4 = 1/t3
      t8 = sqrt((1/t1-t4)*lam1+t4)
      t11 = 1/sqqmin
      t13 = lam1**2
      if (dble(sllmax).le.dble(sqqmin)) then
      sqq = 1/t8
      else
      sqq = 1/((1/sqqmax-t11)*t13+t11)
      end if
      sllup = min(DBLE(sllmax),DBLE(sqq))
      t20 = 1/sllup
      t23 = lam2**2
      if (dble(sllmax).le.dble(sqqmin)) then
      sll = (sllup-sllmin)*lam2+sllmin
      else
      sll = 1/(t20-(t20-1/sllmin)*t23)
      end if
      t27 = sqrt(sll)
      El = t27/2
      t28 = sqrt(sqq)
      Eq = t28/2
      Er = (sqq-sll)/Eq/4
      t32 = Ml**2
      BtQl = sqrt(1-4*t32/sll)
      t37 = 1+BtQl
      t41 = t37/t32*sll/2
      if (ABS(DBLE(Ml**2/sll)).lt.epsMl) then
      spf6 = t41
      else
      spf6 = -t37/(-1+BtQl)
      end if
      if ((ABS(DBLE(spf6)).eq.0d0).or.(CIsBad(spf6))) spf6 = t41
      spf6 = log(spf6)
      AthClc = spf6/2
      BCsThl = tanh(DBLE(AthClc*(2*lam3-1)))
      CsThl = BCsThl/BtQl
      t49 = CsThl**2
      SnThl = sqrt(1-t49)
      BSnThl = BtQl*SnThl
      t51 = PI*lam4
      CsThr = cos(DBLE(t51))
      SnThr = sin(DBLE(t51))
      phil = 2*PI*lam5
      t53 = 1/Scm
      t55 = log(DCMPLX(sqq*t53))
      x1 = exp(lam6*t55)
      x2 = sqq/x1*t53
      BtQlq = sqrt(1-4*t32/sqq)
      BCsThlq = CsThl*BtQlq
      BSnThlq = SnThl*BtQlq
      t65 = cos(DBLE(phil))
      CsThS = CsThl*CsThr-SnThl*SnThr*t65
      BCsThS = BCsThlq*CsThr-BSnThlq*SnThr*t65
      t72 = cos(DBLE(phir))
      t74 = sin(DBLE(phil))
      t76 = sin(DBLE(phir))
      t78 = BSnThl*t65
      phiS = atan2(DBLE(t72*BSnThl*t74+t76*CsThr*t78+t76*SnThr*BCsThl),D
     #BLE(-t76*BSnThl*t74+t72*CsThr*t78+t72*SnThr*BCsThl))

      m2 = sqq
      CALL PDFSeqAll(x1,x2,1d0,1d0,1,1)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C      SUBROUTINE InfoSeq2to2_msv(PARITY,part) ! replaced by InfoSeq2to2(X) in nlocorr_weak_ew.F

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE InfoSeq2to3_msv(PARITY,part)

      IMPLICIT DOUBLE PRECISION (t, s)

      INTEGER PARITY
      DOUBLE PRECISION PARTICLE(3,0:3)

      include 'constcommon_ew.f'
      include 'eventcommon_ew.f'
      include 'particles.f'

C---- For backing up variables when mapping 2->3 event to 2->2 event at soft/collinear limit
      double precision El_bak, Eq_bak, x1_bak, x2_bak
      save El_bak, Eq_bak, x1_bak, x2_bak

      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      IF (PARITY.eq.1) THEN
      El_bak = El
      Eq_bak = Eq
      x1_bak = x1
      x2_bak = x2
      t2 = (x1+x2)*Ep
      t3 = 1/Eq
      t5 = 2*Eq-Er
      t7 = 1/El
      t8 = t3*t5*t7
      t11 = (x1-x2)*Ep
      t12 = t11*t3
      t14 = CsThr*Er*t7
      t18 = (t2*t8/4-t12*t14/4)*El
      t20 = cos(DBLE(phil))
      t22 = SnThr*El*BSnThl*t20
      t24 = t12*t22/2
      t25 = t3*Er
      t26 = t25*t7
      t29 = CsThr*t5*t7
      t34 = (-t2*t26/4+t12*t29/4)*El*BCsThl
      PARTICLE(1,0) = t18-t24+t34
      t35 = sin(DBLE(phir))
      t36 = t35*SnThr
      t37 = t36*Er
      t38 = t37/2
      t39 = cos(DBLE(phir))
      t41 = sin(DBLE(phil))
      t42 = BSnThl*t41
      t43 = t39*El*t42
      t46 = El*BSnThl*t20
      t47 = t35*CsThr*t46
      t48 = t5*BCsThl
      t50 = t36*t48/2
      PARTICLE(1,1) = -t38+t43+t47+t50
      t51 = t39*SnThr
      t52 = t51*Er
      t53 = t52/2
      t55 = t35*El*t42
      t57 = t39*CsThr*t46
      t59 = t51*t48/2
      PARTICLE(1,2) = -t53-t55+t57+t59
      t61 = t2*t3
      t65 = (t11*t8/4-t61*t14/4)*El
      t67 = t61*t22/2
      t73 = (-t11*t26/4+t61*t29/4)*El*BCsThl
      PARTICLE(1,3) = t65-t67+t73
      PARTICLE(2,0) = t18+t24-t34
      PARTICLE(2,1) = -t38-t43-t47-t50
      PARTICLE(2,2) = -t53+t55-t57-t59
      PARTICLE(2,3) = t65+t67-t73
      t76 = t3*CsThr*Er
      PARTICLE(3,0) = t2*t25/2+t11*t76/2
      PARTICLE(3,1) = t37
      PARTICLE(3,2) = t52
      PARTICLE(3,3) = t11*t25/2+t2*t76/2
      ELSE IF (PARITY.eq.-1) THEN
      El_bak = El
      Eq_bak = Eq
      x1_bak = x1
      x2_bak = x2
      t2 = (x1+x2)*Ep
      t3 = 1/Eq
      t5 = 2*Eq-Er
      t7 = 1/El
      t8 = t3*t5*t7
      t11 = (x1-x2)*Ep
      t12 = t11*t3
      t14 = CsThr*Er*t7
      t18 = (t2*t8/4-t12*t14/4)*El
      t20 = cos(DBLE(phil))
      t22 = SnThr*El*BSnThl*t20
      t24 = t12*t22/2
      t25 = t3*Er
      t26 = t25*t7
      t29 = CsThr*t5*t7
      t34 = (-t2*t26/4+t12*t29/4)*El*BCsThl
      PARTICLE(1,0) = t18-t24+t34
      t35 = sin(DBLE(phir))
      t36 = t35*SnThr
      t37 = t36*Er
      t38 = t37/2
      t39 = cos(DBLE(phir))
      t41 = sin(DBLE(phil))
      t42 = BSnThl*t41
      t43 = t39*El*t42
      t46 = El*BSnThl*t20
      t47 = t35*CsThr*t46
      t48 = t5*BCsThl
      t50 = t36*t48/2
      PARTICLE(1,1) = t38-t43-t47-t50
      t51 = t39*SnThr
      t52 = t51*Er
      t53 = t52/2
      t55 = t35*El*t42
      t57 = t39*CsThr*t46
      t59 = t51*t48/2
      PARTICLE(1,2) = t53+t55-t57-t59
      t61 = t2*t3
      t65 = (t11*t8/4-t61*t14/4)*El
      t67 = t61*t22/2
      t73 = (-t11*t26/4+t61*t29/4)*El*BCsThl
      PARTICLE(1,3) = -t65+t67-t73
      PARTICLE(2,0) = t18+t24-t34
      PARTICLE(2,1) = t38+t43+t47+t50
      PARTICLE(2,2) = t53-t55+t57+t59
      PARTICLE(2,3) = -t65-t67+t73
      t76 = t3*CsThr*Er
      PARTICLE(3,0) = t2*t25/2+t11*t76/2
      PARTICLE(3,1) = -t37
      PARTICLE(3,2) = -t52
      PARTICLE(3,3) = -t11*t25/2-t2*t76/2
      ELSE IF (PARITY.eq.0) THEN
      El = El_bak
      Eq = Eq_bak
      x1 = x1_bak
      x2 = x2_bak
      RETURN
      ELSE IF ((PARITY.eq.2).or.(PARITY.eq.-2)) THEN
      s = sqq
      El = Eq
      BtQ = BtQlq
      Exp2Y = x1/x2
      CsTh = CsThS
      BCsTh = BCsThS
      BSnTh = BtQlq*sqrt(1d0-CsThS**2)
      phi = phiS
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      ELSE IF ((PARITY.eq.3).or.(PARITY.eq.-3)) THEN
      s = sqq
      El = Eq
      BtQ = BtQlq
      Exp2Y = x1/x2
      CsTh = CsThl
      BCsTh = BCsThlq
      BSnTh = BSnThlq
      phi = phil+phir
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      ELSE IF ((PARITY.eq.4).or.(PARITY.eq.-4)) THEN
      s = sqq
      El = Eq
      BtQ = BtQlq
      Exp2Y = x1/x2
      CsTh = -CsThl
      BCsTh = -BCsThlq
      BSnTh = BSnThlq
      phi = PI-phil+phir
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      ELSE IF ((PARITY.eq.5).or.(PARITY.eq.-5)) THEN
      s = sll
      Eq = El
      BtQ = BtQl
      t3 = 1/sqq
      Exp2Y = x1/x2*sll*t3
      x1 = x1*sll*t3
      CsTh = CsThl
      BCsTh = BCsThl
      BSnTh = BSnThl
      phi = phil+phir
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      ELSE IF ((PARITY.eq.6).or.(PARITY.eq.-6)) THEN
      s = sll
      Eq = El
      BtQ = BtQl
      Exp2Y = x1/x2*sqq/sll
      x2 = x2*sll/sqq
      CsTh = -CsThl
      BCsTh = -BCsThl
      BSnTh = BSnThl
      phi = PI-phil+phir
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      END IF

      CALL PassVar(PARTICLE,part)

      RETURN

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C      SUBROUTINE PreSeq2to2_msv(CmFct) ! replaced by PreSeq2to2(X) in nlocorr_weak_ew.F

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE PreSeq2to3_msv(CmFct)

      IMPLICIT DOUBLE PRECISION (t, s)

      include 'constcommon_ew.f'
      include 'eventcommon_ew.f'

      DOUBLE PRECISION CmFct

      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)
      t1 = sqq**2
      t3 = sqqmin**2
      t4 = 1/t3
      t5 = sqqmax**2
      t17 = sqrt((sqq-sqqmin)*(sqqmax-sqqmin)*t4/sqqmax)
      if (dble(sllmax).le.dble(sqqmin)) then
      j1 = t1*sqq*(t4-1/t5)/2
      else
      j1 = 4*sqq*Eq*t17
      end if
      t20 = sllup-sllmin
      t24 = sllup**2
      t29 = sqrt((sllup-sll)*t20/t24/sllmin)
      if (dble(sllmax).le.dble(sqqmin)) then
      j2 = t20
      else
      j2 = 4*sll*El*t29
      end if
      t32 = BCsThl**2
      j3 = (1-t32)*AthClc
      j4 = PI*SnThr
      t34 = 1/Scm
      t36 = log(DCMPLX(sqq*t34))
      j6 = -x1*t36
      flxfct = 1/sqq/2
      t39 = PI**2
      PSfct = Er/t39/PI/Eq/256
      j1n6 = 1/x1*t34
      j0 = prej0*flxfct*PSfct*j1n6*j1*j2*j3*j4*j6
      
      CmFct = j0

      RETURN
      END
