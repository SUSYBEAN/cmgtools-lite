      DOUBLE PRECISION FUNCTION NLOEWVAL_mls(X,N,W)
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
      DOUBLE PRECISION CS0z_mls,CS0v_mls,CS0s_mls,CS0c_mls,CS0n_mls,
     &                 CS0n1_mls,CS0n2_mls,CS0n3_mls,CS0n4_mls,CS0n5_mls,
     &                 CS0n6_mls,CS0n7_mls,CS0n8_mls,CS0n9_mls

      NLOEWVAL_mls = 0d0
      Fsum = 0d0

CCC 2->2 process (LO)

      IF ((EWflag.ge.8).or.(porder_ew.le.0)) THEN ! only regenerate events if weak correction is not on, to avoid call PDFs again
      CALL EventSeq2to2(X) ! generate events and prepare pdf functions
      ENDIF
      CALL PreSeq2to2(CmFct) ! generate common factor (contain phase space, flux and jacobian factors)

      CALL InfoSeq2to2(1,part) ! translate to particle 4-momentum information
      Pcut = constraint(part) ! test cuts
      IF (Pcut.ne.0d0) THEN
         CALL CrsecSeqz_mls() ! calculate component of tree amplitude
         Fi = CS0z_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0z_mls(1,iset)
               FVec(iset) = FVec(iset)*CmFct*W
            ENDDO
         ENDIF
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      CALL ReversePart(part) ! convert to parity reversed event for quark and antiquark from the other hadrons
      Mcut = constraint(part) ! test cuts
      IF (Mcut.ne.0d0) THEN
         if (Pcut.eq.0d0) CALL CrsecSeqz_mls() ! calculate component of tree amplitude only if not calculated before
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! calculate it if pp-bar collider (for pp collider it is exact the same result) or not done before
            Fi = CS0z_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0z_mls(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif ! always safe to remove the 'if' and run them every time because it is just a hack for efficiency in the case of pp collider
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

CCC Perturbation order

      NLOEWVAL_mls = Fsum
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
         CALL CrsecSeqv_mls() ! calculate component of virtual amplitude
         CALL CrsecSeqs_mls() ! calculate component of soft amplitude
         CALL CrsecSeqc_mls() ! calculate component of collinear amplitude
         Fi = CS0v_mls(1,0)
         Fi = Fi+CS0s_mls(1,0)
         Fi = Fi+CS0c_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0v_mls(1,iset)+CS0s_mls(1,iset)+CS0c_mls(1,iset)
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
            CALL CrsecSeqv_mls()
            CALL CrsecSeqs_mls()
            CALL CrsecSeqc_mls()
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then
            Fi = CS0v_mls(2,0)
            Fi = Fi+CS0s_mls(2,0)
            Fi = Fi+CS0c_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0v_mls(2,iset)+CS0s_mls(2,iset)+CS0c_mls(2,iset)
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

      CALL EventSeq2to3_mls(X)
      CALL PreSeq2to3_mls(CmFct)
      DoneOnce = .false.

      CALL InfoSeq2to3_mls(1,part)
      CALL doPhtRecb(part,rlflag)
      IF (rlflag.gt.0) then
         Pcut = constraint(part)
      ELSE ! no photon recombination, don't pass it
         Pcut = 0d0
      ENDIF
      IF (Pcut.ne.0d0) THEN
         CALL CrsecSeqn_mls()
         DoneOnce = .true.
         Fi = CS0n_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n_mls(1,iset)
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
            CALL CrsecSeqn_mls() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n_mls(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

      CALL InfoSeq2to3_mls(8,part)
      Pcut = constraint(part)
      CALL InfoSeq2to3_mls(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_mls()
            DoneOnce = .true.
         ENDIF
         Fi = CS0n7_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n7_mls(1,iset)
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
            CALL CrsecSeqn_mls() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n7_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n7_mls(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

      CALL InfoSeq2to3_mls(10,part)
      Pcut = constraint(part)
      CALL InfoSeq2to3_mls(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_mls()
            DoneOnce = .true.
         ENDIF
         Fi = CS0n9_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n9_mls(1,iset)
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
            CALL CrsecSeqn_mls() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n9_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n9_mls(2,iset)
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

      CALL EventSeq2to3_mls(X)
      CALL PreSeq2to3_mls(CmFct)
      DoneOnce = .false.

      CALL InfoSeq2to3_mls(1,part)
      CALL doPhtRecb(part,rlflag)
      Pcut = constraint(part)
      IF (Pcut.ne.0d0) THEN
         CALL CrsecSeqn_mls() ! calculate component of real radiation amplitude
         DoneOnce = .true.
         Fi = CS0n_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n_mls(1,iset)
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
            CALL CrsecSeqn_mls() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n_mls(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

C --- all subtraction terms start here

      IF (dble(sllup).ge.dble(sqq)) THEN ! do it here for efficiency though cuts can do it
      CALL InfoSeq2to3_mls(2,part) 
      Pcut = constraint(part)
      CALL InfoSeq2to3_mls(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_mls()
            DoneOnce = .true.
         ENDIF
         Fi = CS0n1_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n1_mls(1,iset)
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
            CALL CrsecSeqn_mls() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n1_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n1_mls(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      ENDIF

      IF (dble(sllup).ge.dble(sqq)) THEN ! do it here for efficiency though cuts can do it
      CALL InfoSeq2to3_mls(3,part) 
      Pcut = constraint(part)
      CALL InfoSeq2to3_mls(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_mls()
            DoneOnce = .true.
         ENDIF 
         Fi = CS0n2_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n2_mls(1,iset)
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
            CALL CrsecSeqn_mls() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n2_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n2_mls(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      ENDIF

      IF (dble(sllup).ge.dble(sqq)) THEN ! do it here for efficiency though cuts can do it
      CALL InfoSeq2to3_mls(4,part) 
      Pcut = constraint(part)
      CALL InfoSeq2to3_mls(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_mls()
            DoneOnce = .true.
         ENDIF 
         Fi = CS0n3_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n3_mls(1,iset)
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
            CALL CrsecSeqn_mls() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n3_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n3_mls(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      ENDIF

      IF (dble(sllup).ge.dble(sqq)) THEN ! do it here for efficiency though cuts can do it
      CALL InfoSeq2to3_mls(5,part)
      Pcut = constraint(part)
      CALL InfoSeq2to3_mls(0,part) 
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_mls()
            DoneOnce = .true.
         ENDIF
         Fi = CS0n4_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n4_mls(1,iset)
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
            CALL CrsecSeqn_mls() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n4_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n4_mls(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      ENDIF

      IF (dble(sllup).ge.dble(sqq)) THEN ! do it here for efficiency though cuts can do it
      CALL InfoSeq2to3_mls(6,part)
      Pcut = constraint(part)
      CALL InfoSeq2to3_mls(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_mls()
            DoneOnce = .true.
         ENDIF
         Fi = CS0n5_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n5_mls(1,iset)
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
            CALL CrsecSeqn_mls() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n5_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n5_mls(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF
      ENDIF

      CALL InfoSeq2to3_mls(7,part)
      Pcut = constraint(part)
      CALL InfoSeq2to3_mls(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_mls()
            DoneOnce = .true.
         ENDIF
         Fi = CS0n6_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n6_mls(1,iset)
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
            CALL CrsecSeqn_mls() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n6_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n6_mls(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

      CALL InfoSeq2to3_mls(8,part)
      Pcut = constraint(part)
      CALL InfoSeq2to3_mls(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_mls()
            DoneOnce = .true.
         ENDIF
         Fi = CS0n7_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n7_mls(1,iset)
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
            CALL CrsecSeqn_mls() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n7_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n7_mls(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

      CALL InfoSeq2to3_mls(9,part)
      Pcut = constraint(part)
      CALL InfoSeq2to3_mls(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_mls()
            DoneOnce = .true.
         ENDIF
         Fi = CS0n8_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n8_mls(1,iset)
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
            CALL CrsecSeqn_mls() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n8_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n8_mls(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

      CALL InfoSeq2to3_mls(10,part)
      Pcut = constraint(part)
      CALL InfoSeq2to3_mls(0,part)
      IF (Pcut.ne.0d0) THEN
         IF (.not.DoneOnce) THEN
            CALL CrsecSeqn_mls()
            DoneOnce = .true.
         ENDIF
         Fi = CS0n9_mls(1,0)
         Fi = Fi*CmFct
         IF (NPDFerror.gt.0) THEN
            DO iset=1,NPDFerror
               FVec(iset) = CS0n9_mls(1,iset)
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
            CALL CrsecSeqn_mls() ! calculate component of tree amplitude only if not calculated before
            DoneOnce = .true.
         endif
         if ((ptype.ne.1).or.(Pcut.eq.0d0)) then ! only calculated if it is not pp collider (for pp collider it is exact the same result)
            Fi = CS0n9_mls(2,0)
            Fi = Fi*CmFct
            IF (NPDFerror.gt.0) THEN
               DO iset=1,NPDFerror
                  FVec(iset) = CS0n9_mls(2,iset)
                  FVec(iset) = FVec(iset)*CmFct*W
               ENDDO
            ENDIF
         endif
         Fsum = Fsum + Fi
         CALL fill_hist(part,Fi*W,FVec)
      ENDIF

      NLOEWVAL_mls = Fsum

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C      SUBROUTINE EventSeq2to2_mls(X) ! replaced by EventSeq2to2(X) in nlocorr_weak_ew.F

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE EventSeq2to3_mls(X)

      IMPLICIT DOUBLE PRECISION (t, s)

      DOUBLE PRECISION X(6)

      include 'constcommon_ew.f'
      include 'eventcommon_ew.f'

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
      t32 = PI*lam3
      CsThl = cos(DBLE(t32))
      SnThl = sin(DBLE(t32))
      t33 = PI*lam4
      CsThr = cos(DBLE(t33))
      SnThr = sin(DBLE(t33))
      phil = 2*PI*lam5
      t35 = 1/Scm
      t37 = log(DCMPLX(sqq*t35))
      x1 = exp(lam6*t37)
      x2 = sqq/x1*t35
      t43 = cos(DBLE(phil))
      CsThS = CsThl*CsThr-SnThl*SnThr*t43
      t47 = cos(DBLE(phir))
      t49 = sin(DBLE(phil))
      t51 = sin(DBLE(phir))
      t53 = SnThl*t43
      phiS = atan2(DBLE(t47*SnThl*t49+t51*CsThr*t53+t51*SnThr*CsThl),DBL
     #E(-t51*SnThl*t49+t47*CsThr*t53+t47*SnThr*CsThl))

      m2 = sqq
      CALL PDFSeqAll(x1,x2,1d0,1d0,1,1)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C      SUBROUTINE InfoSeq2to2_mls(PARITY,part) ! replaced by InfoSeq2to2(X) in nlocorr_weak_ew.F

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE InfoSeq2to3_mls(PARITY,part)

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
      t22 = SnThr*El*SnThl*t20
      t24 = t12*t22/2
      t25 = t3*Er
      t26 = t25*t7
      t29 = CsThr*t5*t7
      t34 = (-t2*t26/4+t12*t29/4)*El*CsThl
      PARTICLE(1,0) = t18-t24+t34
      t35 = sin(DBLE(phir))
      t36 = t35*SnThr
      t37 = t36*Er
      t38 = t37/2
      t39 = cos(DBLE(phir))
      t41 = sin(DBLE(phil))
      t42 = SnThl*t41
      t43 = t39*El*t42
      t46 = El*SnThl*t20
      t47 = t35*CsThr*t46
      t48 = t5*CsThl
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
      t73 = (-t11*t26/4+t61*t29/4)*El*CsThl
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
      t22 = SnThr*El*SnThl*t20
      t24 = t12*t22/2
      t25 = t3*Er
      t26 = t25*t7
      t29 = CsThr*t5*t7
      t34 = (-t2*t26/4+t12*t29/4)*El*CsThl
      PARTICLE(1,0) = t18-t24+t34
      t35 = sin(DBLE(phir))
      t36 = t35*SnThr
      t37 = t36*Er
      t38 = t37/2
      t39 = cos(DBLE(phir))
      t41 = sin(DBLE(phil))
      t42 = SnThl*t41
      t43 = t39*El*t42
      t46 = El*SnThl*t20
      t47 = t35*CsThr*t46
      t48 = t5*CsThl
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
      t73 = (-t11*t26/4+t61*t29/4)*El*CsThl
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
      Exp2Y = x1/x2
      CsTh = CsThS
      SnTh = sqrt(1d0-CsThS**2)
      phi = phiS
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      ELSE IF ((PARITY.eq.3).or.(PARITY.eq.-3)) THEN
      s = sqq
      El = Eq
      Exp2Y = x1/x2
      CsTh = CsThl
      SnTh = SnThl
      phi = phil+phir
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      ELSE IF ((PARITY.eq.4).or.(PARITY.eq.-4)) THEN
      s = sqq
      El = Eq
      Exp2Y = x1/x2
      CsTh = CsThr
      SnTh = SnThr
      phi = phir
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      ELSE IF ((PARITY.eq.5).or.(PARITY.eq.-5)) THEN
      s = sqq
      El = Eq
      Exp2Y = x1/x2
      CsTh = -CsThl
      SnTh = SnThl
      phi = PI-phil+phir
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      ELSE IF ((PARITY.eq.6).or.(PARITY.eq.-6)) THEN
      s = sqq
      El = Eq
      Exp2Y = x1/x2
      CsTh = -CsThr
      SnTh = SnThr
      phi = phir
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      ELSE IF ((PARITY.eq.7).or.(PARITY.eq.-7)) THEN
      s = sll
      Eq = El
      t3 = 1/sqq
      Exp2Y = x1/x2*sll*t3
      x1 = x1*sll*t3
      CsTh = CsThl
      SnTh = SnThl
      phi = phil+phir
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      ELSE IF ((PARITY.eq.8).or.(PARITY.eq.-8)) THEN
      s = sqq
      El = Eq
      Exp2Y = x1/x2
      CsTh = CsThr
      SnTh = SnThr
      phi = phir
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      ELSE IF ((PARITY.eq.9).or.(PARITY.eq.-9)) THEN
      s = sll
      Eq = El
      Exp2Y = x1/x2*sqq/sll
      x2 = x2*sll/sqq
      CsTh = -CsThl
      SnTh = SnThl
      phi = PI-phil+phir
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      ELSE IF ((PARITY.eq.10).or.(PARITY.eq.-10)) THEN
      s = sqq
      El = Eq
      Exp2Y = x1/x2
      CsTh = -CsThr
      SnTh = SnThr
      phi = phir
      CALL InfoSeq2to2(ISIGN(1,PARITY),part)
      RETURN
      END IF

      CALL PassVar(PARTICLE,part)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C      SUBROUTINE PreSeq2to2_mls(CmFct) ! replaced by PreSeq2to2(X) in nlocorr_weak_ew.F

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE PreSeq2to3_mls(CmFct)

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
      j3 = PI*SnThl
      j4 = PI*SnThr
      t32 = 1/Scm
      t34 = log(DCMPLX(sqq*t32))
      j6 = -x1*t34
      flxfct = 1/sqq/2
      t37 = PI**2
      PSfct = Er/t37/PI/Eq/512
      j1n6 = 1/x1*t32
      j0 = prej0*flxfct*PSfct*j1n6*j1*j2*j3*j4*j6

      CmFct = j0

      RETURN
      END
