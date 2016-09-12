CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC Translate kinematical variables into particle four momentum in lab frame

      subroutine kinvar2partmom(x1,x2,s13,s23,
     & s14,s24,s34,s1v,s2v,s3v,s4v,s1l1,s2l1,s1l2,s2l2,
     & s3l2,s4l2,part)

      implicit none
      double precision pi
      parameter (pi=3.141592653589793d0)

      include 'cutcontrol.f'
      include 'particles.f'

      integer i,j

CCCCC Passed kinematics of the event
      double precision x1,x2,s13,s23,s14,s24,s34,
     & s1v,s2v,s3v,s4v,s1l1,s2l1,s1l2,s2l2,
     & s3l2,s4l2

CCCCC Kinematics of the Z boson and leptons, jets

      integer topflag

      double precision Enw,ptw,pzw,
     & En3,pt3,pz3,
     & En4,pt4,pz4,
     & Enl2,ptl2,pzl2,
     & Enl1,ptl1,pzl1

      double precision cos3,phi3,totphi3,
     & cos4,phi4,totphi4,
     & cosl2,phil2,totphil2,      
     & cosl1,phil1,totphil1

      double precision tmpval

CCCCC Main program parameters
      double precision m2,scm
      COMMON/params/ m2,scm

CCC ---------------------- Definition of particle 4-momenta -----------------
CCC   Particle number:
CCC     1. Z-boson
CCC     2. Lepton
CCC     3. Anti-lepton
CCC     4. Parton 1
CCC     5. Parton 2
CCC     6. Photon
CCC -------------------------------------------------------------------------
      
      do j=1,4
      do i=1,NPART
         part(i,j) =0d0
      enddo
      enddo

CCC -------------------------------------------------------------------------
CCC ----- Fill 1,4 components of 4-momenta

      part(2,1) = -0.5d0*dsqrt(scm)*(x2*s1l2+x1*s2l2)
      part(3,1) = -0.5d0*dsqrt(scm)*(x2*s1l1+x1*s2l1)
      part(4,1) = -0.5d0*dsqrt(scm)*(x2*s13+x1*s23)
      part(5,1) = -0.5d0*dsqrt(scm)*(x2*s14+x1*s24)
      part(1,1) = 0.5d0/x1/x2/dsqrt(scm)*(x2*(m2-x1*x2*scm*s1v)
     & +x1*(m2-x1*x2*scm*s2v))

      part(2,4) = 0.5d0*dsqrt(scm)*(x2*s1l2-x1*s2l2) 
      part(3,4) = 0.5d0*dsqrt(scm)*(x2*s1l1-x1*s2l1)
      part(4,4) = 0.5d0*dsqrt(scm)*(x2*s13-x1*s23) 
      part(5,4) = 0.5d0*dsqrt(scm)*(x2*s14-x1*s24)
      part(1,4) = 0.5d0/x1/x2/dsqrt(scm)*(-x2*(m2-x1*x2*scm*s1v)
     & +x1*(m2-x1*x2*scm*s2v))

CCC ----------------------- Energy and momenta of particles ------------------
CCC   l1=positron, l2=lepton

      Enw = part(1,1)
      pzw = part(1,4)
      Enl2 = part(2,1)
      pzl2 = part(2,4)
      Enl1 = part(3,1)
      pzl1 = part(3,4)
      En3 = part(4,1)
      pz3 = part(4,4)
      En4 = part(5,1)
      pz4 = part(5,4)

      ptl2 = dsqrt(dabs(x1*x2*scm*s1l2*s2l2))
      ptl1 = dsqrt(dabs(x1*x2*scm*s1l1*s2l1))
      pt3 = dsqrt(dabs(x1*x2*scm*s13*s23))
      pt4 = dsqrt(dabs(x1*x2*scm*s14*s24))
      ptw = dsqrt(dabs((m2-x1*x2*scm*s1v)*(m2-x1*x2*scm*s2v)
     & /x1/x2/scm-m2))

CCC -------------------------------------------------------------------------
CCC ---- Get the tranverse plane components of the 4-momenta
 
C     This flag denotes the multiplicity in the final state
C     topflag:
C     0 -> Z (pT=0) + 0 hard partons + 0 photon
C     1 -> Z (pT>0) + 1 hard parton
C     2 -> Z (pT>0) + 2 hard partons
C     3 -> Z (pT=0) + 2 balancing hard partons

      if (((pt3.lt.ptbuf).and.(pt4.lt.ptbuf)).or. 
     &   ((ptw.lt.ptbuf).and.((pt3.lt.ptbuf).or.(pt4.lt.ptbuf)))) then

C ---- If any two pTs are zero, all three must be.  Leptons are back-to-back in 
C ---- tranverse plane, put lepton at 0, positron at pi

         topflag = 0
         part(2,2)=ptl2
         part(3,2)=-ptl1

      elseif ((pt3.gt.ptbuf).and.(pt4.lt.ptbuf)) then

C ---- Z, parton 3 are back-to-back in transverse plane, put Z at 0,
C ---- parton 3 at pi

         topflag = 1
         part(1,2)=ptw
         part(4,2)=-pt3

      elseif ((pt4.gt.ptbuf).and.(pt3.lt.ptbuf)) then

C ---- Z, parton 4 are back-to-back in transverse plane, put Z at 0,
C ---- parton 4 at pi

         topflag = 1
         part(1,2)=ptw
         part(5,2)=-pt4

      elseif (ptw.gt.ptbuf) then

C ---- Z at non-zero pT and two hard partons; Z at phi=0

         topflag = 2

         cos3=(En3*Enw-pz3*pzw-0.5d0*(x1*x2*scm*s3v-m2))/pt3/ptw
         cos4=(En4*Enw-pz4*pzw-0.5d0*(x1*x2*scm*s4v-m2))/pt4/ptw       
         if (cos3.lt.(-1d0)) then
               cos3=-1d0
         else if (cos3.gt.1d0) then
               cos3=1d0
         endif
         if (cos4.lt.(-1d0)) then
               cos4=-1d0
         else if (cos4.gt.1d0) then
               cos4=1d0
         endif
         phi3=dacos(cos3)
         phi4=dacos(cos4)
         totphi3=phi3
         totphi4=2d0*pi-phi4      

         part(1,2)=ptw
         part(4,2)=pt3*dcos(totphi3)
         part(4,3)=pt3*dsin(totphi3)
         part(5,2)=pt4*dcos(totphi4)
         part(5,3)=pt4*dsin(totphi4)

      else

C ---- Z at zero pT, balanced by two jets: Mercedes Benz topology
C ---- Put lepton at phi=0, positron at pi, measure jets from these

         topflag=3

         part(2,2)=ptl2
         part(3,2)=-ptl1
         
         cos3=(En3*Enl2-pz3*pzl2-0.5d0*(x1*x2*scm*s3l2-m2))/pt3/ptl2
         cos4=(En4*Enl2-pz4*pzl2-0.5d0*(x1*x2*scm*s4l2-m2))/pt4/ptl2
         if (cos3.lt.(-1d0)) then
               cos3=-1d0
         else if (cos3.gt.1d0) then
               cos3=1d0
         endif
         if (cos4.lt.(-1d0)) then
               cos4=-1d0
         else if (cos4.gt.1d0) then
               cos4=1d0
         endif
         phi3=dacos(cos3)
         phi4=dacos(cos4)
         totphi3=phi3
         totphi4=2d0*pi-phi4

         part(4,2)=pt3*dcos(totphi3)
         part(4,3)=pt3*dsin(totphi3)
         part(5,2)=pt4*dcos(totphi4)
         part(5,3)=pt4*dsin(totphi4)

      endif

C ---- For the cases with Z pT<>0, use this for leptons

      if ((topflag.eq.1).or.(topflag.eq.2)) then

         cosl2=(Enl2*Enw-pzl2*pzw-0.5d0*m2)/ptw/ptl2
         cosl1=(Enl1*Enw-pzl1*pzw-0.5d0*m2)/ptw/ptl1
         if (cosl2.lt.(-1d0)) then
               cosl2=-1d0
         elseif (cosl2.gt.1d0) then
               cosl2=1d0
         endif
         if (cosl1.lt.(-1d0)) then
               cosl1=-1d0
         elseif (cosl1.gt.1d0) then   
               cosl1=1d0
         endif  
         phil2=dacos(cosl2)
         phil1=dacos(cosl1)
         totphil2=phil2
         totphil1=2d0*pi-phil1

         part(2,2)=ptl2*dcos(totphil2)
         part(2,3)=ptl2*dsin(totphil2)
         part(3,2)=ptl1*dcos(totphil1)
         part(3,3)=ptl1*dsin(totphil1)

      endif

CCC -------------------------------------------------------------------------
CCC ------------------------- Order partons here ----------------------------
CCC ------------------------- redefine parton pT and eta --------------------

      if (pt4.gt.pt3) then

         do j=1,4
            tmpval = part(4,j)
            part(4,j) = part(5,j)
            part(5,j) = tmpval
         enddo

      endif

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC    Jet algorithm; combines jets--should be called before
CCC    constraint or histogramming

      SUBROUTINE doJetAlgo(part,jaflag) ! return a flag for caching particle info
      IMPLICIT NONE
      INCLUDE 'cutcontrol.f'
      INCLUDE 'particles.f'

      INTEGER jaflag

      DOUBLE PRECISION pi
      parameter (pi=3.141592653589793d0)

      double precision jetcom(4)

      DOUBLE PRECISION eta3,eta4,pt3sq,pt4sq
      double precision deta,deta3com,deta4com,
     &     dphi,dphi3com,dphi4com

      INTEGER j

CCC -------------------------------------------------------------------------
CCC ------------------------- Jet Algorithm ---------------------------------

      eta3 = 0.5d0*dlog(dabs((part(4,1)+part(4,4))/
     &                       (part(4,1)-part(4,4))))

      IF (part(5,1).le.0d0) THEN ! no jet algorithm needed
         IF (dabs(eta3).gt.maxetaj4obs) THEN
            jaflag = 4 ! both jets down beam, set to zero
            CALL modifyPart(part,jaflag)
         ELSE
            jaflag = 0
         ENDIF
         RETURN
      ENDIF

      eta4 = 0.5d0*dlog(dabs((part(5,1)+part(5,4))/
     &                       (part(5,1)-part(5,4))))

CCC   Test here for efficiency
      IF ((dabs(eta3).gt.maxetaj4obs).and.(dabs(eta4).gt.maxetaj4obs)) THEN
         jaflag = 4 ! both jets down beam, set to zero
         CALL modifyPart(part,jaflag)
         RETURN
      ELSE
         jaflag = 0
      ENDIF

      IF (maxnumjet.ge.2) THEN ! implement

         pt3sq = part(4,2)**2+part(4,3)**2 ! store squares instead of taking root
         pt4sq = part(5,2)**2+part(5,3)**2

      if (pt3sq.ge.ptbuf**2.and.pt4sq.ge.ptbuf**2) then

         deta = 0.0d0
         dphi = 0.0d0

         deta = eta3-eta4
         dphi = dabs(datan2(part(4,3),part(4,2))-
     &        datan2(part(5,3),part(5,2)))
         if (dphi.gt.pi) dphi = 2d0*pi-dphi

CCC   k_T ALGORITHM
         if (jetalgo.eq.'ktal') then

c            if ((pt4*dsqrt(deta**2+dphi**2)).lt.
c     &           (pt4*jjdRmin)) then
            IF (pt4sq*(deta**2+dphi**2).lt.pt4sq*jjdRmin**2) THEN
               jaflag = 1
            endif

CCC   Anti-k_T ALGORITHM
         elseif (jetalgo.eq.'aktal') then

c            if (pt3.ne.0d0) then ! just to be safe
c               if ((dsqrt(deta**2+dphi**2)/pt3).lt.
c     &              (jjdRmin/pt3)) then ! why are we dividing by pt3 on both sides?
             IF (deta**2+dphi**2.lt.jjdRmin**2) THEN
                  jaflag = 1
             endif
c            endif

CCC   Cambridge-Aachen ALGORITHM
c         elseif (jetalgo.eq.'CAal') then
c
c            if ((dsqrt(deta**2+dphi**2)).lt.
c     &           (jjdRmin)) then
c               jaflag = 1
c            endif

CCC   CONE ALGORITHM
         elseif (jetalgo.eq.'cone') then

            do j = 1,4
               jetcom(j) = part(4,j)+part(5,j)
            enddo

            deta3com = eta3-0.5d0*dlog(dabs((jetcom(1)+jetcom(4))/
     &           (jetcom(1)-jetcom(4))))
            dphi3com = dabs(datan2(part(4,3),part(4,2))-
     &           datan2(jetcom(3),jetcom(2)))
            if (dphi3com.gt.pi) dphi3com = 2d0*pi-dphi3com
         
            deta4com = eta4-0.5d0*dlog(dabs((jetcom(1)+jetcom(4))/
     &           (jetcom(1)-jetcom(4))))
            dphi4com = dabs(datan2(part(5,3),part(5,2))-
     &           datan2(jetcom(3),jetcom(2)))
            if (dphi4com.gt.pi) dphi4com = 2d0*pi-dphi4com

c            if (((dsqrt(deta3com**2+dphi3com**2)).lt.jjdRmin).and.
c     &           ((dsqrt(deta4com**2+dphi4com**2)).lt.jjdRmin).and.
c     &           ((dsqrt(deta**2+dphi**2)).lt.(jjdRmin*Rsep))) then
            IF (deta3com**2+dphi3com**2.lt.jjdRmin**2.and.
     &          deta4com**2+dphi4com**2.lt.jjdRmin**2.and.
     &          deta**2+dphi**2.lt.jjdRmin**2*Rsep**2) THEN
               jaflag = 1
            endif

         endif
      endif
      ENDIF

      IF (jaflag.gt.0) then ! jet algorithm implemented, recalculate rapidity
         eta3 = 0.5d0*dlog(dabs((part(4,1)+part(5,1)+part(4,4)+part(5,4))/
     &                          (part(4,1)+part(5,1)-part(4,4)-part(5,4))))
         IF (dabs(eta3).gt.maxetaj4obs) jaflag = 4 ! set both jets to zero
      ELSE
         IF (dabs(eta3).gt.maxetaj4obs) jaflag = 2 ! 1st jet down beam, set to zero
         IF (dabs(eta4).gt.maxetaj4obs) jaflag = 3 ! 2nd jet down beam, set to zero
      ENDIF

      if (jaflag.gt.0) CALL modifyPart(part,jaflag) ! merges if we should merge

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC    Photon Recombination; combines photon with leptons--should
CCC    be called before constraint or histogramming

      SUBROUTINE doPhtRecb(part,rlflag)
      IMPLICIT NONE
      INCLUDE 'cutcontrol.f'
      INCLUDE 'particles.f'

      INTEGER rlflag

      DOUBLE PRECISION pi
      parameter (pi=3.141592653589793d0)

      DOUBLE PRECISION etal1,etal2,etar
      double precision deta,dphi
      double precision dRsqrl1,dRsqrl2
      double precision tmpdRsq

      INTEGER j

      rlflag = 0

CCC more efficient to see if photon 4-momentum is filled before do photon recombination
      if ((part(6,1).le.0d0).or.(lrdRrecb.eq.0d0)) return

      etar = 0.5d0*dlog(dabs((part(6,1)+part(6,4))/(part(6,1)-part(6,4))))

CCC Whether photon is down beam will be tested again later, do it here before recombination for efficiency
CCC Though the original (anti)lepton may fail the cut, A lepton recombined with photon may pass the it
CCC Since we do not know if photon will be recombined or not, we need to be careful and lepton pseudorapidity cut is also used
      tmpdRsq = max(maxetar4obs**2,
     &              lrdRrecb**2+max(dabs(lepetamax),dabs(alepetamax),dabs(lepetamin),dabs(alepetamin))**2)
      IF (DBLE(etar**2).gt.tmpdRsq) THEN ! photon down beam, set to zero
         DO j=1,4
            part(6,j) = 0d0
         ENDDO
         RETURN
      ENDIF

CCC perform photon recombination
CCC recombination is done regardless of whether photon cuts are turned on or off

      etal1 = 0.5d0*dlog(dabs((part(3,1)+part(3,4))/(part(3,1)-part(3,4))))
      etal2 = 0.5d0*dlog(dabs((part(2,1)+part(2,4))/(part(2,1)-part(2,4))))

C      if (phtrecbon) then

         deta = etal1-etar
         dphi = dabs(datan2(part(3,3),part(3,2))-
     &        datan2(part(6,3),part(6,2)))
         if (dphi.gt.pi) dphi = 2d0*pi-dphi
         dRsqrl1 = dphi**2+deta**2

         deta = etal2-etar
         dphi = dabs(datan2(part(2,3),part(2,2))-
     &        datan2(part(6,3),part(6,2)))
         if (dphi.gt.pi) dphi = 2d0*pi-dphi
         dRsqrl2 = dphi**2+deta**2

         IF (DBLE(dRsqrl2).le.DBLE(dRsqrl1)) THEN ! lepton has the priority of recombining with photon than anti-lepton
            IF (DBLE(dRsqrl2).lt.DBLE(lrdRrecb**2)) THEN
               rlflag = 1
               DO j=1,4
                  part(2,j) = part(2,j)+part(6,j)
                  part(1,j) = part(1,j)+part(6,j) ! reconstructed Z from lepton-pair
                  part(6,j) = 0d0
               ENDDO
            END IF
         ELSE
            IF (DBLE(dRsqrl1).lt.DBLE(lrdRrecb**2)) THEN
               rlflag = 2
               DO j=1,4
                  part(3,j) = part(3,j)+part(6,j)
                  part(1,j) = part(1,j)+part(6,j) ! reconstructed Z from lepton-pair
                  part(6,j) = 0d0
               ENDDO
            END IF
         END IF

C      endif ! phtrecbon

      IF ((part(6,1).ne.0).and.(dabs(etar).gt.maxetar4obs)) THEN ! photon down beam, set to zero
         DO j=1,4
            part(6,j) = 0d0
         ENDDO
      ENDIF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC Determine whether physical cuts are passed

      double precision function constraint(part)

      implicit none
      double precision pi
      parameter (pi=3.141592653589793d0)

      include 'cutcontrol.f'
      include 'particles.f'
      include 'pdfcontrol.f'
CCC      include 'parameters.f'

CCCCC Kinematics of the Z boson and leptons, jets

      double precision minv,tinv

      double precision ptw,Y,
     & pt3,eta3,
     & pt4,eta4,
     & ptl2,etal2,
     & ptl1,etal1,
     & ptr,etar,
     & deta,dphi

CCC   Parameters for ordering the leptons by pT
      double precision hardleppT,softleppT,hardlepeta,softlepeta

      integer numofjets,numofphts

CCC   Coin flip store
      double precision tmprand
      COMMON/randstore/tmprand

CCC   HISTOGRAM ROUTINE
c      external fill_hist

CCC ---------------------- Definition of particle 4-momenta -----------------
CCC   Particle number:
CCC     1. Z-boson
CCC     2. Lepton
CCC     3. Anti-lepton
CCC     4. Parton 1
CCC     5. Parton 2
CCC -------------------------------------------------------------------------

CCC   Note that z should not be used to form constraints; it is used as the 
CCC   convolution variable also, not just the mass

      constraint=1.0d0

CCC   Invariant mass cut

c      minv = dsqrt(m2)
      minv = dsqrt(part(1,1)**2
     &       -part(1,2)**2-part(1,3)**2-part(1,4)**2)

      if ((minv.lt.(invmassmin)).or.(minv.gt.(invmassmax))) then
         constraint = 0.0d0
         RETURN
      endif

CCC   Transverse Mass cut
      tinv = dsqrt(
     &     2d0 * dsqrt( (part(2,2)**2+part(2,3)**2) *
     &     (part(3,2)**2+part(3,3)**2) ) -
     &     2d0 * (part(2,2)*part(3,2) + part(2,3)*part(3,3))
     &     )

      if ((tinv.lt.(tmassmin)).or.(tinv.gt.(tmassmax))) then
         constraint = 0.0d0
         RETURN
      endif

      ptl2 = dsqrt(part(2,2)**2+part(2,3)**2)
      if ((ptl2.lt.leppTmin).or.(ptl2.gt.leppTmax)) then
         constraint = 0.0d0
         RETURN
      endif

      ptl1 = dsqrt(part(3,2)**2+part(3,3)**2)
      if ((ptl1.lt.aleppTmin).or.(ptl1.gt.aleppTmax)) then
         constraint = 0.0d0
         RETURN
      endif
     
      ptw = dsqrt(part(1,2)**2+part(1,3)**2)

CCC   Z pT, Y cuts

      if ((ptw.lt.ZpTmin).or.(ptw.gt.ZpTmax)) then
         constraint = 0.0d0
         RETURN
      endif

      etal2 = 0.5d0*dlog(dabs((part(2,1)+part(2,4))/(part(2,1)-part(2,4))))

CCC   Lepton eta cuts

      if (lepabseta.eq.1) then
         if ((dabs(etal2).lt.lepetamin).or.
     &        (dabs(etal2).gt.lepetamax)) then
            constraint = 0.0d0
            RETURN
         endif
      elseif (lepabseta.eq.0) then
         if ((etal2.lt.lepetamin).or.
     &        (etal2.gt.lepetamax)) then
            constraint = 0.0d0
            RETURN
         endif
      endif

      etal1 = 0.5d0*dlog(dabs((part(3,1)+part(3,4))/(part(3,1)-part(3,4))))

      if (alepabseta.eq.1) then
         if ((dabs(etal1).lt.alepetamin).or.
     &        (dabs(etal1).gt.alepetamax)) then
            constraint = 0.0d0
            RETURN
         endif
      elseif (alepabseta.eq.0) then
         if ((etal1.lt.alepetamin).or.
     &        (etal1.gt.alepetamax)) then
            constraint = 0.0d0
            RETURN
         endif
      endif

      Y = 0.5d0*dlog(dabs((part(1,1)+part(1,4))/(part(1,1)-part(1,4))))

      if ((Y.lt.ZYmin).or.(Y.gt.ZYmax)) then
         constraint = 0.0d0
         RETURN
      endif

      if (ptw.gt.ptbuf) then
         if (ptl1.gt.ptl2) then
            hardleppT = ptl1
            softleppT = ptl2
            hardlepeta = etal1
            softlepeta = etal2
         else
            hardleppT = ptl2
            softleppT = ptl1
            hardlepeta = etal2
            softlepeta = etal1
         endif
      else 
         if (tmprand.gt.5d-1) then
            hardleppT = ptl1
            softleppT = ptl2
            hardlepeta = etal1
            softlepeta = etal2
         else
            hardleppT = ptl2
            softleppT = ptl1
            hardlepeta = etal2
            softlepeta = etal1
         endif
      endif

CCC   Lepton pT cuts

      if ((hardleppT.lt.hardleppTmin).or.(hardleppT.gt.hardleppTmax)) 
     & then
         constraint = 0.0d0
         RETURN
      endif

      if ((softleppT.lt.softleppTmin).or.(softleppT.gt.softleppTmax)) 
     & then
         constraint = 0.0d0
         RETURN
      endif

CCC   Determine Number of Jets in Event

      pt3 = dsqrt(part(4,2)**2+part(4,3)**2)
      pt4 = dsqrt(part(5,2)**2+part(5,3)**2)

      numofjets = 0

      if (pt3.ge.minptj4obs) then
         numofjets = 1
         eta3 = 0.5d0*dlog(dabs((part(4,1)+part(4,4))/
     &        (part(4,1)-part(4,4))))
         if (pt4.ge.minptj4obs) then
            numofjets = 2
            eta4 = 0.5d0*dlog(dabs((part(5,1)+part(5,4))/
     &           (part(5,1)-part(5,4))))
         endif
      endif
         
      if ((numofjets.lt.minnumjet).or.
     &     (numofjets.gt.maxnumjet)) then
         constraint = 0d0
         RETURN
      endif

CCC  Jet pT cuts

cc  Leading jet

      if (pt3.lt.jet1ptmin.and.minnumjet.ge.1) then !.or.(pt3.gt.jet1ptmax)) then
         constraint = 0.0d0
         RETURN
      endif

cc  jet2
c      if ((pt4.lt.jet2ptmin).or.(pt4.gt.jet2ptmax)) then
c         constraint = 0.0d0
c         RETURN
c      endif

CCC   Jet eta cuts--removed and put into jet algorithm

cc   jet1
c      if (numOfJets.ge.1) then
c      if (jet1abseta.eq.1) then
c         if ((dabs(eta3).lt.jet1etamin).or.
c     &        (dabs(eta3).gt.jet1etamax)) then
c            constraint = 0.0d0
c            RETURN
c         endif
c      elseif (jet1abseta.eq.0) then
c         if ((eta3.lt.jet1etamin).or.
c     &        (eta3.gt.jet1etamax)) then
c            constraint = 0.0d0
c            RETURN
c         endif
c      endif
c      endif

cc   jet2
c      if (numOfJets.ge.2) then
c      if (jet2abseta.eq.1) then
c         if ((dabs(eta4).lt.jet2etamin).or.
c     &        (dabs(eta4).gt.jet2etamax)) then
c            constraint = 0.0d0
c            RETURN
c         endif
c      elseif (jet2abseta.eq.0) then
c         if ((eta4.lt.jet2etamin).or.
c     &        (eta4.gt.jet2etamax)) then
c            constraint = 0.0d0
c            RETURN
c         endif
c      endif
c      endif

CCC   Determine Number of Photon in Event

      ptr = dsqrt(part(6,2)**2+part(6,3)**2)

      numofphts = 0
      if (ptr.ge.minptr4obs) then
         numofphts = 1
         etar = 0.5d0*dlog(dabs((part(6,1)+part(6,4))/(part(6,1)-part(6,4))))
      endif

      if ((numofphts.lt.minnumpht).or.
     &     (numofphts.gt.maxnumpht)) then
         constraint = 0d0
         RETURN
      endif

CCC   Only differentiate soft, hard leptons for eta cut if
CCC   the soft lepton pt is less than the hard lepton cut
CCC   This is from the CMS analysis of Stoyan, Hwidong et al.

      if (softleppT.lt.hardleppTmin) then

         if (hardlepabseta.eq.1) then
            if ((dabs(hardlepeta).lt.hardlepetamin).or.
     &           (dabs(hardlepeta).gt.hardlepetamax)) then
               constraint = 0.0d0
               RETURN
            endif
         elseif (hardlepabseta.eq.0) then
            if ((hardlepeta.lt.hardlepetamin).or.
     &           (hardlepeta.gt.hardlepetamax)) then
               constraint = 0.0d0
               RETURN
            endif
         endif
         
         if (softlepabseta.eq.1) then
            if ((dabs(softlepeta).lt.softlepetamin).or.
     &           (dabs(softlepeta).gt.softlepetamax)) then
               constraint = 0.0d0
               RETURN
            endif
         elseif (softlepabseta.eq.0) then
            if ((softlepeta.lt.softlepetamin).or.
     &           (softlepeta.gt.softlepetamax)) then
               constraint = 0.0d0
               RETURN
            endif
         endif
         
      else

CCC   Need to differentiate the cases with absolute-value cuts
CCC   Need to have one of the two leptons pass both eta cuts;
CCC   in this case since both are greater then the hard pT cut,
CCC   it doesn't matter which lepton

         if ((softlepabseta.eq.1).and.(hardlepabseta.eq.1)) then
            if (((dabs(etal1).lt.hardlepetamin).or.(dabs(etal1).gt.
     &           hardlepetamax).or.(dabs(etal2).lt.softlepetamin).or.
     &           (dabs(etal2).gt.softlepetamax)).and.((dabs(etal2).lt.
     &           hardlepetamin).or.(dabs(etal2).gt.hardlepetamax).or.
     &           (dabs(etal1).lt.softlepetamin).or.(dabs(etal1).gt.
     &           softlepetamax))) then
               constraint = 0.0d0
               RETURN
            endif
         elseif ((softlepabseta.eq.0).and.(hardlepabseta.eq.1)) then
            if (((dabs(etal1).lt.hardlepetamin).or.(dabs(etal1).gt.
     &           hardlepetamax).or.(etal2.lt.softlepetamin).or.
     &           (etal2.gt.softlepetamax)).and.((dabs(etal2).lt.
     &           hardlepetamin).or.(dabs(etal2).gt.hardlepetamax).or.
     &           (etal1.lt.softlepetamin).or.(etal1.gt.
     &           softlepetamax))) then
               constraint = 0.0d0
               RETURN
            endif
         elseif ((softlepabseta.eq.1).and.(hardlepabseta.eq.0)) then
            if (((etal1.lt.hardlepetamin).or.(etal1.gt.
     &           hardlepetamax).or.(dabs(etal2).lt.softlepetamin).or.
     &           (dabs(etal2).gt.softlepetamax)).and.((etal2.lt.
     &           hardlepetamin).or.(etal2.gt.hardlepetamax).or.
     &           (dabs(etal1).lt.softlepetamin).or.(dabs(etal1).gt.
     &           softlepetamax))) then
               constraint = 0.0d0
               RETURN
            endif
         elseif ((softlepabseta.eq.0).and.(hardlepabseta.eq.0)) then
            if (((etal1.lt.hardlepetamin).or.(etal1.gt.
     &           hardlepetamax).or.(etal2.lt.softlepetamin).or.
     &           (etal2.gt.softlepetamax)).and.((etal2.lt.
     &           hardlepetamin).or.(etal2.gt.hardlepetamax).or.
     &           (etal1.lt.softlepetamin).or.(etal1.gt.
     &           softlepetamax))) then
               constraint = 0.0d0
               RETURN
            endif
         
         endif
      endif

CCC   Isolation cuts

ccc   lepton-lepton iso cut

      deta = 0.0d0
      dphi = 0.0d0

      deta = etal2-etal1
      dphi = dabs(datan2(part(2,3),part(2,2))-
     &     datan2(part(3,3),part(3,2)))
      if (dphi.gt.pi) then
         dphi = 2d0*pi-dphi
      endif

c      if (dsqrt(deta**2+dphi**2).lt.lldRmin) then
      IF (deta**2+dphi**2.lt.lldRmin**2) THEN
         constraint = 0.0d0
         RETURN
      endif

      if ((dphi.lt.lldPhimin).or.(dphi.gt.lldPhimax)) then
         constraint = 0.0d0
         RETURN
      endif

ccc   lepton-jet iso cut

cc    LEP & JET1
      IF (numOfJets.ge.1) THEN

      deta = 0.0d0
      dphi = 0.0d0

      deta = etal2-eta3
      dphi = dabs(datan2(part(2,3),part(2,2))-
     &     datan2(part(4,3),part(4,2)))
      if (dphi.gt.pi) then
         dphi = 2d0*pi-dphi
      endif

c      if (dsqrt(deta**2+dphi**2).lt.ljdRmin) then
      IF (deta**2+dphi**2.lt.ljdRmin**2) THEN
         constraint = 0.0d0
         RETURN
      endif

cc    aLEP & JET1
      deta = 0.0d0
      dphi = 0.0d0

      deta = etal1-eta3
      dphi = dabs(datan2(part(3,3),part(3,2))-
     &     datan2(part(4,3),part(4,2)))
      if (dphi.gt.pi) then
         dphi = 2d0*pi-dphi
      endif
      
c      if (dsqrt(deta**2+dphi**2).lt.ljdRmin) then
      IF (deta**2+dphi**2.lt.ljdRmin**2) THEN
         constraint = 0.0d0
         RETURN
      endif

      endif

cc    LEP & JET2
      if (numOfJets.ge.2) then

      deta = 0.0d0
      dphi = 0.0d0

      deta = etal2-eta4
      dphi = dabs(datan2(part(2,3),part(2,2))-
     &     datan2(part(5,3),part(5,2)))
      if (dphi.gt.pi) then
         dphi = 2d0*pi-dphi
      endif

c      if (dsqrt(deta**2+dphi**2).lt.ljdRmin) then
      IF (deta**2+dphi**2.lt.ljdRmin**2) THEN
         constraint = 0.0d0
         RETURN
      endif

cc    aLEP & JET2
      deta = 0.0d0
      dphi = 0.0d0

      deta = etal1-eta4
      dphi = dabs(datan2(part(3,3),part(3,2))-
     &     datan2(part(5,3),part(5,2)))
      if (dphi.gt.pi) then
         dphi = 2d0*pi-dphi
      endif
      
c      if (dsqrt(deta**2+dphi**2).lt.ljdRmin) then
      IF (deta**2+dphi**2.lt.ljdRmin**2) THEN
         constraint = 0.0d0
         RETURN
      endif

      endif

cc    LEP & PHT
      IF (numOfphts.ge.1) THEN

      deta = 0.0d0
      dphi = 0.0d0

      deta = etal2-etar
      dphi = dabs(datan2(part(2,3),part(2,2))-
     &     datan2(part(6,3),part(6,2)))
      if (dphi.gt.pi) then
         dphi = 2d0*pi-dphi
      endif

      if (deta**2+dphi**2.lt.lrdRmin**2) then
         constraint = 0.0d0
         RETURN
      endif

cc    aLEP & PHT
      deta = 0.0d0
      dphi = 0.0d0

      deta = etal1-etar
      dphi = dabs(datan2(part(3,3),part(3,2))-
     &     datan2(part(6,3),part(6,2)))
      if (dphi.gt.pi) then
         dphi = 2d0*pi-dphi
      endif

      if (deta**2+dphi**2.lt.lrdRmin**2) then
         constraint = 0.0d0
         RETURN
      endif

      ENDIF

      RETURN
      END

CCC -------------------------------------------------------------------------


