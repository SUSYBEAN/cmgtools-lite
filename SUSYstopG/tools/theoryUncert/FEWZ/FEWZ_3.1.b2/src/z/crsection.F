CCC   Auxiliary subroutine for Z production only
CCC   Define the normalization factor for the photon channel; 
CCC   it comes with a 1/m^4 instead of the Breit-Wigner
CCC   Also do this for the gamma-Z interference

      subroutine getNorm()
      IMPLICIT NONE
      INCLUDE 'parameters.f'
      double complex Ngz2
      COMMON/gzint/ Ngz2
      SAVE /gzint/

      if (photflag.eq.0) then
         Ngg = ((m2-Mz**2)**2+Mz**2*GammaZ**2)/m2**2
         Ngz2 = Ngz*(m2-conjg(Mz2))/m2
      else
         Ngg = 0d0
         Ngz2 = 0d0
      endif

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine LOe0(ndim,var,ncomp,ff)
      IMPLICIT double precision(a-z)
      INTEGER ndim,ncomp
      double precision var(ndim),ff(ncomp)

      INCLUDE 'pdfcontrol.f'

      double precision w1

      include 'parameters.f'

CCC   Need to fill m2, BWjac since no longer done in Fjet

      IF ((porder_ew.le.0).and.(alphaqed0.eq.0d0)) then ! for efficiency b/c EW module doesn't need it
         CALL getM2(Mz,GammaZ,var(11))
         CALL getNorm()
      ENDIF

CCC   Define the weight to be passed internally through function calls
      w1 = norm
      CoeffFuncOrder = 0

      if ((porder_ew.le.0).and.(alphaqed0.eq.0d0)) then
         t1 = norm*qqbLO1e0int(var,ndim,w1)
      else
         t1 = nloewval(var,ndim,1d0) !nloewval(var,ndim,w1/norm)
      endif
      ff(1) = t1

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine NLOe0(ndim,var,ncomp,ff)
      IMPLICIT double precision(a-z)
      INTEGER ndim,ncomp
      double precision var(ndim),ff(ncomp)

      INCLUDE 'pdfcontrol.f'

      double precision w1

      include 'parameters.f'

CCC   Define m2 here so that it can be used in the scale dependence;
CCC   this depends on ln(muf^2/m2)

      CALL getM2(Mz,GammaZ,var(11))
      CALL getNorm()

CCCCC LO and NLO pieces

      w1 = as*norm
      CoeffFuncOrder = 1

      t1 = as*norm*(
     & qNLOR1e0int(var,ndim,w1)+qNLOR2e0int(var,ndim,w1)+
     & qNLOR3e0int(var,ndim,w1)
     & )

      t2 = as*norm*(
     & qNLOV1e0int(var,ndim,w1)+qNLOV2e0int(var,ndim,w1)
     & )

      t3 = as*norm*(
     & qNLOc1e0int(var,ndim,w1)+qNLOc2e0int(var,ndim,w1)
     & )

      w1 = norm
      CoeffFuncOrder = 0

      if ((porder_ew.le.0).and.(alphaqed0.eq.0d0)) then
         t4 = norm*qqbLO1e0int(var,ndim,w1)
      else
         t4 = nloewval(var,ndim,1d0) !nloewval(var,ndim,w1/norm)
      endif

      ff(1) = t1+t2+t3+t4

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine NNLOe0(ndim,var,ncomp,ff)

      IMPLICIT double precision(a-z)
      INTEGER ndim,ncomp
      double precision var(ndim),ff(ncomp)
     
      double precision s0,sinf,y0,yinf,yjac,yvar
      double precision w1,oldw1,oldv1,oldv2,oldv3

      INCLUDE 'pdfcontrol.f'
      include 'sectorcontrol.f'
      include 'parameters.f'

CCC   Define the beta function for use in the muR<>muF scale dependence

      double precision NF,b0,Lrf

      IF ((isect.gt.1).or.((porder_ew.le.0).and.(alphaqed0.eq.0d0))) then ! for efficiency b/c EW module doesn't need it

CCC   Define m2 here so that it can be used in the scale dependence;
CCC   this depends on ln(muf^2/m2)

      CALL getM2(Mz,GammaZ,var(11))
      CALL getNorm()

CCC   Parameters for the muR<>muF scale dependence

      NF = 5d0
      b0 = 11d0/4d0-NF/6d0
      Lrf = dlog(muR**2/muF**2)

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if (isect.eq.1) then
         w1 = norm
         CoeffFuncOrder = 0
         if ((porder_ew.le.0).and.(alphaqed0.eq.0d0)) then
            t1 = norm*qqbLO1e0int(var,ndim,w1)
         else
            t1 = nloewval(var,ndim,1d0) !nloewval(var,ndim,w1/norm)
         endif
         ff(1) = t1
      elseif (isect.eq.2) then
         w1 = as*norm
         CoeffFuncOrder = 1
         t1 = w1*qNLOR1e0int(var,ndim,w1)
         w1 = as**2*norm*b0*Lrf
         CoeffFuncOrder = 2
         t2 = w1*qNLOR1e0int(var,ndim,w1)
         ff(1) = t1+t2
      elseif (isect.eq.3) then
         w1 = as*norm
         CoeffFuncOrder = 1
         t1 = w1*qNLOR2e0int(var,ndim,w1)
         w1 = as**2*norm*b0*Lrf
         CoeffFuncOrder = 2
         t2 = w1*qNLOR2e0int(var,ndim,w1)
         ff(1) = t1+t2
      elseif (isect.eq.4) then
         w1 = as*norm
         CoeffFuncOrder = 1
         t1 = w1*qNLOR3e0int(var,ndim,w1)
         w1 = as**2*norm*b0*Lrf
         CoeffFuncOrder = 2      
         t2 = w1*qNLOR3e0int(var,ndim,w1)
         ff(1) = t1+t2
      elseif (isect.eq.5) then
         w1 = as*norm
         CoeffFuncOrder = 1
         t1 = w1*(qNLOV1e0int(var,ndim,w1)+qNLOV2e0int(var,ndim,w1)
     &           +qNLOc2e0int(var,ndim,w1))
         w1 = as**2*norm*b0*Lrf
         CoeffFuncOrder = 2      
         t2 = w1*(qNLOV1e0int(var,ndim,w1)+qNLOV2e0int(var,ndim,w1)
     &           +qNLOc2e0int(var,ndim,w1))
         ff(1) = t1+t2
      elseif (isect.eq.6) then
         w1 = as*norm
         CoeffFuncOrder = 1
         t1 = w1*qNLOc1e0int(var,ndim,w1)
         w1 = as**2*norm*b0*Lrf
         CoeffFuncOrder = 2      
         t2 = w1*qNLOc1e0int(var,ndim,w1)
         ff(1) = t1+t2
      endif

CCCCC Smooth out variables at endpoints for NNLO sectors

      w1 = as**2*norm*6d0**4*var(1)*(1-var(1))*var(2)*(1-var(2))
     &     *var(3)*(1-var(3))*var(4)*(1-var(4))
      var(1) = 3*var(1)**2-2*var(1)**3
      var(2) = 3*var(2)**2-2*var(2)**3
      var(3) = 3*var(3)**2-2*var(3)**3
      var(4) = 3*var(4)**2-2*var(4)**3

      if ((isect.ge.7.and.isect.le.20).or. ! qiqi,qiqj,qiqbi
     &    (isect.eq.39).or. ! C2
     &    (isect.ge.41.and.isect.le.42).or. ! B2, BC
     &    (isect.ge.123) ! gg
     &   ) then ! do not smooth z->1, terms are linear or higher
      else ! terms logarithmic in z->1, smooth
         w1 = w1*2*(1-var(5))
         var(5) = 2*var(5)-var(5)**2
      endif

      oldw1=w1
      oldv1=var(1)
      oldv2=var(2)
      oldv3=var(3)
      CoeffFuncOrder = 2 ! all other sectors are O(as^2)

CCCCC qiqi RR pieces

      if (isect.eq.7) then
         w1=2*w1
         if (var(2).lt.5d-1) then
         var(2)=2*var(2)
         t1 = w1*qiqirm1e0int(var,ndim,w1)
         else
         var(2)=2*(1-var(2))
         t1 = w1*qiqirm2e0int(var,ndim,w1)
         endif
         ff(1) = t1
      elseif (isect.eq.8) then
         if (var(1).lt.var(2)) then
           var(1) = var(1)/var(2)
           w1=2*w1/var(2)
           if (var(3).lt.5d-1) then
             var(3)=2*var(3)
             t1 = w1*qiqirm3e0int(var,ndim,w1)
           else
             var(3)=2*(1-var(3))
             t1 = w1*qiqirm4e0int(var,ndim,w1)
           endif
         else
           var(2) = var(2)/var(1)
           w1=2*w1/var(1)
           if (var(3).lt.5d-1) then
             var(3)=2*var(3)
             t1 = w1*qiqirm5e0int(var,ndim,w1)
           else
             var(3)=2*(1-var(3))
             t1 = w1*qiqirm6e0int(var,ndim,w1)
           endif
         endif
         ff(1) = t1
      elseif (isect.eq.9) then
         t1 = w1*(qiqirm7e0int(var,ndim,w1)+qiqirm8e0int(var,ndim,w1)) ! 15+16
         ff(1) = t1

CCCCC THIS HAS BEEN CHANGED CCCCCCCC
CCCCC qiqi convols+scale; lump into 1 call
      elseif (isect.eq.10) then
         t1 = w1*(qiqic1e0int(var,ndim,w1)+qiqic2e0int(var,ndim,w1))
         ff(1) = t1
      elseif (isect.eq.11) then
         t1 = w1*qiqic3e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.12) then
         t1 = w1*(qiqisc1e0int(var,ndim,w1)+qiqisc2e0int(var,ndim,w1)
     &           +qiqisc3e0int(var,ndim,w1)+qiqisc4e0int(var,ndim,w1))
         ff(1) = t1
CCCCC qiqj RR piece
      elseif (isect.eq.13) then 
         t1 = w1*(qiqjrm1e0int(var,ndim,w1)+qiqjrm2e0int(var,ndim,w1)) ! 24+25
         ff(1) = t1
CCCCC THIS HAS BEEN CHANGED CCCCCCCC
CCCCC qiqj convols+scale; lump into 1 piece

      elseif (isect.eq.14) then
         t1 = w1*(qiqjc1e0int(var,ndim,w1)+qiqjc2e0int(var,ndim,w1)
     &           +qiqjc3e0int(var,ndim,w1)+qiqjc4e0int(var,ndim,w1)
     &           +qiqjc5e0int(var,ndim,w1)+qiqjc6e0int(var,ndim,w1))
         ff(1) = t1
      elseif (isect.eq.15) then
         t1 = w1*qiqjc7e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.16) then
         t1 = w1*qiqjc8e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.17) then
         t1 = w1*qiqjc9e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.18) then
         t1 = w1*(qiqjsc1e0int(var,ndim,w1)+qiqjsc2e0int(var,ndim,w1)
     &           +qiqjsc3e0int(var,ndim,w1)+qiqjsc10e0int(var,ndim,w1))
         ff(1) = t1
      elseif (isect.eq.19) then
         t1 = w1*(qiqjsc4e0int(var,ndim,w1)+qiqjsc5e0int(var,ndim,w1)
     &           +qiqjsc6e0int(var,ndim,w1)+qiqjsc11e0int(var,ndim,w1))
         ff(1) = t1
      elseif (isect.eq.20) then
         t1 = w1*(qiqjsc7e0int(var,ndim,w1)+qiqjsc8e0int(var,ndim,w1)
     &           +qiqjsc9e0int(var,ndim,w1)+qiqjsc12e0int(var,ndim,w1))
         ff(1) = t1

CCCCC qqb iq7 terms
      elseif (isect.ge.21.and.isect.le.22) then
         w1=w1/2
         t1 = w1*qbiq71e0int(var,ndim,w1)
         ff(1)=t1
      elseif (isect.ge.23.and.isect.le.26) then
         w1=w1/4
         t1 = w1*qbiq72e0int(var,ndim,w1)
         ff(1)=t1
      elseif (isect.eq.27) then
         t1 = w1*qbiq73e0int(var,ndim,w1)
         ff(1)=t1
      elseif (isect.eq.28) then
         t1 = w1*qbiq74e0int(var,ndim,w1)
         ff(1)=t1
      elseif (isect.eq.29) then
         t1 = w1*qbiq75e0int(var,ndim,w1)
         ff(1)=t1
      elseif (isect.eq.30) then
         t1 = w1*qbiq76e0int(var,ndim,w1)
         ff(1)=t1

CCCCC qqb sq3 terms
      elseif (isect.eq.31) then
         if (var(1).gt.var(2)) then
           var(2) = 1-var(2)
           var(1) = (1-var(1))/var(2)
           w1=w1/var(2)
           t1 = w1*sq31e0int(var,ndim,w1)
         else
           var(1) = 1-var(1)
           var(2) = (1-var(2))/var(1)
           w1=w1/var(1)
           t1 = w1*sq32e0int(var,ndim,w1)
         endif
         ff(1) = t1   
      elseif (isect.eq.32) then
c         var(3)=1-var(3)   
         if (var(1).gt.var(2)) then
           var(2) = 1-var(2)
           var(1) = (1-var(1))/var(2)
           w1=w1/var(2)
           t1 = w1*sq33e0int(var,ndim,w1)
         else
           var(1) = 1-var(1)
           var(2) = (1-var(2))/var(1)
           w1=w1/var(1)
           t1 = w1*sq34e0int(var,ndim,w1)
         endif
         ff(1) = t1   
      elseif (isect.eq.33) then
         if (var(1).lt.var(2)) then
           var(1) = var(1)/var(2)
           w1=w1/var(2)
           t1 = w1*sq35e0int(var,ndim,w1)
         else
           var(2) = var(2)/var(1)
           w1=w1/var(1)
           t1 = w1*sq36e0int(var,ndim,w1)
         endif
         ff(1) = t1   
      elseif (isect.eq.34) then
         if (var(1).lt.var(2)) then
           var(1) = var(1)/var(2)
           w1=w1/var(2)
           t1 = w1*sq37e0int(var,ndim,w1)
         else
           var(2) = var(2)/var(1)
           w1=w1/var(1)
           t1 = w1*sq38e0int(var,ndim,w1)
         endif
         ff(1) = t1   
      elseif (isect.eq.35) then ! split from sq311-2 for speed
c         if (var(1).gt.var(2)) then
c           var(2) = 1-var(2)
c           var(1) = (1-var(1))/var(2)
           w1=2*w1!/var(2)
           if (var(3).lt.5d-1) then
             var(3)=2*var(3)
             t1 = w1*sq39e0int(var,ndim,w1)
           else
             var(3)=2*(1-var(3))
             t1 = w1*sq310e0int(var,ndim,w1)
           endif
         ff(1) = t1
      elseif (isect.eq.36) then
c         else
c           var(1) = 1-var(1)
c           var(2) = (1-var(2))/var(1)
           w1=2*w1!/var(1)
           if (var(3).lt.5d-1) then
             var(3)=2*var(3)
             t1 = w1*sq311e0int(var,ndim,w1)
           else
             var(3)=2*(1-var(3))
             t1 = w1*sq312e0int(var,ndim,w1)
           endif
c         endif
         ff(1) = t1   
      elseif (isect.eq.37) then
         if (var(1).lt.var(2)) then
           var(1) = var(1)/var(2)
           w1=2*w1/var(2)
           if (var(3).lt.5d-1) then
             var(3)=2*var(3)
             t1 = w1*sq313e0int(var,ndim,w1)
           else
             var(3)=2*(1-var(3))
             t1 = w1*sq314e0int(var,ndim,w1)
           endif
         else
           var(2) = var(2)/var(1)
           w1=2*w1/var(1)
           if (var(3).lt.5d-1) then
             var(3)=2*var(3)
             t1 = w1*sq315e0int(var,ndim,w1)
           else
             var(3)=2*(1-var(3))
             t1 = w1*sq316e0int(var,ndim,w1)
           endif
         endif
         ff(1) = t1
CCCCC AC, C2, A2
      elseif (isect.eq.38) then
         w1=2*w1
         if (var(3).lt.5d-1) then   
         var(3)=2*var(3)
         t1 = w1*AC1e0int(var,ndim,w1)
         else
         var(3)=2*(1-var(3))
         t1 = w1*AC2e0int(var,ndim,w1)
         endif
         ff(1) = t1
      elseif (isect.eq.39) then
         t1 = w1*(C21e0int(var,ndim,w1)+C22e0int(var,ndim,w1)) ! 70+71
         ff(1) = t1      
      elseif (isect.eq.40) then
         w1=2*w1
       	 if (var(3).lt.5d-1) then
       	 var(3)=2*var(3)
         t1 = w1*(A21e0int(var,ndim,w1)+rm11e0int(var,ndim,w1)) ! 72+223
       	 else
       	 var(3)=2*(1-var(3))
       	 t1 = w1*(A22e0int(var,ndim,w1)+rm12e0int(var,ndim,w1)) ! 73+224
       	 endif

         w1=8*oldw1
         var(1)=oldv1
         var(2)=oldv2
         var(3)=oldv3
         if (var(1).lt.5d-1) then
           var(1)=2*oldv1
           if (var(2).lt.5d-1) then
             var(2)=2*oldv2
             if (var(3).lt.5d-1) then
               var(3)=2*oldv3
               t2 = w1*rm22e0int(var,ndim,w1)
             else
               var(3)=2*(1-oldv3)
               t2 = w1*rm23e0int(var,ndim,w1)
             endif
           else
             var(2)=2*(1-oldv2)
              if (var(3).lt.5d-1) then
               var(3)=2*oldv3
               t2 = w1*rm24e0int(var,ndim,w1)
             else
               var(3)=2*(1-oldv3)
               t2 = w1*rm25e0int(var,ndim,w1)
             endif
           endif
         else
           var(1)=2*(1-oldv1)
           if (var(2).lt.5d-1) then
             var(2)=2*oldv2
             if (var(3).lt.5d-1) then
               var(3)=2*oldv3
               t2 = w1*rm26e0int(var,ndim,w1)
             else
               var(3)=2*(1-oldv3)
               t2 = w1*rm27e0int(var,ndim,w1)
             endif
           else
             var(2)=2*(1-oldv2)
              if (var(3).lt.5d-1) then
               var(3)=2*oldv3
               t2 = w1*rm28e0int(var,ndim,w1)
             else
               var(3)=2*(1-oldv3)
               t2 = w1*rm29e0int(var,ndim,w1)
             endif
           endif
         endif
         ff(1) = t1+t2

CCCCC B2, BC, CD pieces
      elseif (isect.eq.41) then
         t1 = w1*B21e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.42) then
         t1 = w1*(BC1e0int(var,ndim,w1)+BC2e0int(var,ndim,w1))
         ff(1) = t1
      elseif (isect.eq.43) then
         if (var(1).gt.var(2)) then
           var(2) = 1-var(2)
           var(1) = (1-var(1))/var(2)
           w1=2*w1/var(2)
           if (var(3).lt.5d-1) then
             var(3)=2*var(3)
             t1 = w1*CD1e0int(var,ndim,w1)
           else
             var(3)=2*(1-var(3))
             t1 = w1*CD2e0int(var,ndim,w1)
           endif
         else
           var(1) = 1-var(1)
       	   var(2) = (1-var(2))/var(1)
       	   w1=2*w1/var(1)
       	   if (var(3).lt.5d-1) then
       	     var(3)=2*var(3)
             t1 = w1*CD3e0int(var,ndim,w1)
       	   else
       	     var(3)=2*(1-var(3))
       	     t1	= w1*CD4e0int(var,ndim,w1)
       	   endif
         endif
         ff(1) = t1      

CCCCC qqb rm terms
      elseif (isect.eq.44) then
         t1 = w1*rm13e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.45) then
         w1=2*w1
         if (var(2).lt.5d-1) then
           var(2)=2*var(2)
           t1 = w1*rm14e0int(var,ndim,w1)
         else
           var(2)=2*(1-var(2))
           t1 = w1*rm15e0int(var,ndim,w1)
         endif
         ff(1) = t1
      elseif (isect.eq.46) then
c         w1=2*w1
         w1=4*(1-var(5))*w1
         var(5) = 2*var(5)-var(5)**2 ! this sector likes extra z smoothing
       	 if (var(3).lt.5d-1) then
       	   var(3)=2*var(3)
           t1 = w1*rm16e0int(var,ndim,w1)
       	 else
       	   var(3)=2*(1-var(3))
       	   t1 =	w1*rm17e0int(var,ndim,w1)
       	 endif
         ff(1) = t1
      elseif (isect.eq.47) then
         t1 = w1*rm21e0int(var,ndim,w1)
         ff(1) = t1

CCCCC qqb RV terms
      elseif (isect.eq.48) then
c         var(1)=1-var(1)
         t1 = w1*RV1e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.ge.49.and.isect.le.50) then
         if (var(1).lt.5d-1) then   
         var(1)=2*var(1)
         t1 = w1*(RV2e0int(var,ndim,w1)+RV4e0int(var,ndim,w1)) ! 168+170  
         else
         var(1)=2*(1-var(1))
         t1 = w1*(RV3e0int(var,ndim,w1)+RV5e0int(var,ndim,w1)) ! 169+171
         endif
         ff(1) = t1
      elseif (isect.eq.51) then
         w1=2*w1
         if (var(1).lt.5d-1) then
         var(1)=2*var(1)
         t1 = w1*RV6e0int(var,ndim,w1)   
         else
         var(1)=2*(1-var(1))
         t1 = w1*RV7e0int(var,ndim,w1)
         endif
         ff(1) = t1
      elseif (isect.eq.52) then
c         var(1)=1-var(1)
         t1 = w1*RV8e0int(var,ndim,w1)   
         ff(1) = t1
      elseif (isect.eq.53) then
c         var(1)=1-var(1)
         t1 = w1*RV9e0int(var,ndim,w1)   
         ff(1) = t1
      elseif (isect.eq.54) then
         t1 = w1*RV10e0int(var,ndim,w1)   
         ff(1) = t1

CCCCC THIS HAS BEEN CHANGED CCCCCCCC
CCCCC qqb convols+scale; lump into 1 piece
      elseif (isect.eq.55) then
         t1 = w1*(qqbc1e0int(var,ndim,w1)+qqbc7e0int(var,ndim,w1)
     &           +qqbc12e0int(var,ndim,w1)+qqbc13e0int(var,ndim,w1))
         ff(1) = t1
      elseif (isect.eq.56) then
         t1 = w1*(qqbc2e0int(var,ndim,w1)+qqbc5e0int(var,ndim,w1)
     &           +qqbc8e0int(var,ndim,w1)+qqbc11e0int(var,ndim,w1))
         ff(1) = t1
      elseif (isect.eq.57) then
         t1 = w1*qqbc3e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.58) then
         t1 = w1*qqbc4e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.59) then
         t1 = w1*qqbc6e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.60) then
         t1 = w1*qqbc9e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.61) then
         t1 = w1*qqbc10e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.62) then
         t1 = w1*qqbc14e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.63) then
         t1 = w1*qqbc15e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.64) then
         t1 = w1*(qqbc16e0int(var,ndim,w1)+qqbc17e0int(var,ndim,w1)
     &           +qqbc18e0int(var,ndim,w1))
         ff(1) = t1
      elseif (isect.eq.65) then
         t1 = w1*qqbc19e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.66) then
         t1 = w1*qqbc20e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.67) then
         t1 = w1*qqbc21e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.68) then
         t1 = w1*(qqbsc1e0int(var,ndim,w1)+qqbsc2e0int(var,ndim,w1)
     &           +qqbsc3e0int(var,ndim,w1)+qqbsc4e0int(var,ndim,w1)
     &           +qqbsc7e0int(var,ndim,w1)+qqbsc10e0int(var,ndim,w1)
     &           +qqbsc15e0int(var,ndim,w1))
         ff(1) = t1
      elseif (isect.eq.69) then
         t1 = w1*qqbsc5e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.70) then
         t1 = w1*(qqbsc6e0int(var,ndim,w1)+qqbsc8e0int(var,ndim,w1)
     &           +qqbsc14e0int(var,ndim,w1)+qqbsc16e0int(var,ndim,w1)
     &           +qqbsc17e0int(var,ndim,w1))
         ff(1) = t1
      elseif (isect.eq.71) then
         t1 = w1*qqbsc9e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.72) then
         t1 = w1*(qqbsc11e0int(var,ndim,w1)+qqbsc23e0int(var,ndim,w1)) !208+220
         ff(1) = t1
      elseif (isect.eq.73) then
         t1 = w1*(qqbsc12e0int(var,ndim,w1)+qqbsc25e0int(var,ndim,w1)) !209+222
         ff(1) = t1
      elseif (isect.eq.74) then
         t1 = w1*qqbsc13e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.75) then
         t1 = w1*qqbsc18e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.76) then
         t1 = w1*qqbsc19e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.77) then
         t1 = w1*(qqbsc20e0int(var,ndim,w1)
     &           +qqbsc21e0int(var,ndim,w1)+qqbsc22e0int(var,ndim,w1)) !218+219
         ff(1) = t1
      elseif (isect.eq.78) then
         t1 = w1*qqbsc24e0int(var,ndim,w1)
         ff(1) = t1

CCCCC qg iq7 terms
      elseif (isect.eq.79) then
         t1 = w1*qgiq71e0int(var,ndim,w1)
         ff(1)=t1
      elseif (isect.eq.80) then
         t1 = w1*qgiq72e0int(var,ndim,w1)
         ff(1)=t1
      elseif (isect.eq.81) then
         w1=2*w1
         if (var(3).lt.5d-1) then
             var(3)=2*var(3)
             t1 = w1*qgiq73e0int(var,ndim,w1)
         else
             var(3)=2*(1-var(3))
             t1 = w1*qgiq74e0int(var,ndim,w1)
           endif
         ff(1)=t1
      elseif (isect.eq.82) then
         t1 = w1*qgiq75e0int(var,ndim,w1)
         ff(1)=t1
      elseif (isect.eq.83) then
         t1 = w1*qgiq76e0int(var,ndim,w1)
         ff(1)=t1

CCCCC qg sq3 terms
      elseif (isect.eq.84) then ! split from qgsq33-4 for speed
c         if (var(1).lt.var(2)) then
c           var(1) = var(1)/var(2)
           w1=2*w1!/var(2)
           if (var(3).lt.5d-1) then
             var(3)=2*var(3)
             t1 = w1*qgsq31e0int(var,ndim,w1)
           else
             var(3)=2*(1-var(3))
             t1 = w1*qgsq32e0int(var,ndim,w1)
           endif
         ff(1) = t1
c         else
      elseif (isect.eq.85.or.isect.eq.86) then
c           var(2) = var(2)/var(1)
c           w1=2*w1/var(1)
c           w1=2*w1/3
           if (var(3).lt.5d-1) then
             var(3)=2*var(3)
             t1 = w1*qgsq33e0int(var,ndim,w1)
           else
             var(3)=2*(1-var(3))
             t1 = w1*qgsq34e0int(var,ndim,w1)
           endif
c         endif
         ff(1) = t1
      elseif (isect.eq.87) then ! split from qgsq37-8 for speed          
c         if (var(1).gt.var(2)) then
c           var(2) = 1-var(2)
c           var(1) = (1-var(1))/var(2)
           w1=2*w1!/var(2)/2
           if (var(3).lt.5d-1) then
             var(3)=2*var(3)
             t1 = w1*qgsq35e0int(var,ndim,w1)
           else
             var(3)=2*(1-var(3))
             t1 = w1*qgsq36e0int(var,ndim,w1)
           endif
         ff(1) = t1
      elseif (isect.ge.88.and.isect.le.90) then        
c         else
c           var(1) = 1-var(1)
c           var(2) = (1-var(2))/var(1)
c           w1=w1/var(1)/2
           w1=2*w1/3
           if (var(3).lt.5d-1) then
             var(3)=2*var(3)
             t1 = w1*qgsq37e0int(var,ndim,w1)
           else
             var(3)=2*(1-var(3))
             t1 = w1*qgsq38e0int(var,ndim,w1)
           endif
c         endif
         ff(1) = t1
CCCCC qg rm terms
      elseif (isect.eq.91) then            
         t1 = w1*qgrm1e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.ge.92.and.isect.le.94) then  
         w1=w1/3
         t1 = w1*(qgrm2e0int(var,ndim,w1)+qgrm3e0int(var,ndim,w1))
         ff(1) = t1
      elseif (isect.eq.95) then            
         t1 = w1*qgrm4e0int(var,ndim,w1) 
         ff(1) = t1
      elseif (isect.eq.96) then            
         t1 = w1*qgrm5e0int(var,ndim,w1) 
         ff(1) = t1
      elseif (isect.eq.97) then         
         w1=2*w1
         if (var(3).lt.5d-1) then
         var(3)=2*var(3)
         t1 = w1*qgrm6e0int(var,ndim,w1)
         else
         var(3)=2*(1-var(3))
         t1 = w1*qgrm7e0int(var,ndim,w1)
         endif
         ff(1) = t1
      elseif (isect.eq.98) then            
         w1=2*w1
         if (var(2).lt.5d-1) then
         var(2)=2*var(2)
         t1 = w1*qgrm8e0int(var,ndim,w1)
         else
         var(2)=2*(1-var(2))
         t1 = w1*qgrm9e0int(var,ndim,w1)
         endif
         ff(1) = t1
      elseif (isect.eq.99) then            
         t1 = w1*qgrm10e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.100) then            
         w1=2*w1
         if (var(3).lt.5d-1) then
         var(3)=2*var(3)
         t1 = w1*qgrm11e0int(var,ndim,w1)
         else
         var(3)=2*(1-var(3))
         t1 = w1*qgrm12e0int(var,ndim,w1)
         endif
         ff(1) = t1
      elseif (isect.eq.101) then         
         t1 = w1*(qgrm13e0int(var,ndim,w1)+qgrm14e0int(var,ndim,w1)) !151+152
         ff(1) = t1

CCCCC qgRV piece
      elseif (isect.eq.102) then            
c         var(1)=1-var(1)
         t1 = w1*qgRV1e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.103) then            
         t1 = w1*qgRV2e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.104) then            
         t1 = w1*qgRV13e0int(var,ndim,w1)
         w1=2*w1
         if (var(1).lt.5d-1) then
         var(1)=2*var(1)
         t2 = w1*(qgRV3e0int(var,ndim,w1)+qgRV5e0int(var,ndim,w1)) !123+125
         else
         var(1)=2*(1-var(1))
         t2 = w1*(qgRV4e0int(var,ndim,w1)+qgRV6e0int(var,ndim,w1)) ! 124+126
         endif
         ff(1) = t1+t2
      elseif (isect.eq.105) then            
         t1 = w1*qgRV7e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.106) then            
         t1 = w1*qgRV8e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.107) then            
       	 w1=2*w1
       	 if (var(1).lt.5d-1) then
       	 var(1)=2*var(1)
         t1 = w1*(qgRV9e0int(var,ndim,w1)+qgRV11e0int(var,ndim,w1))
         else
       	 var(1)=2*(1-var(1))
         t1 = w1*(qgRV10e0int(var,ndim,w1)+qgRV12e0int(var,ndim,w1))
         endif
         ff(1) = t1
      elseif (isect.eq.108) then            
         t1 = w1*qgRV14e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.109) then            
         t1 = w1*(qgRV15e0int(var,ndim,w1)+qgRV17e0int(var,ndim,w1))
         ff(1) = t1
      elseif (isect.eq.110) then            
         t1 = w1*qgRV16e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.111) then            
         t1 = w1*qgRV18e0int(var,ndim,w1)
         ff(1) = t1
CCCCC THIS HAS BEEN CHANGED CCCCCCCC
CCCCC qg convols+scalel lump in 1 piece
      elseif (isect.eq.112) then            
         t1 = w1*(qgc1e0int(var,ndim,w1)+qgc2e0int(var,ndim,w1)
     &           +qgc3e0int(var,ndim,w1)+qgc4e0int(var,ndim,w1)
     &           +qgc5e0int(var,ndim,w1)+qgc6e0int(var,ndim,w1)
     &           +qgc7e0int(var,ndim,w1)) ! not sure if 7 should be in here
         ff(1) = t1
      elseif (isect.eq.113) then            
         t1 = w1*qgc8e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.114) then            
         t1 = w1*(qgc9e0int(var,ndim,w1)+qgc10e0int(var,ndim,w1)) ! not sure if good combo
         ff(1) = t1
      elseif (isect.eq.115) then            
         t1 = w1*qgc11e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.116) then            
         t1 = w1*qgc12e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.117) then            
         t1 = w1*qgc13e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.118) then            
         t1 = w1*qgc14e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.119) then            
         t1 = w1*qgc15e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.120) then            
         t1 = w1*qgc16e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.121) then            
         w1=2*w1
         if (var(1).lt.5d-1) then
         var(1)=2*var(1)
         t1 = w1*qgc17e0int(var,ndim,w1)
         else
         var(1)=2*(1-var(1))
         t1 = w1*qgc18e0int(var,ndim,w1)
         endif
         ff(1) = t1
      elseif (isect.eq.122) then            
         t1 = w1*(qgsc1e0int(var,ndim,w1)+qgsc2e0int(var,ndim,w1)
     &           +qgsc3e0int(var,ndim,w1)+qgsc4e0int(var,ndim,w1)
     &           +qgsc5e0int(var,ndim,w1)+qgsc6e0int(var,ndim,w1)
     &           +qgsc7e0int(var,ndim,w1)+qgsc8e0int(var,ndim,w1)
     &           +qgsc9e0int(var,ndim,w1)+qgsc10e0int(var,ndim,w1)
     &           +qgsc11e0int(var,ndim,w1)+qgsc12e0int(var,ndim,w1)
     &           +qgsc13e0int(var,ndim,w1)+qgsc14e0int(var,ndim,w1)
     &           +qgsc15e0int(var,ndim,w1)+qgsc16e0int(var,ndim,w1)
     &           +qgsc17e0int(var,ndim,w1)+qgsc18e0int(var,ndim,w1)
     &           +qgsc19e0int(var,ndim,w1)+qgsc20e0int(var,ndim,w1))   
       	 w1=2*w1
       	 if (var(1).lt.5d-1) then
       	 var(1)=2*var(1)
         t2 = w1*qgsc21e0int(var,ndim,w1)
       	 else
       	 var(1)=2*(1-var(1))
       	 t2 = w1*qgsc22e0int(var,ndim,w1)
       	 endif
         ff(1) = t1+t2

CCCCC gg RR piece
      elseif (isect.eq.123) then            
         t1 = w1*ggrm1e0int(var,ndim,w1)
         ff(1) = t1
      elseif (isect.eq.124) then         
         w1=8*w1
         if (var(1).lt.5d-1) then
           var(1)=2*var(1)
           if (var(2).lt.5d-1) then
             var(2)=2*var(2)
             if (var(3).lt.5d-1) then   
               var(3)=2*var(3)
               t1 = w1*ggrm2e0int(var,ndim,w1)  
             else
               var(3)=2*(1-var(3))
               t1 = w1*ggrm3e0int(var,ndim,w1)
             endif
           else
             var(2)=2*(1-var(2))
              if (var(3).lt.5d-1) then
               var(3)=2*var(3)
               t1 = w1*ggrm4e0int(var,ndim,w1)
       	     else
               var(3)=2*(1-var(3))
       	       t1 = w1*ggrm5e0int(var,ndim,w1)
       	     endif
           endif
         else
       	   var(1)=2*(1-var(1))
       	   if (var(2).lt.5d-1) then
       	     var(2)=2*var(2)
       	     if	(var(3).lt.5d-1) then
               var(3)=2*var(3)
               t1 = w1*ggrm6e0int(var,ndim,w1)
       	     else
               var(3)=2*(1-var(3))
       	       t1 = w1*ggrm7e0int(var,ndim,w1)
       	     endif
           else
             var(2)=2*(1-var(2)) 
              if (var(3).lt.5d-1) then
               var(3)=2*var(3)
               t1 = w1*ggrm8e0int(var,ndim,w1)
             else
               var(3)=2*(1-var(3))
               t1 = w1*ggrm9e0int(var,ndim,w1)
             endif
       	   endif
         endif       	 
         ff(1) = t1
      elseif (isect.eq.125) then            
         w1=2*w1
         if (var(2).lt.5d-1) then
         var(2)=2*var(2)
         t1 = w1*(ggrm10e0int(var,ndim,w1)+ggrm14e0int(var,ndim,w1)
     &           +ggrm16e0int(var,ndim,w1))
         else
         var(2)=2*(1-var(2))
         t1 = w1*(ggrm11e0int(var,ndim,w1)+ggrm15e0int(var,ndim,w1)
     &           +ggrm17e0int(var,ndim,w1))
         endif
         ff(1) = t1
      elseif (isect.eq.126) then            
         t1 = w1*(ggrm12e0int(var,ndim,w1)+ggrm13e0int(var,ndim,w1)) !58+59
         ff(1) = t1

CCCCC THIS HAS BEEN CHANGED CCCCCCCC
CCCCC gg convols+scale; lump into 1 piece
      elseif (isect.eq.127) then            
         t1 = w1*(ggc2e0int(var,ndim,w1)+ggsc2e0int(var,ndim,w1) !65+67
     &           +ggsc1e0int(var,ndim,w1))
         ff(1) = t1
      endif

      if ((isect.gt.127).or.(isect.lt.0)) then
         write(*,*) 'SECTOR NUMBER TOO LARGE'
      endif

      RETURN
      END
