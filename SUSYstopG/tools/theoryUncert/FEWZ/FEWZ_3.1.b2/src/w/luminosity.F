
      double precision function lumin(x1,x2,fpflag,lum,iset)
      implicit none

      include 'parameters.f'
      INCLUDE 'pdfcontrol.f'

      double precision x1,x2,fpflag,lum
      integer iset

      double precision tlum
c      double precision asgz, asyl

      double precision dv1,dv2,ds1,ds2,uv1,uv2,us1,us2,sv1,sv2,ss1,ss2
c!      DOUBLE PRECISION str1,str2

      tlum = 0d0

CCC   Call the PDFs

      CALL getPDFs(x1,x2,muf,iset)

CCC   Define the up/down distributions for pp and ppbar colliders
CCC   1 is always the proton
CCC   variable ?v? such as uv1 means quark, variable ?s? such as us1 means anti-quark 

      if ((ptype.eq.1).or.(ptype.eq.3)) then
         dv1 = dval1+dsea1
         dv2 = dval2+dsea2
         ds1 = dsea1
         ds2 = dsea2
         uv1 = uval1+usea1
         uv2 = uval2+usea2
         us1 = usea1
         us2 = usea2
         sv1 = sval1+ssea1
         sv2 = sval2+ssea2
         ss1 = ssea1
         ss2 = ssea2
c      elseif ((ptype.eq.2).or.(ptype.eq.4)) then
       else    
         if (fpflag.lt.(0.0d0)) then
            dv1 = dval1+dsea1
            dv2 = dsea2
            ds1 = dsea1
            ds2 = dval2+dsea2
            uv1 = uval1+usea1
            uv2 = usea2
            us1 = usea1
            us2 = uval2+usea2
            sv1 = sval1+ssea1
            sv2 = ssea2
            ss1 = ssea1
            ss2 = sval2+ssea2
         else
            dv1 = dsea1
            dv2 = dval2+dsea2
            ds1 = dval1+dsea1
            ds2 = dsea2
            uv1 = usea1
            uv2 = uval2+usea2
            us1 = uval1+usea1
            us2 = usea2
            sv1 = ssea1
            sv2 = sval2+ssea2
            ss1 = sval1+ssea1
            ss2 = ssea2
         endif 
      endif

CCC   For now, assume no "valence" strange, override above

c!      str1 = ssea1
c!      str2 = ssea2

CCC   Ignore the override

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCC  W- luminosity  CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if ((ptype.eq.1).or.(ptype.eq.2)) then

CCC   All channels of the form qiqj,qbariqbarj,qiqbarj
CCC   Add in the quark couplings to vector bosons here

CCC   These are the v^2+a^2 pieces, not the axial pieces

CCC   All channels of the form qiqj->qiqi
      if (lum.eq.(1.0d0)) then
c!         tlum = Nww*(V11*uv1*dv2+V12*uv1*str2+V13*uv1*btm2
c!     &     +V12*chm1*dv2+V22*chm1*str2+V23*chm1*btm2+V11*ds1*us2
c!     &     +V12*ds1*chm2+V12*str1*us2+V22*str1*chm2+V13*btm1*us2
c!     &     +V23*btm1*chm2)
         tlum = Nww*(V11*uv1*dv2+V12*uv1*sv2+V13*uv1*btm2
     &     +V12*chm1*dv2+V22*chm1*sv2+V23*chm1*btm2+V11*ds1*us2
     &     +V12*ds1*chm2+V12*ss1*us2+V22*ss1*chm2+V13*btm1*us2
     &     +V23*btm1*chm2)
CCC   All channels of the form qiqi, qbariqbari
      elseif (lum.eq.(2.0d0)) then
c!         tlum = Nww*(V11*dv1*dv2+V12*dv1*dv2+V12*str1*str2
c!     &    +V22*str1*str2+V13*btm1*btm2+V23*btm1*btm2+us1*us2
c!     &    +chm1*chm2)
         tlum = Nww*(V11*dv1*dv2+V12*dv1*dv2+V12*sv1*sv2
     &    +V22*sv1*sv2+V13*btm1*btm2+V23*btm1*btm2+us1*us2
     &    +chm1*chm2)
CCC   All channels of the form qiqbari
      elseif (lum.eq.(3.0d0)) then
c!         tlum = Nww*(V11*dv1*us2+V12*dv1*chm2+V12*str1*us2
c!     &    +V22*str1*chm2+V13*btm1*us2+V23*btm1*chm2+V11*us1*dv2
c!     &    +V12*us1*str2+V13*us1*btm2+V12*chm1*dv2+V22*chm1*str2
c!     &    +V23*chm1*btm2)
         tlum = Nww*(V11*dv1*us2+V12*dv1*chm2+V12*sv1*us2
     &    +V22*sv1*chm2+V13*btm1*us2+V23*btm1*chm2+V11*us1*dv2
     &    +V12*us1*sv2+V13*us1*btm2+V12*chm1*dv2+V22*chm1*sv2
     &    +V23*chm1*btm2)
CCC   The gg channel
      elseif (lum.eq.(4.0d0)) then
         tlum = Nww*glu1*glu2*2.0d0
CCC   The qg channel
      elseif (lum.eq.(5.0d0)) then
c!          tlum = Nww*glu2*(V11*dv1+V12*dv1
c!     &    +V12*str1+V22*str1+V13*btm1+V23*btm1+us1+chm1)
          tlum = Nww*glu2*(V11*dv1+V12*dv1
     &    +V12*sv1+V22*sv1+V13*btm1+V23*btm1+us1+chm1)
CCC   The qi*qbari channel with up quarks; enters BC and a chunk of C^2
      elseif (lum.eq.(11.0d0)) then
         tlum = Nww*(uv1*us2+chm1*chm2)
CCC   The qi*qbari channel with down quarks; enters BC and a chunk of C^2
      elseif (lum.eq.(13.0d0)) then
c!         tlum = Nww*((V11+V12)*dv1*ds2+(V12+V22)*str1*str2
c!     &     +(V13+V23)*btm1*btm2)
         tlum = Nww*((V11+V12)*dv1*ds2+(V12+V22)*sv1*ss2
     &     +(V13+V23)*btm1*btm2)
CCC   The uidj->non-identical quarks channel
      elseif (lum.eq.(15.0d0)) then     
c!         tlum = Nww*((V12*uv1*dv2+V22*uv1*str2+V23*uv1*btm2
c!     &     +V11*chm1*dv2+V12*chm1*str2+V13*chm1*btm2)
c!     &     +((V12+V13)*ds1*us2+(V22+V23)*ds1*chm2+(V11+V13)
c!     &     *str1*us2+(V12+V23)*str1*chm2+(V11+V12)*btm1*us2
c!     &     +(V12+V22)*btm1*chm2))
         tlum = Nww*((V12*uv1*dv2+V22*uv1*sv2+V23*uv1*btm2
     &     +V11*chm1*dv2+V12*chm1*sv2+V13*chm1*btm2)
     &     +((V12+V13)*ds1*us2+(V22+V23)*ds1*chm2+(V11+V13)
     &     *ss1*us2+(V12+V23)*ss1*chm2+(V11+V12)*btm1*us2
     &     +(V12+V22)*btm1*chm2))
CCC   The didj+ubariubarj, i<>j channel.
      elseif (lum.eq.(17.0d0)) then
CCC   didj part
c!         tlum = Nww*((V11*dv1*str2+V12*dv1*str2+V11*dv1*btm2
c!     &     +V12*dv1*btm2+V12*str1*dv2+V22*str1*dv2+V12*str1*btm2
c!     &     +V22*str1*btm2+V13*btm1*dv2+V23*btm1*dv2+V13*btm1
c!     &     *str2+V23*btm1*str2)
         tlum = Nww*((V11*dv1*sv2+V12*dv1*sv2+V11*dv1*btm2
     &     +V12*dv1*btm2+V12*sv1*dv2+V22*sv1*dv2+V12*sv1*btm2
     &     +V22*sv1*btm2+V13*btm1*dv2+V23*btm1*dv2+V13*btm1
     &     *sv2+V23*btm1*sv2)
CCC   ubariubarj part
     &     +(us1*chm2+chm1*us2))
CCC   Terms of the form d*ubar->qi*qbarj, i<>j
      elseif (lum.eq.(19.0d0)) then
C Pieces with emission from down-type
c!         tlum = Nww*(V12*dv1*us2+V22*str1*us2+V23*btm1*us2
c!     &     +V11*dv1*chm2+V12*str1*chm2+V13*btm1*chm2)
         tlum = Nww*(V12*dv1*us2+V22*sv1*us2+V23*btm1*us2
     &     +V11*dv1*chm2+V12*sv1*chm2+V13*btm1*chm2)
      elseif (lum.eq.(25.0d0)) then
C Pieces with emission from anti-up type
c!         tlum = Nww*((V12+V13)*dv1*us2+(V11+V13)*str1*us2
c!     &     +(V11+V12)*btm1*us2+(V22+V23)*chm2*dv1+(V12+V23)*str1*chm2
c!     &     +(V12+V22)*btm1*chm2)
         tlum = Nww*((V12+V13)*dv1*us2+(V11+V13)*sv1*us2
     &     +(V11+V12)*btm1*us2+(V22+V23)*chm2*dv1+(V12+V23)*sv1*chm2
     &     +(V12+V22)*btm1*chm2)
CCC   Terms of the form ui*ubarj, i<>j
      elseif (lum.eq.(21.0d0)) then
         tlum = Nww*(uv1*chm2+chm1*us2)
CCC   Terms of the form di*dbarj, i<>j
      elseif (lum.eq.(23.0d0)) then
c!         tlum = Nww*((V11+V12)*dv1*(str2+btm2)
c!     &     +(V12+V22)*str1*(ds2+btm2)+(V13+V23)*btm1*(ds2+str2))
         tlum = Nww*((V11+V12)*dv1*(ss2+btm2)
     &     +(V12+V22)*sv1*(ds2+btm2)+(V13+V23)*btm1*(ds2+ss2))
CCC   B^2 terms
      elseif (lum.eq.(27.0d0)) then
c!         tlum = Nww*(dv1*ds2+uv1*us2+str1*str2+chm1*chm2+btm1*btm2)
         tlum = Nww*(dv1*ds2+uv1*us2+sv1*ss2+chm1*chm2+btm1*btm2)
CCC   BC terms, ddbar channel
      elseif (lum.eq.(29.0d0)) then
c!         tlum = Nww*((V11+V12)*dv1*ds2+(V12+V22)*str1*str2
c!     &     +(V13+V23)*btm1*btm2)
         tlum = Nww*((V11+V12)*dv1*ds2+(V12+V22)*sv1*ss2
     &     +(V13+V23)*btm1*btm2)
CCC   The BC, uubar channel
      elseif (lum.eq.(31.0d0)) then
         tlum = Nww*(uv1*us2+chm1*chm2)      
      endif  


CCC   These are the asymmetry pieces.  Need to multiply by gvl*gal/(gvl^2+gal^2);
CCC   have absorbed gvl^2+gal^2 into definition of GammaZll
CCC   Also need to multiply by gvq*gaq

CCC   Note that the sign of the terms with p2 as the quark should be minus, 
CCC   becaude of the anti-symemtry of the MEs

      if (lum.eq.(6.0d0)) then
c!         tlum = (0.25d0)*Nww*(V11*dv1*us2+V12*dv1*chm2+V12*str1*us2
c!     &    +V22*str1*chm2+V13*btm1*us2+V23*btm1*chm2-V11*us1*dv2
c!     &    -V12*us1*str2-V13*us1*btm2-V12*chm1*dv2-V22*chm1*str2
c!     &    -V23*chm1*btm2)
         tlum = (0.25d0)*Nww*(V11*dv1*us2+V12*dv1*chm2+V12*sv1*us2
     &    +V22*sv1*chm2+V13*btm1*us2+V23*btm1*chm2-V11*us1*dv2
     &    -V12*us1*sv2-V13*us1*btm2-V12*chm1*dv2-V22*chm1*sv2
     &    -V23*chm1*btm2)
      elseif (lum.eq.(7.0d0)) then
c!         tlum = (0.25d0)*Nww*glu2*(V11*dv1+V12*dv1
c!     &    +V12*str1+V22*str1+V13*btm1+V23*btm1-us1-chm1)
         tlum = (0.25d0)*Nww*glu2*(V11*dv1+V12*dv1
     &    +V12*sv1+V22*sv1+V13*btm1+V23*btm1-us1-chm1)
      elseif (lum.eq.(8.0d0)) then
         tlum = 0.0d0
      elseif (lum.eq.(9.0d0)) then 
c!         tlum = (0.25d0)*Nww*(V11*dv1*dv2+V12*dv1*dv2+V12*str1*str2
c!     &    +V22*str1*str2+V13*btm1*btm2+V23*btm1*btm2-V11*us1*us2
c!     &    -V12*us1*us2-V13*us1*us2-V12*chm1*chm2-V22*chm1*chm2
c!     &    -V23*chm1*chm2)
         tlum = (0.25d0)*Nww*(V11*dv1*dv2+V12*dv1*dv2+V12*sv1*sv2
     &    +V22*sv1*sv2+V13*btm1*btm2+V23*btm1*btm2-V11*us1*us2
     &    -V12*us1*us2-V13*us1*us2-V12*chm1*chm2-V22*chm1*chm2
     &    -V23*chm1*chm2)
      elseif (lum.eq.(10.0d0)) then 
c!         tlum = (0.25d0)*Nww*(V11*uv1*dv2+V12*uv1*str2+V13*uv1*btm2
c!     &     +V12*chm1*dv2+V22*chm1*str2+V23*chm1*btm2-V11*ds1*us2
c!     &     -V12*ds1*chm2-V12*str1*us2-V22*str1*chm2-V13*btm1*us2
c!     &     -V23*btm1*chm2)
         tlum = (0.25d0)*Nww*(V11*uv1*dv2+V12*uv1*sv2+V13*uv1*btm2
     &     +V12*chm1*dv2+V22*chm1*sv2+V23*chm1*btm2-V11*ds1*us2
     &     -V12*ds1*chm2-V12*ss1*us2-V22*ss1*chm2-V13*btm1*us2
     &     -V23*btm1*chm2)
      elseif (lum.eq.(12.0d0)) then
         tlum = (0.25d0)*Nww*(uv1*us2+chm1*chm2)
      elseif (lum.eq.(14.0d0)) then
c!         tlum = (0.25d0)*Nww*((V11+V12)*dv1*ds2+(V12+V22)*str1*str2
c!     &     +(V13+V23)*btm1*btm2)
         tlum = (0.25d0)*Nww*((V11+V12)*dv1*ds2+(V12+V22)*sv1*ss2
     &     +(V13+V23)*btm1*btm2)
      elseif (lum.eq.(16.0d0)) then
c!         tlum = (0.25d0)*Nww*((V12*uv1*dv2+V22*uv1*str2
c!     &     +V23*uv1*btm2+V11*chm1*dv2+V12*chm1*str2+V13*chm1*btm2)
c!     &     -((V12+V13)*ds1*us2+(V22+V23)*ds1*chm2+(V11+V13)*str1*us2
c!     &     +(V12+V23)*str1*chm2+(V11+V12)*btm1*us2+(V12+V22)*btm1*chm2))
         tlum = (0.25d0)*Nww*((V12*uv1*dv2+V22*uv1*sv2
     &     +V23*uv1*btm2+V11*chm1*dv2+V12*chm1*sv2+V13*chm1*btm2)
     &     -((V12+V13)*ds1*us2+(V22+V23)*ds1*chm2+(V11+V13)*ss1*us2
     &     +(V12+V23)*ss1*chm2+(V11+V12)*btm1*us2+(V12+V22)*btm1*chm2))
      elseif (lum.eq.(18.0d0)) then
c!         tlum = (0.25d0)*Nww*((V11*dv1*str2+V12*dv1*str2
c!     &     +V11*dv1*btm2+V12*dv1*btm2+V12*str1*dv2+V22*str1*dv2
c!     &     +V12*str1*btm2+V22*str1*btm2+V13*btm1*dv2+V23*btm1*dv2
c!     &     +V13*btm1*str2+V23*btm1*str2)
         tlum = (0.25d0)*Nww*((V11*dv1*sv2+V12*dv1*sv2
     &     +V11*dv1*btm2+V12*dv1*btm2+V12*sv1*dv2+V22*sv1*dv2
     &     +V12*sv1*btm2+V22*sv1*btm2+V13*btm1*dv2+V23*btm1*dv2
     &     +V13*btm1*sv2+V23*btm1*sv2)
CCC   ubariubarj part
     &    -(us1*chm2+chm1*us2))
      elseif (lum.eq.(20.0d0)) then
c!         tlum = (0.25d0)*Nww*(V12*dv1*us2+V22*str1*us2+V23*btm1*us2
c!     &     +V11*dv1*chm2+V12*str1*chm2+V13*btm1*chm2)
         tlum = (0.25d0)*Nww*(V12*dv1*us2+V22*sv1*us2+V23*btm1*us2
     &     +V11*dv1*chm2+V12*sv1*chm2+V13*btm1*chm2)
      elseif (lum.eq.(26.0d0)) then
c!         tlum = (0.25d0)*Nww*((V12+V13)*dv1*us2+(V11+V13)*str1*us2
c!     &     +(V11+V12)*btm1*us2+(V22+V23)*chm2*dv1+(V12+V23)*str1*chm2
c!     &     +(V12+V22)*btm1*chm2)
         tlum = (0.25d0)*Nww*((V12+V13)*dv1*us2+(V11+V13)*sv1*us2
     &     +(V11+V12)*btm1*us2+(V22+V23)*chm2*dv1+(V12+V23)*sv1*chm2
     &     +(V12+V22)*btm1*chm2)
      elseif (lum.eq.(22.0d0)) then
         tlum = (0.25d0)*Nww*(uv1*chm2+chm1*us2)
CCC   Terms of the form di*dbarj, i<>j
      elseif (lum.eq.(24.0d0)) then
c!         tlum = (0.25d0)*Nww*((V11+V12)*dv1*(str2+btm2)
c!     &     +(V12+V22)*str1*(ds2+btm2)+(V13+V23)*btm1*(ds2+str2))
         tlum = (0.25d0)*Nww*((V11+V12)*dv1*(ss2+btm2)
     &     +(V12+V22)*sv1*(ds2+btm2)+(V13+V23)*btm1*(ds2+ss2))
CCC  B2 channel
      elseif (lum.eq.(28.0d0)) then
c!         tlum = (0.25d0)*Nww*(dv1*ds2+uv1*us2+str1*str2
c!     &     +chm1*chm2+btm1*btm2)
         tlum = (0.25d0)*Nww*(dv1*ds2+uv1*us2+sv1*ss2
     &     +chm1*chm2+btm1*btm2)
CCC   The BC, ddbar channel   
      elseif (lum.eq.(30.0d0)) then
c!         tlum = (0.25d0)*Nww*((V11+V12)*dv1*ds2
c!     &     +(V12+V22)*str1*str2+(V13+V23)*btm1*btm2)
         tlum = (0.25d0)*Nww*((V11+V12)*dv1*ds2
     &     +(V12+V22)*sv1*ss2+(V13+V23)*btm1*btm2)
CCC   The BC, uubar channel
      elseif (lum.eq.(32.0d0)) then   
         tlum = (0.25d0)*Nww*(uv1*us2+chm1*chm2)
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCC  W+ luminosity  CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      elseif ((ptype.eq.3).or.(ptype.eq.4)) then

CCC   All channels of the form qiqj->qiqi; same as for W-
CCC   Note that the down-type should be made 1 and the up-type 2, 
CCC   to account for emission off different lines in the MEs
      if (lum.eq.(1.0d0)) then
c!          tlum = Nww*(V11*uv2*dv1+V12*uv2*str1+V13*uv2*btm1
c!     &     +V12*chm2*dv1+V22*chm2*str1+V23*chm2*btm1+V11*ds2*us1
c!     &     +V12*ds2*chm1+V12*str2*us1+V22*str2*chm1+V13*btm2*us1
c!     &     +V23*btm2*chm1)
          tlum = Nww*(V11*uv2*dv1+V12*uv2*sv1+V13*uv2*btm1
     &     +V12*chm2*dv1+V22*chm2*sv1+V23*chm2*btm1+V11*ds2*us1
     &     +V12*ds2*chm1+V12*ss2*us1+V22*ss2*chm1+V13*btm2*us1
     &     +V23*btm2*chm1)
CCC   All channels of the form qiqi, qbariqbari
      elseif (lum.eq.(2.0d0)) then
c!         tlum = Nww*(uv1*uv2+chm1*chm2+(V11+V12)*ds1*ds2
c!     &    +(V12+V22)*str1*str2+(V13+V23)*btm1*btm2)
         tlum = Nww*(uv1*uv2+chm1*chm2+(V11+V12)*ds1*ds2
     &    +(V12+V22)*ss1*ss2+(V13+V23)*btm1*btm2)
      elseif (lum.eq.(3.0d0)) then
c!         tlum = Nww*(V11*uv1*ds2+V12*uv1*str2+V13*uv1*btm2
c!     &    +V12*chm1*ds2+V22*chm1*str2+V23*chm1*btm2+V11*ds1*uv2
c!     &    +V12*ds1*chm2+V12*str1*uv2+V22*str1*chm2+V13*btm1*uv2
c!     &    +V23*btm1*chm2)
         tlum = Nww*(V11*uv1*ds2+V12*uv1*ss2+V13*uv1*btm2
     &    +V12*chm1*ds2+V22*chm1*ss2+V23*chm1*btm2+V11*ds1*uv2
     &    +V12*ds1*chm2+V12*ss1*uv2+V22*ss1*chm2+V13*btm1*uv2
     &    +V23*btm1*chm2)
CCC   The gg channel
      elseif (lum.eq.(4.0d0)) then
         tlum = Nww*glu1*glu2*2.0d0
CCC   The qg channel
      elseif (lum.eq.(5.0d0)) then
c!         tlum = Nww*glu2*(uv1+chm1
c!     &    +(V11+V12)*ds1+(V12+V22)*str1+(V13+V23)*btm1)
         tlum = Nww*glu2*(uv1+chm1
     &    +(V11+V12)*ds1+(V12+V22)*ss1+(V13+V23)*btm1)
CCC   The qi*qbari channel with up quarks; enters BC and a chunk of C^2
CCC   Note that for W+, this is 13, not 11 like for W-; needed since 
CCC   the matrix element structure changes
      elseif (lum.eq.(13.0d0)) then
         tlum = Nww*(uv1*us2+chm1*chm2)
CCC   The qi*qbari channel with down quarks; enters BC and a chunk of C^2
CCC   Have switched spots with uubar for W+
      elseif (lum.eq.(11.0d0)) then
c!         tlum = Nww*((V11+V12)*dv1*ds2+(V12+V22)*str1*str2
c!     &     +(V13+V23)*btm1*btm2)
         tlum = Nww*((V11+V12)*dv1*ds2+(V12+V22)*sv1*ss2
     &     +(V13+V23)*btm1*btm2)
CCC   The uidj->non-identical quarks channel
      elseif (lum.eq.(15.0d0)) then
c!         tlum = Nww*((V12+V13)*dv1*uv2+(V11+V13)*str1*uv2
c!     &     +(V11+V12)*btm1*uv2+(V22+V23)*dv1*chm2+(V12+V23)*str1*chm2
c!     &     +(V12+V22)*btm1*chm2+V12*us1*ds2+V11*chm1*ds2+V22*us1*str2
c!     &     +V12*chm1*str2+V23*us1*btm2+V13*chm1*btm2)
         tlum = Nww*((V12+V13)*dv1*uv2+(V11+V13)*sv1*uv2
     &     +(V11+V12)*btm1*uv2+(V22+V23)*dv1*chm2+(V12+V23)*sv1*chm2
     &     +(V12+V22)*btm1*chm2+V12*us1*ds2+V11*chm1*ds2+V22*us1*ss2
     &     +V12*chm1*ss2+V23*us1*btm2+V13*chm1*btm2)
CCC   The uiuj+dbaridbarj, i<>j channel.
      elseif (lum.eq.(17.0d0)) then
CCC   uiuj part
         tlum = Nww*(uv1*chm2+chm1*uv2
CCC   dbaridbarj part
c!     &    +((V11+V12)*ds1*(str2+btm2)+(V12+V22)*str1*
c!     &    (ds2+btm2)+(V13+V23)*btm1*(ds2+str2)))
     &    +((V11+V12)*ds1*(ss2+btm2)+(V12+V22)*ss1*
     &    (ds2+btm2)+(V13+V23)*btm1*(ds2+ss2)))
CCC   Terms of the form u*dbar->qi*qbarj, i<>j
      elseif (lum.eq.(19.0d0)) then
C Pieces with emission from up-type
c!         tlum = Nww*((V12+V13)*uv1*ds2+(V22+V23)*chm1*ds2
c!     &     +(V11+V13)*uv1*str2+(V12+V23)*chm1*str2+(V11+V12)*uv1*btm2
c!     &     +(V12+V22)*chm1*btm2)
         tlum = Nww*((V12+V13)*uv1*ds2+(V22+V23)*chm1*ds2
     &     +(V11+V13)*uv1*ss2+(V12+V23)*chm1*ss2+(V11+V12)*uv1*btm2
     &     +(V12+V22)*chm1*btm2)
      elseif (lum.eq.(25.0d0)) then
C Pieces with emission from anti-down type
c!         tlum = Nww*(V12*uv1*ds2+V11*chm1*ds2+V22*uv1*str2
c!     &     +V12*chm1*str2+V23*uv1*btm2+V13*chm1*btm2)
         tlum = Nww*(V12*uv1*ds2+V11*chm1*ds2+V22*uv1*ss2
     &     +V12*chm1*ss2+V23*uv1*btm2+V13*chm1*btm2)
CCC   Terms of the form di*dbarj, i<>j
CCC   Note that this is 21, not 23 like for W-; need to account 
CCC   emission from dbar, not d
      elseif (lum.eq.(21.0d0)) then
c!         tlum = Nww*((V11+V12)*(str1+btm1)*ds2
c!     &     +(V12+V22)*(dv1+btm1)*str2+(V13+V23)*(dv1+str1)*btm2)
         tlum = Nww*((V11+V12)*(sv1+btm1)*ds2
     &     +(V12+V22)*(dv1+btm1)*ss2+(V13+V23)*(dv1+sv1)*btm2)
CCC   Terms of the form ui*ubarj, i<>j
CCC   Note that this is 23, not 21 like for W-; need to account
CCC   emission from u, not ubar
      elseif (lum.eq.(23.0d0)) then
         tlum = Nww*(uv1*chm2+chm1*us2)
CCC   B^2 terms
      elseif (lum.eq.(27.0d0)) then
c!         tlum = Nww*(dv1*ds2+uv1*us2+str1*str2+chm1*chm2
c!     &     +btm1*btm2)
         tlum = Nww*(dv1*ds2+uv1*us2+sv1*ss2+chm1*chm2
     &     +btm1*btm2)
CCC   The BC, uubar channel
CCC   Note that this is 29, not 31 as for W-, because the structure changes 
      elseif (lum.eq.(29.0d0)) then
         tlum = Nww*(uv1*us2+chm1*chm2)
CCC   BC terms, ddbar channel
CCC   Note that this is 29, not 31 as for W-, because the structure changes
      elseif (lum.eq.(31.0d0)) then
c!         tlum = Nww*((V11+V12)*dv1*ds2+(V12+V22)*str1*str2
c!     &    +(V13+V23)*btm1*btm2)
         tlum = Nww*((V11+V12)*dv1*ds2+(V12+V22)*sv1*ss2
     &    +(V13+V23)*btm1*btm2)
      endif


      if (lum.eq.(6.0d0)) then
c!         tlum = (0.25d0)*Nww*(V11*uv1*ds2+V12*uv1*str2+V13*uv1*btm2
c!     &    +V12*chm1*ds2+V22*chm1*str2+V23*chm1*btm2-V11*ds1*uv2
c!     &    -V12*ds1*chm2-V12*str1*uv2-V22*str1*chm2-V13*btm1*uv2
c!     &    -V23*btm1*chm2)
         tlum = (0.25d0)*Nww*(V11*uv1*ds2+V12*uv1*ss2+V13*uv1*btm2
     &    +V12*chm1*ds2+V22*chm1*ss2+V23*chm1*btm2-V11*ds1*uv2
     &    -V12*ds1*chm2-V12*ss1*uv2-V22*ss1*chm2-V13*btm1*uv2
     &    -V23*btm1*chm2)
      elseif (lum.eq.(7.0d0)) then
c!         tlum = (0.25d0)*Nww*glu2*(uv1+chm1
c!     &    -(V11+V12)*ds1-(V12+V22)*str1-(V13+V23)*btm1)
         tlum = (0.25d0)*Nww*glu2*(uv1+chm1
     &    -(V11+V12)*ds1-(V12+V22)*ss1-(V13+V23)*btm1)
CCC   The gg channel
      elseif (lum.eq.(8.0d0)) then 
         tlum = 0.0d0
CCC   All channels of the form qiqi, qbariqbari
      elseif (lum.eq.(9.0d0)) then
c!         tlum = (0.25d0)*Nww*(uv1*uv2+chm1*chm2-(V11+V12)*ds1*ds2
c!     &    -(V12+V22)*str1*str2-(V13+V23)*btm1*btm2)
         tlum = (0.25d0)*Nww*(uv1*uv2+chm1*chm2-(V11+V12)*ds1*ds2
     &    -(V12+V22)*ss1*ss2-(V13+V23)*btm1*btm2)
CCC   Channels of the form qiqj->identical
CCC   Note that the down-type should be 1, up-type 2
      elseif (lum.eq.(10.0d0)) then
c!         tlum = (0.25d0)*Nww*(V11*uv2*dv1+V12*uv2*str1+V13*uv2*btm1
c!     &     +V12*chm2*dv1+V22*chm2*str1+V23*chm2*btm1-V11*ds2*us1
c!     &     -V12*ds2*chm1-V12*str2*us1-V22*str2*chm1-V13*btm2*us1
c!     &     -V23*btm2*chm1)
         tlum = (0.25d0)*Nww*(V11*uv2*dv1+V12*uv2*sv1+V13*uv2*btm1
     &     +V12*chm2*dv1+V22*chm2*sv1+V23*chm2*btm1-V11*ds2*us1
     &     -V12*ds2*chm1-V12*ss2*us1-V22*ss2*chm1-V13*btm2*us1
     &     -V23*btm2*chm1)
CCC   Channels of the form uubar
CCC   Note that for W+, this is 14, not 12 like for W-; needed since
CCC   the matrix element structure changes
      elseif (lum.eq.(14.0d0)) then
         tlum = (0.25d0)*Nww*(uv1*us2+chm1*chm2)
CCC   Channels of the form ddbar   
CCC   Note that for W+, this is 12, not 13 like for W-
      elseif (lum.eq.(12.0d0)) then
c!         tlum = (0.25d0)*Nww*((V11+V12)*dv1*ds2+(V12+V22)*str1*str2
c!     &     +(V13+V23)*btm1*btm2)
         tlum = (0.25d0)*Nww*((V11+V12)*dv1*ds2+(V12+V22)*sv1*ss2
     &     +(V13+V23)*btm1*btm2)
CCC   The uidj->non-identical quarks channel
      elseif (lum.eq.(16.0d0)) then
c!         tlum = (0.25d0)*Nww*((V12+V13)*dv1*uv2+(V11+V13)*str1*uv2
c!     &     +(V11+V12)*btm1*uv2+(V22+V23)*dv1*chm2+(V12+V23)*str1*chm2
c!     &     +(V12+V22)*btm1*chm2-V12*us1*ds2-V11*chm1*ds2-V22*us1*str2
c!     &     -V12*chm1*str2-V23*us1*btm2-V13*chm1*btm2)
         tlum = (0.25d0)*Nww*((V12+V13)*dv1*uv2+(V11+V13)*sv1*uv2
     &     +(V11+V12)*btm1*uv2+(V22+V23)*dv1*chm2+(V12+V23)*sv1*chm2
     &     +(V12+V22)*btm1*chm2-V12*us1*ds2-V11*chm1*ds2-V22*us1*ss2
     &     -V12*chm1*ss2-V23*us1*btm2-V13*chm1*btm2)
CCC   The uiuj+dbaridbarj, i<>j channel.
      elseif (lum.eq.(18.0d0)) then
CCC   uiuj part
         tlum = (0.25d0)*Nww*(uv1*chm2+chm1*uv2
CCC   dbaridbarj part
c!     &    -((V11+V12)*ds1*(str2+btm2)+(V12+V22)*str1*
c!     &    (ds2+btm2)+(V13+V23)*btm1*(ds2+str2)))
     &    -((V11+V12)*ds1*(ss2+btm2)+(V12+V22)*ss1*
     &    (ds2+btm2)+(V13+V23)*btm1*(ds2+ss2)))
CCC   Terms of the form u*dbar->qi*qbarj, i<>j
      elseif (lum.eq.(20.0d0)) then
C Pieces with emission from up-type
c!         tlum = (0.25d0)*Nww*((V12+V13)*uv1*ds2+(V22+V23)*chm1*ds2
c!     &     +(V11+V13)*uv1*str2+(V12+V23)*chm1*str2+(V11+V12)*uv1*btm2
c!     &     +(V12+V22)*chm1*btm2)
         tlum = (0.25d0)*Nww*((V12+V13)*uv1*ds2+(V22+V23)*chm1*ds2
     &     +(V11+V13)*uv1*ss2+(V12+V23)*chm1*ss2+(V11+V12)*uv1*btm2
     &     +(V12+V22)*chm1*btm2)
      elseif (lum.eq.(26.0d0)) then
C Pieces with emission from anti-down type
c!         tlum = (0.25d0)*Nww*(V12*uv1*ds2+V11*chm1*ds2+V22*uv1*str2
c!     &     +V12*chm1*str2+V23*uv1*btm2+V13*chm1*btm2)
         tlum = (0.25d0)*Nww*(V12*uv1*ds2+V11*chm1*ds2+V22*uv1*ss2
     &     +V12*chm1*ss2+V23*uv1*btm2+V13*chm1*btm2)
CCC   Terms of the form di*dbarj, i<>j
CCC   Note that this is 22, not 24 like for W-; need to account
CCC   emission from dbar, not d
      elseif (lum.eq.(22.0d0)) then
c!         tlum = (0.25d0)*Nww*((V11+V12)*(str1+btm1)*ds2
c!     &     +(V12+V22)*(dv1+btm1)*str2+(V13+V23)*(dv1+str1)*btm2)
         tlum = (0.25d0)*Nww*((V11+V12)*(sv1+btm1)*ds2
     &     +(V12+V22)*(dv1+btm1)*ss2+(V13+V23)*(dv1+sv1)*btm2)
CCC   Terms of the form ui*ubarj, i<>j
CCC   Note that this is 24, not 22 like for W-; need to account
CCC   emission from u, not ubar
      elseif (lum.eq.(24.0d0)) then
         tlum = (0.25d0)*Nww*(uv1*chm2+chm1*us2)
CCC  B2 channel
      elseif (lum.eq.(28.0d0)) then
c!         tlum = (0.25d0)*Nww*(dv1*ds2+uv1*us2+str1*str2
c!     &     +chm1*chm2+btm1*btm2)
         tlum = (0.25d0)*Nww*(dv1*ds2+uv1*us2+sv1*ss2
     &     +chm1*chm2+btm1*btm2)
CCC   The BC, uubar channel
CCC   Note that this is 30, not 32 as for W-
      elseif (lum.eq.(30.0d0)) then
         tlum = (0.25d0)*Nww*(uv1*us2+chm1*chm2)
CCC   The BC, ddbar channel
CCC   Note that this is 32, not 30 as for W-
      elseif (lum.eq.(32.0d0)) then
c!         tlum = (0.25d0)*Nww*((V11+V12)*dv1*ds2
c!     &     +(V12+V22)*str1*str2+(V13+V23)*btm1*btm2)
         tlum = (0.25d0)*Nww*((V11+V12)*dv1*ds2
     &     +(V12+V22)*sv1*ss2+(V13+V23)*btm1*btm2)
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      endif

      lumin = tlum

      END
