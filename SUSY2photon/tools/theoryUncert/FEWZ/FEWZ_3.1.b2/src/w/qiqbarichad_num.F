CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqbicnum1e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t11 = z**2
      taux1 = -8*z*(1+z)*t5*t6+4.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux2 = 32*z*(1+z)*t5*t6-16.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoeff
     #*tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tcoe
     #ff*tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoeff
     #*tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tcoe
     #ff*tcoff4)

CCC Sum

      qiqbicnum1e0 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tc
     #off4
      RETURN
      END

      double precision function qiqbicnum1e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t7 = z**2
      taux1 = -8*z*(1+z)*t3+4.D0/3.D0*(z-1)*(4*t7+7*z+4)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux2 = -96*z*(1+z)*t5*t6+16*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoeff
     #*tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tcoe
     #ff*tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoeff
     #*tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tcoe
     #ff*tcoff4)

CCC Sum

      qiqbicnum1e1 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tc
     #off4
      RETURN
      END

      double precision function qiqbicnum1e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux1 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux2 = 64*z*(1+z)*t5*t6-32.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoeff
     #*tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tcoe
     #ff*tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoeff
     #*tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tcoe
     #ff*tcoff4)

CCC Sum

      qiqbicnum1e2 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tc
     #off4
      RETURN
      END

      double precision function qiqbicnum1e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux1 = 0
      taux2 = 0

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoeff
     #*tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tcoe
     #ff*tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoeff
     #*tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tcoe
     #ff*tcoff4)

CCC Sum

      qiqbicnum1e3 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tc
     #off4
      RETURN
      END

      double precision function qiqbicnum2e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t11 = z**2
      taux3 = -8*z*(1+z)*t5*t6+4.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux4 = 32*z*(1+z)*t5*t6-16.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,z,1.D0,Tcoeff
     #*tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,z,1.D0,Tcoe
     #ff*tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,z,1.D0,Tcoeff
     #*tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,z,1.D0,Tcoe
     #ff*tcoff8)

CCC Sum

      qiqbicnum2e0 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tc
     #off8
      RETURN
      END

      double precision function qiqbicnum2e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t7 = z**2
      taux3 = -8*z*(1+z)*t3+4.D0/3.D0*(z-1)*(4*t7+7*z+4)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux4 = -96*z*(1+z)*t5*t6+16*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,z,1.D0,Tcoeff
     #*tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,z,1.D0,Tcoe
     #ff*tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,z,1.D0,Tcoeff
     #*tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,z,1.D0,Tcoe
     #ff*tcoff8)

CCC Sum

      qiqbicnum2e1 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tc
     #off8
      RETURN
      END

      double precision function qiqbicnum2e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux3 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux4 = 64*z*(1+z)*t5*t6-32.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,z,1.D0,Tcoeff
     #*tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,z,1.D0,Tcoe
     #ff*tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,z,1.D0,Tcoeff
     #*tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,z,1.D0,Tcoe
     #ff*tcoff8)

CCC Sum

      qiqbicnum2e2 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tc
     #off8
      RETURN
      END

      double precision function qiqbicnum2e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux3 = 0
      taux4 = 0

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,z,1.D0,Tcoeff
     #*tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,z,1.D0,Tcoe
     #ff*tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,z,1.D0,Tcoeff
     #*tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,z,1.D0,Tcoe
     #ff*tcoff8)

CCC Sum

      qiqbicnum2e3 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tc
     #off8
      RETURN
      END

      double precision function qiqbicnum3e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux5 = 4*z*(1+z)*t5*t7-4.D0/3.D0*z*(3+15*z+8*t12)*t5*t6+8.D0/9.D0
     #*(z-1)*(28*t12+z+10)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux6 = -16*z*(1+z)*t5*t7+16.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-32.D0/
     #9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux5
      tcoff11 = taux6
      tcoff12 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoeff
     #*tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tco
     #eff*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoef
     #f*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tco
     #eff*tcoff12)

CCC Sum

      qiqbicnum3e0 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet
     #12*tcoff12
      RETURN
      END

      double precision function qiqbicnum3e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t4 = t3**2
      t8 = z**2
      taux5 = 4*z*(1+z)*t4-4.D0/3.D0*z*(3+15*z+8*t8)*t3+8.D0/9.D0*(z-1)*
     #(28*t8+z+10)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux6 = 48*z*(1+z)*t5*t7-16*z*(3+15*z+8*t12)*t5*t6+32.D0/3.D0*(z-1
     #)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux5
      tcoff11 = taux6
      tcoff12 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoeff
     #*tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tco
     #eff*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoef
     #f*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tco
     #eff*tcoff12)

CCC Sum

      qiqbicnum3e1 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet
     #12*tcoff12
      RETURN
      END

      double precision function qiqbicnum3e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux5 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux6 = -32*z*(1+z)*t5*t7+32.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-64.D0/
     #9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux5
      tcoff11 = taux6
      tcoff12 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoeff
     #*tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tco
     #eff*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoef
     #f*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tco
     #eff*tcoff12)

CCC Sum

      qiqbicnum3e2 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet
     #12*tcoff12
      RETURN
      END

      double precision function qiqbicnum3e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux5 = 0
      taux6 = 0

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux5
      tcoff11 = taux6
      tcoff12 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoeff
     #*tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tco
     #eff*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoef
     #f*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tco
     #eff*tcoff12)

CCC Sum

      qiqbicnum3e3 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet
     #12*tcoff12
      RETURN
      END

      double precision function qiqbicnum4e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux7 = 4*z*(1+z)*t5*t7-4.D0/3.D0*z*(3+15*z+8*t12)*t5*t6+8.D0/9.D0
     #*(z-1)*(28*t12+z+10)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux8 = -16*z*(1+z)*t5*t7+16.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-32.D0/
     #9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff13 = taux7
      tcoff14 = taux7
      tcoff15 = taux8
      tcoff16 = taux8

CCC Jet functions

      tFjet13 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,z,1.D0,Tcoef
     #f*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,z,1.D0,Tco
     #eff*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,z,1.D0,Tcoef
     #f*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,z,1.D0,Tco
     #eff*tcoff16)

CCC Sum

      qiqbicnum4e0 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFj
     #et16*tcoff16
      RETURN
      END

      double precision function qiqbicnum4e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t4 = t3**2
      t8 = z**2
      taux7 = 4*z*(1+z)*t4-4.D0/3.D0*z*(3+15*z+8*t8)*t3+8.D0/9.D0*(z-1)*
     #(28*t8+z+10)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux8 = 48*z*(1+z)*t5*t7-16*z*(3+15*z+8*t12)*t5*t6+32.D0/3.D0*(z-1
     #)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff13 = taux7
      tcoff14 = taux7
      tcoff15 = taux8
      tcoff16 = taux8

CCC Jet functions

      tFjet13 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,z,1.D0,Tcoef
     #f*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,z,1.D0,Tco
     #eff*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,z,1.D0,Tcoef
     #f*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,z,1.D0,Tco
     #eff*tcoff16)

CCC Sum

      qiqbicnum4e1 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFj
     #et16*tcoff16
      RETURN
      END

      double precision function qiqbicnum4e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux7 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux8 = -32*z*(1+z)*t5*t7+32.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-64.D0/
     #9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff13 = taux7
      tcoff14 = taux7
      tcoff15 = taux8
      tcoff16 = taux8

CCC Jet functions

      tFjet13 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,z,1.D0,Tcoef
     #f*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,z,1.D0,Tco
     #eff*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,z,1.D0,Tcoef
     #f*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,z,1.D0,Tco
     #eff*tcoff16)

CCC Sum

      qiqbicnum4e2 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFj
     #et16*tcoff16
      RETURN
      END

      double precision function qiqbicnum4e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux7 = 0
      taux8 = 0

CCC Coefficients

      tcoff13 = taux7
      tcoff14 = taux7
      tcoff15 = taux8
      tcoff16 = taux8

CCC Jet functions

      tFjet13 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,z,1.D0,Tcoef
     #f*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,z,1.D0,Tco
     #eff*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,z,1.D0,Tcoef
     #f*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,z,1.D0,Tco
     #eff*tcoff16)

CCC Sum

      qiqbicnum4e3 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFj
     #et16*tcoff16
      RETURN
      END

      double precision function qiqbicnum5e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t11 = z**2
      taux9 = -8*z*(1+z)*t5*t6+4.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux10 = 32*z*(1+z)*t5*t6-16.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff17 = taux9
      tcoff18 = taux9
      tcoff19 = taux10
      tcoff20 = taux10

CCC Jet functions

      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,15.D0,z,1.D0,Tcoef
     #f*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,15.D0,z,1.D0,Tco
     #eff*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,16.D0,z,1.D0,Tcoef
     #f*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,16.D0,z,1.D0,Tco
     #eff*tcoff20)

CCC Sum

      qiqbicnum5e0 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFj
     #et20*tcoff20
      RETURN
      END

      double precision function qiqbicnum5e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t7 = z**2
      taux9 = -8*z*(1+z)*t3+4.D0/3.D0*(z-1)*(4*t7+7*z+4)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux10 = -96*z*(1+z)*t5*t6+16*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff17 = taux9
      tcoff18 = taux9
      tcoff19 = taux10
      tcoff20 = taux10

CCC Jet functions

      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,15.D0,z,1.D0,Tcoef
     #f*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,15.D0,z,1.D0,Tco
     #eff*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,16.D0,z,1.D0,Tcoef
     #f*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,16.D0,z,1.D0,Tco
     #eff*tcoff20)

CCC Sum

      qiqbicnum5e1 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFj
     #et20*tcoff20
      RETURN
      END

      double precision function qiqbicnum5e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux9 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux10 = 64*z*(1+z)*t5*t6-32.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff17 = taux9
      tcoff18 = taux9
      tcoff19 = taux10
      tcoff20 = taux10

CCC Jet functions

      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,15.D0,z,1.D0,Tcoef
     #f*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,15.D0,z,1.D0,Tco
     #eff*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,16.D0,z,1.D0,Tcoef
     #f*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,16.D0,z,1.D0,Tco
     #eff*tcoff20)

CCC Sum

      qiqbicnum5e2 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFj
     #et20*tcoff20
      RETURN
      END

      double precision function qiqbicnum5e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux9 = 0
      taux10 = 0

CCC Coefficients

      tcoff17 = taux9
      tcoff18 = taux9
      tcoff19 = taux10
      tcoff20 = taux10

CCC Jet functions

      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,15.D0,z,1.D0,Tcoef
     #f*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,15.D0,z,1.D0,Tco
     #eff*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,16.D0,z,1.D0,Tcoef
     #f*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,16.D0,z,1.D0,Tco
     #eff*tcoff20)

CCC Sum

      qiqbicnum5e3 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFj
     #et20*tcoff20
      RETURN
      END

      double precision function qiqbicnum6e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux11 = 4*z*(1+z)*t5*t7-4.D0/3.D0*z*(3+15*z+8*t12)*t5*t6+8.D0/9.D
     #0*(z-1)*(28*t12+z+10)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux12 = -16*z*(1+z)*t5*t7+16.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-32.D0
     #/9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff21 = taux11
      tcoff22 = taux11
      tcoff23 = taux12
      tcoff24 = taux12

CCC Jet functions

      tFjet21 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,15.D0,z,1.D0,Tcoef
     #f*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,15.D0,z,1.D0,Tco
     #eff*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,16.D0,z,1.D0,Tcoef
     #f*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,16.D0,z,1.D0,Tco
     #eff*tcoff24)

CCC Sum

      qiqbicnum6e0 = tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFj
     #et24*tcoff24
      RETURN
      END

      double precision function qiqbicnum6e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t4 = t3**2
      t8 = z**2
      taux11 = 4*z*(1+z)*t4-4.D0/3.D0*z*(3+15*z+8*t8)*t3+8.D0/9.D0*(z-1)
     #*(28*t8+z+10)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux12 = 48*z*(1+z)*t5*t7-16*z*(3+15*z+8*t12)*t5*t6+32.D0/3.D0*(z-
     #1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff21 = taux11
      tcoff22 = taux11
      tcoff23 = taux12
      tcoff24 = taux12

CCC Jet functions

      tFjet21 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,15.D0,z,1.D0,Tcoef
     #f*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,15.D0,z,1.D0,Tco
     #eff*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,16.D0,z,1.D0,Tcoef
     #f*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,16.D0,z,1.D0,Tco
     #eff*tcoff24)

CCC Sum

      qiqbicnum6e1 = tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFj
     #et24*tcoff24
      RETURN
      END

      double precision function qiqbicnum6e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux11 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux12 = -32*z*(1+z)*t5*t7+32.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-64.D0
     #/9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff21 = taux11
      tcoff22 = taux11
      tcoff23 = taux12
      tcoff24 = taux12

CCC Jet functions

      tFjet21 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,15.D0,z,1.D0,Tcoef
     #f*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,15.D0,z,1.D0,Tco
     #eff*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,16.D0,z,1.D0,Tcoef
     #f*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,16.D0,z,1.D0,Tco
     #eff*tcoff24)

CCC Sum

      qiqbicnum6e2 = tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFj
     #et24*tcoff24
      RETURN
      END

      double precision function qiqbicnum6e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux11 = 0
      taux12 = 0

CCC Coefficients

      tcoff21 = taux11
      tcoff22 = taux11
      tcoff23 = taux12
      tcoff24 = taux12

CCC Jet functions

      tFjet21 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,15.D0,z,1.D0,Tcoef
     #f*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,15.D0,z,1.D0,Tco
     #eff*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,16.D0,z,1.D0,Tcoef
     #f*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,16.D0,z,1.D0,Tco
     #eff*tcoff24)

CCC Sum

      qiqbicnum6e3 = tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFj
     #et24*tcoff24
      RETURN
      END

      double precision function qiqbicnum7e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t11 = z**2
      taux13 = -8*z*(1+z)*t5*t6+4.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux14 = 32*z*(1+z)*t5*t6-16.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff25 = taux13
      tcoff26 = taux13
      tcoff27 = taux14
      tcoff28 = taux14

CCC Jet functions

      tFjet25 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,17.D0,1.D0,z,Tco
     #eff*tcoff25)
      tFjet26 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,17.D0,1.D0,z,Tcoef
     #f*tcoff26)
      tFjet27 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,18.D0,1.D0,z,Tco
     #eff*tcoff27)
      tFjet28 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,18.D0,1.D0,z,Tcoef
     #f*tcoff28)

CCC Sum

      qiqbicnum7e0 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFj
     #et28*tcoff28
      RETURN
      END

      double precision function qiqbicnum7e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t7 = z**2
      taux13 = -8*z*(1+z)*t3+4.D0/3.D0*(z-1)*(4*t7+7*z+4)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux14 = -96*z*(1+z)*t5*t6+16*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff25 = taux13
      tcoff26 = taux13
      tcoff27 = taux14
      tcoff28 = taux14

CCC Jet functions

      tFjet25 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,17.D0,1.D0,z,Tco
     #eff*tcoff25)
      tFjet26 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,17.D0,1.D0,z,Tcoef
     #f*tcoff26)
      tFjet27 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,18.D0,1.D0,z,Tco
     #eff*tcoff27)
      tFjet28 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,18.D0,1.D0,z,Tcoef
     #f*tcoff28)

CCC Sum

      qiqbicnum7e1 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFj
     #et28*tcoff28
      RETURN
      END

      double precision function qiqbicnum7e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux13 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux14 = 64*z*(1+z)*t5*t6-32.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff25 = taux13
      tcoff26 = taux13
      tcoff27 = taux14
      tcoff28 = taux14

CCC Jet functions

      tFjet25 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,17.D0,1.D0,z,Tco
     #eff*tcoff25)
      tFjet26 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,17.D0,1.D0,z,Tcoef
     #f*tcoff26)
      tFjet27 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,18.D0,1.D0,z,Tco
     #eff*tcoff27)
      tFjet28 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,18.D0,1.D0,z,Tcoef
     #f*tcoff28)

CCC Sum

      qiqbicnum7e2 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFj
     #et28*tcoff28
      RETURN
      END

      double precision function qiqbicnum7e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux13 = 0
      taux14 = 0

CCC Coefficients

      tcoff25 = taux13
      tcoff26 = taux13
      tcoff27 = taux14
      tcoff28 = taux14

CCC Jet functions

      tFjet25 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,17.D0,1.D0,z,Tco
     #eff*tcoff25)
      tFjet26 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,17.D0,1.D0,z,Tcoef
     #f*tcoff26)
      tFjet27 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,18.D0,1.D0,z,Tco
     #eff*tcoff27)
      tFjet28 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,18.D0,1.D0,z,Tcoef
     #f*tcoff28)

CCC Sum

      qiqbicnum7e3 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFj
     #et28*tcoff28
      RETURN
      END

      double precision function qiqbicnum8e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux15 = 4*z*(1+z)*t5*t7-4.D0/3.D0*z*(3+15*z+8*t12)*t5*t6+8.D0/9.D
     #0*(z-1)*(28*t12+z+10)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux16 = -16*z*(1+z)*t5*t7+16.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-32.D0
     #/9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff29 = taux15
      tcoff30 = taux15
      tcoff31 = taux16
      tcoff32 = taux16

CCC Jet functions

      tFjet29 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,17.D0,1.D0,z,Tco
     #eff*tcoff29)
      tFjet30 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,17.D0,1.D0,z,Tcoef
     #f*tcoff30)
      tFjet31 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,18.D0,1.D0,z,Tco
     #eff*tcoff31)
      tFjet32 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,18.D0,1.D0,z,Tcoef
     #f*tcoff32)

CCC Sum

      qiqbicnum8e0 = tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFj
     #et32*tcoff32
      RETURN
      END

      double precision function qiqbicnum8e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t4 = t3**2
      t8 = z**2
      taux15 = 4*z*(1+z)*t4-4.D0/3.D0*z*(3+15*z+8*t8)*t3+8.D0/9.D0*(z-1)
     #*(28*t8+z+10)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux16 = 48*z*(1+z)*t5*t7-16*z*(3+15*z+8*t12)*t5*t6+32.D0/3.D0*(z-
     #1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff29 = taux15
      tcoff30 = taux15
      tcoff31 = taux16
      tcoff32 = taux16

CCC Jet functions

      tFjet29 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,17.D0,1.D0,z,Tco
     #eff*tcoff29)
      tFjet30 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,17.D0,1.D0,z,Tcoef
     #f*tcoff30)
      tFjet31 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,18.D0,1.D0,z,Tco
     #eff*tcoff31)
      tFjet32 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,18.D0,1.D0,z,Tcoef
     #f*tcoff32)

CCC Sum

      qiqbicnum8e1 = tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFj
     #et32*tcoff32
      RETURN
      END

      double precision function qiqbicnum8e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux15 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux16 = -32*z*(1+z)*t5*t7+32.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-64.D0
     #/9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff29 = taux15
      tcoff30 = taux15
      tcoff31 = taux16
      tcoff32 = taux16

CCC Jet functions

      tFjet29 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,17.D0,1.D0,z,Tco
     #eff*tcoff29)
      tFjet30 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,17.D0,1.D0,z,Tcoef
     #f*tcoff30)
      tFjet31 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,18.D0,1.D0,z,Tco
     #eff*tcoff31)
      tFjet32 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,18.D0,1.D0,z,Tcoef
     #f*tcoff32)

CCC Sum

      qiqbicnum8e2 = tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFj
     #et32*tcoff32
      RETURN
      END

      double precision function qiqbicnum8e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux15 = 0
      taux16 = 0

CCC Coefficients

      tcoff29 = taux15
      tcoff30 = taux15
      tcoff31 = taux16
      tcoff32 = taux16

CCC Jet functions

      tFjet29 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,17.D0,1.D0,z,Tco
     #eff*tcoff29)
      tFjet30 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,17.D0,1.D0,z,Tcoef
     #f*tcoff30)
      tFjet31 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,18.D0,1.D0,z,Tco
     #eff*tcoff31)
      tFjet32 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,18.D0,1.D0,z,Tcoef
     #f*tcoff32)

CCC Sum

      qiqbicnum8e3 = tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFj
     #et32*tcoff32
      RETURN
      END

      double precision function qiqbicnum9e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t11 = z**2
      taux17 = -8*z*(1+z)*t5*t6+4.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux18 = 32*z*(1+z)*t5*t6-16.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff33 = taux17
      tcoff34 = taux17
      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet33 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,19.D0,1.D0,z,Tco
     #eff*tcoff33)
      tFjet34 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,19.D0,1.D0,z,Tcoef
     #f*tcoff34)
      tFjet35 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,20.D0,1.D0,z,Tco
     #eff*tcoff35)
      tFjet36 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,20.D0,1.D0,z,Tcoef
     #f*tcoff36)

CCC Sum

      qiqbicnum9e0 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFj
     #et36*tcoff36
      RETURN
      END

      double precision function qiqbicnum9e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t7 = z**2
      taux17 = -8*z*(1+z)*t3+4.D0/3.D0*(z-1)*(4*t7+7*z+4)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux18 = -96*z*(1+z)*t5*t6+16*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff33 = taux17
      tcoff34 = taux17
      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet33 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,19.D0,1.D0,z,Tco
     #eff*tcoff33)
      tFjet34 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,19.D0,1.D0,z,Tcoef
     #f*tcoff34)
      tFjet35 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,20.D0,1.D0,z,Tco
     #eff*tcoff35)
      tFjet36 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,20.D0,1.D0,z,Tcoef
     #f*tcoff36)

CCC Sum

      qiqbicnum9e1 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFj
     #et36*tcoff36
      RETURN
      END

      double precision function qiqbicnum9e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux17 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux18 = 64*z*(1+z)*t5*t6-32.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff33 = taux17
      tcoff34 = taux17
      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet33 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,19.D0,1.D0,z,Tco
     #eff*tcoff33)
      tFjet34 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,19.D0,1.D0,z,Tcoef
     #f*tcoff34)
      tFjet35 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,20.D0,1.D0,z,Tco
     #eff*tcoff35)
      tFjet36 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,20.D0,1.D0,z,Tcoef
     #f*tcoff36)

CCC Sum

      qiqbicnum9e2 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFj
     #et36*tcoff36
      RETURN
      END

      double precision function qiqbicnum9e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux17 = 0
      taux18 = 0

CCC Coefficients

      tcoff33 = taux17
      tcoff34 = taux17
      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet33 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,19.D0,1.D0,z,Tco
     #eff*tcoff33)
      tFjet34 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,19.D0,1.D0,z,Tcoef
     #f*tcoff34)
      tFjet35 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,20.D0,1.D0,z,Tco
     #eff*tcoff35)
      tFjet36 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,20.D0,1.D0,z,Tcoef
     #f*tcoff36)

CCC Sum

      qiqbicnum9e3 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFj
     #et36*tcoff36
      RETURN
      END

      double precision function qiqbicnum10e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t11 = z**2
      taux19 = -8*z*(1+z)*t5*t6+4.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux20 = 32*z*(1+z)*t5*t6-16.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff37 = taux19
      tcoff38 = taux19
      tcoff39 = taux20
      tcoff40 = taux20

CCC Jet functions

      tFjet37 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,25.D0,z,1.D0,Tcoef
     #f*tcoff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,25.D0,z,1.D0,Tco
     #eff*tcoff38)
      tFjet39 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,26.D0,z,1.D0,Tcoef
     #f*tcoff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,26.D0,z,1.D0,Tco
     #eff*tcoff40)

CCC Sum

      qiqbicnum10e0 = tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tF
     #jet40*tcoff40
      RETURN
      END

      double precision function qiqbicnum10e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t7 = z**2
      taux19 = -8*z*(1+z)*t3+4.D0/3.D0*(z-1)*(4*t7+7*z+4)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux20 = -96*z*(1+z)*t5*t6+16*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff37 = taux19
      tcoff38 = taux19
      tcoff39 = taux20
      tcoff40 = taux20

CCC Jet functions

      tFjet37 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,25.D0,z,1.D0,Tcoef
     #f*tcoff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,25.D0,z,1.D0,Tco
     #eff*tcoff38)
      tFjet39 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,26.D0,z,1.D0,Tcoef
     #f*tcoff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,26.D0,z,1.D0,Tco
     #eff*tcoff40)

CCC Sum

      qiqbicnum10e1 = tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tF
     #jet40*tcoff40
      RETURN
      END

      double precision function qiqbicnum10e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux19 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux20 = 64*z*(1+z)*t5*t6-32.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff37 = taux19
      tcoff38 = taux19
      tcoff39 = taux20
      tcoff40 = taux20

CCC Jet functions

      tFjet37 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,25.D0,z,1.D0,Tcoef
     #f*tcoff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,25.D0,z,1.D0,Tco
     #eff*tcoff38)
      tFjet39 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,26.D0,z,1.D0,Tcoef
     #f*tcoff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,26.D0,z,1.D0,Tco
     #eff*tcoff40)

CCC Sum

      qiqbicnum10e2 = tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tF
     #jet40*tcoff40
      RETURN
      END

      double precision function qiqbicnum10e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux19 = 0
      taux20 = 0

CCC Coefficients

      tcoff37 = taux19
      tcoff38 = taux19
      tcoff39 = taux20
      tcoff40 = taux20

CCC Jet functions

      tFjet37 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,25.D0,z,1.D0,Tcoef
     #f*tcoff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,25.D0,z,1.D0,Tco
     #eff*tcoff38)
      tFjet39 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,26.D0,z,1.D0,Tcoef
     #f*tcoff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,26.D0,z,1.D0,Tco
     #eff*tcoff40)

CCC Sum

      qiqbicnum10e3 = tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tF
     #jet40*tcoff40
      RETURN
      END

      double precision function qiqbicnum11e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux21 = 4*z*(1+z)*t5*t7-4.D0/3.D0*z*(3+15*z+8*t12)*t5*t6+8.D0/9.D
     #0*(z-1)*(28*t12+z+10)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux22 = -16*z*(1+z)*t5*t7+16.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-32.D0
     #/9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff41 = taux21
      tcoff42 = taux21
      tcoff43 = taux22
      tcoff44 = taux22

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,19.D0,1.D0,z,Tco
     #eff*tcoff41)
      tFjet42 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,19.D0,1.D0,z,Tcoef
     #f*tcoff42)
      tFjet43 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,20.D0,1.D0,z,Tco
     #eff*tcoff43)
      tFjet44 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,20.D0,1.D0,z,Tcoef
     #f*tcoff44)

CCC Sum

      qiqbicnum11e0 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tF
     #jet44*tcoff44
      RETURN
      END

      double precision function qiqbicnum11e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t4 = t3**2
      t8 = z**2
      taux21 = 4*z*(1+z)*t4-4.D0/3.D0*z*(3+15*z+8*t8)*t3+8.D0/9.D0*(z-1)
     #*(28*t8+z+10)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux22 = 48*z*(1+z)*t5*t7-16*z*(3+15*z+8*t12)*t5*t6+32.D0/3.D0*(z-
     #1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff41 = taux21
      tcoff42 = taux21
      tcoff43 = taux22
      tcoff44 = taux22

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,19.D0,1.D0,z,Tco
     #eff*tcoff41)
      tFjet42 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,19.D0,1.D0,z,Tcoef
     #f*tcoff42)
      tFjet43 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,20.D0,1.D0,z,Tco
     #eff*tcoff43)
      tFjet44 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,20.D0,1.D0,z,Tcoef
     #f*tcoff44)

CCC Sum

      qiqbicnum11e1 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tF
     #jet44*tcoff44
      RETURN
      END

      double precision function qiqbicnum11e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux21 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux22 = -32*z*(1+z)*t5*t7+32.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-64.D0
     #/9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff41 = taux21
      tcoff42 = taux21
      tcoff43 = taux22
      tcoff44 = taux22

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,19.D0,1.D0,z,Tco
     #eff*tcoff41)
      tFjet42 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,19.D0,1.D0,z,Tcoef
     #f*tcoff42)
      tFjet43 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,20.D0,1.D0,z,Tco
     #eff*tcoff43)
      tFjet44 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,20.D0,1.D0,z,Tcoef
     #f*tcoff44)

CCC Sum

      qiqbicnum11e2 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tF
     #jet44*tcoff44
      RETURN
      END

      double precision function qiqbicnum11e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux21 = 0
      taux22 = 0

CCC Coefficients

      tcoff41 = taux21
      tcoff42 = taux21
      tcoff43 = taux22
      tcoff44 = taux22

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,19.D0,1.D0,z,Tco
     #eff*tcoff41)
      tFjet42 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,19.D0,1.D0,z,Tcoef
     #f*tcoff42)
      tFjet43 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,20.D0,1.D0,z,Tco
     #eff*tcoff43)
      tFjet44 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,20.D0,1.D0,z,Tcoef
     #f*tcoff44)

CCC Sum

      qiqbicnum11e3 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tF
     #jet44*tcoff44
      RETURN
      END

      double precision function qiqbicnum12e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux23 = 4*z*(1+z)*t5*t7-4.D0/3.D0*z*(3+15*z+8*t12)*t5*t6+8.D0/9.D
     #0*(z-1)*(28*t12+z+10)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux24 = -16*z*(1+z)*t5*t7+16.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-32.D0
     #/9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff45 = taux23
      tcoff46 = taux23
      tcoff47 = taux24
      tcoff48 = taux24

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,25.D0,z,1.D0,Tcoef
     #f*tcoff45)
      tFjet46 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,25.D0,z,1.D0,Tco
     #eff*tcoff46)
      tFjet47 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,26.D0,z,1.D0,Tcoef
     #f*tcoff47)
      tFjet48 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,26.D0,z,1.D0,Tco
     #eff*tcoff48)

CCC Sum

      qiqbicnum12e0 = tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tF
     #jet48*tcoff48
      RETURN
      END

      double precision function qiqbicnum12e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t4 = t3**2
      t8 = z**2
      taux23 = 4*z*(1+z)*t4-4.D0/3.D0*z*(3+15*z+8*t8)*t3+8.D0/9.D0*(z-1)
     #*(28*t8+z+10)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux24 = 48*z*(1+z)*t5*t7-16*z*(3+15*z+8*t12)*t5*t6+32.D0/3.D0*(z-
     #1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff45 = taux23
      tcoff46 = taux23
      tcoff47 = taux24
      tcoff48 = taux24

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,25.D0,z,1.D0,Tcoef
     #f*tcoff45)
      tFjet46 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,25.D0,z,1.D0,Tco
     #eff*tcoff46)
      tFjet47 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,26.D0,z,1.D0,Tcoef
     #f*tcoff47)
      tFjet48 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,26.D0,z,1.D0,Tco
     #eff*tcoff48)

CCC Sum

      qiqbicnum12e1 = tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tF
     #jet48*tcoff48
      RETURN
      END

      double precision function qiqbicnum12e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux23 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux24 = -32*z*(1+z)*t5*t7+32.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-64.D0
     #/9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff45 = taux23
      tcoff46 = taux23
      tcoff47 = taux24
      tcoff48 = taux24

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,25.D0,z,1.D0,Tcoef
     #f*tcoff45)
      tFjet46 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,25.D0,z,1.D0,Tco
     #eff*tcoff46)
      tFjet47 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,26.D0,z,1.D0,Tcoef
     #f*tcoff47)
      tFjet48 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,26.D0,z,1.D0,Tco
     #eff*tcoff48)

CCC Sum

      qiqbicnum12e2 = tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tF
     #jet48*tcoff48
      RETURN
      END

      double precision function qiqbicnum12e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux23 = 0
      taux24 = 0

CCC Coefficients

      tcoff45 = taux23
      tcoff46 = taux23
      tcoff47 = taux24
      tcoff48 = taux24

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,25.D0,z,1.D0,Tcoef
     #f*tcoff45)
      tFjet46 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,25.D0,z,1.D0,Tco
     #eff*tcoff46)
      tFjet47 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,26.D0,z,1.D0,Tcoef
     #f*tcoff47)
      tFjet48 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,26.D0,z,1.D0,Tco
     #eff*tcoff48)

CCC Sum

      qiqbicnum12e3 = tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tF
     #jet48*tcoff48
      RETURN
      END

      double precision function qiqbicnum13e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t11 = z**2
      taux25 = -8*z*(1+z)*t5*t6+4.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux26 = 32*z*(1+z)*t5*t6-16.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff49 = taux25
      tcoff50 = taux25
      tcoff51 = taux26
      tcoff52 = taux26

CCC Jet functions

      tFjet49 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,23.D0,1.D0,z,Tcoef
     #f*tcoff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,23.D0,1.D0,z,Tco
     #eff*tcoff50)
      tFjet51 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,24.D0,1.D0,z,Tcoef
     #f*tcoff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,24.D0,1.D0,z,Tco
     #eff*tcoff52)

CCC Sum

      qiqbicnum13e0 = tFjet49*tcoff49+tFjet50*tcoff50+tFjet51*tcoff51+tF
     #jet52*tcoff52
      RETURN
      END

      double precision function qiqbicnum13e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t7 = z**2
      taux25 = -8*z*(1+z)*t3+4.D0/3.D0*(z-1)*(4*t7+7*z+4)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux26 = -96*z*(1+z)*t5*t6+16*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff49 = taux25
      tcoff50 = taux25
      tcoff51 = taux26
      tcoff52 = taux26

CCC Jet functions

      tFjet49 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,23.D0,1.D0,z,Tcoef
     #f*tcoff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,23.D0,1.D0,z,Tco
     #eff*tcoff50)
      tFjet51 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,24.D0,1.D0,z,Tcoef
     #f*tcoff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,24.D0,1.D0,z,Tco
     #eff*tcoff52)

CCC Sum

      qiqbicnum13e1 = tFjet49*tcoff49+tFjet50*tcoff50+tFjet51*tcoff51+tF
     #jet52*tcoff52
      RETURN
      END

      double precision function qiqbicnum13e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux25 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux26 = 64*z*(1+z)*t5*t6-32.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff49 = taux25
      tcoff50 = taux25
      tcoff51 = taux26
      tcoff52 = taux26

CCC Jet functions

      tFjet49 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,23.D0,1.D0,z,Tcoef
     #f*tcoff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,23.D0,1.D0,z,Tco
     #eff*tcoff50)
      tFjet51 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,24.D0,1.D0,z,Tcoef
     #f*tcoff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,24.D0,1.D0,z,Tco
     #eff*tcoff52)

CCC Sum

      qiqbicnum13e2 = tFjet49*tcoff49+tFjet50*tcoff50+tFjet51*tcoff51+tF
     #jet52*tcoff52
      RETURN
      END

      double precision function qiqbicnum13e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux25 = 0
      taux26 = 0

CCC Coefficients

      tcoff49 = taux25
      tcoff50 = taux25
      tcoff51 = taux26
      tcoff52 = taux26

CCC Jet functions

      tFjet49 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,23.D0,1.D0,z,Tcoef
     #f*tcoff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,23.D0,1.D0,z,Tco
     #eff*tcoff50)
      tFjet51 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,24.D0,1.D0,z,Tcoef
     #f*tcoff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,24.D0,1.D0,z,Tco
     #eff*tcoff52)

CCC Sum

      qiqbicnum13e3 = tFjet49*tcoff49+tFjet50*tcoff50+tFjet51*tcoff51+tF
     #jet52*tcoff52
      RETURN
      END

      double precision function qiqbicnum14e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t11 = z**2
      taux27 = -8*z*(1+z)*t5*t6+4.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux28 = 32*z*(1+z)*t5*t6-16.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff53 = taux27
      tcoff54 = taux27
      tcoff55 = taux28
      tcoff56 = taux28

CCC Jet functions

      tFjet53 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,21.D0,z,1.D0,Tcoef
     #f*tcoff53)
      tFjet54 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,21.D0,z,1.D0,Tco
     #eff*tcoff54)
      tFjet55 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,22.D0,z,1.D0,Tcoef
     #f*tcoff55)
      tFjet56 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,22.D0,z,1.D0,Tco
     #eff*tcoff56)

CCC Sum

      qiqbicnum14e0 = tFjet53*tcoff53+tFjet54*tcoff54+tFjet55*tcoff55+tF
     #jet56*tcoff56
      RETURN
      END

      double precision function qiqbicnum14e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t7 = z**2
      taux27 = -8*z*(1+z)*t3+4.D0/3.D0*(z-1)*(4*t7+7*z+4)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux28 = -96*z*(1+z)*t5*t6+16*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff53 = taux27
      tcoff54 = taux27
      tcoff55 = taux28
      tcoff56 = taux28

CCC Jet functions

      tFjet53 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,21.D0,z,1.D0,Tcoef
     #f*tcoff53)
      tFjet54 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,21.D0,z,1.D0,Tco
     #eff*tcoff54)
      tFjet55 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,22.D0,z,1.D0,Tcoef
     #f*tcoff55)
      tFjet56 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,22.D0,z,1.D0,Tco
     #eff*tcoff56)

CCC Sum

      qiqbicnum14e1 = tFjet53*tcoff53+tFjet54*tcoff54+tFjet55*tcoff55+tF
     #jet56*tcoff56
      RETURN
      END

      double precision function qiqbicnum14e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux27 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux28 = 64*z*(1+z)*t5*t6-32.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff53 = taux27
      tcoff54 = taux27
      tcoff55 = taux28
      tcoff56 = taux28

CCC Jet functions

      tFjet53 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,21.D0,z,1.D0,Tcoef
     #f*tcoff53)
      tFjet54 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,21.D0,z,1.D0,Tco
     #eff*tcoff54)
      tFjet55 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,22.D0,z,1.D0,Tcoef
     #f*tcoff55)
      tFjet56 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,22.D0,z,1.D0,Tco
     #eff*tcoff56)

CCC Sum

      qiqbicnum14e2 = tFjet53*tcoff53+tFjet54*tcoff54+tFjet55*tcoff55+tF
     #jet56*tcoff56
      RETURN
      END

      double precision function qiqbicnum14e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux27 = 0
      taux28 = 0

CCC Coefficients

      tcoff53 = taux27
      tcoff54 = taux27
      tcoff55 = taux28
      tcoff56 = taux28

CCC Jet functions

      tFjet53 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,21.D0,z,1.D0,Tcoef
     #f*tcoff53)
      tFjet54 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,21.D0,z,1.D0,Tco
     #eff*tcoff54)
      tFjet55 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,22.D0,z,1.D0,Tcoef
     #f*tcoff55)
      tFjet56 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,22.D0,z,1.D0,Tco
     #eff*tcoff56)

CCC Sum

      qiqbicnum14e3 = tFjet53*tcoff53+tFjet54*tcoff54+tFjet55*tcoff55+tF
     #jet56*tcoff56
      RETURN
      END

      double precision function qiqbicnum15e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux29 = 4*z*(1+z)*t5*t7-4.D0/3.D0*z*(3+15*z+8*t12)*t5*t6+8.D0/9.D
     #0*(z-1)*(28*t12+z+10)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux30 = -16*z*(1+z)*t5*t7+16.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-32.D0
     #/9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff57 = taux29
      tcoff58 = taux29
      tcoff59 = taux30
      tcoff60 = taux30

CCC Jet functions

      tFjet57 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,23.D0,1.D0,z,Tcoef
     #f*tcoff57)
      tFjet58 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,23.D0,1.D0,z,Tco
     #eff*tcoff58)
      tFjet59 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,24.D0,1.D0,z,Tcoef
     #f*tcoff59)
      tFjet60 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,24.D0,1.D0,z,Tco
     #eff*tcoff60)

CCC Sum

      qiqbicnum15e0 = tFjet57*tcoff57+tFjet58*tcoff58+tFjet59*tcoff59+tF
     #jet60*tcoff60
      RETURN
      END

      double precision function qiqbicnum15e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t4 = t3**2
      t8 = z**2
      taux29 = 4*z*(1+z)*t4-4.D0/3.D0*z*(3+15*z+8*t8)*t3+8.D0/9.D0*(z-1)
     #*(28*t8+z+10)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux30 = 48*z*(1+z)*t5*t7-16*z*(3+15*z+8*t12)*t5*t6+32.D0/3.D0*(z-
     #1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff57 = taux29
      tcoff58 = taux29
      tcoff59 = taux30
      tcoff60 = taux30

CCC Jet functions

      tFjet57 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,23.D0,1.D0,z,Tcoef
     #f*tcoff57)
      tFjet58 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,23.D0,1.D0,z,Tco
     #eff*tcoff58)
      tFjet59 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,24.D0,1.D0,z,Tcoef
     #f*tcoff59)
      tFjet60 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,24.D0,1.D0,z,Tco
     #eff*tcoff60)

CCC Sum

      qiqbicnum15e1 = tFjet57*tcoff57+tFjet58*tcoff58+tFjet59*tcoff59+tF
     #jet60*tcoff60
      RETURN
      END

      double precision function qiqbicnum15e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux29 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux30 = -32*z*(1+z)*t5*t7+32.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-64.D0
     #/9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff57 = taux29
      tcoff58 = taux29
      tcoff59 = taux30
      tcoff60 = taux30

CCC Jet functions

      tFjet57 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,23.D0,1.D0,z,Tcoef
     #f*tcoff57)
      tFjet58 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,23.D0,1.D0,z,Tco
     #eff*tcoff58)
      tFjet59 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,24.D0,1.D0,z,Tcoef
     #f*tcoff59)
      tFjet60 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,24.D0,1.D0,z,Tco
     #eff*tcoff60)

CCC Sum

      qiqbicnum15e2 = tFjet57*tcoff57+tFjet58*tcoff58+tFjet59*tcoff59+tF
     #jet60*tcoff60
      RETURN
      END

      double precision function qiqbicnum15e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux29 = 0
      taux30 = 0

CCC Coefficients

      tcoff57 = taux29
      tcoff58 = taux29
      tcoff59 = taux30
      tcoff60 = taux30

CCC Jet functions

      tFjet57 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,23.D0,1.D0,z,Tcoef
     #f*tcoff57)
      tFjet58 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,23.D0,1.D0,z,Tco
     #eff*tcoff58)
      tFjet59 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,24.D0,1.D0,z,Tcoef
     #f*tcoff59)
      tFjet60 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,24.D0,1.D0,z,Tco
     #eff*tcoff60)

CCC Sum

      qiqbicnum15e3 = tFjet57*tcoff57+tFjet58*tcoff58+tFjet59*tcoff59+tF
     #jet60*tcoff60
      RETURN
      END

      double precision function qiqbicnum16e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux31 = 4*z*(1+z)*t5*t7-4.D0/3.D0*z*(3+15*z+8*t12)*t5*t6+8.D0/9.D
     #0*(z-1)*(28*t12+z+10)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux32 = -16*z*(1+z)*t5*t7+16.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-32.D0
     #/9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff61 = taux31
      tcoff62 = taux31
      tcoff63 = taux32
      tcoff64 = taux32

CCC Jet functions

      tFjet61 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,21.D0,z,1.D0,Tcoef
     #f*tcoff61)
      tFjet62 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,21.D0,z,1.D0,Tco
     #eff*tcoff62)
      tFjet63 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,22.D0,z,1.D0,Tcoef
     #f*tcoff63)
      tFjet64 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,22.D0,z,1.D0,Tco
     #eff*tcoff64)

CCC Sum

      qiqbicnum16e0 = tFjet61*tcoff61+tFjet62*tcoff62+tFjet63*tcoff63+tF
     #jet64*tcoff64
      RETURN
      END

      double precision function qiqbicnum16e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t4 = t3**2
      t8 = z**2
      taux31 = 4*z*(1+z)*t4-4.D0/3.D0*z*(3+15*z+8*t8)*t3+8.D0/9.D0*(z-1)
     #*(28*t8+z+10)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux32 = 48*z*(1+z)*t5*t7-16*z*(3+15*z+8*t12)*t5*t6+32.D0/3.D0*(z-
     #1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff61 = taux31
      tcoff62 = taux31
      tcoff63 = taux32
      tcoff64 = taux32

CCC Jet functions

      tFjet61 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,21.D0,z,1.D0,Tcoef
     #f*tcoff61)
      tFjet62 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,21.D0,z,1.D0,Tco
     #eff*tcoff62)
      tFjet63 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,22.D0,z,1.D0,Tcoef
     #f*tcoff63)
      tFjet64 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,22.D0,z,1.D0,Tco
     #eff*tcoff64)

CCC Sum

      qiqbicnum16e1 = tFjet61*tcoff61+tFjet62*tcoff62+tFjet63*tcoff63+tF
     #jet64*tcoff64
      RETURN
      END

      double precision function qiqbicnum16e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux31 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux32 = -32*z*(1+z)*t5*t7+32.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-64.D0
     #/9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff61 = taux31
      tcoff62 = taux31
      tcoff63 = taux32
      tcoff64 = taux32

CCC Jet functions

      tFjet61 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,21.D0,z,1.D0,Tcoef
     #f*tcoff61)
      tFjet62 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,21.D0,z,1.D0,Tco
     #eff*tcoff62)
      tFjet63 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,22.D0,z,1.D0,Tcoef
     #f*tcoff63)
      tFjet64 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,22.D0,z,1.D0,Tco
     #eff*tcoff64)

CCC Sum

      qiqbicnum16e2 = tFjet61*tcoff61+tFjet62*tcoff62+tFjet63*tcoff63+tF
     #jet64*tcoff64
      RETURN
      END

      double precision function qiqbicnum16e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux31 = 0
      taux32 = 0

CCC Coefficients

      tcoff61 = taux31
      tcoff62 = taux31
      tcoff63 = taux32
      tcoff64 = taux32

CCC Jet functions

      tFjet61 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,21.D0,z,1.D0,Tcoef
     #f*tcoff61)
      tFjet62 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,21.D0,z,1.D0,Tco
     #eff*tcoff62)
      tFjet63 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,22.D0,z,1.D0,Tcoef
     #f*tcoff63)
      tFjet64 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,22.D0,z,1.D0,Tco
     #eff*tcoff64)

CCC Sum

      qiqbicnum16e3 = tFjet61*tcoff61+tFjet62*tcoff62+tFjet63*tcoff63+tF
     #jet64*tcoff64
      RETURN
      END

      double precision function qiqbicnum17e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp-1
      t10 = zp*l1
      t12 = -64*l1-64*zp+64*t10
      t17 = -64*l1+64*t10+64-128*zp
      t18 = -s1l2-s2l2-zp
      t23 = -s1l1-s2l1-zp
      t28 = -64*t6
      taux33 = 2.D0/3.D0*(1+t2)/yp*((128*t6*l1*s1l2+t12*s2l2+t17*t18)*s1
     #l1+(t12*s2l1+t17*t23)*s1l2+t28*t18*s2l1+t28*t23*s2l2+128*t23*t18)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = -l1-zp+t6
      t10 = -l1+t6+1-2*zp
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = -zp+1
      taux34 = 2.D0/3.D0*(1+t2)/yp*((256*t7*s2l2+256*t10*t11)*s1l1+(-256
     #*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2)
      t2 = (1-yp)**2
      t6 = zp-1
      t10 = zp*l1
      t12 = -64*l1-64*zp+64*t10
      t17 = -64*l1+64*t10+64-128*zp
      t18 = -s1l2-s2l2-zp
      t23 = -s1l1-s2l1-zp
      t28 = -64*t6
      taux35 = 2.D0/3.D0*(1+t2)/yp*((128*t6*l1*s1l2+t12*s2l2+t17*t18)*s1
     #l1+(t12*s2l1+t17*t23)*s1l2+t28*t18*s2l1+t28*t23*s2l2+128*t23*t18)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = -l1-zp+t6
      t10 = -l1+t6+1-2*zp
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = -zp+1
      taux36 = 2.D0/3.D0*(1+t2)/yp*((256*t7*s2l2+256*t10*t11)*s1l1+(-256
     #*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2)

CCC Coefficients

      tcoff65 = taux33
      tcoff66 = taux33
      tcoff67 = taux34
      tcoff68 = taux34
      tcoff69 = taux35
      tcoff70 = taux35
      tcoff71 = -taux36
      tcoff72 = -taux36

CCC Jet functions

      tFjet65 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*tco
     #ff65)
      tFjet66 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*t
     #coff66)
      tFjet67 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*tco
     #ff67)
      tFjet68 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*t
     #coff68)
      tFjet69 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*t
     #coff69)
      tFjet70 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*tco
     #ff70)
      tFjet71 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*t
     #coff71)
      tFjet72 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*tco
     #ff72)

CCC Sum

      qiqbicnum17e0 = tFjet65*tcoff65+tFjet66*tcoff66+tFjet67*tcoff67+tF
     #jet68*tcoff68+tFjet69*tcoff69+tFjet70*tcoff70+tFjet71*tcoff71+tFje
     #t72*tcoff72
      RETURN
      END

      double precision function qiqbicnum17e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp*l1
      t8 = 64-64*l1+64*t6
      t25 = zp**2
      t28 = l1**2
      taux33 = 2.D0/3.D0*(1+t2)/yp*(t8*s2l2*s1l1+t8*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t8*(-s1l2-s2l2-zp))*s2l1+t8*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(2*t25*l1-t25*t28-2*t25+2*zp*t28-1+2*zp-2*t6-t28))
      t2 = (1-yp)**2
      t6 = 3*zp
      t7 = 3*l1
      t8 = zp*l1
      t9 = 3*t8
      t10 = 1+t6+t7-t9
      t13 = t7-t9-3+6*zp
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = t6+t8-l1-3
      taux34 = 2.D0/3.D0*(1+t2)/yp*((256*t10*s2l2+256*t13*t14)*s1l1+(-25
     #6*t10*s2l1-256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)
      t2 = (1-yp)**2
      t6 = zp*l1
      t8 = 64-64*l1+64*t6
      t25 = zp**2
      t28 = l1**2
      taux35 = 2.D0/3.D0*(1+t2)/yp*(t8*s2l2*s1l1+t8*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t8*(-s1l2-s2l2-zp))*s2l1+t8*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(2*t25*l1-t25*t28-2*t25+2*zp*t28-1+2*zp-2*t6-t28))
      t2 = (1-yp)**2
      t6 = 3*zp
      t7 = 3*l1
      t8 = zp*l1
      t9 = 3*t8
      t10 = 1+t6+t7-t9
      t13 = t7-t9-3+6*zp
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = t6+t8-l1-3
      taux36 = 2.D0/3.D0*(1+t2)/yp*((256*t10*s2l2+256*t13*t14)*s1l1+(-25
     #6*t10*s2l1-256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)

CCC Coefficients

      tcoff65 = taux33
      tcoff66 = taux33
      tcoff67 = taux34
      tcoff68 = taux34
      tcoff69 = taux35
      tcoff70 = taux35
      tcoff71 = -taux36
      tcoff72 = -taux36

CCC Jet functions

      tFjet65 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*tco
     #ff65)
      tFjet66 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*t
     #coff66)
      tFjet67 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*tco
     #ff67)
      tFjet68 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*t
     #coff68)
      tFjet69 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*t
     #coff69)
      tFjet70 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*tco
     #ff70)
      tFjet71 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*t
     #coff71)
      tFjet72 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*tco
     #ff72)

CCC Sum

      qiqbicnum17e1 = tFjet65*tcoff65+tFjet66*tcoff66+tFjet67*tcoff67+tF
     #jet68*tcoff68+tFjet69*tcoff69+tFjet70*tcoff70+tFjet71*tcoff71+tFje
     #t72*tcoff72
      RETURN
      END

      double precision function qiqbicnum17e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp**2
      t7 = l1**2
      taux33 = -128.D0/3.D0*(1+t2)/yp*zp*(t6*t7+2*zp*l1-2*zp*t7+1+t7-2*l
     #1)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 2*zp
      t8 = t6-3-l1-t7
      t11 = -l1+t6+1-t7
      t12 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t26 = -3*t6+1+3*l1-t7
      taux34 = 2.D0/3.D0*(1+t2)/yp*((256*t8*s2l2+512*t11*t12)*s1l1+(-256
     #*t8*s2l1-512*t11*t19)*s1l2+256*t26*t12*s2l1-256*t26*t19*s2l2)
      t2 = (1-yp)**2
      t6 = zp**2
      t7 = l1**2
      taux35 = -128.D0/3.D0*(1+t2)/yp*zp*(t6*t7+2*zp*l1-2*zp*t7+1+t7-2*l
     #1)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 2*zp
      t8 = t6-3-l1-t7
      t11 = -l1+t6+1-t7
      t12 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t26 = -3*t6+1+3*l1-t7
      taux36 = 2.D0/3.D0*(1+t2)/yp*((256*t8*s2l2+512*t11*t12)*s1l1+(-256
     #*t8*s2l1-512*t11*t19)*s1l2+256*t26*t12*s2l1-256*t26*t19*s2l2)

CCC Coefficients

      tcoff65 = taux33
      tcoff66 = taux33
      tcoff67 = taux34
      tcoff68 = taux34
      tcoff69 = taux35
      tcoff70 = taux35
      tcoff71 = -taux36
      tcoff72 = -taux36

CCC Jet functions

      tFjet65 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*tco
     #ff65)
      tFjet66 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*t
     #coff66)
      tFjet67 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*tco
     #ff67)
      tFjet68 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*t
     #coff68)
      tFjet69 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*t
     #coff69)
      tFjet70 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*tco
     #ff70)
      tFjet71 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*t
     #coff71)
      tFjet72 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*tco
     #ff72)

CCC Sum

      qiqbicnum17e2 = tFjet65*tcoff65+tFjet66*tcoff66+tFjet67*tcoff67+tF
     #jet68*tcoff68+tFjet69*tcoff69+tFjet70*tcoff70+tFjet71*tcoff71+tFje
     #t72*tcoff72
      RETURN
      END

      double precision function qiqbicnum17e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux33 = 0
      t2 = (1-yp)**2
      t7 = zp*l1+1-l1
      taux34 = 2.D0/3.D0*(1+t2)/yp*(512*t7*s2l2*s1l1-512*t7*s2l1*s1l2+51
     #2*t7*(-s1l2-s2l2-zp)*s2l1-512*t7*(-s1l1-s2l1-zp)*s2l2)
      taux35 = 0
      t2 = (1-yp)**2
      t7 = zp*l1+1-l1
      taux36 = 2.D0/3.D0*(1+t2)/yp*(512*t7*s2l2*s1l1-512*t7*s2l1*s1l2+51
     #2*t7*(-s1l2-s2l2-zp)*s2l1-512*t7*(-s1l1-s2l1-zp)*s2l2)

CCC Coefficients

      tcoff65 = taux33
      tcoff66 = taux33
      tcoff67 = taux34
      tcoff68 = taux34
      tcoff69 = taux35
      tcoff70 = taux35
      tcoff71 = -taux36
      tcoff72 = -taux36

CCC Jet functions

      tFjet65 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*tco
     #ff65)
      tFjet66 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*t
     #coff66)
      tFjet67 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*tco
     #ff67)
      tFjet68 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*t
     #coff68)
      tFjet69 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*t
     #coff69)
      tFjet70 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*tco
     #ff70)
      tFjet71 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*t
     #coff71)
      tFjet72 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*tco
     #ff72)

CCC Sum

      qiqbicnum17e3 = tFjet65*tcoff65+tFjet66*tcoff66+tFjet67*tcoff67+tF
     #jet68*tcoff68+tFjet69*tcoff69+tFjet70*tcoff70+tFjet71*tcoff71+tFje
     #t72*tcoff72
      RETURN
      END

      double precision function qiqbicnum18e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp-1
      t10 = zp*l1
      t12 = -64*l1-64*zp+64*t10
      t17 = -64*l1+64*t10+64-128*zp
      t18 = -s1l2-s2l2-zp
      t23 = -s1l1-s2l1-zp
      t28 = -64*t6
      taux37 = 2.D0/3.D0*(1+t2)/yp*((128*t6*l1*s1l2+t12*s2l2+t17*t18)*s1
     #l1+(t12*s2l1+t17*t23)*s1l2+t28*t18*s2l1+t28*t23*s2l2+128*t23*t18)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = -l1-zp+t6
      t10 = -l1+t6+1-2*zp
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = -zp+1
      taux38 = 2.D0/3.D0*(1+t2)/yp*((256*t7*s2l2+256*t10*t11)*s1l1+(-256
     #*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2)

CCC Coefficients

      tcoff73 = taux37
      tcoff74 = taux37
      tcoff75 = taux38
      tcoff76 = taux38

CCC Jet functions

      tFjet73 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,15.D0,yp,1.D0,Tcoeff*tco
     #ff73)
      tFjet74 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,15.D0,yp,1.D0,Tcoeff*t
     #coff74)
      tFjet75 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,16.D0,yp,1.D0,Tcoeff*tco
     #ff75)
      tFjet76 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,16.D0,yp,1.D0,Tcoeff*t
     #coff76)

CCC Sum

      qiqbicnum18e0 = tFjet73*tcoff73+tFjet74*tcoff74+tFjet75*tcoff75+tF
     #jet76*tcoff76
      RETURN
      END

      double precision function qiqbicnum18e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp*l1
      t8 = 64-64*l1+64*t6
      t25 = zp**2
      t28 = l1**2
      taux37 = 2.D0/3.D0*(1+t2)/yp*(t8*s2l2*s1l1+t8*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t8*(-s1l2-s2l2-zp))*s2l1+t8*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(2*t25*l1-t25*t28-2*t25+2*zp*t28-1+2*zp-2*t6-t28))
      t2 = (1-yp)**2
      t6 = 3*zp
      t7 = 3*l1
      t8 = zp*l1
      t9 = 3*t8
      t10 = 1+t6+t7-t9
      t13 = t7-t9-3+6*zp
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = t6+t8-l1-3
      taux38 = 2.D0/3.D0*(1+t2)/yp*((256*t10*s2l2+256*t13*t14)*s1l1+(-25
     #6*t10*s2l1-256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)

CCC Coefficients

      tcoff73 = taux37
      tcoff74 = taux37
      tcoff75 = taux38
      tcoff76 = taux38

CCC Jet functions

      tFjet73 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,15.D0,yp,1.D0,Tcoeff*tco
     #ff73)
      tFjet74 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,15.D0,yp,1.D0,Tcoeff*t
     #coff74)
      tFjet75 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,16.D0,yp,1.D0,Tcoeff*tco
     #ff75)
      tFjet76 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,16.D0,yp,1.D0,Tcoeff*t
     #coff76)

CCC Sum

      qiqbicnum18e1 = tFjet73*tcoff73+tFjet74*tcoff74+tFjet75*tcoff75+tF
     #jet76*tcoff76
      RETURN
      END

      double precision function qiqbicnum18e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp**2
      t7 = l1**2
      taux37 = -128.D0/3.D0*(1+t2)/yp*zp*(t6*t7+2*zp*l1+1-2*zp*t7+t7-2*l
     #1)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 2*zp
      t8 = t6-3-l1-t7
      t11 = -l1+t6+1-t7
      t12 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t26 = -3*t6+1+3*l1-t7
      taux38 = 2.D0/3.D0*(1+t2)/yp*((256*t8*s2l2+512*t11*t12)*s1l1+(-256
     #*t8*s2l1-512*t11*t19)*s1l2+256*t26*t12*s2l1-256*t26*t19*s2l2)

CCC Coefficients

      tcoff73 = taux37
      tcoff74 = taux37
      tcoff75 = taux38
      tcoff76 = taux38

CCC Jet functions

      tFjet73 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,15.D0,yp,1.D0,Tcoeff*tco
     #ff73)
      tFjet74 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,15.D0,yp,1.D0,Tcoeff*t
     #coff74)
      tFjet75 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,16.D0,yp,1.D0,Tcoeff*tco
     #ff75)
      tFjet76 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,16.D0,yp,1.D0,Tcoeff*t
     #coff76)

CCC Sum

      qiqbicnum18e2 = tFjet73*tcoff73+tFjet74*tcoff74+tFjet75*tcoff75+tF
     #jet76*tcoff76
      RETURN
      END

      double precision function qiqbicnum18e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux37 = 0
      t2 = (1-yp)**2
      t7 = zp*l1+1-l1
      taux38 = 2.D0/3.D0*(1+t2)/yp*(512*t7*s2l2*s1l1-512*t7*s2l1*s1l2+51
     #2*t7*(-s1l2-s2l2-zp)*s2l1-512*t7*(-s1l1-s2l1-zp)*s2l2)

CCC Coefficients

      tcoff73 = taux37
      tcoff74 = taux37
      tcoff75 = taux38
      tcoff76 = taux38

CCC Jet functions

      tFjet73 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,15.D0,yp,1.D0,Tcoeff*tco
     #ff73)
      tFjet74 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,15.D0,yp,1.D0,Tcoeff*t
     #coff74)
      tFjet75 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,16.D0,yp,1.D0,Tcoeff*tco
     #ff75)
      tFjet76 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,16.D0,yp,1.D0,Tcoeff*t
     #coff76)

CCC Sum

      qiqbicnum18e3 = tFjet73*tcoff73+tFjet74*tcoff74+tFjet75*tcoff75+tF
     #jet76*tcoff76
      RETURN
      END

      double precision function qiqbicnum19e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp-1
      t10 = zp*l1
      t12 = -64*l1-64*zp+64*t10
      t17 = -64*l1+64*t10+64-128*zp
      t18 = -s1l2-s2l2-zp
      t23 = -s1l1-s2l1-zp
      t28 = -64*t6
      taux39 = 2.D0/3.D0*(1+t2)/yp*((128*t6*l1*s1l2+t12*s2l2+t17*t18)*s1
     #l1+(t12*s2l1+t17*t23)*s1l2+t28*t18*s2l1+t28*t23*s2l2+128*t23*t18)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = -l1-zp+t6
      t10 = -l1+t6+1-2*zp
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = -zp+1
      taux40 = 2.D0/3.D0*(1+t2)/yp*((256*t7*s2l2+256*t10*t11)*s1l1+(-256
     #*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2)

CCC Coefficients

      tcoff77 = taux39
      tcoff78 = taux39
      tcoff79 = taux40
      tcoff80 = taux40

CCC Jet functions

      tFjet77 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,17.D0,1.D0,yp,Tcoeff*t
     #coff77)
      tFjet78 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,17.D0,1.D0,yp,Tcoeff*tco
     #ff78)
      tFjet79 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,18.D0,1.D0,yp,Tcoeff*t
     #coff79)
      tFjet80 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,18.D0,1.D0,yp,Tcoeff*tco
     #ff80)

CCC Sum

      qiqbicnum19e0 = tFjet77*tcoff77+tFjet78*tcoff78+tFjet79*tcoff79+tF
     #jet80*tcoff80
      RETURN
      END

      double precision function qiqbicnum19e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp*l1
      t8 = 64-64*l1+64*t6
      t25 = zp**2
      t28 = l1**2
      taux39 = 2.D0/3.D0*(1+t2)/yp*(t8*s2l2*s1l1+t8*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t8*(-s1l2-s2l2-zp))*s2l1+t8*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(2*t25*l1-t25*t28-2*t25+2*zp*t28-1+2*zp-2*t6-t28))
      t2 = (1-yp)**2
      t6 = 3*zp
      t7 = 3*l1
      t8 = zp*l1
      t9 = 3*t8
      t10 = 1+t6+t7-t9
      t13 = t7-t9-3+6*zp
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = t6+t8-l1-3
      taux40 = 2.D0/3.D0*(1+t2)/yp*((256*t10*s2l2+256*t13*t14)*s1l1+(-25
     #6*t10*s2l1-256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)

CCC Coefficients

      tcoff77 = taux39
      tcoff78 = taux39
      tcoff79 = taux40
      tcoff80 = taux40

CCC Jet functions

      tFjet77 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,17.D0,1.D0,yp,Tcoeff*t
     #coff77)
      tFjet78 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,17.D0,1.D0,yp,Tcoeff*tco
     #ff78)
      tFjet79 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,18.D0,1.D0,yp,Tcoeff*t
     #coff79)
      tFjet80 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,18.D0,1.D0,yp,Tcoeff*tco
     #ff80)

CCC Sum

      qiqbicnum19e1 = tFjet77*tcoff77+tFjet78*tcoff78+tFjet79*tcoff79+tF
     #jet80*tcoff80
      RETURN
      END

      double precision function qiqbicnum19e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp**2
      t7 = l1**2
      taux39 = -128.D0/3.D0*(1+t2)/yp*zp*(t6*t7+2*zp*l1+1-2*zp*t7+t7-2*l
     #1)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 2*zp
      t8 = t6-3-l1-t7
      t11 = -l1+t6+1-t7
      t12 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t26 = -3*t6+1+3*l1-t7
      taux40 = 2.D0/3.D0*(1+t2)/yp*((256*t8*s2l2+512*t11*t12)*s1l1+(-256
     #*t8*s2l1-512*t11*t19)*s1l2+256*t26*t12*s2l1-256*t26*t19*s2l2)

CCC Coefficients

      tcoff77 = taux39
      tcoff78 = taux39
      tcoff79 = taux40
      tcoff80 = taux40

CCC Jet functions

      tFjet77 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,17.D0,1.D0,yp,Tcoeff*t
     #coff77)
      tFjet78 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,17.D0,1.D0,yp,Tcoeff*tco
     #ff78)
      tFjet79 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,18.D0,1.D0,yp,Tcoeff*t
     #coff79)
      tFjet80 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,18.D0,1.D0,yp,Tcoeff*tco
     #ff80)

CCC Sum

      qiqbicnum19e2 = tFjet77*tcoff77+tFjet78*tcoff78+tFjet79*tcoff79+tF
     #jet80*tcoff80
      RETURN
      END

      double precision function qiqbicnum19e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux39 = 0
      t2 = (1-yp)**2
      t7 = zp*l1+1-l1
      taux40 = 2.D0/3.D0*(1+t2)/yp*(512*t7*s2l2*s1l1-512*t7*s2l1*s1l2+51
     #2*t7*(-s1l2-s2l2-zp)*s2l1-512*t7*(-s1l1-s2l1-zp)*s2l2)

CCC Coefficients

      tcoff77 = taux39
      tcoff78 = taux39
      tcoff79 = taux40
      tcoff80 = taux40

CCC Jet functions

      tFjet77 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,17.D0,1.D0,yp,Tcoeff*t
     #coff77)
      tFjet78 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,17.D0,1.D0,yp,Tcoeff*tco
     #ff78)
      tFjet79 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,18.D0,1.D0,yp,Tcoeff*t
     #coff79)
      tFjet80 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,18.D0,1.D0,yp,Tcoeff*tco
     #ff80)

CCC Sum

      qiqbicnum19e3 = tFjet77*tcoff77+tFjet78*tcoff78+tFjet79*tcoff79+tF
     #jet80*tcoff80
      RETURN
      END

      double precision function qiqbicnum20e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp-1
      t10 = zp*l1
      t12 = -64*l1-64*zp+64*t10
      t17 = -64*l1+64*t10+64-128*zp
      t18 = -s1l2-s2l2-zp
      t23 = -s1l1-s2l1-zp
      t28 = -64*t6
      taux41 = 2.D0/3.D0*(1+t2)/yp*((128*t6*l1*s1l2+t12*s2l2+t17*t18)*s1
     #l1+(t12*s2l1+t17*t23)*s1l2+t28*t18*s2l1+t28*t23*s2l2+128*t23*t18)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = -l1-zp+t6
      t10 = -l1+t6+1-2*zp
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = -zp+1
      taux42 = 2.D0/3.D0*(1+t2)/yp*((256*t7*s2l2+256*t10*t11)*s1l1+(-256
     #*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2)
      t2 = (1-yp)**2
      t6 = zp-1
      t10 = zp*l1
      t12 = -64*l1-64*zp+64*t10
      t17 = -64*l1+64*t10+64-128*zp
      t18 = -s1l2-s2l2-zp
      t23 = -s1l1-s2l1-zp
      t28 = -64*t6
      taux43 = 2.D0/3.D0*(1+t2)/yp*((128*t6*l1*s1l2+t12*s2l2+t17*t18)*s1
     #l1+(t12*s2l1+t17*t23)*s1l2+t28*t18*s2l1+t28*t23*s2l2+128*t23*t18)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = -l1-zp+t6
      t10 = -l1+t6+1-2*zp
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = -zp+1
      taux44 = 2.D0/3.D0*(1+t2)/yp*((256*t7*s2l2+256*t10*t11)*s1l1+(-256
     #*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2)

CCC Coefficients

      tcoff81 = taux41
      tcoff82 = taux41
      tcoff83 = taux42
      tcoff84 = taux42
      tcoff85 = taux43
      tcoff86 = taux43
      tcoff87 = -taux44
      tcoff88 = -taux44

CCC Jet functions

      tFjet81 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,19.D0,1.D0,yp,Tcoeff*tco
     #ff81)
      tFjet82 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,19.D0,1.D0,yp,Tcoeff*t
     #coff82)
      tFjet83 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,20.D0,1.D0,yp,Tcoeff*tco
     #ff83)
      tFjet84 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,20.D0,1.D0,yp,Tcoeff*t
     #coff84)
      tFjet85 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,25.D0,yp,1.D0,Tcoeff*tco
     #ff85)
      tFjet86 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,25.D0,yp,1.D0,Tcoeff*t
     #coff86)
      tFjet87 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,26.D0,yp,1.D0,Tcoeff*tco
     #ff87)
      tFjet88 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,26.D0,yp,1.D0,Tcoeff*t
     #coff88)

CCC Sum

      qiqbicnum20e0 = tFjet81*tcoff81+tFjet82*tcoff82+tFjet83*tcoff83+tF
     #jet84*tcoff84+tFjet85*tcoff85+tFjet86*tcoff86+tFjet87*tcoff87+tFje
     #t88*tcoff88
      RETURN
      END

      double precision function qiqbicnum20e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp*l1
      t8 = 64-64*l1+64*t6
      t25 = zp**2
      t28 = l1**2
      taux41 = 2.D0/3.D0*(1+t2)/yp*(t8*s2l2*s1l1+t8*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t8*(-s1l2-s2l2-zp))*s2l1+t8*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(2*t25*l1-t25*t28-2*t25+2*zp*t28-1+2*zp-2*t6-t28))
      t2 = (1-yp)**2
      t6 = 3*zp
      t7 = 3*l1
      t8 = zp*l1
      t9 = 3*t8
      t10 = 1+t6+t7-t9
      t13 = t7-t9-3+6*zp
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = t6+t8-l1-3
      taux42 = 2.D0/3.D0*(1+t2)/yp*((256*t10*s2l2+256*t13*t14)*s1l1+(-25
     #6*t10*s2l1-256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)
      t2 = (1-yp)**2
      t6 = zp*l1
      t8 = 64-64*l1+64*t6
      t25 = zp**2
      t28 = l1**2
      taux43 = 2.D0/3.D0*(1+t2)/yp*(t8*s2l2*s1l1+t8*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t8*(-s1l2-s2l2-zp))*s2l1+t8*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(2*t25*l1-t25*t28-2*t25+2*zp*t28-1+2*zp-2*t6-t28))
      t2 = (1-yp)**2
      t6 = 3*zp
      t7 = 3*l1
      t8 = zp*l1
      t9 = 3*t8
      t10 = 1+t6+t7-t9
      t13 = t7-t9-3+6*zp
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = t6+t8-l1-3
      taux44 = 2.D0/3.D0*(1+t2)/yp*((256*t10*s2l2+256*t13*t14)*s1l1+(-25
     #6*t10*s2l1-256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)

CCC Coefficients

      tcoff81 = taux41
      tcoff82 = taux41
      tcoff83 = taux42
      tcoff84 = taux42
      tcoff85 = taux43
      tcoff86 = taux43
      tcoff87 = -taux44
      tcoff88 = -taux44

CCC Jet functions

      tFjet81 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,19.D0,1.D0,yp,Tcoeff*tco
     #ff81)
      tFjet82 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,19.D0,1.D0,yp,Tcoeff*t
     #coff82)
      tFjet83 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,20.D0,1.D0,yp,Tcoeff*tco
     #ff83)
      tFjet84 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,20.D0,1.D0,yp,Tcoeff*t
     #coff84)
      tFjet85 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,25.D0,yp,1.D0,Tcoeff*tco
     #ff85)
      tFjet86 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,25.D0,yp,1.D0,Tcoeff*t
     #coff86)
      tFjet87 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,26.D0,yp,1.D0,Tcoeff*tco
     #ff87)
      tFjet88 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,26.D0,yp,1.D0,Tcoeff*t
     #coff88)

CCC Sum

      qiqbicnum20e1 = tFjet81*tcoff81+tFjet82*tcoff82+tFjet83*tcoff83+tF
     #jet84*tcoff84+tFjet85*tcoff85+tFjet86*tcoff86+tFjet87*tcoff87+tFje
     #t88*tcoff88
      RETURN
      END

      double precision function qiqbicnum20e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp**2
      t7 = l1**2
      taux41 = -128.D0/3.D0*(1+t2)/yp*zp*(t6*t7+2*zp*l1+1-2*zp*t7+t7-2*l
     #1)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 2*zp
      t8 = t6-3-l1-t7
      t11 = -l1+t6+1-t7
      t12 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t26 = -3*t6+1+3*l1-t7
      taux42 = 2.D0/3.D0*(1+t2)/yp*((256*t8*s2l2+512*t11*t12)*s1l1+(-256
     #*t8*s2l1-512*t11*t19)*s1l2+256*t26*t12*s2l1-256*t26*t19*s2l2)
      t2 = (1-yp)**2
      t6 = zp**2
      t7 = l1**2
      taux43 = -128.D0/3.D0*(1+t2)/yp*zp*(t6*t7+2*zp*l1+1-2*zp*t7+t7-2*l
     #1)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 2*zp
      t8 = t6-3-l1-t7
      t11 = -l1+t6+1-t7
      t12 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t26 = -3*t6+1+3*l1-t7
      taux44 = 2.D0/3.D0*(1+t2)/yp*((256*t8*s2l2+512*t11*t12)*s1l1+(-256
     #*t8*s2l1-512*t11*t19)*s1l2+256*t26*t12*s2l1-256*t26*t19*s2l2)

CCC Coefficients

      tcoff81 = taux41
      tcoff82 = taux41
      tcoff83 = taux42
      tcoff84 = taux42
      tcoff85 = taux43
      tcoff86 = taux43
      tcoff87 = -taux44
      tcoff88 = -taux44

CCC Jet functions

      tFjet81 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,19.D0,1.D0,yp,Tcoeff*tco
     #ff81)
      tFjet82 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,19.D0,1.D0,yp,Tcoeff*t
     #coff82)
      tFjet83 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,20.D0,1.D0,yp,Tcoeff*tco
     #ff83)
      tFjet84 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,20.D0,1.D0,yp,Tcoeff*t
     #coff84)
      tFjet85 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,25.D0,yp,1.D0,Tcoeff*tco
     #ff85)
      tFjet86 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,25.D0,yp,1.D0,Tcoeff*t
     #coff86)
      tFjet87 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,26.D0,yp,1.D0,Tcoeff*tco
     #ff87)
      tFjet88 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,26.D0,yp,1.D0,Tcoeff*t
     #coff88)

CCC Sum

      qiqbicnum20e2 = tFjet81*tcoff81+tFjet82*tcoff82+tFjet83*tcoff83+tF
     #jet84*tcoff84+tFjet85*tcoff85+tFjet86*tcoff86+tFjet87*tcoff87+tFje
     #t88*tcoff88
      RETURN
      END

      double precision function qiqbicnum20e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux41 = 0
      t2 = (1-yp)**2
      t7 = 1-l1+zp*l1
      taux42 = 2.D0/3.D0*(1+t2)/yp*(512*t7*s2l2*s1l1-512*t7*s2l1*s1l2+51
     #2*t7*(-s1l2-s2l2-zp)*s2l1-512*t7*(-s1l1-s2l1-zp)*s2l2)
      taux43 = 0
      t2 = (1-yp)**2
      t7 = 1-l1+zp*l1
      taux44 = 2.D0/3.D0*(1+t2)/yp*(512*t7*s2l2*s1l1-512*t7*s2l1*s1l2+51
     #2*t7*(-s1l2-s2l2-zp)*s2l1-512*t7*(-s1l1-s2l1-zp)*s2l2)

CCC Coefficients

      tcoff81 = taux41
      tcoff82 = taux41
      tcoff83 = taux42
      tcoff84 = taux42
      tcoff85 = taux43
      tcoff86 = taux43
      tcoff87 = -taux44
      tcoff88 = -taux44

CCC Jet functions

      tFjet81 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,19.D0,1.D0,yp,Tcoeff*tco
     #ff81)
      tFjet82 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,19.D0,1.D0,yp,Tcoeff*t
     #coff82)
      tFjet83 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,20.D0,1.D0,yp,Tcoeff*tco
     #ff83)
      tFjet84 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,20.D0,1.D0,yp,Tcoeff*t
     #coff84)
      tFjet85 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,25.D0,yp,1.D0,Tcoeff*tco
     #ff85)
      tFjet86 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,25.D0,yp,1.D0,Tcoeff*t
     #coff86)
      tFjet87 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,26.D0,yp,1.D0,Tcoeff*tco
     #ff87)
      tFjet88 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,26.D0,yp,1.D0,Tcoeff*t
     #coff88)

CCC Sum

      qiqbicnum20e3 = tFjet81*tcoff81+tFjet82*tcoff82+tFjet83*tcoff83+tF
     #jet84*tcoff84+tFjet85*tcoff85+tFjet86*tcoff86+tFjet87*tcoff87+tFje
     #t88*tcoff88
      RETURN
      END

      double precision function qiqbicnum21e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp-1
      t10 = zp*l1
      t12 = -64*l1-64*zp+64*t10
      t17 = -64*l1+64*t10+64-128*zp
      t18 = -s1l2-s2l2-zp
      t23 = -s1l1-s2l1-zp
      t28 = -64*t6
      taux45 = 2.D0/3.D0*(1+t2)/yp*((128*t6*l1*s1l2+t12*s2l2+t17*t18)*s1
     #l1+(t12*s2l1+t17*t23)*s1l2+t28*t18*s2l1+t28*t23*s2l2+128*t23*t18)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = -l1-zp+t6
      t10 = -l1+t6+1-2*zp
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = -zp+1
      taux46 = 2.D0/3.D0*(1+t2)/yp*((256*t7*s2l2+256*t10*t11)*s1l1+(-256
     #*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2)
      t2 = (1-yp)**2
      t6 = zp-1
      t10 = zp*l1
      t12 = -64*l1-64*zp+64*t10
      t17 = -64*l1+64*t10+64-128*zp
      t18 = -s1l2-s2l2-zp
      t23 = -s1l1-s2l1-zp
      t28 = -64*t6
      taux47 = 2.D0/3.D0*(1+t2)/yp*((128*t6*l1*s1l2+t12*s2l2+t17*t18)*s1
     #l1+(t12*s2l1+t17*t23)*s1l2+t28*t18*s2l1+t28*t23*s2l2+128*t23*t18)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = -l1-zp+t6
      t10 = -l1+t6+1-2*zp
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = -zp+1
      taux48 = 2.D0/3.D0*(1+t2)/yp*((256*t7*s2l2+256*t10*t11)*s1l1+(-256
     #*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2)

CCC Coefficients

      tcoff89 = taux45
      tcoff90 = taux45
      tcoff91 = taux46
      tcoff92 = taux46
      tcoff93 = taux47
      tcoff94 = taux47
      tcoff95 = -taux48
      tcoff96 = -taux48

CCC Jet functions

      tFjet89 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,23.D0,1.D0,yp,Tcoeff*tco
     #ff89)
      tFjet90 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,23.D0,1.D0,yp,Tcoeff*t
     #coff90)
      tFjet91 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,24.D0,1.D0,yp,Tcoeff*tco
     #ff91)
      tFjet92 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,24.D0,1.D0,yp,Tcoeff*t
     #coff92)
      tFjet93 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,21.D0,yp,1.D0,Tcoeff*t
     #coff93)
      tFjet94 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,21.D0,yp,1.D0,Tcoeff*tco
     #ff94)
      tFjet95 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,22.D0,yp,1.D0,Tcoeff*t
     #coff95)
      tFjet96 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,22.D0,yp,1.D0,Tcoeff*tco
     #ff96)

CCC Sum

      qiqbicnum21e0 = tFjet89*tcoff89+tFjet90*tcoff90+tFjet91*tcoff91+tF
     #jet92*tcoff92+tFjet93*tcoff93+tFjet94*tcoff94+tFjet95*tcoff95+tFje
     #t96*tcoff96
      RETURN
      END

      double precision function qiqbicnum21e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp*l1
      t8 = 64-64*l1+64*t6
      t25 = zp**2
      t28 = l1**2
      taux45 = 2.D0/3.D0*(1+t2)/yp*(t8*s2l2*s1l1+t8*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t8*(-s1l2-s2l2-zp))*s2l1+t8*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(2*t25*l1-t25*t28-2*t25+2*zp*t28-1+2*zp-2*t6-t28))
      t2 = (1-yp)**2
      t6 = 3*zp
      t7 = 3*l1
      t8 = zp*l1
      t9 = 3*t8
      t10 = 1+t6+t7-t9
      t13 = t7-t9-3+6*zp
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = t6+t8-l1-3
      taux46 = 2.D0/3.D0*(1+t2)/yp*((256*t10*s2l2+256*t13*t14)*s1l1+(-25
     #6*t10*s2l1-256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)
      t2 = (1-yp)**2
      t6 = zp*l1
      t8 = 64-64*l1+64*t6
      t25 = zp**2
      t28 = l1**2
      taux47 = 2.D0/3.D0*(1+t2)/yp*(t8*s2l2*s1l1+t8*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t8*(-s1l2-s2l2-zp))*s2l1+t8*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(2*t25*l1-t25*t28-2*t25+2*zp*t28-1+2*zp-2*t6-t28))
      t2 = (1-yp)**2
      t6 = 3*zp
      t7 = 3*l1
      t8 = zp*l1
      t9 = 3*t8
      t10 = 1+t6+t7-t9
      t13 = t7-t9-3+6*zp
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = t6+t8-l1-3
      taux48 = 2.D0/3.D0*(1+t2)/yp*((256*t10*s2l2+256*t13*t14)*s1l1+(-25
     #6*t10*s2l1-256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)

CCC Coefficients

      tcoff89 = taux45
      tcoff90 = taux45
      tcoff91 = taux46
      tcoff92 = taux46
      tcoff93 = taux47
      tcoff94 = taux47
      tcoff95 = -taux48
      tcoff96 = -taux48

CCC Jet functions

      tFjet89 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,23.D0,1.D0,yp,Tcoeff*tco
     #ff89)
      tFjet90 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,23.D0,1.D0,yp,Tcoeff*t
     #coff90)
      tFjet91 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,24.D0,1.D0,yp,Tcoeff*tco
     #ff91)
      tFjet92 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,24.D0,1.D0,yp,Tcoeff*t
     #coff92)
      tFjet93 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,21.D0,yp,1.D0,Tcoeff*t
     #coff93)
      tFjet94 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,21.D0,yp,1.D0,Tcoeff*tco
     #ff94)
      tFjet95 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,22.D0,yp,1.D0,Tcoeff*t
     #coff95)
      tFjet96 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,22.D0,yp,1.D0,Tcoeff*tco
     #ff96)

CCC Sum

      qiqbicnum21e1 = tFjet89*tcoff89+tFjet90*tcoff90+tFjet91*tcoff91+tF
     #jet92*tcoff92+tFjet93*tcoff93+tFjet94*tcoff94+tFjet95*tcoff95+tFje
     #t96*tcoff96
      RETURN
      END

      double precision function qiqbicnum21e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp**2
      t7 = l1**2
      taux45 = -128.D0/3.D0*(1+t2)/yp*zp*(t6*t7+2*zp*l1+1-2*zp*t7+t7-2*l
     #1)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 2*zp
      t8 = t6-3-l1-t7
      t11 = -l1+t6+1-t7
      t12 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t26 = -3*t6+1+3*l1-t7
      taux46 = 2.D0/3.D0*(1+t2)/yp*((256*t8*s2l2+512*t11*t12)*s1l1+(-256
     #*t8*s2l1-512*t11*t19)*s1l2+256*t26*t12*s2l1-256*t26*t19*s2l2)
      t2 = (1-yp)**2
      t6 = zp**2
      t7 = l1**2
      taux47 = -128.D0/3.D0*(1+t2)/yp*zp*(t6*t7+2*zp*l1+1-2*zp*t7+t7-2*l
     #1)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 2*zp
      t8 = t6-3-l1-t7
      t11 = -l1+t6+1-t7
      t12 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t26 = -3*t6+1+3*l1-t7
      taux48 = 2.D0/3.D0*(1+t2)/yp*((256*t8*s2l2+512*t11*t12)*s1l1+(-256
     #*t8*s2l1-512*t11*t19)*s1l2+256*t26*t12*s2l1-256*t26*t19*s2l2)

CCC Coefficients

      tcoff89 = taux45
      tcoff90 = taux45
      tcoff91 = taux46
      tcoff92 = taux46
      tcoff93 = taux47
      tcoff94 = taux47
      tcoff95 = -taux48
      tcoff96 = -taux48

CCC Jet functions

      tFjet89 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,23.D0,1.D0,yp,Tcoeff*tco
     #ff89)
      tFjet90 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,23.D0,1.D0,yp,Tcoeff*t
     #coff90)
      tFjet91 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,24.D0,1.D0,yp,Tcoeff*tco
     #ff91)
      tFjet92 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,24.D0,1.D0,yp,Tcoeff*t
     #coff92)
      tFjet93 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,21.D0,yp,1.D0,Tcoeff*t
     #coff93)
      tFjet94 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,21.D0,yp,1.D0,Tcoeff*tco
     #ff94)
      tFjet95 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,22.D0,yp,1.D0,Tcoeff*t
     #coff95)
      tFjet96 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,22.D0,yp,1.D0,Tcoeff*tco
     #ff96)

CCC Sum

      qiqbicnum21e2 = tFjet89*tcoff89+tFjet90*tcoff90+tFjet91*tcoff91+tF
     #jet92*tcoff92+tFjet93*tcoff93+tFjet94*tcoff94+tFjet95*tcoff95+tFje
     #t96*tcoff96
      RETURN
      END

      double precision function qiqbicnum21e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux45 = 0
      t2 = (1-yp)**2
      t7 = 1-l1+zp*l1
      taux46 = 2.D0/3.D0*(1+t2)/yp*(512*t7*s2l2*s1l1-512*t7*s2l1*s1l2+51
     #2*t7*(-s1l2-s2l2-zp)*s2l1-512*t7*(-s1l1-s2l1-zp)*s2l2)
      taux47 = 0
      t2 = (1-yp)**2
      t7 = 1-l1+zp*l1
      taux48 = 2.D0/3.D0*(1+t2)/yp*(512*t7*s2l2*s1l1-512*t7*s2l1*s1l2+51
     #2*t7*(-s1l2-s2l2-zp)*s2l1-512*t7*(-s1l1-s2l1-zp)*s2l2)

CCC Coefficients

      tcoff89 = taux45
      tcoff90 = taux45
      tcoff91 = taux46
      tcoff92 = taux46
      tcoff93 = taux47
      tcoff94 = taux47
      tcoff95 = -taux48
      tcoff96 = -taux48

CCC Jet functions

      tFjet89 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,23.D0,1.D0,yp,Tcoeff*tco
     #ff89)
      tFjet90 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,23.D0,1.D0,yp,Tcoeff*t
     #coff90)
      tFjet91 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,24.D0,1.D0,yp,Tcoeff*tco
     #ff91)
      tFjet92 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,24.D0,1.D0,yp,Tcoeff*t
     #coff92)
      tFjet93 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,21.D0,yp,1.D0,Tcoeff*t
     #coff93)
      tFjet94 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,21.D0,yp,1.D0,Tcoeff*tco
     #ff94)
      tFjet95 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,22.D0,yp,1.D0,Tcoeff*t
     #coff95)
      tFjet96 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,22.D0,yp,1.D0,Tcoeff*tco
     #ff96)

CCC Sum

      qiqbicnum21e3 = tFjet89*tcoff89+tFjet90*tcoff90+tFjet91*tcoff91+tF
     #jet92*tcoff92+tFjet93*tcoff93+tFjet94*tcoff94+tFjet95*tcoff95+tFje
     #t96*tcoff96
      RETURN
      END

