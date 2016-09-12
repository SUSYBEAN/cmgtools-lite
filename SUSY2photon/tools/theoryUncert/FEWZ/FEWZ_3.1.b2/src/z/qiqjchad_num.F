CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqjcnum1e0(var,ndim,z,ux,l1,
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
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,1.D0,z,1.D0,Tcoeff*
     #tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,1.D0,z,1.D0,Tcoef
     #f*tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,10.D0,z,1.D0,Tcoeff
     #*tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,10.D0,z,1.D0,Tcoe
     #ff*tcoff4)

CCC Sum

      qiqjcnum1e0 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tco
     #ff4
      RETURN
      END

      double precision function qiqjcnum1e1(var,ndim,z,ux,l1,
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
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,1.D0,z,1.D0,Tcoeff*
     #tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,1.D0,z,1.D0,Tcoef
     #f*tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,10.D0,z,1.D0,Tcoeff
     #*tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,10.D0,z,1.D0,Tcoe
     #ff*tcoff4)

CCC Sum

      qiqjcnum1e1 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tco
     #ff4
      RETURN
      END

      double precision function qiqjcnum1e2(var,ndim,z,ux,l1,
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
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,1.D0,z,1.D0,Tcoeff*
     #tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,1.D0,z,1.D0,Tcoef
     #f*tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,10.D0,z,1.D0,Tcoeff
     #*tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,10.D0,z,1.D0,Tcoe
     #ff*tcoff4)

CCC Sum

      qiqjcnum1e2 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tco
     #ff4
      RETURN
      END

      double precision function qiqjcnum1e3(var,ndim,z,ux,l1,
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
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,1.D0,z,1.D0,Tcoeff*
     #tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,1.D0,z,1.D0,Tcoef
     #f*tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,10.D0,z,1.D0,Tcoeff
     #*tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,10.D0,z,1.D0,Tcoe
     #ff*tcoff4)

CCC Sum

      qiqjcnum1e3 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tco
     #ff4
      RETURN
      END

      double precision function qiqjcnum2e0(var,ndim,z,ux,l1,
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
      taux3 = 4*z*(1+z)*t5*t7-4.D0/3.D0*z*(3+15*z+8*t12)*t5*t6+8.D0/9.D0
     #*(z-1)*(28*t12+z+10)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux4 = -16*z*(1+z)*t5*t7+16.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-32.D0/
     #9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,1.D0,z,1.D0,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,1.D0,z,1.D0,Tcoef
     #f*tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,10.D0,z,1.D0,Tcoeff
     #*tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,10.D0,z,1.D0,Tcoe
     #ff*tcoff8)

CCC Sum

      qiqjcnum2e0 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tco
     #ff8
      RETURN
      END

      double precision function qiqjcnum2e1(var,ndim,z,ux,l1,
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
      taux3 = 4*z*(1+z)*t4-4.D0/3.D0*z*(3+15*z+8*t8)*t3+8.D0/9.D0*(z-1)*
     #(28*t8+z+10)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t7 = t6**2
      t12 = z**2
      taux4 = 48*z*(1+z)*t5*t7-16*z*(3+15*z+8*t12)*t5*t6+32.D0/3.D0*(z-1
     #)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,1.D0,z,1.D0,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,1.D0,z,1.D0,Tcoef
     #f*tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,10.D0,z,1.D0,Tcoeff
     #*tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,10.D0,z,1.D0,Tcoe
     #ff*tcoff8)

CCC Sum

      qiqjcnum2e1 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tco
     #ff8
      RETURN
      END

      double precision function qiqjcnum2e2(var,ndim,z,ux,l1,
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
      t7 = t6**2
      t12 = z**2
      taux4 = -32*z*(1+z)*t5*t7+32.D0/3.D0*z*(3+15*z+8*t12)*t5*t6-64.D0/
     #9.D0*(z-1)*(28*t12+z+10)*t5

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,1.D0,z,1.D0,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,1.D0,z,1.D0,Tcoef
     #f*tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,10.D0,z,1.D0,Tcoeff
     #*tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,10.D0,z,1.D0,Tcoe
     #ff*tcoff8)

CCC Sum

      qiqjcnum2e2 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tco
     #ff8
      RETURN
      END

      double precision function qiqjcnum2e3(var,ndim,z,ux,l1,
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
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,1.D0,z,1.D0,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,1.D0,z,1.D0,Tcoef
     #f*tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,10.D0,z,1.D0,Tcoeff
     #*tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,10.D0,z,1.D0,Tcoe
     #ff*tcoff8)

CCC Sum

      qiqjcnum2e3 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tco
     #ff8
      RETURN
      END

      double precision function qiqjcnum3e0(var,ndim,z,ux,l1,
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
      taux5 = -8*z*(1+z)*t5*t6+4.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux6 = 32*z*(1+z)*t5*t6-16.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux5
      tcoff11 = taux6
      tcoff12 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,1.D0,z,Tcoeff
     #*tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,1.D0,z,Tco
     #eff*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,1.D0,z,Tcoef
     #f*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,1.D0,z,Tco
     #eff*tcoff12)

CCC Sum

      qiqjcnum3e0 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet1
     #2*tcoff12
      RETURN
      END

      double precision function qiqjcnum3e1(var,ndim,z,ux,l1,
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
      taux5 = -8*z*(1+z)*t3+4.D0/3.D0*(z-1)*(4*t7+7*z+4)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux6 = -96*z*(1+z)*t5*t6+16*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux5
      tcoff11 = taux6
      tcoff12 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,1.D0,z,Tcoeff
     #*tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,1.D0,z,Tco
     #eff*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,1.D0,z,Tcoef
     #f*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,1.D0,z,Tco
     #eff*tcoff12)

CCC Sum

      qiqjcnum3e1 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet1
     #2*tcoff12
      RETURN
      END

      double precision function qiqjcnum3e2(var,ndim,z,ux,l1,
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
      t11 = z**2
      taux6 = 64*z*(1+z)*t5*t6-32.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux5
      tcoff11 = taux6
      tcoff12 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,1.D0,z,Tcoeff
     #*tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,1.D0,z,Tco
     #eff*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,1.D0,z,Tcoef
     #f*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,1.D0,z,Tco
     #eff*tcoff12)

CCC Sum

      qiqjcnum3e2 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet1
     #2*tcoff12
      RETURN
      END

      double precision function qiqjcnum3e3(var,ndim,z,ux,l1,
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
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,1.D0,z,Tcoeff
     #*tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,1.D0,z,Tco
     #eff*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,1.D0,z,Tcoef
     #f*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,1.D0,z,Tco
     #eff*tcoff12)

CCC Sum

      qiqjcnum3e3 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet1
     #2*tcoff12
      RETURN
      END

      double precision function qiqjcnum4e0(var,ndim,z,ux,l1,
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
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,1.D0,z,Tcoef
     #f*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,1.D0,z,Tco
     #eff*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,1.D0,z,Tcoef
     #f*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,1.D0,z,Tco
     #eff*tcoff16)

CCC Sum

      qiqjcnum4e0 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFje
     #t16*tcoff16
      RETURN
      END

      double precision function qiqjcnum4e1(var,ndim,z,ux,l1,
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
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,1.D0,z,Tcoef
     #f*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,1.D0,z,Tco
     #eff*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,1.D0,z,Tcoef
     #f*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,1.D0,z,Tco
     #eff*tcoff16)

CCC Sum

      qiqjcnum4e1 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFje
     #t16*tcoff16
      RETURN
      END

      double precision function qiqjcnum4e2(var,ndim,z,ux,l1,
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
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,1.D0,z,Tcoef
     #f*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,1.D0,z,Tco
     #eff*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,1.D0,z,Tcoef
     #f*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,1.D0,z,Tco
     #eff*tcoff16)

CCC Sum

      qiqjcnum4e2 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFje
     #t16*tcoff16
      RETURN
      END

      double precision function qiqjcnum4e3(var,ndim,z,ux,l1,
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
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,11.D0,1.D0,z,Tcoef
     #f*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,11.D0,1.D0,z,Tco
     #eff*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,12.D0,1.D0,z,Tcoef
     #f*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,12.D0,1.D0,z,Tco
     #eff*tcoff16)

CCC Sum

      qiqjcnum4e3 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFje
     #t16*tcoff16
      RETURN
      END

      double precision function qiqjcnum5e0(var,ndim,z,ux,l1,
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
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoef
     #f*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tco
     #eff*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoef
     #f*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tco
     #eff*tcoff20)

CCC Sum

      qiqjcnum5e0 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFje
     #t20*tcoff20
      RETURN
      END

      double precision function qiqjcnum5e1(var,ndim,z,ux,l1,
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
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoef
     #f*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tco
     #eff*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoef
     #f*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tco
     #eff*tcoff20)

CCC Sum

      qiqjcnum5e1 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFje
     #t20*tcoff20
      RETURN
      END

      double precision function qiqjcnum5e2(var,ndim,z,ux,l1,
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
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoef
     #f*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tco
     #eff*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoef
     #f*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tco
     #eff*tcoff20)

CCC Sum

      qiqjcnum5e2 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFje
     #t20*tcoff20
      RETURN
      END

      double precision function qiqjcnum5e3(var,ndim,z,ux,l1,
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
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoef
     #f*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tco
     #eff*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoef
     #f*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tco
     #eff*tcoff20)

CCC Sum

      qiqjcnum5e3 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFje
     #t20*tcoff20
      RETURN
      END

      double precision function qiqjcnum6e0(var,ndim,z,ux,l1,
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
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoef
     #f*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tco
     #eff*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoef
     #f*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tco
     #eff*tcoff24)

CCC Sum

      qiqjcnum6e0 = tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFje
     #t24*tcoff24
      RETURN
      END

      double precision function qiqjcnum6e1(var,ndim,z,ux,l1,
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
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoef
     #f*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tco
     #eff*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoef
     #f*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tco
     #eff*tcoff24)

CCC Sum

      qiqjcnum6e1 = tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFje
     #t24*tcoff24
      RETURN
      END

      double precision function qiqjcnum6e2(var,ndim,z,ux,l1,
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
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoef
     #f*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tco
     #eff*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoef
     #f*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tco
     #eff*tcoff24)

CCC Sum

      qiqjcnum6e2 = tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFje
     #t24*tcoff24
      RETURN
      END

      double precision function qiqjcnum6e3(var,ndim,z,ux,l1,
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
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,13.D0,1.D0,z,Tcoef
     #f*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,13.D0,1.D0,z,Tco
     #eff*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,14.D0,1.D0,z,Tcoef
     #f*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,14.D0,1.D0,z,Tco
     #eff*tcoff24)

CCC Sum

      qiqjcnum6e3 = tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFje
     #t24*tcoff24
      RETURN
      END

      double precision function qiqjcnum7e0(var,ndim,z,ux,l1,
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
      t17 = 64-64*l1-128*zp+64*t10
      t18 = -s1l2-s2l2-zp
      t23 = -s1l1-s2l1-zp
      t28 = -64*t6
      taux13 = 2.D0/3.D0*(1+t2)/yp*((128*t6*l1*s1l2+t12*s2l2+t17*t18)*s1
     #l1+(t12*s2l1+t17*t23)*s1l2+t28*t18*s2l1+t28*t23*s2l2+128*t23*t18)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = -zp+t6-l1
      t10 = 1-l1-2*zp+t6
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = -zp+1
      taux14 = 2.D0/3.D0*(1+t2)/yp*((256*t7*s2l2+256*t10*t11)*s1l1+(-256
     #*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2)

CCC Coefficients

      tcoff25 = taux13
      tcoff26 = taux13
      tcoff27 = taux14
      tcoff28 = taux14

CCC Jet functions

      tFjet25 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,1.D0,yp,1.D0,Tcoeff*tcof
     #f25)
      tFjet26 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,1.D0,yp,1.D0,Tcoeff*tc
     #off26)
      tFjet27 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,10.D0,yp,1.D0,Tcoeff*tco
     #ff27)
      tFjet28 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,10.D0,yp,1.D0,Tcoeff*t
     #coff28)

CCC Sum

      qiqjcnum7e0 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFje
     #t28*tcoff28
      RETURN
      END

      double precision function qiqjcnum7e1(var,ndim,z,ux,l1,
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
      t8 = -64*l1+64*t6+64
      t25 = zp**2
      t30 = l1**2
      taux13 = 2.D0/3.D0*(1+t2)/yp*(t8*s2l2*s1l1+t8*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t8*(-s1l2-s2l2-zp))*s2l1+t8*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(2*t25*l1-2*t25+2*zp-t25*t30+2*zp*t30-2*t6-t30-1))
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 3*t6
      t8 = 3*l1
      t9 = 3*zp
      t10 = -t7+t8+t9+1
      t13 = -3+t8+6*zp-t7
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = -3-l1+t6+t9
      taux14 = 2.D0/3.D0*(1+t2)/yp*((256*t10*s2l2+256*t13*t14)*s1l1+(-25
     #6*t10*s2l1-256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)

CCC Coefficients

      tcoff25 = taux13
      tcoff26 = taux13
      tcoff27 = taux14
      tcoff28 = taux14

CCC Jet functions

      tFjet25 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,1.D0,yp,1.D0,Tcoeff*tcof
     #f25)
      tFjet26 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,1.D0,yp,1.D0,Tcoeff*tc
     #off26)
      tFjet27 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,10.D0,yp,1.D0,Tcoeff*tco
     #ff27)
      tFjet28 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,10.D0,yp,1.D0,Tcoeff*t
     #coff28)

CCC Sum

      qiqjcnum7e1 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFje
     #t28*tcoff28
      RETURN
      END

      double precision function qiqjcnum7e2(var,ndim,z,ux,l1,
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
      t8 = zp**2
      t9 = l1**2
      taux13 = -128.D0/3.D0*(1+t2)/yp*zp*(2*zp*l1+t8*t9-2*l1+1-2*zp*t9+t
     #9)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 2*zp
      t8 = -3-l1+t6-t7
      t11 = 1-l1-t7+t6
      t12 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t26 = 1+3*l1-3*t6-t7
      taux14 = 2.D0/3.D0*(1+t2)/yp*((256*t8*s2l2+512*t11*t12)*s1l1+(-256
     #*t8*s2l1-512*t11*t19)*s1l2+256*t26*t12*s2l1-256*t26*t19*s2l2)

CCC Coefficients

      tcoff25 = taux13
      tcoff26 = taux13
      tcoff27 = taux14
      tcoff28 = taux14

CCC Jet functions

      tFjet25 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,1.D0,yp,1.D0,Tcoeff*tcof
     #f25)
      tFjet26 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,1.D0,yp,1.D0,Tcoeff*tc
     #off26)
      tFjet27 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,10.D0,yp,1.D0,Tcoeff*tco
     #ff27)
      tFjet28 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,10.D0,yp,1.D0,Tcoeff*t
     #coff28)

CCC Sum

      qiqjcnum7e2 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFje
     #t28*tcoff28
      RETURN
      END

      double precision function qiqjcnum7e3(var,ndim,z,ux,l1,
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

      taux13 = 0
      t2 = (1-yp)**2
      t7 = -l1+zp*l1+1
      taux14 = 2.D0/3.D0*(1+t2)/yp*(512*t7*s2l2*s1l1-512*t7*s2l1*s1l2+51
     #2*t7*(-s1l2-s2l2-zp)*s2l1-512*t7*(-s1l1-s2l1-zp)*s2l2)

CCC Coefficients

      tcoff25 = taux13
      tcoff26 = taux13
      tcoff27 = taux14
      tcoff28 = taux14

CCC Jet functions

      tFjet25 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,1.D0,yp,1.D0,Tcoeff*tcof
     #f25)
      tFjet26 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,1.D0,yp,1.D0,Tcoeff*tc
     #off26)
      tFjet27 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,10.D0,yp,1.D0,Tcoeff*tco
     #ff27)
      tFjet28 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,10.D0,yp,1.D0,Tcoeff*t
     #coff28)

CCC Sum

      qiqjcnum7e3 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFje
     #t28*tcoff28
      RETURN
      END

      double precision function qiqjcnum8e0(var,ndim,z,ux,l1,
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
      t17 = 64-64*l1-128*zp+64*t10
      t18 = -s1l2-s2l2-zp
      t23 = -s1l1-s2l1-zp
      t28 = -64*t6
      taux15 = 2.D0/3.D0*(1+t2)/yp*((128*t6*l1*s1l2+t12*s2l2+t17*t18)*s1
     #l1+(t12*s2l1+t17*t23)*s1l2+t28*t18*s2l1+t28*t23*s2l2+128*t23*t18)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = -zp+t6-l1
      t10 = 1-l1-2*zp+t6
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = -zp+1
      taux16 = 2.D0/3.D0*(1+t2)/yp*((256*t7*s2l2+256*t10*t11)*s1l1+(-256
     #*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2)

CCC Coefficients

      tcoff29 = taux15
      tcoff30 = taux15
      tcoff31 = -taux16
      tcoff32 = -taux16

CCC Jet functions

      tFjet29 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*tco
     #ff29)
      tFjet30 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*t
     #coff30)
      tFjet31 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*tco
     #ff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*t
     #coff32)

CCC Sum

      qiqjcnum8e0 = tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFje
     #t32*tcoff32
      RETURN
      END

      double precision function qiqjcnum8e1(var,ndim,z,ux,l1,
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
      t8 = -64*l1+64*t6+64
      t25 = zp**2
      t30 = l1**2
      taux15 = 2.D0/3.D0*(1+t2)/yp*(t8*s2l2*s1l1+t8*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t8*(-s1l2-s2l2-zp))*s2l1+t8*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(2*t25*l1-2*t25+2*zp-t25*t30+2*zp*t30-2*t6-t30-1))
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 3*t6
      t8 = 3*l1
      t9 = 3*zp
      t10 = -t7+t8+t9+1
      t13 = -3+t8+6*zp-t7
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = -3-l1+t6+t9
      taux16 = 2.D0/3.D0*(1+t2)/yp*((256*t10*s2l2+256*t13*t14)*s1l1+(-25
     #6*t10*s2l1-256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)

CCC Coefficients

      tcoff29 = taux15
      tcoff30 = taux15
      tcoff31 = -taux16
      tcoff32 = -taux16

CCC Jet functions

      tFjet29 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*tco
     #ff29)
      tFjet30 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*t
     #coff30)
      tFjet31 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*tco
     #ff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*t
     #coff32)

CCC Sum

      qiqjcnum8e1 = tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFje
     #t32*tcoff32
      RETURN
      END

      double precision function qiqjcnum8e2(var,ndim,z,ux,l1,
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
      t8 = zp**2
      t9 = l1**2
      taux15 = -128.D0/3.D0*(1+t2)/yp*zp*(2*zp*l1+t8*t9-2*l1+1-2*zp*t9+t
     #9)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 2*zp
      t8 = -l1+t6-3-t7
      t11 = 1-l1-t7+t6
      t12 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t26 = 3*l1-3*t6+1-t7
      taux16 = 2.D0/3.D0*(1+t2)/yp*((256*t8*s2l2+512*t11*t12)*s1l1+(-256
     #*t8*s2l1-512*t11*t19)*s1l2+256*t26*t12*s2l1-256*t26*t19*s2l2)

CCC Coefficients

      tcoff29 = taux15
      tcoff30 = taux15
      tcoff31 = -taux16
      tcoff32 = -taux16

CCC Jet functions

      tFjet29 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*tco
     #ff29)
      tFjet30 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*t
     #coff30)
      tFjet31 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*tco
     #ff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*t
     #coff32)

CCC Sum

      qiqjcnum8e2 = tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFje
     #t32*tcoff32
      RETURN
      END

      double precision function qiqjcnum8e3(var,ndim,z,ux,l1,
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

      taux15 = 0
      t2 = (1-yp)**2
      t7 = -l1+zp*l1+1
      taux16 = 2.D0/3.D0*(1+t2)/yp*(512*t7*s2l2*s1l1-512*t7*s2l1*s1l2+51
     #2*t7*(-s1l2-s2l2-zp)*s2l1-512*t7*(-s1l1-s2l1-zp)*s2l2)

CCC Coefficients

      tcoff29 = taux15
      tcoff30 = taux15
      tcoff31 = -taux16
      tcoff32 = -taux16

CCC Jet functions

      tFjet29 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*tco
     #ff29)
      tFjet30 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,11.D0,yp,1.D0,Tcoeff*t
     #coff30)
      tFjet31 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*tco
     #ff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0
     #,zp,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,12.D0,yp,1.D0,Tcoeff*t
     #coff32)

CCC Sum

      qiqjcnum8e3 = tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFje
     #t32*tcoff32
      RETURN
      END

      double precision function qiqjcnum9e0(var,ndim,z,ux,l1,
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
      t17 = 64-64*l1-128*zp+64*t10
      t18 = -s1l2-s2l2-zp
      t23 = -s1l1-s2l1-zp
      t28 = -64*t6
      taux17 = 2.D0/3.D0*(1+t2)/yp*((128*t6*l1*s1l2+t12*s2l2+t17*t18)*s1
     #l1+(t12*s2l1+t17*t23)*s1l2+t28*t18*s2l1+t28*t23*s2l2+128*t23*t18)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = -zp+t6-l1
      t10 = 1-l1-2*zp+t6
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = -zp+1
      taux18 = 2.D0/3.D0*(1+t2)/yp*((256*t7*s2l2+256*t10*t11)*s1l1+(-256
     #*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2)

CCC Coefficients

      tcoff33 = taux17
      tcoff34 = taux17
      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet33 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*t
     #coff33)
      tFjet34 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*tco
     #ff34)
      tFjet35 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*t
     #coff35)
      tFjet36 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*tco
     #ff36)

CCC Sum

      qiqjcnum9e0 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFje
     #t36*tcoff36
      RETURN
      END

      double precision function qiqjcnum9e1(var,ndim,z,ux,l1,
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
      t8 = -64*l1+64*t6+64
      t25 = zp**2
      t30 = l1**2
      taux17 = 2.D0/3.D0*(1+t2)/yp*(t8*s2l2*s1l1+t8*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t8*(-s1l2-s2l2-zp))*s2l1+t8*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(2*t25*l1-2*t25+2*zp-t25*t30+2*zp*t30-2*t6-t30-1))
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 3*t6
      t8 = 3*l1
      t9 = 3*zp
      t10 = -t7+t8+t9+1
      t13 = -3+t8+6*zp-t7
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = -3-l1+t6+t9
      taux18 = 2.D0/3.D0*(1+t2)/yp*((256*t10*s2l2+256*t13*t14)*s1l1+(-25
     #6*t10*s2l1-256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)

CCC Coefficients

      tcoff33 = taux17
      tcoff34 = taux17
      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet33 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*t
     #coff33)
      tFjet34 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*tco
     #ff34)
      tFjet35 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*t
     #coff35)
      tFjet36 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*tco
     #ff36)

CCC Sum

      qiqjcnum9e1 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFje
     #t36*tcoff36
      RETURN
      END

      double precision function qiqjcnum9e2(var,ndim,z,ux,l1,
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
      t8 = zp**2
      t9 = l1**2
      taux17 = -128.D0/3.D0*(1+t2)/yp*zp*(2*zp*l1+t8*t9-2*l1+1-2*zp*t9+t
     #9)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 2*zp
      t8 = -l1+t6-3-t7
      t11 = 1-l1-t7+t6
      t12 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t26 = 3*l1-3*t6+1-t7
      taux18 = 2.D0/3.D0*(1+t2)/yp*((256*t8*s2l2+512*t11*t12)*s1l1+(-256
     #*t8*s2l1-512*t11*t19)*s1l2+256*t26*t12*s2l1-256*t26*t19*s2l2)

CCC Coefficients

      tcoff33 = taux17
      tcoff34 = taux17
      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet33 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*t
     #coff33)
      tFjet34 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*tco
     #ff34)
      tFjet35 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*t
     #coff35)
      tFjet36 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*tco
     #ff36)

CCC Sum

      qiqjcnum9e2 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFje
     #t36*tcoff36
      RETURN
      END

      double precision function qiqjcnum9e3(var,ndim,z,ux,l1,
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

      taux17 = 0
      t2 = (1-yp)**2
      t7 = -l1+zp*l1+1
      taux18 = 2.D0/3.D0*(1+t2)/yp*(512*t7*s2l2*s1l1-512*t7*s2l1*s1l2+51
     #2*t7*(-s1l2-s2l2-zp)*s2l1-512*t7*(-s1l1-s2l1-zp)*s2l2)

CCC Coefficients

      tcoff33 = taux17
      tcoff34 = taux17
      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet33 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*t
     #coff33)
      tFjet34 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,13.D0,1.D0,yp,Tcoeff*tco
     #ff34)
      tFjet35 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*t
     #coff35)
      tFjet36 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,14.D0,1.D0,yp,Tcoeff*tco
     #ff36)

CCC Sum

      qiqjcnum9e3 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFje
     #t36*tcoff36
      RETURN
      END

