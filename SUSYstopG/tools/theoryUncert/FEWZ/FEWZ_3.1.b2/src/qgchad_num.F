CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgcnum1e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t10 = s1l1*s2l2
      taux1 = -48*t1*t3-24*t1*z+48*t1*t2+48*t10*t3+24*t10*z-48*t10*t2

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff*
     #tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoef
     #f*tcoff2)

CCC Sum

      qgcnum1e0 = tFjet1*tcoff1+tFjet2*tcoff2
      RETURN
      END

      double precision function qgcnum1e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t10 = s1l1*s2l2
      taux1 = 144*t1*t3+72*t1*z-144*t1*t2-144*t10*t3-72*t10*z+144*t10*t2

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff*
     #tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoef
     #f*tcoff2)

CCC Sum

      qgcnum1e1 = tFjet1*tcoff1+tFjet2*tcoff2
      RETURN
      END

      double precision function qgcnum1e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t10 = s1l1*s2l2
      taux1 = -96*t1*t3-48*t1*z+96*t1*t2+96*t10*t3+48*t10*z-96*t10*t2

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff*
     #tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoef
     #f*tcoff2)

CCC Sum

      qgcnum1e2 = tFjet1*tcoff1+tFjet2*tcoff2
      RETURN
      END

      double precision function qgcnum1e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff*
     #tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoef
     #f*tcoff2)

CCC Sum

      qgcnum1e3 = tFjet1*tcoff1+tFjet2*tcoff2
      RETURN
      END

      double precision function qgcnum2e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t6 = -s1l2*s2l1-s1l1*s2l2
      t7 = dlog(z)
      t11 = z**2
      t17 = dlog(1-z)
      t21 = t11*z
      taux2 = -18*z*(4*z+1)*t6*t7-18*z*(2*t11-2*z+1)*t6*t17+t6*(-24+4*NF
     #*t21-4*NF*t11+2*NF*z-78*t11-51*z+120*t21)/2
      t6 = s1l1*s2l2-s1l2*s2l1
      t7 = dlog(z)
      t11 = z**2
      t17 = dlog(1-z)
      t21 = t11*z
      taux3 = 72*z*(4*z+1)*t6*t7+72*z*(2*t11-2*z+1)*t6*t17-2*t6*(-24+4*N
     #F*t21-4*NF*t11+2*NF*z-78*t11-51*z+120*t21)

CCC Coefficients

      tcoff3 = taux2
      tcoff4 = taux2
      tcoff5 = taux3
      tcoff6 = taux3

CCC Jet functions

      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff*
     #tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoef
     #f*tcoff4)
      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoef
     #f*tcoff6)

CCC Sum

      qgcnum2e0 = tFjet3*tcoff3+tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tcoff
     #6
      RETURN
      END

      double precision function qgcnum2e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t4 = dlog(z)
      t7 = z**2
      t13 = dlog(1-z)
      t16 = t7*z
      taux2 = -18*z*(4*z+1)*t4-18*z*(2*t7-2*z+1)*t13-12+2*NF*t16-2*NF*t7
     #+NF*z-39*t7-51.D0/2.D0*z+60*t16
      t6 = s1l1*s2l2-s1l2*s2l1
      t7 = dlog(z)
      t11 = z**2
      t17 = dlog(1-z)
      t23 = t11*z
      taux3 = -216*z*(4*z+1)*t6*t7-216*z*(2*t11-2*z+1)*t6*t17+6*t6*(-24+
     #4*NF*t23-4*NF*t11+2*NF*z-78*t11-51*z+120*t23)

CCC Coefficients

      tcoff3 = taux2
      tcoff4 = taux2
      tcoff5 = taux3
      tcoff6 = taux3

CCC Jet functions

      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff*
     #tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoef
     #f*tcoff4)
      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoef
     #f*tcoff6)

CCC Sum

      qgcnum2e1 = tFjet3*tcoff3+tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tcoff
     #6
      RETURN
      END

      double precision function qgcnum2e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux2 = 0
      t6 = s1l1*s2l2-s1l2*s2l1
      t7 = dlog(z)
      t11 = z**2
      t17 = dlog(1-z)
      t21 = t11*z
      taux3 = 144*z*(4*z+1)*t6*t7+144*z*(2*t11-2*z+1)*t6*t17-4*t6*(-24+4
     #*NF*t21-4*NF*t11+2*NF*z-78*t11-51*z+120*t21)

CCC Coefficients

      tcoff3 = taux2
      tcoff4 = taux2
      tcoff5 = taux3
      tcoff6 = taux3

CCC Jet functions

      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff*
     #tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoef
     #f*tcoff4)
      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoef
     #f*tcoff6)

CCC Sum

      qgcnum2e2 = tFjet3*tcoff3+tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tcoff
     #6
      RETURN
      END

      double precision function qgcnum2e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux2 = 0
      taux3 = 0

CCC Coefficients

      tcoff3 = taux2
      tcoff4 = taux2
      tcoff5 = taux3
      tcoff6 = taux3

CCC Jet functions

      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff*
     #tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoef
     #f*tcoff4)
      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoef
     #f*tcoff6)

CCC Sum

      qgcnum2e3 = tFjet3*tcoff3+tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tcoff
     #6
      RETURN
      END

      double precision function qgcnum3e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = 2*z
      t2 = z**2
      t8 = -s1l2*s2l1-s1l1*s2l2
      t9 = dlog(z)
      t17 = dlog(1-z)
      taux4 = 4*z*(1-t1+4*t2)*t8*t9-8*z*(2*t2-t1+1)*t8*t17-2*z*(-1+4*z)*
     #t8
      t1 = 2*z
      t2 = z**2
      t8 = s1l1*s2l2-s1l2*s2l1
      t9 = dlog(z)
      t17 = dlog(1-z)
      taux5 = -16*z*(1-t1+4*t2)*t8*t9+32*z*(2*t2-t1+1)*t8*t17+8*z*(-1+4*
     #z)*t8

CCC Coefficients

      tcoff7 = taux4
      tcoff8 = taux4
      tcoff9 = taux5
      tcoff10 = taux5

CCC Jet functions

      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff*
     #tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoef
     #f*tcoff8)
      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff*
     #tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff10)

CCC Sum

      qgcnum3e0 = tFjet7*tcoff7+tFjet8*tcoff8+tFjet9*tcoff9+tFjet10*tcof
     #f10
      RETURN
      END

      double precision function qgcnum3e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = 2*z
      t2 = z**2
      t6 = dlog(z)
      t13 = dlog(1-z)
      taux4 = 4*z*(1-t1+4*t2)*t6-8*z*(2*t2-t1+1)*t13-2*z*(-1+4*z)
      t1 = 2*z
      t2 = z**2
      t8 = s1l1*s2l2-s1l2*s2l1
      t9 = dlog(z)
      t17 = dlog(1-z)
      taux5 = 48*z*(1-t1+4*t2)*t8*t9-96*z*(2*t2-t1+1)*t8*t17-24*z*(-1+4*
     #z)*t8

CCC Coefficients

      tcoff7 = taux4
      tcoff8 = taux4
      tcoff9 = taux5
      tcoff10 = taux5

CCC Jet functions

      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff*
     #tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoef
     #f*tcoff8)
      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff*
     #tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff10)

CCC Sum

      qgcnum3e1 = tFjet7*tcoff7+tFjet8*tcoff8+tFjet9*tcoff9+tFjet10*tcof
     #f10
      RETURN
      END

      double precision function qgcnum3e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux4 = 0
      t1 = 2*z
      t2 = z**2
      t8 = s1l1*s2l2-s1l2*s2l1
      t9 = dlog(z)
      t17 = dlog(1-z)
      taux5 = -32*z*(1-t1+4*t2)*t8*t9+64*z*(2*t2-t1+1)*t8*t17+16*z*(-1+4
     #*z)*t8

CCC Coefficients

      tcoff7 = taux4
      tcoff8 = taux4
      tcoff9 = taux5
      tcoff10 = taux5

CCC Jet functions

      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff*
     #tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoef
     #f*tcoff8)
      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff*
     #tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff10)

CCC Sum

      qgcnum3e2 = tFjet7*tcoff7+tFjet8*tcoff8+tFjet9*tcoff9+tFjet10*tcof
     #f10
      RETURN
      END

      double precision function qgcnum3e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux4 = 0
      taux5 = 0

CCC Coefficients

      tcoff7 = taux4
      tcoff8 = taux4
      tcoff9 = taux5
      tcoff10 = taux5

CCC Jet functions

      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff*
     #tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoef
     #f*tcoff8)
      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff*
     #tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff10)

CCC Sum

      qgcnum3e3 = tFjet7*tcoff7+tFjet8*tcoff8+tFjet9*tcoff9+tFjet10*tcof
     #f10
      RETURN
      END

      double precision function qgcnum4e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t10 = s1l1*s2l2
      taux6 = 12*t1*t3+6*t1*z-12*t1*t2+12*t10*t3+6*t10*z-12*t10*t2
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t10 = s1l1*s2l2
      taux7 = -48*t1*t3-24*t1*z+48*t1*t2+48*t10*t3+24*t10*z-48*t10*t2

CCC Coefficients

      tcoff11 = taux6
      tcoff12 = taux6
      tcoff13 = taux7
      tcoff14 = taux7

CCC Jet functions

      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff12)
      tFjet13 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff14)

CCC Sum

      qgcnum4e0 = tFjet11*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet1
     #4*tcoff14
      RETURN
      END

      double precision function qgcnum4e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z**2
      taux6 = -12*t1*z-6*z+12*t1
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t10 = s1l1*s2l2
      taux7 = 144*t1*t3+72*t1*z-144*t1*t2-144*t10*t3-72*t10*z+144*t10*t2

CCC Coefficients

      tcoff11 = taux6
      tcoff12 = taux6
      tcoff13 = taux7
      tcoff14 = taux7

CCC Jet functions

      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff12)
      tFjet13 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff14)

CCC Sum

      qgcnum4e1 = tFjet11*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet1
     #4*tcoff14
      RETURN
      END

      double precision function qgcnum4e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux6 = 0
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t10 = s1l1*s2l2
      taux7 = -96*t1*t3-48*t1*z+96*t1*t2+96*t10*t3+48*t10*z-96*t10*t2

CCC Coefficients

      tcoff11 = taux6
      tcoff12 = taux6
      tcoff13 = taux7
      tcoff14 = taux7

CCC Jet functions

      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff12)
      tFjet13 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff14)

CCC Sum

      qgcnum4e2 = tFjet11*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet1
     #4*tcoff14
      RETURN
      END

      double precision function qgcnum4e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux6 = 0
      taux7 = 0

CCC Coefficients

      tcoff11 = taux6
      tcoff12 = taux6
      tcoff13 = taux7
      tcoff14 = taux7

CCC Jet functions

      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff12)
      tFjet13 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff14)

CCC Sum

      qgcnum4e3 = tFjet11*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet1
     #4*tcoff14
      RETURN
      END

      double precision function qgcnum5e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = s1l1*s2l2
      t2 = z**2
      t3 = t2*z
      t4 = dlog(z)
      t5 = t3*t4
      t7 = dlog(1+z)
      t8 = t5*t7
      t11 = s1l2*s2l1
      t12 = Li2(-z)
      t13 = z*t12
      t16 = t3*t12
      t19 = 1-z
      t22 = dlog(t19/z)
      t23 = t22**2
      t24 = t23*z
      t27 = t23*t2
      t30 = dlog(t19)
      t31 = t30**2
      t32 = t31*t2
      t35 = t4**2
      t36 = t35*z
      t39 = t22*z
      t42 = t31*z
      t48 = t35*t2
      t57 = -36*t1*t8-18*t11*t13-36*t11*t16+4*t1*t24-8*t1*t27+18*t11*t32
     #-11*t11*t36-8*t11*t39-9*t11*t42-8*t11*t27+20*t11-14*t1*t48+4*t11*t
     #24-14*t11*t48-8*t1*t39
      t58 = t22*t2
      t67 = t31*t3
      t76 = z*t30
      t79 = t2*t4
      t82 = t22*t3
      t85 = t23*t3
      t88 = t2*t30
      t91 = t3*t30
      t96 = 16*t1*t58-9*t1*t42+18*t1*t32-18*t1*t13-18*t11*t67-18*t1*t67-
     #11*t1*t36+16*t11*t58+8*t11*t76+80*t1*t79-16*t1*t82+8*t1*t85-36*t11
     #*t88+36*t11*t91-36*t1*t88
      t100 = z*t4
      t115 = 0.3141592653589793D1**2
      t116 = t115*t3
      t130 = 36*t1*t91+7*t1*t100+80*t11*t79+7*t11*t100-16*t11*t82+8*t11*
     #t85+132*t1*t5-36*t1*t16-8.D0/3.D0*t1*t116+132*t11*t5+8*t1*t76+20*t
     #1+10*t1*z+167*t1*t2-178*t1*t3
      t139 = t79*t7
      t142 = t115*z
      t145 = t115*t2
      t148 = t2*t12
      t159 = t100*t7
      t166 = -178*t11*t3+10*t11*z+167*t11*t2-36*t11*t8-36*t11*t139-4.D0/
     #3.D0*t11*t142-10.D0/3.D0*t1*t145-36*t1*t148-4.D0/3.D0*t1*t142-36*t
     #11*t148-8.D0/3.D0*t11*t116-10.D0/3.D0*t11*t145-18*t11*t159-36*t1*t
     #139-18*t1*t159
      taux8 = t57+t96+t130+t166
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = dlog(z)
      t4 = t2*t3
      t6 = dlog(1+z)
      t7 = t4*t6
      t14 = z*t3
      t15 = t14*t6
      t18 = s1l1*s2l2
      t24 = 1-z
      t27 = dlog(t24/z)
      t28 = t27*t2
      t31 = dlog(t24)
      t32 = t31**2
      t33 = t32*z
      t36 = t32*t2
      t39 = Li2(-z)
      t40 = z*t39
      t43 = t2*z
      t44 = t32*t43
      t47 = t43*t3
      t48 = t47*t6
      t51 = 0.3141592653589793D1**2
      t52 = t51*t2
      t55 = t2*t39
      t58 = 144*t1*t7-40*t1*z-668*t1*t2+72*t1*t15-144*t18*t7-80*t1-72*t1
     #8*t15+64*t18*t28-36*t18*t33+72*t18*t36-72*t18*t40+72*t1*t44-144*t1
     #8*t48-40.D0/3.D0*t18*t52-144*t18*t55
      t59 = t51*z
      t64 = t3**2
      t65 = t64*z
      t70 = t31*z
      t77 = t51*t43
      t88 = t27*z
      t93 = t27*t43
      t96 = -16.D0/3.D0*t18*t59-72*t18*t44-44*t18*t65-64*t1*t28-32*t1*t7
     #0+320*t18*t4+144*t1*t55+32.D0/3.D0*t1*t77+40.D0/3.D0*t1*t52+16.D0/
     #3.D0*t1*t59-72*t1*t36+44*t1*t65+32*t1*t88+36*t1*t33-64*t18*t93
      t98 = t27**2
      t99 = t98*t43
      t102 = t31*t2
      t105 = t43*t31
      t122 = t43*t39
      t125 = t98*z
      t128 = t98*t2
      t135 = 32*t18*t99+144*t1*t102-144*t1*t105-144*t18*t102+144*t18*t10
     #5+28*t18*t14-320*t1*t4-28*t1*t14+64*t1*t93+72*t1*t40+144*t1*t122+1
     #6*t18*t125-32*t18*t128-32*t1*t99+528*t18*t47
      t157 = t64*t2
      t166 = -144*t18*t122-32.D0/3.D0*t18*t77-528*t1*t47+40*t18*z+668*t1
     #8*t2-712*t18*t43+712*t1*t43+32*t18*t70+80*t18+144*t1*t48+32*t1*t12
     #8-56*t18*t157-16*t1*t125+56*t1*t157-32*t18*t88
      taux9 = t58+t96+t135+t166

CCC Coefficients

      tcoff15 = taux8
      tcoff16 = taux8
      tcoff17 = taux9
      tcoff18 = taux9

CCC Jet functions

      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff16)
      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff18)

CCC Sum

      qgcnum5e0 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet1
     #8*tcoff18
      RETURN
      END

      double precision function qgcnum5e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z**2
      t2 = t1*z
      t5 = Li2(-z)
      t9 = 1-z
      t10 = dlog(t9)
      t15 = dlog(z)
      t16 = t15**2
      t19 = z*t15
      t23 = t2*t15
      t25 = t1*t15
      t29 = dlog(t9/z)
      t30 = t29**2
      t37 = -20+178*t2-167*t1+36*t2*t5-10*z+36*t1*t10-8*z*t10+14*t1*t16-
     #7*t19-36*t2*t10-132*t23-80*t25+8*t30*t1-4*z*t30-16*t29*t1
      t46 = t10**2
      t53 = 0.3141592653589793D1**2
      t65 = dlog(1+z)
      t72 = 11*t16*z+8*t29*z-8*t2*t30+16*t29*t2+18*t46*t2-18*t46*t1+9*t4
     #6*z+10.D0/3.D0*t53*t1+4.D0/3.D0*t53*z+36*t1*t5+18*z*t5+8.D0/3.D0*t
     #53*t2+36*t23*t65+36*t25*t65+18*t19*t65
      taux8 = t37+t72
      t1 = s1l2*s2l1
      t2 = 1-z
      t5 = dlog(t2/z)
      t6 = t5**2
      t7 = t6*z
      t10 = dlog(z)
      t11 = t10**2
      t12 = z**2
      t13 = t11*t12
      t16 = s1l1*s2l2
      t17 = z*t5
      t20 = t5*t12
      t23 = dlog(t2)
      t24 = t23**2
      t25 = t24*z
      t28 = t24*t12
      t31 = Li2(-z)
      t32 = t31*z
      t35 = z*t12
      t36 = t24*t35
      t41 = t11*z
      t46 = t23*z
      t49 = t12*t10
      t52 = t5*t35
      t55 = t6*t35
      t58 = 48*t1*t7-168*t1*t13+96*t16*t17-192*t16*t20+108*t16*t25-216*t
     #16*t28+216*t16*t32-216*t1*t36+216*t16*t36+132*t16*t41+192*t1*t20+9
     #6*t46*t1-960*t16*t49+192*t16*t52-96*t16*t55
      t59 = t12*t23
      t62 = t35*t23
      t69 = z*t10
      t72 = 0.3141592653589793D1**2
      t73 = t72*t35
      t76 = t12*t31
      t79 = t72*z
      t82 = t72*t12
      t97 = -432*t1*t59+432*t1*t62+432*t16*t59-432*t16*t62-84*t16*t69-32
     #*t1*t73-432*t1*t76+16*t16*t79+40*t16*t82+432*t16*t76-120*t16*z-200
     #4*t16*t12+2136*t16*t35-2136*t1*t35+120*t1*z
      t102 = dlog(1+z)
      t103 = t69*t102
      t106 = t49*t102
      t111 = t35*t10
      t112 = t111*t102
      t121 = t35*t31
      t126 = t6*t12
      t137 = 2004*t1*t12-216*t1*t103+432*t16*t106+216*t16*t103+432*t16*t
     #112-40*t1*t82-16*t1*t79-216*t1*t32-432*t1*t121-48*t16*t7+96*t16*t1
     #26+216*t1*t28-132*t1*t41-96*t1*t17-108*t1*t25
      t166 = -96*t1*t126+168*t16*t13-432*t1*t112-432*t1*t106+240*t1-240*
     #t16+960*t1*t49+84*t1*t69-192*t1*t52+96*t1*t55-1584*t16*t111+432*t1
     #6*t121+32*t16*t73+1584*t1*t111-96*t16*t46
      taux9 = t58+t97+t137+t166

CCC Coefficients

      tcoff15 = taux8
      tcoff16 = taux8
      tcoff17 = taux9
      tcoff18 = taux9

CCC Jet functions

      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff16)
      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff18)

CCC Sum

      qgcnum5e1 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet1
     #8*tcoff18
      RETURN
      END

      double precision function qgcnum5e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux8 = 0
      t1 = s1l1*s2l2
      t2 = dlog(z)
      t3 = z*t2
      t5 = dlog(1+z)
      t6 = t3*t5
      t9 = z**2
      t10 = t9*z
      t11 = t10*t2
      t12 = t11*t5
      t21 = s1l2*s2l1
      t30 = 0.3141592653589793D1**2
      t31 = t30*t9
      t34 = Li2(-z)
      t35 = t9*t34
      t38 = t30*z
      t43 = t30*t10
      t48 = -144*t1*t6-288*t1*t12+80*t1*z+1336*t1*t9-1424*t1*t10+1424*t2
     #1*t10-80*t21*z-1336*t21*t9+56*t1*t3-80.D0/3.D0*t1*t31-288*t1*t35-3
     #2.D0/3.D0*t1*t38+288*t21*t35+64.D0/3.D0*t21*t43+80.D0/3.D0*t21*t31
      t51 = z*t34
      t54 = t10*t34
      t57 = 1-z
      t60 = dlog(t57/z)
      t61 = t60**2
      t62 = t61*z
      t65 = t61*t9
      t68 = dlog(t57)
      t69 = t68**2
      t70 = t69*t9
      t73 = t2**2
      t74 = t73*z
      t77 = t60*z
      t80 = t69*z
      t85 = t73*t9
      t94 = t60*t9
      t97 = 32.D0/3.D0*t21*t38+144*t21*t51+288*t21*t54+32*t1*t62-64*t1*t
     #65-144*t21*t70+88*t21*t74+64*t21*t77+72*t21*t80+64*t21*t65-112*t1*
     #t85-32*t21*t62+112*t21*t85-64*t1*t77+128*t1*t94
      t105 = t69*t10
      t114 = z*t68
      t117 = t9*t2
      t120 = t60*t10
      t123 = t61*t10
      t126 = t9*t68
      t129 = t10*t68
      t136 = -72*t1*t80+144*t1*t70-144*t1*t51+144*t21*t105-144*t1*t105-8
     #8*t1*t74-128*t21*t94-64*t21*t114+640*t1*t117-128*t1*t120+64*t1*t12
     #3+288*t21*t126-288*t21*t129-288*t1*t126-640*t21*t117
      t155 = t117*t5
      t166 = -56*t21*t3+128*t21*t120-64*t21*t123+1056*t1*t11-288*t1*t54-
     #64.D0/3.D0*t1*t43-1056*t21*t11+64*t1*t114+288*t21*t12+288*t21*t155
     #-160*t21+160*t1+288*t1*t129+144*t21*t6-288*t1*t155
      taux9 = t48+t97+t136+t166

CCC Coefficients

      tcoff15 = taux8
      tcoff16 = taux8
      tcoff17 = taux9
      tcoff18 = taux9

CCC Jet functions

      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff16)
      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff18)

CCC Sum

      qgcnum5e2 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet1
     #8*tcoff18
      RETURN
      END

      double precision function qgcnum5e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux8 = 0
      taux9 = 0

CCC Coefficients

      tcoff15 = taux8
      tcoff16 = taux8
      tcoff17 = taux9
      tcoff18 = taux9

CCC Jet functions

      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff16)
      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff18)

CCC Sum

      qgcnum5e3 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet1
     #8*tcoff18
      RETURN
      END

      double precision function qgcnum6e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = s1l2*s2l1
      t2 = z**2
      t3 = z*t2
      t5 = 12*t1*t3
      t7 = 6*t1*z
      t9 = 12*t1*t2
      t10 = s1l1*s2l2
      t12 = 12*t10*t3
      t14 = 6*t10*z
      t16 = 12*t10*t2
      t18 = dlog(1-z)
      taux10 = t5+t7-t9+t12+t14-t16+4.D0/3.D0*t18*(t5+t7-t9+t12+t14-t16)
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = z*t2
      t5 = 48*t1*t3
      t7 = 24*t1*z
      t9 = 48*t1*t2
      t10 = s1l1*s2l2
      t12 = 48*t10*t3
      t14 = 24*t10*z
      t16 = 48*t10*t2
      t18 = dlog(1-z)
      taux11 = -t5-t7+t9+t12+t14-t16+4.D0/3.D0*t18*(-t5-t7+t9+t12+t14-t1
     #6)

CCC Coefficients

      tcoff19 = taux10
      tcoff20 = taux10
      tcoff21 = taux11
      tcoff22 = taux11

CCC Jet functions

      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff20)
      tFjet21 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff22)

CCC Sum

      qgcnum6e0 = tFjet19*tcoff19+tFjet20*tcoff20+tFjet21*tcoff21+tFjet2
     #2*tcoff22
      RETURN
      END

      double precision function qgcnum6e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z**2
      t3 = 12*t1*z
      t4 = 6*z
      t5 = 12*t1
      t7 = dlog(1-z)
      taux10 = -t3-t4+t5+4.D0/3.D0*t7*(-t3-t4+t5)
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = z*t2
      t5 = 144*t1*t3
      t7 = 72*t1*z
      t9 = 144*t1*t2
      t10 = s1l1*s2l2
      t12 = 144*t10*t3
      t14 = 72*t10*z
      t16 = 144*t10*t2
      t18 = dlog(1-z)
      taux11 = t5+t7-t9-t12-t14+t16+4.D0/3.D0*t18*(t5+t7-t9-t12-t14+t16)

CCC Coefficients

      tcoff19 = taux10
      tcoff20 = taux10
      tcoff21 = taux11
      tcoff22 = taux11

CCC Jet functions

      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff20)
      tFjet21 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff22)

CCC Sum

      qgcnum6e1 = tFjet19*tcoff19+tFjet20*tcoff20+tFjet21*tcoff21+tFjet2
     #2*tcoff22
      RETURN
      END

      double precision function qgcnum6e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux10 = 0
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = z*t2
      t5 = 96*t1*t3
      t7 = 48*t1*z
      t9 = 96*t1*t2
      t10 = s1l1*s2l2
      t12 = 96*t10*t3
      t14 = 48*t10*z
      t16 = 96*t10*t2
      t18 = dlog(1-z)
      taux11 = -t5-t7+t9+t12+t14-t16+4.D0/3.D0*t18*(-t5-t7+t9+t12+t14-t1
     #6)

CCC Coefficients

      tcoff19 = taux10
      tcoff20 = taux10
      tcoff21 = taux11
      tcoff22 = taux11

CCC Jet functions

      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff20)
      tFjet21 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff22)

CCC Sum

      qgcnum6e2 = tFjet19*tcoff19+tFjet20*tcoff20+tFjet21*tcoff21+tFjet2
     #2*tcoff22
      RETURN
      END

      double precision function qgcnum6e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux10 = 0
      taux11 = 0

CCC Coefficients

      tcoff19 = taux10
      tcoff20 = taux10
      tcoff21 = taux11
      tcoff22 = taux11

CCC Jet functions

      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff20)
      tFjet21 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff22)

CCC Sum

      qgcnum6e3 = tFjet19*tcoff19+tFjet20*tcoff20+tFjet21*tcoff21+tFjet2
     #2*tcoff22
      RETURN
      END

      double precision function qgcnum7e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t6 = z**2
      t9 = (1-z)*l2+z
      t10 = t9**2
      t13 = 1/t9
      t16 = (1-z*t13)**2
      taux12 = (24*s1l2*s2l1+24*s1l1*s2l2)*z*(t6/t10/4+t16/4)*t13
      t6 = z**2
      t9 = (1-z)*l2+z
      t10 = t9**2
      t13 = 1/t9
      t16 = (1-z*t13)**2
      taux13 = (-96*s1l2*s2l1+96*s1l1*s2l2)*z*(t6/t10/4+t16/4)*t13

CCC Coefficients

      tcoff23 = taux12
      tcoff24 = taux12
      tcoff25 = taux13
      tcoff26 = taux13

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,t3,z/t3,Tcoef
     #f*tcoff23)
      t4 = (1-z)*l2+z
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,t4,z/t4,Tco
     #eff*tcoff24)
      t3 = (1-z)*l2+z
      tFjet25 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,t3,z/t3,Tcoef
     #f*tcoff25)
      t4 = (1-z)*l2+z
      tFjet26 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,t4,z/t4,Tco
     #eff*tcoff26)

CCC Sum

      qgcnum7e0 = tFjet23*tcoff23+tFjet24*tcoff24+tFjet25*tcoff25+tFjet2
     #6*tcoff26
      RETURN
      END

      double precision function qgcnum7e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z**2
      t4 = (1-z)*l2+z
      t5 = t4**2
      t8 = 1/t4
      t11 = (1-z*t8)**2
      taux12 = -24*z*(t1/t5/4+t11/4)*t8
      t6 = z**2
      t9 = (1-z)*l2+z
      t10 = t9**2
      t13 = 1/t9
      t16 = (1-z*t13)**2
      taux13 = (288*s1l2*s2l1-288*s1l1*s2l2)*z*(t6/t10/4+t16/4)*t13

CCC Coefficients

      tcoff23 = taux12
      tcoff24 = taux12
      tcoff25 = taux13
      tcoff26 = taux13

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,t3,z/t3,Tcoef
     #f*tcoff23)
      t4 = (1-z)*l2+z
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,t4,z/t4,Tco
     #eff*tcoff24)
      t3 = (1-z)*l2+z
      tFjet25 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,t3,z/t3,Tcoef
     #f*tcoff25)
      t4 = (1-z)*l2+z
      tFjet26 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,t4,z/t4,Tco
     #eff*tcoff26)

CCC Sum

      qgcnum7e1 = tFjet23*tcoff23+tFjet24*tcoff24+tFjet25*tcoff25+tFjet2
     #6*tcoff26
      RETURN
      END

      double precision function qgcnum7e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux12 = 0
      t6 = z**2
      t9 = (1-z)*l2+z
      t10 = t9**2
      t13 = 1/t9
      t16 = (1-z*t13)**2
      taux13 = (-192*s1l2*s2l1+192*s1l1*s2l2)*z*(t6/t10/4+t16/4)*t13

CCC Coefficients

      tcoff23 = taux12
      tcoff24 = taux12
      tcoff25 = taux13
      tcoff26 = taux13

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,t3,z/t3,Tcoef
     #f*tcoff23)
      t4 = (1-z)*l2+z
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,t4,z/t4,Tco
     #eff*tcoff24)
      t3 = (1-z)*l2+z
      tFjet25 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,t3,z/t3,Tcoef
     #f*tcoff25)
      t4 = (1-z)*l2+z
      tFjet26 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,t4,z/t4,Tco
     #eff*tcoff26)

CCC Sum

      qgcnum7e2 = tFjet23*tcoff23+tFjet24*tcoff24+tFjet25*tcoff25+tFjet2
     #6*tcoff26
      RETURN
      END

      double precision function qgcnum7e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux12 = 0
      taux13 = 0

CCC Coefficients

      tcoff23 = taux12
      tcoff24 = taux12
      tcoff25 = taux13
      tcoff26 = taux13

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,t3,z/t3,Tcoef
     #f*tcoff23)
      t4 = (1-z)*l2+z
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,t4,z/t4,Tco
     #eff*tcoff24)
      t3 = (1-z)*l2+z
      tFjet25 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,t3,z/t3,Tcoef
     #f*tcoff25)
      t4 = (1-z)*l2+z
      tFjet26 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,t4,z/t4,Tco
     #eff*tcoff26)

CCC Sum

      qgcnum7e3 = tFjet23*tcoff23+tFjet24*tcoff24+tFjet25*tcoff25+tFjet2
     #6*tcoff26
      RETURN
      END

      double precision function qgcnum8e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z**2
      t3 = (1-z)**2
      taux14 = z*(t1/4+t3/4)*(128*s1l1*s2l1+64*s1l1+64*s2l1)
      t1 = z**2
      t3 = (1-z)**2
      taux15 = z*(t1/4+t3/4)*(-256*s2l1+256*s1l1)

CCC Coefficients

      tcoff27 = taux14
      tcoff28 = taux14
      tcoff29 = taux15
      tcoff30 = taux15

CCC Jet functions

      tFjet27 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff27)
      tFjet28 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff28)
      tFjet29 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff29)
      tFjet30 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff30)

CCC Sum

      qgcnum8e0 = tFjet27*tcoff27+tFjet28*tcoff28+tFjet29*tcoff29+tFjet3
     #0*tcoff30
      RETURN
      END

      double precision function qgcnum8e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z**2
      t3 = (1-z)**2
      taux14 = z*(t1/4+t3/4)*(-96*s2l1-96*s1l1-64*s1l1*s2l1)
      t1 = z**2
      t3 = (1-z)**2
      taux15 = z*(t1/4+t3/4)*(896*s2l1-896*s1l1)

CCC Coefficients

      tcoff27 = taux14
      tcoff28 = taux14
      tcoff29 = taux15
      tcoff30 = taux15

CCC Jet functions

      tFjet27 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff27)
      tFjet28 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff28)
      tFjet29 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff29)
      tFjet30 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff30)

CCC Sum

      qgcnum8e1 = tFjet27*tcoff27+tFjet28*tcoff28+tFjet29*tcoff29+tFjet3
     #0*tcoff30
      RETURN
      END

      double precision function qgcnum8e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z**2
      t3 = (1-z)**2
      taux14 = z*(t1/4+t3/4)*(96*s1l1+96*s2l1+128*s1l1*s2l1)
      t1 = z**2
      t3 = (1-z)**2
      taux15 = z*(t1/4+t3/4)*(-896*s2l1+896*s1l1)

CCC Coefficients

      tcoff27 = taux14
      tcoff28 = taux14
      tcoff29 = taux15
      tcoff30 = taux15

CCC Jet functions

      tFjet27 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff27)
      tFjet28 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff28)
      tFjet29 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff29)
      tFjet30 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff30)

CCC Sum

      qgcnum8e2 = tFjet27*tcoff27+tFjet28*tcoff28+tFjet29*tcoff29+tFjet3
     #0*tcoff30
      RETURN
      END

      double precision function qgcnum8e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z**2
      t3 = (1-z)**2
      taux14 = z*(t1/4+t3/4)*(-64*s1l1-64*s2l1)
      t1 = z**2
      t3 = (1-z)**2
      taux15 = z*(t1/4+t3/4)*(384*s2l1-384*s1l1)

CCC Coefficients

      tcoff27 = taux14
      tcoff28 = taux14
      tcoff29 = taux15
      tcoff30 = taux15

CCC Jet functions

      tFjet27 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff27)
      tFjet28 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff28)
      tFjet29 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff29)
      tFjet30 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff30)

CCC Sum

      qgcnum8e3 = tFjet27*tcoff27+tFjet28*tcoff28+tFjet29*tcoff29+tFjet3
     #0*tcoff30
      RETURN
      END

      double precision function qgcnum9e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z-1
      t5 = z*l1
      t7 = -64*l1-64*z+64*t5
      t12 = 64-64*l1-128*z+64*t5
      t13 = -s1l2-s2l2-zp
      t18 = -s1l1-s2l1-zp
      t23 = -64*t1
      taux16 = (128*t1*l1*s1l2+t7*s2l2+t12*t13)*s1l1+(t7*s2l1+t12*t18)*s
     #1l2+t23*t13*s2l1+t23*t18*s2l2+128*t18*t13
      t1 = z*l1
      t2 = -z+t1-l1
      t5 = 1-l1-2*z+t1
      t6 = -s1l2-s2l2-zp
      t12 = -s1l1-s2l1-zp
      t17 = 1-z
      taux17 = (256*t2*s2l2+256*t5*t6)*s1l1+(-256*t2*s2l1-256*t5*t12)*s1
     #l2+256*t17*t6*s2l1-256*t17*t12*s2l2

CCC Coefficients

      tcoff31 = taux16
      tcoff32 = taux16
      tcoff33 = taux17
      tcoff34 = taux17

CCC Jet functions

      tFjet31 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff32)
      tFjet33 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff33)
      tFjet34 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff34)

CCC Sum

      qgcnum9e0 = tFjet31*tcoff31+tFjet32*tcoff32+tFjet33*tcoff33+tFjet3
     #4*tcoff34
      RETURN
      END

      double precision function qgcnum9e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z*l1
      t3 = 64-64*l1+64*t1
      t20 = z**2
      t23 = l1**2
      taux16 = t3*s2l2*s1l1+t3*s2l1*s1l2+(128*(z-1)*(-1+l1)*s2l2+t3*(-s1
     #l2-s2l2-zp))*s2l1+t3*(-s1l1-s2l1-zp)*s2l2-64*z*(2*t20*l1-t20*t23+2
     #*z-2*t20-2*t1+2*z*t23-1-t23)
      t1 = 3*z
      t2 = 3*l1
      t3 = z*l1
      t4 = 3*t3
      t5 = t1+1+t2-t4
      t8 = -3+t2+6*z-t4
      t9 = -s1l2-s2l2-zp
      t15 = -s1l1-s2l1-zp
      t20 = -3+t3+t1-l1
      taux17 = (256*t5*s2l2+256*t8*t9)*s1l1+(-256*t5*s2l1-256*t8*t15)*s1
     #l2+256*t20*t9*s2l1-256*t20*t15*s2l2

CCC Coefficients

      tcoff31 = taux16
      tcoff32 = taux16
      tcoff33 = taux17
      tcoff34 = taux17

CCC Jet functions

      tFjet31 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff32)
      tFjet33 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff33)
      tFjet34 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff34)

CCC Sum

      qgcnum9e1 = tFjet31*tcoff31+tFjet32*tcoff32+tFjet33*tcoff33+tFjet3
     #4*tcoff34
      RETURN
      END

      double precision function qgcnum9e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z**2
      t2 = l1**2
      taux16 = -64*z*(t1*t2+1+t2+2*z*l1-2*z*t2-2*l1)
      t1 = z*l1
      t2 = 2*z
      t3 = t1-3-l1-t2
      t6 = 1-l1-t2+t1
      t7 = -s1l2-s2l2-zp
      t14 = -s1l1-s2l1-zp
      t21 = -3*t1+1+3*l1-t2
      taux17 = (256*t3*s2l2+512*t6*t7)*s1l1+(-256*t3*s2l1-512*t6*t14)*s1
     #l2+256*t21*t7*s2l1-256*t21*t14*s2l2

CCC Coefficients

      tcoff31 = taux16
      tcoff32 = taux16
      tcoff33 = taux17
      tcoff34 = taux17

CCC Jet functions

      tFjet31 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff32)
      tFjet33 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff33)
      tFjet34 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff34)

CCC Sum

      qgcnum9e2 = tFjet31*tcoff31+tFjet32*tcoff32+tFjet33*tcoff33+tFjet3
     #4*tcoff34
      RETURN
      END

      double precision function qgcnum9e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux16 = 0
      t2 = z*l1+1-l1
      taux17 = 512*t2*s2l2*s1l1-512*t2*s2l1*s1l2+512*t2*(-s1l2-s2l2-zp)*
     #s2l1-512*t2*(-s1l1-s2l1-zp)*s2l2

CCC Coefficients

      tcoff31 = taux16
      tcoff32 = taux16
      tcoff33 = taux17
      tcoff34 = taux17

CCC Jet functions

      tFjet31 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff32)
      tFjet33 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff33)
      tFjet34 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff34)

CCC Sum

      qgcnum9e3 = tFjet31*tcoff31+tFjet32*tcoff32+tFjet33*tcoff33+tFjet3
     #4*tcoff34
      RETURN
      END

      double precision function qgcnum10e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z*l1
      t2 = -z+t1-l1
      t5 = 1-l1-2*z+t1
      t6 = -s1l2-s2l2-zp
      t12 = -s1l1-s2l1-zp
      t17 = 1-z
      taux18 = (256*t2*s2l2+256*t5*t6)*s1l1+(-256*t2*s2l1-256*t5*t12)*s1
     #l2+256*t17*t6*s2l1-256*t17*t12*s2l2

CCC Coefficients

      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet35 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff35)
      tFjet36 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff36)

CCC Sum

      qgcnum10e0 = tFjet35*tcoff35+tFjet36*tcoff36
      RETURN
      END

      double precision function qgcnum10e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = 3*z
      t2 = 3*l1
      t3 = z*l1
      t4 = 3*t3
      t5 = t1+1+t2-t4
      t8 = -3+t2+6*z-t4
      t9 = -s1l2-s2l2-zp
      t15 = -s1l1-s2l1-zp
      t20 = -3+t3+t1-l1
      taux18 = (256*t5*s2l2+256*t8*t9)*s1l1+(-256*t5*s2l1-256*t8*t15)*s1
     #l2+256*t20*t9*s2l1-256*t20*t15*s2l2

CCC Coefficients

      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet35 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff35)
      tFjet36 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff36)

CCC Sum

      qgcnum10e1 = tFjet35*tcoff35+tFjet36*tcoff36
      RETURN
      END

      double precision function qgcnum10e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z*l1
      t2 = 2*z
      t3 = t1-3-l1-t2
      t6 = 1-l1-t2+t1
      t7 = -s1l2-s2l2-zp
      t14 = -s1l1-s2l1-zp
      t21 = -3*t1+1+3*l1-t2
      taux18 = (256*t3*s2l2+512*t6*t7)*s1l1+(-256*t3*s2l1-512*t6*t14)*s1
     #l2+256*t21*t7*s2l1-256*t21*t14*s2l2

CCC Coefficients

      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet35 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff35)
      tFjet36 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff36)

CCC Sum

      qgcnum10e2 = tFjet35*tcoff35+tFjet36*tcoff36
      RETURN
      END

      double precision function qgcnum10e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t2 = z*l1+1-l1
      taux18 = 512*t2*s2l2*s1l1-512*t2*s2l1*s1l2+512*t2*(-s1l2-s2l2-zp)*
     #s2l1-512*t2*(-s1l1-s2l1-zp)*s2l2

CCC Coefficients

      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet35 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff35)
      tFjet36 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff36)

CCC Sum

      qgcnum10e3 = tFjet35*tcoff35+tFjet36*tcoff36
      RETURN
      END

      double precision function qgcnum11e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z-1
      t3 = -t1
      t4 = dlog(t3)
      t6 = 3+4*t4
      t10 = z*l1
      t12 = (-z+t10-l1)*t6
      t17 = (1-l1-2*z+t10)*t6
      t18 = -s1l2-s2l2-zp
      t24 = -s1l1-s2l1-zp
      t29 = t3*t6
      taux19 = (128.D0/3.D0*t1*l1*t6*s1l2+64.D0/3.D0*t12*s2l2+64.D0/3.D0
     #*t17*t18)*s1l1+(64.D0/3.D0*t12*s2l1+64.D0/3.D0*t17*t24)*s1l2+64.D0
     #/3.D0*t29*t18*s2l1+64.D0/3.D0*t29*t24*s2l2+(128+512.D0/3.D0*t4)*t1
     #8*t24
      t1 = z*l1
      t3 = 1-z
      t4 = dlog(t3)
      t6 = 3+4*t4
      t7 = (-z+t1-l1)*t6
      t11 = (1-l1-2*z+t1)*t6
      t12 = -s1l2-s2l2-zp
      t18 = -s1l1-s2l1-zp
      t23 = t3*t6
      taux20 = (256.D0/3.D0*t7*s2l2+256.D0/3.D0*t11*t12)*s1l1+(-256.D0/3
     #.D0*t7*s2l1-256.D0/3.D0*t11*t18)*s1l2+256.D0/3.D0*t23*t12*s2l1-256
     #.D0/3.D0*t23*t18*s2l2

CCC Coefficients

      tcoff37 = taux19
      tcoff38 = taux19
      tcoff39 = taux20
      tcoff40 = taux20

CCC Jet functions

      tFjet37 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff38)
      tFjet39 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff40)

CCC Sum

      qgcnum11e0 = tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tFjet
     #40*tcoff40
      RETURN
      END

      double precision function qgcnum11e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z*l1
      t3 = 1-z
      t4 = dlog(t3)
      t6 = 3+4*t4
      t7 = (t1+1-l1)*t6
      t29 = z**2
      t32 = l1**2
      taux19 = 64.D0/3.D0*t7*s2l2*s1l1+64.D0/3.D0*t7*s2l1*s1l2+(-128.D0/
     #3.D0*t3*(-1+l1)*t6*s2l2+64.D0/3.D0*t7*(-s1l2-s2l2-zp))*s2l1+64.D0/
     #3.D0*t7*(-s1l1-s2l1-zp)*s2l2-64.D0/3.D0*z*(2*t29*l1-t29*t32+2*z-2*
     #t29-2*t1+2*z*t32-1-t32)*t6
      t1 = 3*z
      t2 = 3*l1
      t3 = z*l1
      t4 = 3*t3
      t7 = dlog(1-z)
      t9 = 3+4*t7
      t10 = (t1+1+t2-t4)*t9
      t14 = (-3+t2+6*z-t4)*t9
      t15 = -s1l2-s2l2-zp
      t21 = -s1l1-s2l1-zp
      t27 = (-3+t3+t1-l1)*t9
      taux20 = (256.D0/3.D0*t10*s2l2+256.D0/3.D0*t14*t15)*s1l1+(-256.D0/
     #3.D0*t10*s2l1-256.D0/3.D0*t14*t21)*s1l2+256.D0/3.D0*t27*t15*s2l1-2
     #56.D0/3.D0*t27*t21*s2l2

CCC Coefficients

      tcoff37 = taux19
      tcoff38 = taux19
      tcoff39 = taux20
      tcoff40 = taux20

CCC Jet functions

      tFjet37 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff38)
      tFjet39 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff40)

CCC Sum

      qgcnum11e1 = tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tFjet
     #40*tcoff40
      RETURN
      END

      double precision function qgcnum11e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z**2
      t2 = l1**2
      t12 = dlog(1-z)
      taux19 = -64.D0/3.D0*z*(t1*t2-2*z*t2+2*z*l1+1+t2-2*l1)*(3+4*t12)
      t1 = z*l1
      t2 = 2*z
      t5 = dlog(1-z)
      t7 = 3+4*t5
      t8 = (t1-3-l1-t2)*t7
      t12 = (1-l1-t2+t1)*t7
      t13 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t28 = (-3*t1+1+3*l1-t2)*t7
      taux20 = (256.D0/3.D0*t8*s2l2+512.D0/3.D0*t12*t13)*s1l1+(-256.D0/3
     #.D0*t8*s2l1-512.D0/3.D0*t12*t20)*s1l2+256.D0/3.D0*t28*t13*s2l1-256
     #.D0/3.D0*t28*t20*s2l2

CCC Coefficients

      tcoff37 = taux19
      tcoff38 = taux19
      tcoff39 = taux20
      tcoff40 = taux20

CCC Jet functions

      tFjet37 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff38)
      tFjet39 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff40)

CCC Sum

      qgcnum11e2 = tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tFjet
     #40*tcoff40
      RETURN
      END

      double precision function qgcnum11e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux19 = 0
      t4 = dlog(1-z)
      t7 = (z*l1+1-l1)*(3+4*t4)
      taux20 = 512.D0/3.D0*t7*s2l2*s1l1-512.D0/3.D0*t7*s2l1*s1l2+512.D0/
     #3.D0*t7*(-s1l2-s2l2-zp)*s2l1-512.D0/3.D0*t7*(-s1l1-s2l1-zp)*s2l2

CCC Coefficients

      tcoff37 = taux19
      tcoff38 = taux19
      tcoff39 = taux20
      tcoff40 = taux20

CCC Jet functions

      tFjet37 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff38)
      tFjet39 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff40)

CCC Sum

      qgcnum11e3 = tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tFjet
     #40*tcoff40
      RETURN
      END

      double precision function qgcnum12e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t3 = zp-1
      t7 = zp*l1
      t9 = -64*l1-64*zp+64*t7
      t14 = 64-64*l1-128*zp+64*t7
      t15 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = -64*t3
      taux21 = (-2.D0/3.D0-2.D0/3.D0*yp)*((128*t3*l1*s1l2+t9*s2l2+t14*t1
     #5)*s1l1+(t9*s2l1+t14*t20)*s1l2+t25*t15*s2l1+t25*t20*s2l2+128*t20*t
     #15)
      t3 = zp*l1
      t4 = -zp+t3-l1
      t7 = 1-l1-2*zp+t3
      t8 = -s1l2-s2l2-zp
      t14 = -s1l1-s2l1-zp
      t19 = -zp+1
      taux22 = (-2.D0/3.D0-2.D0/3.D0*yp)*((256*t4*s2l2+256*t7*t8)*s1l1+(
     #-256*t4*s2l1-256*t7*t14)*s1l2+256*t19*t8*s2l1-256*t19*t14*s2l2)

CCC Coefficients

      tcoff41 = taux21
      tcoff42 = taux21
      tcoff43 = taux22
      tcoff44 = taux22

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f41)
      tFjet42 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off42)
      tFjet43 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f43)
      tFjet44 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off44)

CCC Sum

      qgcnum12e0 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFjet
     #44*tcoff44
      RETURN
      END

      double precision function qgcnum12e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t3 = zp*l1
      t5 = 64-64*l1+64*t3
      t22 = zp**2
      t25 = l1**2
      taux21 = (-2.D0/3.D0-2.D0/3.D0*yp)*(t5*s2l2*s1l1+t5*s2l1*s1l2+(128
     #*(zp-1)*(-1+l1)*s2l2+t5*(-s1l2-s2l2-zp))*s2l1+t5*(-s1l1-s2l1-zp)*s
     #2l2-64*zp*(2*t22*l1-t22*t25+2*zp-2*t22-2*t3+2*zp*t25-1-t25))
      t3 = 3*zp
      t4 = 3*l1
      t5 = zp*l1
      t6 = 3*t5
      t7 = t3+1+t4-t6
      t10 = -3+t4+6*zp-t6
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = -3+t5+t3-l1
      taux22 = (-2.D0/3.D0-2.D0/3.D0*yp)*((256*t7*s2l2+256*t10*t11)*s1l1
     #+(-256*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2
     #)

CCC Coefficients

      tcoff41 = taux21
      tcoff42 = taux21
      tcoff43 = taux22
      tcoff44 = taux22

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f41)
      tFjet42 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off42)
      tFjet43 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f43)
      tFjet44 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off44)

CCC Sum

      qgcnum12e1 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFjet
     #44*tcoff44
      RETURN
      END

      double precision function qgcnum12e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t4 = zp**2
      t5 = l1**2
      taux21 = -64*(-2.D0/3.D0-2.D0/3.D0*yp)*zp*(t4*t5-2*zp*t5+2*zp*l1+1
     #+t5-2*l1)
      t3 = zp*l1
      t4 = 2*zp
      t5 = t3-3-l1-t4
      t8 = 1-l1-t4+t3
      t9 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t23 = -3*t3+1+3*l1-t4
      taux22 = (-2.D0/3.D0-2.D0/3.D0*yp)*((256*t5*s2l2+512*t8*t9)*s1l1+(
     #-256*t5*s2l1-512*t8*t16)*s1l2+256*t23*t9*s2l1-256*t23*t16*s2l2)

CCC Coefficients

      tcoff41 = taux21
      tcoff42 = taux21
      tcoff43 = taux22
      tcoff44 = taux22

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f41)
      tFjet42 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off42)
      tFjet43 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f43)
      tFjet44 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off44)

CCC Sum

      qgcnum12e2 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFjet
     #44*tcoff44
      RETURN
      END

      double precision function qgcnum12e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux21 = 0
      t4 = zp*l1+1-l1
      taux22 = (-2.D0/3.D0-2.D0/3.D0*yp)*(512*t4*s2l2*s1l1-512*t4*s2l1*s
     #1l2+512*t4*(-s1l2-s2l2-zp)*s2l1-512*t4*(-s1l1-s2l1-zp)*s2l2)

CCC Coefficients

      tcoff41 = taux21
      tcoff42 = taux21
      tcoff43 = taux22
      tcoff44 = taux22

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f41)
      tFjet42 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off42)
      tFjet43 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f43)
      tFjet44 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off44)

CCC Sum

      qgcnum12e3 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFjet
     #44*tcoff44
      RETURN
      END

      double precision function qgcnum13e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = zp-1
      t5 = zp*l1
      t7 = -64*l1-64*zp+64*t5
      t12 = 64-64*l1-128*zp+64*t5
      t13 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t25 = -64*t1
      taux23 = 4.D0/3.D0*(128*t1*l1*s1l2+t7*s2l2+t12*t13)*s1l1+4.D0/3.D0
     #*(t7*s2l1+t12*t19)*s1l2+4.D0/3.D0*t25*t13*s2l1+4.D0/3.D0*t25*t19*s
     #2l2+512.D0/3.D0*t19*t13
      t1 = zp*l1
      t2 = -zp+t1-l1
      t5 = 1-l1-2*zp+t1
      t6 = -s1l2-s2l2-zp
      t13 = -s1l1-s2l1-zp
      t19 = -zp+1
      taux24 = 4.D0/3.D0*(256*t2*s2l2+256*t5*t6)*s1l1+4.D0/3.D0*(-256*t2
     #*s2l1-256*t5*t13)*s1l2+1024.D0/3.D0*t19*t6*s2l1-1024.D0/3.D0*t19*t
     #13*s2l2

CCC Coefficients

      tcoff45 = taux23
      tcoff46 = taux23
      tcoff47 = taux24
      tcoff48 = taux24

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f45)
      tFjet46 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off46)
      tFjet47 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f47)
      tFjet48 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off48)

CCC Sum

      qgcnum13e0 = tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tFjet
     #48*tcoff48
      RETURN
      END

      double precision function qgcnum13e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = zp*l1
      t3 = 64-64*l1+64*t1
      t24 = zp**2
      t27 = l1**2
      taux23 = 4.D0/3.D0*t3*s2l2*s1l1+4.D0/3.D0*t3*s2l1*s1l2+4.D0/3.D0*(
     #128*(zp-1)*(-1+l1)*s2l2+t3*(-s1l2-s2l2-zp))*s2l1+4.D0/3.D0*t3*(-s1
     #l1-s2l1-zp)*s2l2-256.D0/3.D0*zp*(2*t24*l1-t24*t27+2*zp-2*t24-2*t1+
     #2*zp*t27-1-t27)
      t1 = 3*zp
      t2 = 3*l1
      t3 = zp*l1
      t4 = 3*t3
      t5 = t1+1+t2-t4
      t8 = -3+t2+6*zp-t4
      t9 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t22 = -3+t3+t1-l1
      taux24 = 4.D0/3.D0*(256*t5*s2l2+256*t8*t9)*s1l1+4.D0/3.D0*(-256*t5
     #*s2l1-256*t8*t16)*s1l2+1024.D0/3.D0*t22*t9*s2l1-1024.D0/3.D0*t22*t
     #16*s2l2

CCC Coefficients

      tcoff45 = taux23
      tcoff46 = taux23
      tcoff47 = taux24
      tcoff48 = taux24

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f45)
      tFjet46 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off46)
      tFjet47 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f47)
      tFjet48 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off48)

CCC Sum

      qgcnum13e1 = tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tFjet
     #48*tcoff48
      RETURN
      END

      double precision function qgcnum13e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = zp**2
      t2 = l1**2
      taux23 = -256.D0/3.D0*zp*(t1*t2-2*zp*t2+2*zp*l1+1+t2-2*l1)
      t1 = zp*l1
      t2 = 2*zp
      t3 = t1-3-l1-t2
      t6 = 1-l1-t2+t1
      t7 = -s1l2-s2l2-zp
      t15 = -s1l1-s2l1-zp
      t23 = -3*t1+1+3*l1-t2
      taux24 = 4.D0/3.D0*(256*t3*s2l2+512*t6*t7)*s1l1+4.D0/3.D0*(-256*t3
     #*s2l1-512*t6*t15)*s1l2+1024.D0/3.D0*t23*t7*s2l1-1024.D0/3.D0*t23*t
     #15*s2l2

CCC Coefficients

      tcoff45 = taux23
      tcoff46 = taux23
      tcoff47 = taux24
      tcoff48 = taux24

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f45)
      tFjet46 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off46)
      tFjet47 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f47)
      tFjet48 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off48)

CCC Sum

      qgcnum13e2 = tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tFjet
     #48*tcoff48
      RETURN
      END

      double precision function qgcnum13e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux23 = 0
      t2 = zp*l1+1-l1
      taux24 = 2048.D0/3.D0*t2*s2l2*s1l1-2048.D0/3.D0*t2*s2l1*s1l2+2048.
     #D0/3.D0*t2*(-s1l2-s2l2-zp)*s2l1-2048.D0/3.D0*t2*(-s1l1-s2l1-zp)*s2
     #l2

CCC Coefficients

      tcoff45 = taux23
      tcoff46 = taux23
      tcoff47 = taux24
      tcoff48 = taux24

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f45)
      tFjet46 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off46)
      tFjet47 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f47)
      tFjet48 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off48)

CCC Sum

      qgcnum13e3 = tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tFjet
     #48*tcoff48
      RETURN
      END

      double precision function qgcnum14e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z-1
      t4 = -t1
      t5 = dlog(t4)
      t7 = -33+2*NF-36*t5
      t11 = z*l1
      t13 = (-l1-z+t11)*t7
      t18 = (1-l1-2*z+t11)*t7
      t19 = -s1l2-s2l2-zp
      t25 = -s1l1-s2l1-zp
      t30 = t4*t7
      taux25 = (-32.D0/3.D0*t1*l1*t7*s1l2-16.D0/3.D0*t13*s2l2-16.D0/3.D0
     #*t18*t19)*s1l1+(-16.D0/3.D0*t13*s2l1-16.D0/3.D0*t18*t25)*s1l2-16.D
     #0/3.D0*t30*t19*s2l1-16.D0/3.D0*t30*t25*s2l2+(352-64.D0/3.D0*NF+384
     #*t5)*t19*t25
      t1 = z*l1
      t4 = 1-z
      t5 = dlog(t4)
      t7 = -33+2*NF-36*t5
      t8 = (-l1-z+t1)*t7
      t12 = (1-l1-2*z+t1)*t7
      t13 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t24 = t4*t7
      taux26 = (-64.D0/3.D0*t8*s2l2-64.D0/3.D0*t12*t13)*s1l1+(64.D0/3.D0
     #*t8*s2l1+64.D0/3.D0*t12*t19)*s1l2-64.D0/3.D0*t24*t13*s2l1+64.D0/3.
     #D0*t24*t19*s2l2

CCC Coefficients

      tcoff49 = taux25
      tcoff50 = taux25
      tcoff51 = taux26
      tcoff52 = taux26

CCC Jet functions

      tFjet49 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff50)
      tFjet51 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff52)

CCC Sum

      qgcnum14e0 = tFjet49*tcoff49+tFjet50*tcoff50+tFjet51*tcoff51+tFjet
     #52*tcoff52
      RETURN
      END

      double precision function qgcnum14e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z*l1
      t4 = 1-z
      t5 = dlog(t4)
      t7 = -33+2*NF-36*t5
      t8 = (t1+1-l1)*t7
      t30 = z**2
      t33 = l1**2
      taux25 = -16.D0/3.D0*t8*s2l2*s1l1-16.D0/3.D0*t8*s2l1*s1l2+(32.D0/3
     #.D0*t4*(-1+l1)*t7*s2l2-16.D0/3.D0*t8*(-s1l2-s2l2-zp))*s2l1-16.D0/3
     #.D0*t8*(-s1l1-s2l1-zp)*s2l2+16.D0/3.D0*z*(2*t30*l1-t30*t33+2*z-2*t
     #30-2*t1+2*z*t33-1-t33)*t7
      t1 = 3*z
      t2 = 3*l1
      t3 = z*l1
      t4 = 3*t3
      t8 = dlog(1-z)
      t10 = -33+2*NF-36*t8
      t11 = (t1+1+t2-t4)*t10
      t15 = (-3+t2+6*z-t4)*t10
      t16 = -s1l2-s2l2-zp
      t22 = -s1l1-s2l1-zp
      t28 = (-3+t3+t1-l1)*t10
      taux26 = (-64.D0/3.D0*t11*s2l2-64.D0/3.D0*t15*t16)*s1l1+(64.D0/3.D
     #0*t11*s2l1+64.D0/3.D0*t15*t22)*s1l2-64.D0/3.D0*t28*t16*s2l1+64.D0/
     #3.D0*t28*t22*s2l2

CCC Coefficients

      tcoff49 = taux25
      tcoff50 = taux25
      tcoff51 = taux26
      tcoff52 = taux26

CCC Jet functions

      tFjet49 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff50)
      tFjet51 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff52)

CCC Sum

      qgcnum14e1 = tFjet49*tcoff49+tFjet50*tcoff50+tFjet51*tcoff51+tFjet
     #52*tcoff52
      RETURN
      END

      double precision function qgcnum14e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = z**2
      t2 = l1**2
      t13 = dlog(1-z)
      taux25 = 16.D0/3.D0*z*(t1*t2+1+t2+2*z*l1-2*z*t2-2*l1)*(-33+2*NF-36
     #*t13)
      t1 = z*l1
      t2 = 2*z
      t6 = dlog(1-z)
      t8 = -33+2*NF-36*t6
      t9 = (t1-3-l1-t2)*t8
      t13 = (1-l1-t2+t1)*t8
      t14 = -s1l2-s2l2-zp
      t21 = -s1l1-s2l1-zp
      t29 = (-3*t1+1+3*l1-t2)*t8
      taux26 = (-64.D0/3.D0*t9*s2l2-128.D0/3.D0*t13*t14)*s1l1+(64.D0/3.D
     #0*t9*s2l1+128.D0/3.D0*t13*t21)*s1l2-64.D0/3.D0*t29*t14*s2l1+64.D0/
     #3.D0*t29*t21*s2l2

CCC Coefficients

      tcoff49 = taux25
      tcoff50 = taux25
      tcoff51 = taux26
      tcoff52 = taux26

CCC Jet functions

      tFjet49 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff50)
      tFjet51 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff52)

CCC Sum

      qgcnum14e2 = tFjet49*tcoff49+tFjet50*tcoff50+tFjet51*tcoff51+tFjet
     #52*tcoff52
      RETURN
      END

      double precision function qgcnum14e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux25 = 0
      t5 = dlog(1-z)
      t8 = (z*l1+1-l1)*(-33+2*NF-36*t5)
      taux26 = -128.D0/3.D0*t8*s2l2*s1l1+128.D0/3.D0*t8*s2l1*s1l2-128.D0
     #/3.D0*t8*(-s1l2-s2l2-zp)*s2l1+128.D0/3.D0*t8*(-s1l1-s2l1-zp)*s2l2

CCC Coefficients

      tcoff49 = taux25
      tcoff50 = taux25
      tcoff51 = taux26
      tcoff52 = taux26

CCC Jet functions

      tFjet49 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff50)
      tFjet51 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff52)

CCC Sum

      qgcnum14e3 = tFjet49*tcoff49+tFjet50*tcoff50+tFjet51*tcoff51+tFjet
     #52*tcoff52
      RETURN
      END

      double precision function qgcnum15e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = 1-yp
      t7 = zp-1
      t11 = zp*l1
      t13 = -64*l1-64*zp+64*t11
      t18 = 64-64*l1-128*zp+64*t11
      t19 = -s1l2-s2l2-zp
      t24 = -s1l1-s2l1-zp
      t29 = -64*t7
      taux27 = (3*yp*t1+3*t1/yp-3)*((128*t7*l1*s1l2+t13*s2l2+t18*t19)*s1
     #l1+(t13*s2l1+t18*t24)*s1l2+t29*t19*s2l1+t29*t24*s2l2+128*t24*t19)
      t1 = 1-yp
      t7 = zp*l1
      t8 = -zp+t7-l1
      t11 = 1-l1-2*zp+t7
      t12 = -s1l2-s2l2-zp
      t18 = -s1l1-s2l1-zp
      t23 = -zp+1
      taux28 = (3*yp*t1+3*t1/yp-3)*((256*t8*s2l2+256*t11*t12)*s1l1+(-256
     #*t8*s2l1-256*t11*t18)*s1l2+256*t23*t12*s2l1-256*t23*t18*s2l2)

CCC Coefficients

      tcoff53 = taux27
      tcoff54 = taux27
      tcoff55 = taux28
      tcoff56 = taux28

CCC Jet functions

      tFjet53 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f53)
      tFjet54 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off54)
      tFjet55 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f55)
      tFjet56 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off56)

CCC Sum

      qgcnum15e0 = tFjet53*tcoff53+tFjet54*tcoff54+tFjet55*tcoff55+tFjet
     #56*tcoff56
      RETURN
      END

      double precision function qgcnum15e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = 1-yp
      t7 = zp*l1
      t9 = 64-64*l1+64*t7
      t26 = zp**2
      t29 = l1**2
      taux27 = (3*yp*t1+3*t1/yp-3)*(t9*s2l2*s1l1+t9*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t9*(-s1l2-s2l2-zp))*s2l1+t9*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(2*t26*l1-t26*t29+2*zp-2*t26-2*t7+2*zp*t29-1-t29))
      t1 = 1-yp
      t7 = 3*zp
      t8 = 3*l1
      t9 = zp*l1
      t10 = 3*t9
      t11 = t7+1+t8-t10
      t14 = -3+t8+6*zp-t10
      t15 = -s1l2-s2l2-zp
      t21 = -s1l1-s2l1-zp
      t26 = -3+t9+t7-l1
      taux28 = (3*yp*t1+3*t1/yp-3)*((256*t11*s2l2+256*t14*t15)*s1l1+(-25
     #6*t11*s2l1-256*t14*t21)*s1l2+256*t26*t15*s2l1-256*t26*t21*s2l2)

CCC Coefficients

      tcoff53 = taux27
      tcoff54 = taux27
      tcoff55 = taux28
      tcoff56 = taux28

CCC Jet functions

      tFjet53 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f53)
      tFjet54 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off54)
      tFjet55 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f55)
      tFjet56 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off56)

CCC Sum

      qgcnum15e1 = tFjet53*tcoff53+tFjet54*tcoff54+tFjet55*tcoff55+tFjet
     #56*tcoff56
      RETURN
      END

      double precision function qgcnum15e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = 1-yp
      t8 = zp**2
      t9 = l1**2
      taux27 = -64*(3*yp*t1+3*t1/yp-3)*zp*(t8*t9+1+t9+2*zp*l1-2*zp*t9-2*
     #l1)
      t1 = 1-yp
      t7 = zp*l1
      t8 = 2*zp
      t9 = t7-3-l1-t8
      t12 = 1-l1-t8+t7
      t13 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t27 = -3*t7+1+3*l1-t8
      taux28 = (3*yp*t1+3*t1/yp-3)*((256*t9*s2l2+512*t12*t13)*s1l1+(-256
     #*t9*s2l1-512*t12*t20)*s1l2+256*t27*t13*s2l1-256*t27*t20*s2l2)

CCC Coefficients

      tcoff53 = taux27
      tcoff54 = taux27
      tcoff55 = taux28
      tcoff56 = taux28

CCC Jet functions

      tFjet53 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f53)
      tFjet54 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off54)
      tFjet55 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f55)
      tFjet56 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off56)

CCC Sum

      qgcnum15e2 = tFjet53*tcoff53+tFjet54*tcoff54+tFjet55*tcoff55+tFjet
     #56*tcoff56
      RETURN
      END

      double precision function qgcnum15e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux27 = 0
      t1 = 1-yp
      t8 = 1-l1+zp*l1
      taux28 = (3*yp*t1+3*t1/yp-3)*(512*t8*s2l2*s1l1-512*t8*s2l1*s1l2+51
     #2*t8*(-s1l2-s2l2-zp)*s2l1-512*t8*(-s1l1-s2l1-zp)*s2l2)

CCC Coefficients

      tcoff53 = taux27
      tcoff54 = taux27
      tcoff55 = taux28
      tcoff56 = taux28

CCC Jet functions

      tFjet53 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f53)
      tFjet54 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off54)
      tFjet55 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f55)
      tFjet56 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off56)

CCC Sum

      qgcnum15e3 = tFjet53*tcoff53+tFjet54*tcoff54+tFjet55*tcoff55+tFjet
     #56*tcoff56
      RETURN
      END

      double precision function qgcnum16e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = zp-1
      t5 = zp*l1
      t7 = -64*l1-64*zp+64*t5
      t12 = 64-64*l1-128*zp+64*t5
      t13 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t25 = -64*t1
      taux29 = 3*(128*t1*l1*s1l2+t7*s2l2+t12*t13)*s1l1+3*(t7*s2l1+t12*t1
     #9)*s1l2+3*t25*t13*s2l1+3*t25*t19*s2l2+384*t19*t13
      t1 = zp*l1
      t2 = -zp+t1-l1
      t5 = 1-l1-2*zp+t1
      t6 = -s1l2-s2l2-zp
      t13 = -s1l1-s2l1-zp
      t19 = -zp+1
      taux30 = 3*(256*t2*s2l2+256*t5*t6)*s1l1+3*(-256*t2*s2l1-256*t5*t13
     #)*s1l2+768*t19*t6*s2l1-768*t19*t13*s2l2

CCC Coefficients

      tcoff57 = taux29
      tcoff58 = taux29
      tcoff59 = taux30
      tcoff60 = taux30

CCC Jet functions

      tFjet57 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f57)
      tFjet58 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off58)
      tFjet59 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f59)
      tFjet60 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off60)

CCC Sum

      qgcnum16e0 = tFjet57*tcoff57+tFjet58*tcoff58+tFjet59*tcoff59+tFjet
     #60*tcoff60
      RETURN
      END

      double precision function qgcnum16e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = zp*l1
      t3 = 64-64*l1+64*t1
      t24 = zp**2
      t27 = l1**2
      taux29 = 3*t3*s2l2*s1l1+3*t3*s2l1*s1l2+3*(128*(zp-1)*(-1+l1)*s2l2+
     #t3*(-s1l2-s2l2-zp))*s2l1+3*t3*(-s1l1-s2l1-zp)*s2l2-192*zp*(2*t24*l
     #1-t24*t27+2*zp-2*t24-2*t1+2*zp*t27-1-t27)
      t1 = 3*zp
      t2 = 3*l1
      t3 = zp*l1
      t4 = 3*t3
      t5 = t1+1+t2-t4
      t8 = -3+t2+6*zp-t4
      t9 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t22 = -3+t3+t1-l1
      taux30 = 3*(256*t5*s2l2+256*t8*t9)*s1l1+3*(-256*t5*s2l1-256*t8*t16
     #)*s1l2+768*t22*t9*s2l1-768*t22*t16*s2l2

CCC Coefficients

      tcoff57 = taux29
      tcoff58 = taux29
      tcoff59 = taux30
      tcoff60 = taux30

CCC Jet functions

      tFjet57 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f57)
      tFjet58 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off58)
      tFjet59 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f59)
      tFjet60 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off60)

CCC Sum

      qgcnum16e1 = tFjet57*tcoff57+tFjet58*tcoff58+tFjet59*tcoff59+tFjet
     #60*tcoff60
      RETURN
      END

      double precision function qgcnum16e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = zp**2
      t2 = l1**2
      taux29 = -192*zp*(t1*t2+1+t2+2*zp*l1-2*zp*t2-2*l1)
      t1 = zp*l1
      t2 = 2*zp
      t3 = t1-3-l1-t2
      t6 = 1-l1-t2+t1
      t7 = -s1l2-s2l2-zp
      t15 = -s1l1-s2l1-zp
      t23 = -3*t1+1+3*l1-t2
      taux30 = 3*(256*t3*s2l2+512*t6*t7)*s1l1+3*(-256*t3*s2l1-512*t6*t15
     #)*s1l2+768*t23*t7*s2l1-768*t23*t15*s2l2

CCC Coefficients

      tcoff57 = taux29
      tcoff58 = taux29
      tcoff59 = taux30
      tcoff60 = taux30

CCC Jet functions

      tFjet57 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f57)
      tFjet58 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off58)
      tFjet59 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f59)
      tFjet60 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off60)

CCC Sum

      qgcnum16e2 = tFjet57*tcoff57+tFjet58*tcoff58+tFjet59*tcoff59+tFjet
     #60*tcoff60
      RETURN
      END

      double precision function qgcnum16e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux29 = 0
      t2 = 1-l1+zp*l1
      taux30 = 1536*t2*s2l2*s1l1-1536*t2*s2l1*s1l2+1536*t2*(-s1l2-s2l2-z
     #p)*s2l1-1536*t2*(-s1l1-s2l1-zp)*s2l2

CCC Coefficients

      tcoff57 = taux29
      tcoff58 = taux29
      tcoff59 = taux30
      tcoff60 = taux30

CCC Jet functions

      tFjet57 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f57)
      tFjet58 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off58)
      tFjet59 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f59)
      tFjet60 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off60)

CCC Sum

      qgcnum16e3 = tFjet57*tcoff57+tFjet58*tcoff58+tFjet59*tcoff59+tFjet
     #60*tcoff60
      RETURN
      END

      double precision function qgcnum17e0(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = yp**2
      t3 = (1-yp)**2
      t6 = zp-1
      t11 = 64+64*zp
      t13 = zp*l1
      t15 = 64*l1+64*zp-64*t13
      t16 = -s1l2-s2l2-zp
      t21 = -s1l1-s2l1-zp
      t30 = -64*l1+64*t13+64
      taux31 = (t1/4+t3/4)*((-128*t6*l1*s1l2+t11*s2l2+t15*t16)*s1l1+(t11
     #*s2l1+t15*t21)*s1l2+(128*t6*(-1+l1)*s2l2+t30*t16)*s2l1+t30*t21*s2l
     #2)
      t1 = yp**2
      t3 = (1-yp)**2
      t6 = 1+zp
      t8 = zp*l1
      t9 = t8-l1-zp
      t10 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t21 = -1+l1-t8
      taux32 = (t1/4+t3/4)*((256*t6*s2l2-256*t9*t10)*s1l1+(-256*t6*s2l1+
     #256*t9*t16)*s1l2+256*t21*t10*s2l1-256*t21*t16*s2l2)

CCC Coefficients

      tcoff61 = taux31
      tcoff62 = taux31
      tcoff63 = taux32
      tcoff64 = taux32

CCC Jet functions

      tFjet61 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f61)
      tFjet62 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off62)
      tFjet63 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f63)
      tFjet64 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off64)

CCC Sum

      qgcnum17e0 = tFjet61*tcoff61+tFjet62*tcoff62+tFjet63*tcoff63+tFjet
     #64*tcoff64
      RETURN
      END

      double precision function qgcnum17e1(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = yp**2
      t3 = (1-yp)**2
      t7 = 64-64*zp
      t8 = -s1l2-s2l2-zp
      t9 = t7*t8
      t11 = -s1l1-s2l1-zp
      t12 = t7*t11
      t18 = zp**2
      t19 = l1**2
      taux31 = (t1/4+t3/4)*(t9*s1l1+t12*s1l2+t9*s2l1+t12*s2l2+128*t11*t8
     #+64*zp*(-2*t18*t19-1-t18+2*t18*l1-4*zp*l1+4*zp*t19-2*t19+2*l1))
      t1 = yp**2
      t3 = (1-yp)**2
      t7 = -3-3*zp
      t10 = 4*zp*l1
      t12 = 4*l1
      t13 = -1-t10+4*zp+t12
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = t10+3-t12
      taux32 = (t1/4+t3/4)*((256*t7*s2l2-256*t13*t14)*s1l1+(-256*t7*s2l1
     #+256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)

CCC Coefficients

      tcoff61 = taux31
      tcoff62 = taux31
      tcoff63 = taux32
      tcoff64 = taux32

CCC Jet functions

      tFjet61 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f61)
      tFjet62 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off62)
      tFjet63 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f63)
      tFjet64 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off64)

CCC Sum

      qgcnum17e1 = tFjet61*tcoff61+tFjet62*tcoff62+tFjet63*tcoff63+tFjet
     #64*tcoff64
      RETURN
      END

      double precision function qgcnum17e2(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      t1 = yp**2
      t3 = (1-yp)**2
      t7 = zp**2
      taux31 = 64*(t1/4+t3/4)*zp*(t7+1-2*zp)
      t1 = yp**2
      t3 = (1-yp)**2
      t6 = 1+zp
      t11 = 4*zp*l1
      t12 = 4*l1
      t13 = -5*zp+3+t11-t12
      t14 = -s1l2-s2l2-zp
      t21 = -s1l1-s2l1-zp
      t26 = -zp-1-t11+t12
      taux32 = (t1/4+t3/4)*((512*t6*s2l2-256*t13*t14)*s1l1+(-512*t6*s2l1
     #+256*t13*t21)*s1l2+256*t26*t14*s2l1-256*t26*t21*s2l2)

CCC Coefficients

      tcoff61 = taux31
      tcoff62 = taux31
      tcoff63 = taux32
      tcoff64 = taux32

CCC Jet functions

      tFjet61 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f61)
      tFjet62 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off62)
      tFjet63 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f63)
      tFjet64 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off64)

CCC Sum

      qgcnum17e2 = tFjet61*tcoff61+tFjet62*tcoff62+tFjet63*tcoff63+tFjet
     #64*tcoff64
      RETURN
      END

      double precision function qgcnum17e3(var,ndim,z,ux,l1,
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
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
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

      taux31 = 0
      t1 = yp**2
      t3 = (1-yp)**2
      t6 = zp-1
      t8 = t6*(-s1l2-s2l2-zp)
      t11 = t6*(-s1l1-s2l1-zp)
      taux32 = (t1/4+t3/4)*(-512*t8*s1l1+512*t11*s1l2+512*t8*s2l1-512*t1
     #1*s2l2)

CCC Coefficients

      tcoff61 = taux31
      tcoff62 = taux31
      tcoff63 = taux32
      tcoff64 = taux32

CCC Jet functions

      tFjet61 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f61)
      tFjet62 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off62)
      tFjet63 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f63)
      tFjet64 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off64)

CCC Sum

      qgcnum17e3 = tFjet61*tcoff61+tFjet62*tcoff62+tFjet63*tcoff63+tFjet
     #64*tcoff64
      RETURN
      END

