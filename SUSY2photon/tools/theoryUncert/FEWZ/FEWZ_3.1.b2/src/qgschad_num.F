CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgscnum1e0(var,ndim,z,ux,l1,
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

      qgscnum1e0 = tFjet1*tcoff1+tFjet2*tcoff2
      RETURN
      END

      double precision function qgscnum1e1(var,ndim,z,ux,l1,
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

      qgscnum1e1 = tFjet1*tcoff1+tFjet2*tcoff2
      RETURN
      END

      double precision function qgscnum1e2(var,ndim,z,ux,l1,
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

      qgscnum1e2 = tFjet1*tcoff1+tFjet2*tcoff2
      RETURN
      END

      double precision function qgscnum1e3(var,ndim,z,ux,l1,
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

      qgscnum1e3 = tFjet1*tcoff1+tFjet2*tcoff2
      RETURN
      END

      double precision function qgscnum2e0(var,ndim,z,ux,l1,
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
      taux2 = -18*z*(1+4*z)*t6*t7-18*z*(2*t11-2*z+1)*t6*t17+t6*(120*t21+
     #2*NF*z+4*NF*t21-24-51*z-4*NF*t11-78*t11)/2
      t6 = -s1l2*s2l1+s1l1*s2l2
      t7 = dlog(z)
      t11 = z**2
      t17 = dlog(1-z)
      t21 = t11*z
      taux3 = 72*z*(1+4*z)*t6*t7+72*z*(2*t11-2*z+1)*t6*t17-2*t6*(120*t21
     #+2*NF*z+4*NF*t21-24-51*z-4*NF*t11-78*t11)

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

      qgscnum2e0 = tFjet3*tcoff3+tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tcof
     #f6
      RETURN
      END

      double precision function qgscnum2e1(var,ndim,z,ux,l1,
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
      taux2 = -18*z*(1+4*z)*t4-18*z*(2*t7-2*z+1)*t13+60*t16+NF*z+2*NF*t1
     #6-12-51.D0/2.D0*z-2*NF*t7-39*t7
      t6 = -s1l2*s2l1+s1l1*s2l2
      t7 = dlog(z)
      t11 = z**2
      t17 = dlog(1-z)
      t23 = t11*z
      taux3 = -216*z*(1+4*z)*t6*t7-216*z*(2*t11-2*z+1)*t6*t17+6*t6*(120*
     #t23+2*NF*z+4*NF*t23-24-51*z-4*NF*t11-78*t11)

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

      qgscnum2e1 = tFjet3*tcoff3+tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tcof
     #f6
      RETURN
      END

      double precision function qgscnum2e2(var,ndim,z,ux,l1,
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
      t6 = -s1l2*s2l1+s1l1*s2l2
      t7 = dlog(z)
      t11 = z**2
      t17 = dlog(1-z)
      t21 = t11*z
      taux3 = 144*z*(1+4*z)*t6*t7+144*z*(2*t11-2*z+1)*t6*t17-4*t6*(120*t
     #21+2*NF*z+4*NF*t21-24-51*z-4*NF*t11-78*t11)

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

      qgscnum2e2 = tFjet3*tcoff3+tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tcof
     #f6
      RETURN
      END

      double precision function qgscnum2e3(var,ndim,z,ux,l1,
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

      qgscnum2e3 = tFjet3*tcoff3+tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tcof
     #f6
      RETURN
      END

      double precision function qgscnum3e0(var,ndim,z,ux,l1,
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
      taux4 = 4*z*(1-t1+4*t2)*t8*t9-8*z*(2*t2-t1+1)*t8*t17-2*z*(4*z-1)*t
     #8
      t1 = 2*z
      t2 = z**2
      t8 = -s1l2*s2l1+s1l1*s2l2
      t9 = dlog(z)
      t17 = dlog(1-z)
      taux5 = -16*z*(1-t1+4*t2)*t8*t9+32*z*(2*t2-t1+1)*t8*t17+8*z*(4*z-1
     #)*t8

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

      qgscnum3e0 = tFjet7*tcoff7+tFjet8*tcoff8+tFjet9*tcoff9+tFjet10*tco
     #ff10
      RETURN
      END

      double precision function qgscnum3e1(var,ndim,z,ux,l1,
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
      taux4 = 4*z*(1-t1+4*t2)*t6-8*z*(2*t2-t1+1)*t13-2*z*(4*z-1)
      t1 = 2*z
      t2 = z**2
      t8 = -s1l2*s2l1+s1l1*s2l2
      t9 = dlog(z)
      t17 = dlog(1-z)
      taux5 = 48*z*(1-t1+4*t2)*t8*t9-96*z*(2*t2-t1+1)*t8*t17-24*z*(4*z-1
     #)*t8

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

      qgscnum3e1 = tFjet7*tcoff7+tFjet8*tcoff8+tFjet9*tcoff9+tFjet10*tco
     #ff10
      RETURN
      END

      double precision function qgscnum3e2(var,ndim,z,ux,l1,
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
      t8 = -s1l2*s2l1+s1l1*s2l2
      t9 = dlog(z)
      t17 = dlog(1-z)
      taux5 = -32*z*(1-t1+4*t2)*t8*t9+64*z*(2*t2-t1+1)*t8*t17+16*z*(4*z-
     #1)*t8

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

      qgscnum3e2 = tFjet7*tcoff7+tFjet8*tcoff8+tFjet9*tcoff9+tFjet10*tco
     #ff10
      RETURN
      END

      double precision function qgscnum3e3(var,ndim,z,ux,l1,
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

      qgscnum3e3 = tFjet7*tcoff7+tFjet8*tcoff8+tFjet9*tcoff9+tFjet10*tco
     #ff10
      RETURN
      END

      double precision function qgscnum4e0(var,ndim,z,ux,l1,
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
      t5 = 12*t1*t3
      t7 = 6*t1*z
      t9 = 12*t1*t2
      t10 = s1l1*s2l2
      t12 = 12*t10*t3
      t14 = 6*t10*z
      t16 = 12*t10*t2
      t18 = dlog(1-z)
      taux6 = t5+t7-t9+t12+t14-t16+4.D0/3.D0*t18*(t5+t7-t9+t12+t14-t16)
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t5 = 48*t1*t3
      t7 = 24*t1*z
      t9 = 48*t1*t2
      t10 = s1l1*s2l2
      t12 = 48*t10*t3
      t14 = 24*t10*z
      t16 = 48*t10*t2
      t18 = dlog(1-z)
      taux7 = -t5-t7+t9+t12+t14-t16+4.D0/3.D0*t18*(-t5-t7+t9+t12+t14-t16
     #)

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

      qgscnum4e0 = tFjet11*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet
     #14*tcoff14
      RETURN
      END

      double precision function qgscnum4e1(var,ndim,z,ux,l1,
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
      taux6 = -t3-t4+t5+4.D0/3.D0*t7*(-t3-t4+t5)
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t5 = 144*t1*t3
      t7 = 72*t1*z
      t9 = 144*t1*t2
      t10 = s1l1*s2l2
      t12 = 144*t10*t3
      t14 = 72*t10*z
      t16 = 144*t10*t2
      t18 = dlog(1-z)
      taux7 = t5+t7-t9-t12-t14+t16+4.D0/3.D0*t18*(t5+t7-t9-t12-t14+t16)

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

      qgscnum4e1 = tFjet11*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet
     #14*tcoff14
      RETURN
      END

      double precision function qgscnum4e2(var,ndim,z,ux,l1,
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
      t5 = 96*t1*t3
      t7 = 48*t1*z
      t9 = 96*t1*t2
      t10 = s1l1*s2l2
      t12 = 96*t10*t3
      t14 = 48*t10*z
      t16 = 96*t10*t2
      t18 = dlog(1-z)
      taux7 = -t5-t7+t9+t12+t14-t16+4.D0/3.D0*t18*(-t5-t7+t9+t12+t14-t16
     #)

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

      qgscnum4e2 = tFjet11*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet
     #14*tcoff14
      RETURN
      END

      double precision function qgscnum4e3(var,ndim,z,ux,l1,
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

      qgscnum4e3 = tFjet11*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet
     #14*tcoff14
      RETURN
      END

      double precision function qgscnum5e0(var,ndim,z,ux,l1,
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
      taux8 = (24*s1l2*s2l1+24*s1l1*s2l2)*z*(t6/t10/4+t16/4)*t13
      t6 = z**2
      t9 = (1-z)*l2+z
      t10 = t9**2
      t13 = 1/t9
      t16 = (1-z*t13)**2
      taux9 = (-96*s1l2*s2l1+96*s1l1*s2l2)*z*(t6/t10/4+t16/4)*t13

CCC Coefficients

      tcoff15 = taux8
      tcoff16 = taux8
      tcoff17 = taux9
      tcoff18 = taux9

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,t3,z/t3,Tcoef
     #f*tcoff15)
      t4 = (1-z)*l2+z
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,t4,z/t4,Tco
     #eff*tcoff16)
      t3 = (1-z)*l2+z
      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,t3,z/t3,Tcoef
     #f*tcoff17)
      t4 = (1-z)*l2+z
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,t4,z/t4,Tco
     #eff*tcoff18)

CCC Sum

      qgscnum5e0 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet
     #18*tcoff18
      RETURN
      END

      double precision function qgscnum5e1(var,ndim,z,ux,l1,
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
      taux8 = -24*z*(t1/t5/4+t11/4)*t8
      t6 = z**2
      t9 = (1-z)*l2+z
      t10 = t9**2
      t13 = 1/t9
      t16 = (1-z*t13)**2
      taux9 = (288*s1l2*s2l1-288*s1l1*s2l2)*z*(t6/t10/4+t16/4)*t13

CCC Coefficients

      tcoff15 = taux8
      tcoff16 = taux8
      tcoff17 = taux9
      tcoff18 = taux9

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,t3,z/t3,Tcoef
     #f*tcoff15)
      t4 = (1-z)*l2+z
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,t4,z/t4,Tco
     #eff*tcoff16)
      t3 = (1-z)*l2+z
      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,t3,z/t3,Tcoef
     #f*tcoff17)
      t4 = (1-z)*l2+z
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,t4,z/t4,Tco
     #eff*tcoff18)

CCC Sum

      qgscnum5e1 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet
     #18*tcoff18
      RETURN
      END

      double precision function qgscnum5e2(var,ndim,z,ux,l1,
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
      t6 = z**2
      t9 = (1-z)*l2+z
      t10 = t9**2
      t13 = 1/t9
      t16 = (1-z*t13)**2
      taux9 = (-192*s1l2*s2l1+192*s1l1*s2l2)*z*(t6/t10/4+t16/4)*t13

CCC Coefficients

      tcoff15 = taux8
      tcoff16 = taux8
      tcoff17 = taux9
      tcoff18 = taux9

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,t3,z/t3,Tcoef
     #f*tcoff15)
      t4 = (1-z)*l2+z
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,t4,z/t4,Tco
     #eff*tcoff16)
      t3 = (1-z)*l2+z
      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,t3,z/t3,Tcoef
     #f*tcoff17)
      t4 = (1-z)*l2+z
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,t4,z/t4,Tco
     #eff*tcoff18)

CCC Sum

      qgscnum5e2 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet
     #18*tcoff18
      RETURN
      END

      double precision function qgscnum5e3(var,ndim,z,ux,l1,
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

      t3 = (1-z)*l2+z
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,t3,z/t3,Tcoef
     #f*tcoff15)
      t4 = (1-z)*l2+z
      tFjet16 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,t4,z/t4,Tco
     #eff*tcoff16)
      t3 = (1-z)*l2+z
      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,t3,z/t3,Tcoef
     #f*tcoff17)
      t4 = (1-z)*l2+z
      tFjet18 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,t4,z/t4,Tco
     #eff*tcoff18)

CCC Sum

      qgscnum5e3 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet
     #18*tcoff18
      RETURN
      END

      double precision function qgscnum6e0(var,ndim,z,ux,l1,
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
      taux10 = 12*t1*t3+6*t1*z-12*t1*t2+12*t10*t3+6*t10*z-12*t10*t2
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t10 = s1l1*s2l2
      taux11 = -48*t1*t3-24*t1*z+48*t1*t2+48*t10*t3+24*t10*z-48*t10*t2

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

      qgscnum6e0 = tFjet19*tcoff19+tFjet20*tcoff20+tFjet21*tcoff21+tFjet
     #22*tcoff22
      RETURN
      END

      double precision function qgscnum6e1(var,ndim,z,ux,l1,
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
      taux10 = -12*t1*z-6*z+12*t1
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t10 = s1l1*s2l2
      taux11 = 144*t1*t3+72*t1*z-144*t1*t2-144*t10*t3-72*t10*z+144*t10*t
     #2

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

      qgscnum6e1 = tFjet19*tcoff19+tFjet20*tcoff20+tFjet21*tcoff21+tFjet
     #22*tcoff22
      RETURN
      END

      double precision function qgscnum6e2(var,ndim,z,ux,l1,
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
      t3 = t2*z
      t10 = s1l1*s2l2
      taux11 = -96*t1*t3-48*t1*z+96*t1*t2+96*t10*t3+48*t10*z-96*t10*t2

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

      qgscnum6e2 = tFjet19*tcoff19+tFjet20*tcoff20+tFjet21*tcoff21+tFjet
     #22*tcoff22
      RETURN
      END

      double precision function qgscnum6e3(var,ndim,z,ux,l1,
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

      qgscnum6e3 = tFjet19*tcoff19+tFjet20*tcoff20+tFjet21*tcoff21+tFjet
     #22*tcoff22
      RETURN
      END

      double precision function qgscnum7e0(var,ndim,z,ux,l1,
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
      taux12 = -18*z*(1+4*z)*t6*t7-18*z*(2*t11-2*z+1)*t6*t17+t6*(120*t21
     #+2*NF*z+4*NF*t21-24-51*z-4*NF*t11-78*t11)/2
      t6 = -s1l2*s2l1+s1l1*s2l2
      t7 = dlog(z)
      t11 = z**2
      t17 = dlog(1-z)
      t21 = t11*z
      taux13 = 72*z*(1+4*z)*t6*t7+72*z*(2*t11-2*z+1)*t6*t17-2*t6*(120*t2
     #1+2*NF*z+4*NF*t21-24-51*z-4*NF*t11-78*t11)

CCC Coefficients

      tcoff23 = taux12
      tcoff24 = taux12
      tcoff25 = taux13
      tcoff26 = taux13

CCC Jet functions

      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff24)
      tFjet25 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff25)
      tFjet26 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff26)

CCC Sum

      qgscnum7e0 = tFjet23*tcoff23+tFjet24*tcoff24+tFjet25*tcoff25+tFjet
     #26*tcoff26
      RETURN
      END

      double precision function qgscnum7e1(var,ndim,z,ux,l1,
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
      taux12 = -18*z*(1+4*z)*t4-18*z*(2*t7-2*z+1)*t13+60*t16+NF*z+2*NF*t
     #16-12-51.D0/2.D0*z-2*NF*t7-39*t7
      t6 = -s1l2*s2l1+s1l1*s2l2
      t7 = dlog(z)
      t11 = z**2
      t17 = dlog(1-z)
      t23 = t11*z
      taux13 = -216*z*(1+4*z)*t6*t7-216*z*(2*t11-2*z+1)*t6*t17+6*t6*(120
     #*t23+2*NF*z+4*NF*t23-24-51*z-4*NF*t11-78*t11)

CCC Coefficients

      tcoff23 = taux12
      tcoff24 = taux12
      tcoff25 = taux13
      tcoff26 = taux13

CCC Jet functions

      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff24)
      tFjet25 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff25)
      tFjet26 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff26)

CCC Sum

      qgscnum7e1 = tFjet23*tcoff23+tFjet24*tcoff24+tFjet25*tcoff25+tFjet
     #26*tcoff26
      RETURN
      END

      double precision function qgscnum7e2(var,ndim,z,ux,l1,
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
      t6 = -s1l2*s2l1+s1l1*s2l2
      t7 = dlog(z)
      t11 = z**2
      t17 = dlog(1-z)
      t21 = t11*z
      taux13 = 144*z*(1+4*z)*t6*t7+144*z*(2*t11-2*z+1)*t6*t17-4*t6*(120*
     #t21+2*NF*z+4*NF*t21-24-51*z-4*NF*t11-78*t11)

CCC Coefficients

      tcoff23 = taux12
      tcoff24 = taux12
      tcoff25 = taux13
      tcoff26 = taux13

CCC Jet functions

      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff24)
      tFjet25 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff25)
      tFjet26 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff26)

CCC Sum

      qgscnum7e2 = tFjet23*tcoff23+tFjet24*tcoff24+tFjet25*tcoff25+tFjet
     #26*tcoff26
      RETURN
      END

      double precision function qgscnum7e3(var,ndim,z,ux,l1,
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

      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff24)
      tFjet25 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff25)
      tFjet26 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff26)

CCC Sum

      qgscnum7e3 = tFjet23*tcoff23+tFjet24*tcoff24+tFjet25*tcoff25+tFjet
     #26*tcoff26
      RETURN
      END

      double precision function qgscnum8e0(var,ndim,z,ux,l1,
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
      taux14 = 4*z*(1-t1+4*t2)*t8*t9-8*z*(2*t2-t1+1)*t8*t17-2*z*(4*z-1)*
     #t8
      t1 = 2*z
      t2 = z**2
      t8 = -s1l2*s2l1+s1l1*s2l2
      t9 = dlog(z)
      t17 = dlog(1-z)
      taux15 = -16*z*(1-t1+4*t2)*t8*t9+32*z*(2*t2-t1+1)*t8*t17+8*z*(4*z-
     #1)*t8

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

      qgscnum8e0 = tFjet27*tcoff27+tFjet28*tcoff28+tFjet29*tcoff29+tFjet
     #30*tcoff30
      RETURN
      END

      double precision function qgscnum8e1(var,ndim,z,ux,l1,
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
      taux14 = 4*z*(1-t1+4*t2)*t6-8*z*(2*t2-t1+1)*t13-2*z*(4*z-1)
      t1 = 2*z
      t2 = z**2
      t8 = -s1l2*s2l1+s1l1*s2l2
      t9 = dlog(z)
      t17 = dlog(1-z)
      taux15 = 48*z*(1-t1+4*t2)*t8*t9-96*z*(2*t2-t1+1)*t8*t17-24*z*(4*z-
     #1)*t8

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

      qgscnum8e1 = tFjet27*tcoff27+tFjet28*tcoff28+tFjet29*tcoff29+tFjet
     #30*tcoff30
      RETURN
      END

      double precision function qgscnum8e2(var,ndim,z,ux,l1,
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

      taux14 = 0
      t1 = 2*z
      t2 = z**2
      t8 = -s1l2*s2l1+s1l1*s2l2
      t9 = dlog(z)
      t17 = dlog(1-z)
      taux15 = -32*z*(1-t1+4*t2)*t8*t9+64*z*(2*t2-t1+1)*t8*t17+16*z*(4*z
     #-1)*t8

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

      qgscnum8e2 = tFjet27*tcoff27+tFjet28*tcoff28+tFjet29*tcoff29+tFjet
     #30*tcoff30
      RETURN
      END

      double precision function qgscnum8e3(var,ndim,z,ux,l1,
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

      taux14 = 0
      taux15 = 0

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

      qgscnum8e3 = tFjet27*tcoff27+tFjet28*tcoff28+tFjet29*tcoff29+tFjet
     #30*tcoff30
      RETURN
      END

      double precision function qgscnum9e0(var,ndim,z,ux,l1,
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
      taux16 = 12*t1*t3+6*t1*z-12*t1*t2+12*t10*t3+6*t10*z-12*t10*t2
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t10 = s1l1*s2l2
      taux17 = -48*t1*t3-24*t1*z+48*t1*t2+48*t10*t3+24*t10*z-48*t10*t2

CCC Coefficients

      tcoff31 = taux16
      tcoff32 = taux16
      tcoff33 = taux17
      tcoff34 = taux17

CCC Jet functions

      tFjet31 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff32)
      tFjet33 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff33)
      tFjet34 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff34)

CCC Sum

      qgscnum9e0 = tFjet31*tcoff31+tFjet32*tcoff32+tFjet33*tcoff33+tFjet
     #34*tcoff34
      RETURN
      END

      double precision function qgscnum9e1(var,ndim,z,ux,l1,
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
      taux16 = -12*t1*z-6*z+12*t1
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t10 = s1l1*s2l2
      taux17 = 144*t1*t3+72*t1*z-144*t1*t2-144*t10*t3-72*t10*z+144*t10*t
     #2

CCC Coefficients

      tcoff31 = taux16
      tcoff32 = taux16
      tcoff33 = taux17
      tcoff34 = taux17

CCC Jet functions

      tFjet31 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff32)
      tFjet33 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff33)
      tFjet34 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff34)

CCC Sum

      qgscnum9e1 = tFjet31*tcoff31+tFjet32*tcoff32+tFjet33*tcoff33+tFjet
     #34*tcoff34
      RETURN
      END

      double precision function qgscnum9e2(var,ndim,z,ux,l1,
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

      taux16 = 0
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t10 = s1l1*s2l2
      taux17 = -96*t1*t3-48*t1*z+96*t1*t2+96*t10*t3+48*t10*z-96*t10*t2

CCC Coefficients

      tcoff31 = taux16
      tcoff32 = taux16
      tcoff33 = taux17
      tcoff34 = taux17

CCC Jet functions

      tFjet31 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff32)
      tFjet33 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff33)
      tFjet34 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff34)

CCC Sum

      qgscnum9e2 = tFjet31*tcoff31+tFjet32*tcoff32+tFjet33*tcoff33+tFjet
     #34*tcoff34
      RETURN
      END

      double precision function qgscnum9e3(var,ndim,z,ux,l1,
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

      taux16 = 0
      taux17 = 0

CCC Coefficients

      tcoff31 = taux16
      tcoff32 = taux16
      tcoff33 = taux17
      tcoff34 = taux17

CCC Jet functions

      tFjet31 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff32)
      tFjet33 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff33)
      tFjet34 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff34)

CCC Sum

      qgscnum9e3 = tFjet31*tcoff31+tFjet32*tcoff32+tFjet33*tcoff33+tFjet
     #34*tcoff34
      RETURN
      END

      double precision function qgscnum10e0(var,ndim,z,ux,l1,
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
      t2 = dlog(z)
      t3 = t2**2
      t4 = z**2
      t5 = t3*t4
      t8 = s1l1*s2l2
      t9 = 1-z
      t10 = dlog(t9)
      t11 = t10**2
      t12 = t4*z
      t13 = t11*t12
      t16 = Li2(-z)
      t17 = t12*t16
      t20 = 0.3141592653589793D1**2
      t21 = t20*t12
      t26 = dlog(t9/z)
      t27 = t26*t12
      t30 = t12*t10
      t33 = z*t2
      t36 = t4*t2
      t39 = t12*t2
      t46 = t26*t4
      t49 = t20*z
      t52 = t20*t4
      t55 = t26*z
      t58 = -14*t1*t5-18*t8*t13-36*t8*t17-8.D0/3.D0*t8*t21-16*t8*t27+36*
     #t1*t30+7*t1*t33+80*t1*t36+132*t8*t39+36*t8*t30-18*t1*t13+16*t1*t46
     #-4.D0/3.D0*t1*t49-10.D0/3.D0*t1*t52-8*t1*t55
      t59 = z*t16
      t76 = z*t10
      t79 = t4*t10
      t82 = t26**2
      t83 = t82*t4
      t88 = t4*t16
      t95 = -18*t1*t59+10*t8*z+167*t8*t4-178*t1*t12+10*t1*z+167*t1*t4+7*
     #t8*t33+80*t8*t36+8*t1*t76-36*t1*t79-8*t1*t83-178*t8*t12-36*t8*t88-
     #4.D0/3.D0*t8*t49-10.D0/3.D0*t8*t52
      t97 = t11*t4
      t104 = t3*z
      t107 = t11*z
      t111 = dlog(1+z)
      t112 = t39*t111
      t115 = t82*t12
      t132 = t82*z
      t135 = 18*t8*t97-36*t1*t88-8*t8*t83-11*t8*t104-9*t1*t107-36*t1*t11
     #2+8*t8*t115-14*t8*t5+8*t1*t115-36*t8*t112-11*t1*t104+132*t1*t39+8*
     #t8*t76-36*t8*t79+4*t1*t132
      t136 = t36*t111
      t139 = t33*t111
      t166 = -36*t8*t136-18*t8*t139-36*t1*t136-18*t1*t139-18*t8*t59-9*t8
     #*t107-36*t1*t17-8.D0/3.D0*t1*t21-16*t1*t27+20*t1+20*t8+18*t1*t97-8
     #*t8*t55+16*t8*t46+4*t8*t132
      taux18 = t58+t95+t135+t166
      t1 = s1l1*s2l2
      t2 = dlog(z)
      t3 = z*t2
      t5 = dlog(1+z)
      t6 = t3*t5
      t9 = s1l2*s2l1
      t10 = z**2
      t11 = t10*t2
      t12 = t11*t5
      t19 = 1-z
      t20 = dlog(t19)
      t21 = t20**2
      t22 = t21*z
      t25 = z*t10
      t26 = Li2(-z)
      t27 = t25*t26
      t30 = 0.3141592653589793D1**2
      t31 = t30*t25
      t36 = dlog(t19/z)
      t37 = t36*t25
      t40 = t21*t25
      t43 = t36**2
      t44 = t43*t25
      t47 = t2**2
      t48 = t47*t10
      t57 = -72*t1*t6+144*t9*t12+72*t9*t6+80*t1-80*t9-36*t1*t22+144*t9*t
     #27+32.D0/3.D0*t9*t31+64*t9*t37+72*t9*t40+32*t1*t44-56*t1*t48-32*t9
     #*t44+56*t9*t48-72*t1*t40
      t64 = t25*t2
      t67 = t25*t20
      t74 = t20*z
      t77 = t10*t20
      t80 = t43*t10
      t83 = t47*z
      t92 = t43*z
      t95 = -144*t1*t27-32.D0/3.D0*t1*t31-64*t1*t37+528*t1*t64+144*t1*t6
     #7+28*t1*t3+320*t1*t11-32*t9*t74+144*t9*t77+32*t9*t80+44*t9*t83-528
     #*t9*t64+32*t1*t74-144*t1*t77-16*t9*t92
      t97 = t36*t10
      t100 = t30*z
      t103 = t30*t10
      t120 = t64*t5
      t125 = t36*z
      t128 = t26*z
      t131 = t10*t26
      t134 = -64*t9*t97+16.D0/3.D0*t9*t100+40.D0/3.D0*t9*t103-320*t9*t11
     #-712*t1*t25+40*t1*z+668*t1*t10+712*t9*t25-40*t9*z-668*t9*t10+144*t
     #9*t120-144*t1*t120+32*t9*t125+72*t9*t128+144*t9*t131
      t141 = t21*t10
      t166 = -32*t1*t80-44*t1*t83+36*t9*t22-72*t9*t141-32*t1*t125+64*t1*
     #t97+16*t1*t92-144*t1*t131-16.D0/3.D0*t1*t100-40.D0/3.D0*t1*t103+72
     #*t1*t141-72*t1*t128-28*t9*t3-144*t9*t67-144*t1*t12
      taux19 = t57+t95+t134+t166

CCC Coefficients

      tcoff35 = taux18
      tcoff36 = taux18
      tcoff37 = taux19
      tcoff38 = taux19

CCC Jet functions

      tFjet35 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff35)
      tFjet36 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff36)
      tFjet37 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff38)

CCC Sum

      qgscnum10e0 = tFjet35*tcoff35+tFjet36*tcoff36+tFjet37*tcoff37+tFje
     #t38*tcoff38
      RETURN
      END

      double precision function qgscnum10e1(var,ndim,z,ux,l1,
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
      t2 = dlog(z)
      t3 = t1*t2
      t5 = dlog(1+z)
      t8 = z*t2
      t11 = t1*z
      t12 = t11*t2
      t15 = 1-z
      t16 = dlog(t15)
      t21 = dlog(t15/z)
      t22 = t21**2
      t25 = t16**2
      t28 = t2**2
      t36 = Li2(-z)
      t42 = -20+36*t3*t5+18*t8*t5+36*t12*t5-36*t11*t16+8*t1*t22+9*t25*z+
     #11*t28*z-80*t3+178*t11-10*z+36*t1*t16+36*t11*t36+18*t25*t11-132*t1
     #2
      t46 = 0.3141592653589793D1**2
      t72 = -7*t8-8*z*t16+10.D0/3.D0*t1*t46+4.D0/3.D0*t46*z+36*t1*t36-4*
     #t22*z-16*t1*t21+8*t21*z+18*t36*z-18*t1*t25+8.D0/3.D0*t46*t11+14*t2
     #8*t1-8*t22*t11+16*t21*t11-167*t1
      taux18 = t42+t72
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = dlog(z)
      t4 = t2*t3
      t6 = dlog(1+z)
      t7 = t4*t6
      t10 = s1l1*s2l2
      t13 = t3*z
      t14 = t13*t6
      t21 = 0.3141592653589793D1**2
      t22 = t21*t2
      t25 = 1-z
      t28 = dlog(t25/z)
      t29 = t28*z
      t32 = Li2(-z)
      t33 = z*t32
      t36 = t2*t32
      t39 = t28**2
      t40 = t39*t2
      t43 = t3**2
      t44 = t43*z
      t47 = dlog(t25)
      t48 = t47**2
      t49 = t48*z
      t52 = t48*t2
      t57 = -432*t1*t7+432*t10*t7+216*t10*t14-216*t1*t14+240*t1-240*t10-
     #40*t1*t22-96*t1*t29-216*t1*t33-432*t1*t36+96*t10*t40+132*t10*t44-1
     #08*t1*t49+216*t1*t52+96*t10*t29
      t58 = t28*t2
      t61 = t39*z
      t66 = t21*z
      t77 = z*t2
      t78 = t77*t32
      t81 = t77*t47
      t84 = t77*t3
      t85 = t84*t6
      t96 = -192*t10*t58-48*t10*t61+432*t36*t10+16*t10*t66+40*t10*t22-21
     #6*t10*t52+216*t10*t33+108*t10*t49-432*t1*t78+432*t1*t81-432*t1*t85
     #+432*t10*t85+84*t1*t13+960*t1*t4-1584*t10*t84
      t104 = t47*z
      t107 = t2*t47
      t130 = -432*t10*t81-84*t10*t13-960*t4*t10+96*t1*t104-432*t1*t107-9
     #6*t1*t40-132*t1*t44+1584*t1*t84-96*t10*t104+432*t10*t107+48*t1*t61
     #+192*t1*t58-16*t1*t66+2136*t10*t77-120*z*t10
      t139 = t21*t77
      t142 = t28*t77
      t145 = t48*t77
      t148 = t39*t77
      t151 = t43*t2
      t166 = -2004*t10*t2-2136*t1*t77+120*t1*z+2004*t1*t2-32*t1*t139-192
     #*t1*t142-216*t1*t145-96*t10*t148+168*t10*t151+96*t1*t148-168*t1*t1
     #51+216*t10*t145+432*t10*t78+32*t10*t139+192*t10*t142
      taux19 = t57+t96+t130+t166

CCC Coefficients

      tcoff35 = taux18
      tcoff36 = taux18
      tcoff37 = taux19
      tcoff38 = taux19

CCC Jet functions

      tFjet35 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff35)
      tFjet36 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff36)
      tFjet37 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff38)

CCC Sum

      qgscnum10e1 = tFjet35*tcoff35+tFjet36*tcoff36+tFjet37*tcoff37+tFje
     #t38*tcoff38
      RETURN
      END

      double precision function qgscnum10e2(var,ndim,z,ux,l1,
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

      taux18 = 0
      t1 = s1l2*s2l1
      t2 = 0.3141592653589793D1**2
      t3 = t2*z
      t6 = z**2
      t7 = t2*t6
      t10 = 1-z
      t13 = dlog(t10/z)
      t14 = t13*z
      t17 = Li2(-z)
      t18 = z*t17
      t21 = t6*t17
      t24 = s1l1*s2l2
      t25 = t13**2
      t26 = t25*t6
      t31 = t6*z
      t32 = t25*t31
      t35 = dlog(z)
      t36 = t35**2
      t37 = t36*t6
      t44 = dlog(t10)
      t45 = t44**2
      t46 = t45*t31
      t49 = t31*t17
      t56 = 32.D0/3.D0*t1*t3+80.D0/3.D0*t1*t7+64*t1*t14+144*t1*t18+288*t
     #1*t21-64*t24*t26-1336*t1*t6+64*t24*t32-112*t24*t37-64*t1*t32+112*t
     #1*t37-144*t24*t46-288*t24*t49-1424*t24*t31+80*t24*z
      t63 = t31*t35
      t66 = z*t44
      t69 = t6*t44
      t72 = t25*z
      t75 = t13*t6
      t78 = t2*t31
      t81 = t13*t31
      t85 = t31*t44
      t88 = z*t35
      t92 = dlog(1+z)
      t93 = t63*t92
      t96 = t36*z
      t99 = 1336*t24*t6+1424*t1*t31-80*t1*z-1056*t1*t63+64*t24*t66-288*t
     #24*t69-32*t1*t72-128*t1*t75-64.D0/3.D0*t24*t78-128*t24*t81-160*t1-
     #288*t1*t85-56*t1*t88+288*t1*t93-88*t24*t96
      t101 = t45*z
      t104 = t45*t6
      t123 = t6*t35
      t124 = t123*t92
      t127 = t88*t92
      t135 = 72*t1*t101-144*t1*t104-64*t24*t14+128*t24*t75+32*t24*t72-28
     #8*t24*t93+288*t1*t49+64.D0/3.D0*t1*t78+128*t1*t81+144*t1*t46-288*t
     #24*t124-144*t24*t127+160*t24+288*t1*t124+144*t1*t127
      t166 = -640*t1*t123+1056*t24*t63+288*t24*t85+56*t24*t88+640*t24*t1
     #23-64*t1*t66+288*t1*t69+64*t1*t26+88*t1*t96-288*t24*t21-32.D0/3.D0
     #*t24*t3-80.D0/3.D0*t24*t7+144*t24*t104-144*t24*t18-72*t24*t101
      taux19 = t56+t99+t135+t166

CCC Coefficients

      tcoff35 = taux18
      tcoff36 = taux18
      tcoff37 = taux19
      tcoff38 = taux19

CCC Jet functions

      tFjet35 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff35)
      tFjet36 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff36)
      tFjet37 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff38)

CCC Sum

      qgscnum10e2 = tFjet35*tcoff35+tFjet36*tcoff36+tFjet37*tcoff37+tFje
     #t38*tcoff38
      RETURN
      END

      double precision function qgscnum10e3(var,ndim,z,ux,l1,
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

      taux18 = 0
      taux19 = 0

CCC Coefficients

      tcoff35 = taux18
      tcoff36 = taux18
      tcoff37 = taux19
      tcoff38 = taux19

CCC Jet functions

      tFjet35 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff35)
      tFjet36 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff36)
      tFjet37 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff38)

CCC Sum

      qgscnum10e3 = tFjet35*tcoff35+tFjet36*tcoff36+tFjet37*tcoff37+tFje
     #t38*tcoff38
      RETURN
      END

      double precision function qgscnum11e0(var,ndim,z,ux,l1,
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
      t5 = 12*t1*t3
      t7 = 6*t1*z
      t9 = 12*t1*t2
      t10 = s1l1*s2l2
      t12 = 12*t10*t3
      t14 = 6*t10*z
      t16 = 12*t10*t2
      t18 = dlog(1-z)
      taux20 = t5+t7-t9+t12+t14-t16+4.D0/3.D0*t18*(t5+t7-t9+t12+t14-t16)
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t5 = 48*t1*t3
      t7 = 24*t1*z
      t9 = 48*t1*t2
      t10 = s1l1*s2l2
      t12 = 48*t10*t3
      t14 = 24*t10*z
      t16 = 48*t10*t2
      t18 = dlog(1-z)
      taux21 = -t5-t7+t9+t12+t14-t16+4.D0/3.D0*t18*(-t5-t7+t9+t12+t14-t1
     #6)

CCC Coefficients

      tcoff39 = taux20
      tcoff40 = taux20
      tcoff41 = taux21
      tcoff42 = taux21

CCC Jet functions

      tFjet39 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff40)
      tFjet41 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff41)
      tFjet42 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff42)

CCC Sum

      qgscnum11e0 = tFjet39*tcoff39+tFjet40*tcoff40+tFjet41*tcoff41+tFje
     #t42*tcoff42
      RETURN
      END

      double precision function qgscnum11e1(var,ndim,z,ux,l1,
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
      taux20 = -t3-t4+t5+4.D0/3.D0*t7*(-t3-t4+t5)
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t5 = 144*t1*t3
      t7 = 72*t1*z
      t9 = 144*t1*t2
      t10 = s1l1*s2l2
      t12 = 144*t10*t3
      t14 = 72*t10*z
      t16 = 144*t10*t2
      t18 = dlog(1-z)
      taux21 = t5+t7-t9-t12-t14+t16+4.D0/3.D0*t18*(t5+t7-t9-t12-t14+t16)

CCC Coefficients

      tcoff39 = taux20
      tcoff40 = taux20
      tcoff41 = taux21
      tcoff42 = taux21

CCC Jet functions

      tFjet39 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff40)
      tFjet41 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff41)
      tFjet42 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff42)

CCC Sum

      qgscnum11e1 = tFjet39*tcoff39+tFjet40*tcoff40+tFjet41*tcoff41+tFje
     #t42*tcoff42
      RETURN
      END

      double precision function qgscnum11e2(var,ndim,z,ux,l1,
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

      taux20 = 0
      t1 = s1l2*s2l1
      t2 = z**2
      t3 = t2*z
      t5 = 96*t1*t3
      t7 = 48*t1*z
      t9 = 96*t1*t2
      t10 = s1l1*s2l2
      t12 = 96*t10*t3
      t14 = 48*t10*z
      t16 = 96*t10*t2
      t18 = dlog(1-z)
      taux21 = -t5-t7+t9+t12+t14-t16+4.D0/3.D0*t18*(-t5-t7+t9+t12+t14-t1
     #6)

CCC Coefficients

      tcoff39 = taux20
      tcoff40 = taux20
      tcoff41 = taux21
      tcoff42 = taux21

CCC Jet functions

      tFjet39 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff40)
      tFjet41 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff41)
      tFjet42 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff42)

CCC Sum

      qgscnum11e2 = tFjet39*tcoff39+tFjet40*tcoff40+tFjet41*tcoff41+tFje
     #t42*tcoff42
      RETURN
      END

      double precision function qgscnum11e3(var,ndim,z,ux,l1,
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

      taux20 = 0
      taux21 = 0

CCC Coefficients

      tcoff39 = taux20
      tcoff40 = taux20
      tcoff41 = taux21
      tcoff42 = taux21

CCC Jet functions

      tFjet39 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff40)
      tFjet41 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff41)
      tFjet42 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff42)

CCC Sum

      qgscnum11e3 = tFjet39*tcoff39+tFjet40*tcoff40+tFjet41*tcoff41+tFje
     #t42*tcoff42
      RETURN
      END

      double precision function qgscnum12e0(var,ndim,z,ux,l1,
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
      taux22 = (24*s1l2*s2l1+24*s1l1*s2l2)*z*(t6/t10/4+t16/4)*t13
      t6 = z**2
      t9 = (1-z)*l2+z
      t10 = t9**2
      t13 = 1/t9
      t16 = (1-z*t13)**2
      taux23 = (-96*s1l2*s2l1+96*s1l1*s2l2)*z*(t6/t10/4+t16/4)*t13

CCC Coefficients

      tcoff43 = taux22
      tcoff44 = taux22
      tcoff45 = taux23
      tcoff46 = taux23

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet43 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,t3,z/t3,Tcoef
     #f*tcoff43)
      t4 = (1-z)*l2+z
      tFjet44 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,t4,z/t4,Tco
     #eff*tcoff44)
      t3 = (1-z)*l2+z
      tFjet45 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,t3,z/t3,Tcoef
     #f*tcoff45)
      t4 = (1-z)*l2+z
      tFjet46 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,t4,z/t4,Tco
     #eff*tcoff46)

CCC Sum

      qgscnum12e0 = tFjet43*tcoff43+tFjet44*tcoff44+tFjet45*tcoff45+tFje
     #t46*tcoff46
      RETURN
      END

      double precision function qgscnum12e1(var,ndim,z,ux,l1,
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
      taux22 = -24*z*(t1/t5/4+t11/4)*t8
      t6 = z**2
      t9 = (1-z)*l2+z
      t10 = t9**2
      t13 = 1/t9
      t16 = (1-z*t13)**2
      taux23 = (288*s1l2*s2l1-288*s1l1*s2l2)*z*(t6/t10/4+t16/4)*t13

CCC Coefficients

      tcoff43 = taux22
      tcoff44 = taux22
      tcoff45 = taux23
      tcoff46 = taux23

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet43 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,t3,z/t3,Tcoef
     #f*tcoff43)
      t4 = (1-z)*l2+z
      tFjet44 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,t4,z/t4,Tco
     #eff*tcoff44)
      t3 = (1-z)*l2+z
      tFjet45 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,t3,z/t3,Tcoef
     #f*tcoff45)
      t4 = (1-z)*l2+z
      tFjet46 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,t4,z/t4,Tco
     #eff*tcoff46)

CCC Sum

      qgscnum12e1 = tFjet43*tcoff43+tFjet44*tcoff44+tFjet45*tcoff45+tFje
     #t46*tcoff46
      RETURN
      END

      double precision function qgscnum12e2(var,ndim,z,ux,l1,
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

      taux22 = 0
      t6 = z**2
      t9 = (1-z)*l2+z
      t10 = t9**2
      t13 = 1/t9
      t16 = (1-z*t13)**2
      taux23 = (-192*s1l2*s2l1+192*s1l1*s2l2)*z*(t6/t10/4+t16/4)*t13

CCC Coefficients

      tcoff43 = taux22
      tcoff44 = taux22
      tcoff45 = taux23
      tcoff46 = taux23

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet43 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,t3,z/t3,Tcoef
     #f*tcoff43)
      t4 = (1-z)*l2+z
      tFjet44 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,t4,z/t4,Tco
     #eff*tcoff44)
      t3 = (1-z)*l2+z
      tFjet45 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,t3,z/t3,Tcoef
     #f*tcoff45)
      t4 = (1-z)*l2+z
      tFjet46 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,t4,z/t4,Tco
     #eff*tcoff46)

CCC Sum

      qgscnum12e2 = tFjet43*tcoff43+tFjet44*tcoff44+tFjet45*tcoff45+tFje
     #t46*tcoff46
      RETURN
      END

      double precision function qgscnum12e3(var,ndim,z,ux,l1,
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

      taux22 = 0
      taux23 = 0

CCC Coefficients

      tcoff43 = taux22
      tcoff44 = taux22
      tcoff45 = taux23
      tcoff46 = taux23

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet43 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,t3,z/t3,Tcoef
     #f*tcoff43)
      t4 = (1-z)*l2+z
      tFjet44 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,t4,z/t4,Tco
     #eff*tcoff44)
      t3 = (1-z)*l2+z
      tFjet45 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,t3,z/t3,Tcoef
     #f*tcoff45)
      t4 = (1-z)*l2+z
      tFjet46 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,t4,z/t4,Tco
     #eff*tcoff46)

CCC Sum

      qgscnum12e3 = tFjet43*tcoff43+tFjet44*tcoff44+tFjet45*tcoff45+tFje
     #t46*tcoff46
      RETURN
      END

      double precision function qgscnum13e0(var,ndim,z,ux,l1,
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
      taux24 = z*(t1/4+t3/4)*(64*s1l1+128*s1l1*s2l1+64*s2l1)
      t1 = z**2
      t3 = (1-z)**2
      taux25 = z*(t1/4+t3/4)*(-256*s2l1+256*s1l1)

CCC Coefficients

      tcoff47 = taux24
      tcoff48 = taux24
      tcoff49 = taux25
      tcoff50 = taux25

CCC Jet functions

      tFjet47 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff47)
      tFjet48 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff48)
      tFjet49 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff50)

CCC Sum

      qgscnum13e0 = tFjet47*tcoff47+tFjet48*tcoff48+tFjet49*tcoff49+tFje
     #t50*tcoff50
      RETURN
      END

      double precision function qgscnum13e1(var,ndim,z,ux,l1,
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
      taux24 = z*(t1/4+t3/4)*(-96*s2l1-96*s1l1-64*s1l1*s2l1)
      t1 = z**2
      t3 = (1-z)**2
      taux25 = z*(t1/4+t3/4)*(896*s2l1-896*s1l1)

CCC Coefficients

      tcoff47 = taux24
      tcoff48 = taux24
      tcoff49 = taux25
      tcoff50 = taux25

CCC Jet functions

      tFjet47 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff47)
      tFjet48 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff48)
      tFjet49 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff50)

CCC Sum

      qgscnum13e1 = tFjet47*tcoff47+tFjet48*tcoff48+tFjet49*tcoff49+tFje
     #t50*tcoff50
      RETURN
      END

      double precision function qgscnum13e2(var,ndim,z,ux,l1,
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
      taux24 = z*(t1/4+t3/4)*(96*s2l1+128*s1l1*s2l1+96*s1l1)
      t1 = z**2
      t3 = (1-z)**2
      taux25 = z*(t1/4+t3/4)*(-896*s2l1+896*s1l1)

CCC Coefficients

      tcoff47 = taux24
      tcoff48 = taux24
      tcoff49 = taux25
      tcoff50 = taux25

CCC Jet functions

      tFjet47 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff47)
      tFjet48 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff48)
      tFjet49 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff50)

CCC Sum

      qgscnum13e2 = tFjet47*tcoff47+tFjet48*tcoff48+tFjet49*tcoff49+tFje
     #t50*tcoff50
      RETURN
      END

      double precision function qgscnum13e3(var,ndim,z,ux,l1,
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
      taux24 = z*(t1/4+t3/4)*(-64*s2l1-64*s1l1)
      t1 = z**2
      t3 = (1-z)**2
      taux25 = z*(t1/4+t3/4)*(384*s2l1-384*s1l1)

CCC Coefficients

      tcoff47 = taux24
      tcoff48 = taux24
      tcoff49 = taux25
      tcoff50 = taux25

CCC Jet functions

      tFjet47 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,5.D0,1.D0,z,Tcoeff
     #*tcoff47)
      tFjet48 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,5.D0,1.D0,z,Tcoe
     #ff*tcoff48)
      tFjet49 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,7.D0,1.D0,z,Tcoeff
     #*tcoff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,7.D0,1.D0,z,Tcoe
     #ff*tcoff50)

CCC Sum

      qgscnum13e3 = tFjet47*tcoff47+tFjet48*tcoff48+tFjet49*tcoff49+tFje
     #t50*tcoff50
      RETURN
      END

      double precision function qgscnum14e0(var,ndim,z,ux,l1,
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
      taux26 = (128*t1*l1*s1l2+t7*s2l2+t12*t13)*s1l1+(t7*s2l1+t12*t18)*s
     #1l2+t23*t13*s2l1+t23*t18*s2l2+128*t18*t13
      t1 = z*l1
      t2 = t1-l1-z
      t5 = 1-l1-2*z+t1
      t6 = -s1l2-s2l2-zp
      t12 = -s1l1-s2l1-zp
      t17 = 1-z
      taux27 = (256*t2*s2l2+256*t5*t6)*s1l1+(-256*t2*s2l1-256*t5*t12)*s1
     #l2+256*t17*t6*s2l1-256*t17*t12*s2l2

CCC Coefficients

      tcoff51 = taux26
      tcoff52 = taux26
      tcoff53 = taux27
      tcoff54 = taux27

CCC Jet functions

      tFjet51 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff52)
      tFjet53 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff53)
      tFjet54 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff54)

CCC Sum

      qgscnum14e0 = tFjet51*tcoff51+tFjet52*tcoff52+tFjet53*tcoff53+tFje
     #t54*tcoff54
      RETURN
      END

      double precision function qgscnum14e1(var,ndim,z,ux,l1,
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
      taux26 = t3*s2l2*s1l1+t3*s2l1*s1l2+(128*(z-1)*(-1+l1)*s2l2+t3*(-s1
     #l2-s2l2-zp))*s2l1+t3*(-s1l1-s2l1-zp)*s2l2-64*z*(2*t20*l1-t20*t23-2
     #*t1-2*t20+2*z+2*z*t23-t23-1)
      t1 = 3*z
      t2 = 3*l1
      t3 = z*l1
      t4 = 3*t3
      t5 = t1+t2-t4+1
      t8 = -3+t2+6*z-t4
      t9 = -s1l2-s2l2-zp
      t15 = -s1l1-s2l1-zp
      t20 = -l1+t1+t3-3
      taux27 = (256*t5*s2l2+256*t8*t9)*s1l1+(-256*t5*s2l1-256*t8*t15)*s1
     #l2+256*t20*t9*s2l1-256*t20*t15*s2l2

CCC Coefficients

      tcoff51 = taux26
      tcoff52 = taux26
      tcoff53 = taux27
      tcoff54 = taux27

CCC Jet functions

      tFjet51 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff52)
      tFjet53 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff53)
      tFjet54 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff54)

CCC Sum

      qgscnum14e1 = tFjet51*tcoff51+tFjet52*tcoff52+tFjet53*tcoff53+tFje
     #t54*tcoff54
      RETURN
      END

      double precision function qgscnum14e2(var,ndim,z,ux,l1,
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

      t3 = z**2
      t4 = l1**2
      taux26 = -64*z*(2*z*l1+t3*t4-2*z*t4+1-2*l1+t4)
      t1 = z*l1
      t2 = 2*z
      t3 = t1-3-l1-t2
      t6 = 1-l1-t2+t1
      t7 = -s1l2-s2l2-zp
      t14 = -s1l1-s2l1-zp
      t21 = -3*t1+1+3*l1-t2
      taux27 = (256*t3*s2l2+512*t6*t7)*s1l1+(-256*t3*s2l1-512*t6*t14)*s1
     #l2+256*t21*t7*s2l1-256*t21*t14*s2l2

CCC Coefficients

      tcoff51 = taux26
      tcoff52 = taux26
      tcoff53 = taux27
      tcoff54 = taux27

CCC Jet functions

      tFjet51 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff52)
      tFjet53 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff53)
      tFjet54 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff54)

CCC Sum

      qgscnum14e2 = tFjet51*tcoff51+tFjet52*tcoff52+tFjet53*tcoff53+tFje
     #t54*tcoff54
      RETURN
      END

      double precision function qgscnum14e3(var,ndim,z,ux,l1,
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

      taux26 = 0
      t2 = 1-l1+z*l1
      taux27 = 512*t2*s2l2*s1l1-512*t2*s2l1*s1l2+512*t2*(-s1l2-s2l2-zp)*
     #s2l1-512*t2*(-s1l1-s2l1-zp)*s2l2

CCC Coefficients

      tcoff51 = taux26
      tcoff52 = taux26
      tcoff53 = taux27
      tcoff54 = taux27

CCC Jet functions

      tFjet51 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff52)
      tFjet53 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff53)
      tFjet54 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff54)

CCC Sum

      qgscnum14e3 = tFjet51*tcoff51+tFjet52*tcoff52+tFjet53*tcoff53+tFje
     #t54*tcoff54
      RETURN
      END

      double precision function qgscnum15e0(var,ndim,z,ux,l1,
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
      t12 = (-l1-z+t10)*t6
      t17 = (1-l1-2*z+t10)*t6
      t18 = -s1l2-s2l2-zp
      t24 = -s1l1-s2l1-zp
      t29 = t3*t6
      taux28 = (128.D0/3.D0*t1*l1*t6*s1l2+64.D0/3.D0*t12*s2l2+64.D0/3.D0
     #*t17*t18)*s1l1+(64.D0/3.D0*t12*s2l1+64.D0/3.D0*t17*t24)*s1l2+64.D0
     #/3.D0*t29*t18*s2l1+64.D0/3.D0*t29*t24*s2l2+(128+512.D0/3.D0*t4)*t1
     #8*t24
      t1 = z*l1
      t3 = 1-z
      t4 = dlog(t3)
      t6 = 3+4*t4
      t7 = (-l1-z+t1)*t6
      t11 = (1-l1-2*z+t1)*t6
      t12 = -s1l2-s2l2-zp
      t18 = -s1l1-s2l1-zp
      t23 = t3*t6
      taux29 = (256.D0/3.D0*t7*s2l2+256.D0/3.D0*t11*t12)*s1l1+(-256.D0/3
     #.D0*t7*s2l1-256.D0/3.D0*t11*t18)*s1l2+256.D0/3.D0*t23*t12*s2l1-256
     #.D0/3.D0*t23*t18*s2l2

CCC Coefficients

      tcoff55 = taux28
      tcoff56 = taux28
      tcoff57 = taux29
      tcoff58 = taux29

CCC Jet functions

      tFjet55 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff55)
      tFjet56 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff56)
      tFjet57 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff57)
      tFjet58 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff58)

CCC Sum

      qgscnum15e0 = tFjet55*tcoff55+tFjet56*tcoff56+tFjet57*tcoff57+tFje
     #t58*tcoff58
      RETURN
      END

      double precision function qgscnum15e1(var,ndim,z,ux,l1,
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
      t7 = (1-l1+t1)*t6
      t29 = z**2
      t32 = l1**2
      taux28 = 64.D0/3.D0*t7*s2l2*s1l1+64.D0/3.D0*t7*s2l1*s1l2+(-128.D0/
     #3.D0*t3*(-1+l1)*t6*s2l2+64.D0/3.D0*t7*(-s1l2-s2l2-zp))*s2l1+64.D0/
     #3.D0*t7*(-s1l1-s2l1-zp)*s2l2-64.D0/3.D0*z*(2*t29*l1-t29*t32-2*t1-2
     #*t29+2*z+2*z*t32-t32-1)*t6
      t1 = 3*z
      t2 = 3*l1
      t3 = z*l1
      t4 = 3*t3
      t7 = dlog(1-z)
      t9 = 3+4*t7
      t10 = (t1+t2-t4+1)*t9
      t14 = (-3+t2+6*z-t4)*t9
      t15 = -s1l2-s2l2-zp
      t21 = -s1l1-s2l1-zp
      t27 = (-l1+t1+t3-3)*t9
      taux29 = (256.D0/3.D0*t10*s2l2+256.D0/3.D0*t14*t15)*s1l1+(-256.D0/
     #3.D0*t10*s2l1-256.D0/3.D0*t14*t21)*s1l2+256.D0/3.D0*t27*t15*s2l1-2
     #56.D0/3.D0*t27*t21*s2l2

CCC Coefficients

      tcoff55 = taux28
      tcoff56 = taux28
      tcoff57 = taux29
      tcoff58 = taux29

CCC Jet functions

      tFjet55 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff55)
      tFjet56 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff56)
      tFjet57 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff57)
      tFjet58 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff58)

CCC Sum

      qgscnum15e1 = tFjet55*tcoff55+tFjet56*tcoff56+tFjet57*tcoff57+tFje
     #t58*tcoff58
      RETURN
      END

      double precision function qgscnum15e2(var,ndim,z,ux,l1,
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

      t3 = z**2
      t4 = l1**2
      t12 = dlog(1-z)
      taux28 = -64.D0/3.D0*z*(2*z*l1+t3*t4-2*z*t4+1-2*l1+t4)*(3+4*t12)
      t1 = z*l1
      t2 = 2*z
      t5 = dlog(1-z)
      t7 = 3+4*t5
      t8 = (-3-l1+t1-t2)*t7
      t12 = (1-l1-t2+t1)*t7
      t13 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t28 = (1+3*l1-3*t1-t2)*t7
      taux29 = (256.D0/3.D0*t8*s2l2+512.D0/3.D0*t12*t13)*s1l1+(-256.D0/3
     #.D0*t8*s2l1-512.D0/3.D0*t12*t20)*s1l2+256.D0/3.D0*t28*t13*s2l1-256
     #.D0/3.D0*t28*t20*s2l2

CCC Coefficients

      tcoff55 = taux28
      tcoff56 = taux28
      tcoff57 = taux29
      tcoff58 = taux29

CCC Jet functions

      tFjet55 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff55)
      tFjet56 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff56)
      tFjet57 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff57)
      tFjet58 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff58)

CCC Sum

      qgscnum15e2 = tFjet55*tcoff55+tFjet56*tcoff56+tFjet57*tcoff57+tFje
     #t58*tcoff58
      RETURN
      END

      double precision function qgscnum15e3(var,ndim,z,ux,l1,
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

      taux28 = 0
      t4 = dlog(1-z)
      t7 = (1-l1+z*l1)*(3+4*t4)
      taux29 = 512.D0/3.D0*t7*s2l2*s1l1-512.D0/3.D0*t7*s2l1*s1l2+512.D0/
     #3.D0*t7*(-s1l2-s2l2-zp)*s2l1-512.D0/3.D0*t7*(-s1l1-s2l1-zp)*s2l2

CCC Coefficients

      tcoff55 = taux28
      tcoff56 = taux28
      tcoff57 = taux29
      tcoff58 = taux29

CCC Jet functions

      tFjet55 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff55)
      tFjet56 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff56)
      tFjet57 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff57)
      tFjet58 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff58)

CCC Sum

      qgscnum15e3 = tFjet55*tcoff55+tFjet56*tcoff56+tFjet57*tcoff57+tFje
     #t58*tcoff58
      RETURN
      END

      double precision function qgscnum16e0(var,ndim,z,ux,l1,
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
      taux30 = (-2.D0/3.D0-2.D0/3.D0*yp)*((128*t3*l1*s1l2+t9*s2l2+t14*t1
     #5)*s1l1+(t9*s2l1+t14*t20)*s1l2+t25*t15*s2l1+t25*t20*s2l2+128*t20*t
     #15)
      t3 = zp*l1
      t4 = t3-l1-zp
      t7 = 1-l1-2*zp+t3
      t8 = -s1l2-s2l2-zp
      t14 = -s1l1-s2l1-zp
      t19 = -zp+1
      taux31 = (-2.D0/3.D0-2.D0/3.D0*yp)*((256*t4*s2l2+256*t7*t8)*s1l1+(
     #-256*t4*s2l1-256*t7*t14)*s1l2+256*t19*t8*s2l1-256*t19*t14*s2l2)

CCC Coefficients

      tcoff59 = taux30
      tcoff60 = taux30
      tcoff61 = taux31
      tcoff62 = taux31

CCC Jet functions

      tFjet59 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f59)
      tFjet60 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off60)
      tFjet61 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f61)
      tFjet62 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off62)

CCC Sum

      qgscnum16e0 = tFjet59*tcoff59+tFjet60*tcoff60+tFjet61*tcoff61+tFje
     #t62*tcoff62
      RETURN
      END

      double precision function qgscnum16e1(var,ndim,z,ux,l1,
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
      taux30 = (-2.D0/3.D0-2.D0/3.D0*yp)*(t5*s2l2*s1l1+t5*s2l1*s1l2+(128
     #*(zp-1)*(-1+l1)*s2l2+t5*(-s1l2-s2l2-zp))*s2l1+t5*(-s1l1-s2l1-zp)*s
     #2l2-64*zp*(2*t22*l1-t22*t25-2*t3-2*t22+2*zp+2*zp*t25-t25-1))
      t3 = 3*zp
      t4 = 3*l1
      t5 = zp*l1
      t6 = 3*t5
      t7 = t3+t4-t6+1
      t10 = -3+t4+6*zp-t6
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = -l1+t3+t5-3
      taux31 = (-2.D0/3.D0-2.D0/3.D0*yp)*((256*t7*s2l2+256*t10*t11)*s1l1
     #+(-256*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2
     #)

CCC Coefficients

      tcoff59 = taux30
      tcoff60 = taux30
      tcoff61 = taux31
      tcoff62 = taux31

CCC Jet functions

      tFjet59 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f59)
      tFjet60 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off60)
      tFjet61 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f61)
      tFjet62 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off62)

CCC Sum

      qgscnum16e1 = tFjet59*tcoff59+tFjet60*tcoff60+tFjet61*tcoff61+tFje
     #t62*tcoff62
      RETURN
      END

      double precision function qgscnum16e2(var,ndim,z,ux,l1,
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

      t6 = zp**2
      t7 = l1**2
      taux30 = -64*(-2.D0/3.D0-2.D0/3.D0*yp)*zp*(2*zp*l1+t6*t7-2*zp*t7+1
     #-2*l1+t7)
      t3 = zp*l1
      t4 = 2*zp
      t5 = t3-3-l1-t4
      t8 = 1-l1-t4+t3
      t9 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t23 = -3*t3+1+3*l1-t4
      taux31 = (-2.D0/3.D0-2.D0/3.D0*yp)*((256*t5*s2l2+512*t8*t9)*s1l1+(
     #-256*t5*s2l1-512*t8*t16)*s1l2+256*t23*t9*s2l1-256*t23*t16*s2l2)

CCC Coefficients

      tcoff59 = taux30
      tcoff60 = taux30
      tcoff61 = taux31
      tcoff62 = taux31

CCC Jet functions

      tFjet59 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f59)
      tFjet60 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off60)
      tFjet61 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f61)
      tFjet62 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off62)

CCC Sum

      qgscnum16e2 = tFjet59*tcoff59+tFjet60*tcoff60+tFjet61*tcoff61+tFje
     #t62*tcoff62
      RETURN
      END

      double precision function qgscnum16e3(var,ndim,z,ux,l1,
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

      taux30 = 0
      t4 = 1-l1+zp*l1
      taux31 = (-2.D0/3.D0-2.D0/3.D0*yp)*(512*t4*s2l2*s1l1-512*t4*s2l1*s
     #1l2+512*t4*(-s1l2-s2l2-zp)*s2l1-512*t4*(-s1l1-s2l1-zp)*s2l2)

CCC Coefficients

      tcoff59 = taux30
      tcoff60 = taux30
      tcoff61 = taux31
      tcoff62 = taux31

CCC Jet functions

      tFjet59 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f59)
      tFjet60 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off60)
      tFjet61 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f61)
      tFjet62 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off62)

CCC Sum

      qgscnum16e3 = tFjet59*tcoff59+tFjet60*tcoff60+tFjet61*tcoff61+tFje
     #t62*tcoff62
      RETURN
      END

      double precision function qgscnum17e0(var,ndim,z,ux,l1,
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
      taux32 = 4.D0/3.D0*(128*t1*l1*s1l2+t7*s2l2+t12*t13)*s1l1+4.D0/3.D0
     #*(t7*s2l1+t12*t19)*s1l2+4.D0/3.D0*t25*t13*s2l1+4.D0/3.D0*t25*t19*s
     #2l2+512.D0/3.D0*t19*t13
      t1 = zp*l1
      t2 = t1-l1-zp
      t5 = 1-l1-2*zp+t1
      t6 = -s1l2-s2l2-zp
      t13 = -s1l1-s2l1-zp
      t19 = -zp+1
      taux33 = 4.D0/3.D0*(256*t2*s2l2+256*t5*t6)*s1l1+4.D0/3.D0*(-256*t2
     #*s2l1-256*t5*t13)*s1l2+1024.D0/3.D0*t19*t6*s2l1-1024.D0/3.D0*t19*t
     #13*s2l2

CCC Coefficients

      tcoff63 = taux32
      tcoff64 = taux32
      tcoff65 = taux33
      tcoff66 = taux33

CCC Jet functions

      tFjet63 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f63)
      tFjet64 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off64)
      tFjet65 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f65)
      tFjet66 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off66)

CCC Sum

      qgscnum17e0 = tFjet63*tcoff63+tFjet64*tcoff64+tFjet65*tcoff65+tFje
     #t66*tcoff66
      RETURN
      END

      double precision function qgscnum17e1(var,ndim,z,ux,l1,
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
      taux32 = 4.D0/3.D0*t3*s2l2*s1l1+4.D0/3.D0*t3*s2l1*s1l2+4.D0/3.D0*(
     #128*(zp-1)*(-1+l1)*s2l2+t3*(-s1l2-s2l2-zp))*s2l1+4.D0/3.D0*t3*(-s1
     #l1-s2l1-zp)*s2l2-256.D0/3.D0*zp*(2*t24*l1-t24*t27-2*t1-2*t24+2*zp+
     #2*zp*t27-t27-1)
      t1 = 3*zp
      t2 = 3*l1
      t3 = zp*l1
      t4 = 3*t3
      t5 = t1+t2-t4+1
      t8 = -3+t2+6*zp-t4
      t9 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t22 = -l1+t1+t3-3
      taux33 = 4.D0/3.D0*(256*t5*s2l2+256*t8*t9)*s1l1+4.D0/3.D0*(-256*t5
     #*s2l1-256*t8*t16)*s1l2+1024.D0/3.D0*t22*t9*s2l1-1024.D0/3.D0*t22*t
     #16*s2l2

CCC Coefficients

      tcoff63 = taux32
      tcoff64 = taux32
      tcoff65 = taux33
      tcoff66 = taux33

CCC Jet functions

      tFjet63 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f63)
      tFjet64 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off64)
      tFjet65 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f65)
      tFjet66 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off66)

CCC Sum

      qgscnum17e1 = tFjet63*tcoff63+tFjet64*tcoff64+tFjet65*tcoff65+tFje
     #t66*tcoff66
      RETURN
      END

      double precision function qgscnum17e2(var,ndim,z,ux,l1,
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

      t3 = zp**2
      t4 = l1**2
      taux32 = -256.D0/3.D0*zp*(2*zp*l1+t3*t4-2*zp*t4+1-2*l1+t4)
      t1 = zp*l1
      t2 = 2*zp
      t3 = -3-l1+t1-t2
      t6 = 1-l1-t2+t1
      t7 = -s1l2-s2l2-zp
      t15 = -s1l1-s2l1-zp
      t23 = 1+3*l1-3*t1-t2
      taux33 = 4.D0/3.D0*(256*t3*s2l2+512*t6*t7)*s1l1+4.D0/3.D0*(-256*t3
     #*s2l1-512*t6*t15)*s1l2+1024.D0/3.D0*t23*t7*s2l1-1024.D0/3.D0*t23*t
     #15*s2l2

CCC Coefficients

      tcoff63 = taux32
      tcoff64 = taux32
      tcoff65 = taux33
      tcoff66 = taux33

CCC Jet functions

      tFjet63 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f63)
      tFjet64 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off64)
      tFjet65 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f65)
      tFjet66 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off66)

CCC Sum

      qgscnum17e2 = tFjet63*tcoff63+tFjet64*tcoff64+tFjet65*tcoff65+tFje
     #t66*tcoff66
      RETURN
      END

      double precision function qgscnum17e3(var,ndim,z,ux,l1,
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

      taux32 = 0
      t2 = 1-l1+zp*l1
      taux33 = 2048.D0/3.D0*t2*s2l2*s1l1-2048.D0/3.D0*t2*s2l1*s1l2+2048.
     #D0/3.D0*t2*(-s1l2-s2l2-zp)*s2l1-2048.D0/3.D0*t2*(-s1l1-s2l1-zp)*s2
     #l2

CCC Coefficients

      tcoff63 = taux32
      tcoff64 = taux32
      tcoff65 = taux33
      tcoff66 = taux33

CCC Jet functions

      tFjet63 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tcof
     #f63)
      tFjet64 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,yp,1.D0,Tcoeff*tc
     #off64)
      tFjet65 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tcof
     #f65)
      tFjet66 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,yp,1.D0,Tcoeff*tc
     #off66)

CCC Sum

      qgscnum17e3 = tFjet63*tcoff63+tFjet64*tcoff64+tFjet65*tcoff65+tFje
     #t66*tcoff66
      RETURN
      END

      double precision function qgscnum18e0(var,ndim,z,ux,l1,
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
      taux34 = (-32.D0/3.D0*t1*l1*t7*s1l2-16.D0/3.D0*t13*s2l2-16.D0/3.D0
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
      taux35 = (-64.D0/3.D0*t8*s2l2-64.D0/3.D0*t12*t13)*s1l1+(64.D0/3.D0
     #*t8*s2l1+64.D0/3.D0*t12*t19)*s1l2-64.D0/3.D0*t24*t13*s2l1+64.D0/3.
     #D0*t24*t19*s2l2

CCC Coefficients

      tcoff67 = taux34
      tcoff68 = taux34
      tcoff69 = taux35
      tcoff70 = taux35

CCC Jet functions

      tFjet67 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff67)
      tFjet68 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff68)
      tFjet69 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff69)
      tFjet70 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff70)

CCC Sum

      qgscnum18e0 = tFjet67*tcoff67+tFjet68*tcoff68+tFjet69*tcoff69+tFje
     #t70*tcoff70
      RETURN
      END

      double precision function qgscnum18e1(var,ndim,z,ux,l1,
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
      taux34 = -16.D0/3.D0*t8*s2l2*s1l1-16.D0/3.D0*t8*s2l1*s1l2+(32.D0/3
     #.D0*t4*(-1+l1)*t7*s2l2-16.D0/3.D0*t8*(-s1l2-s2l2-zp))*s2l1-16.D0/3
     #.D0*t8*(-s1l1-s2l1-zp)*s2l2+16.D0/3.D0*z*(2*t30*l1-t30*t33-2*t1-2*
     #t30+2*z+2*z*t33-t33-1)*t7
      t1 = 3*z
      t2 = 3*l1
      t3 = z*l1
      t4 = 3*t3
      t8 = dlog(1-z)
      t10 = -33+2*NF-36*t8
      t11 = (t1+t2-t4+1)*t10
      t15 = (-3+t2+6*z-t4)*t10
      t16 = -s1l2-s2l2-zp
      t22 = -s1l1-s2l1-zp
      t28 = (-l1+t1+t3-3)*t10
      taux35 = (-64.D0/3.D0*t11*s2l2-64.D0/3.D0*t15*t16)*s1l1+(64.D0/3.D
     #0*t11*s2l1+64.D0/3.D0*t15*t22)*s1l2-64.D0/3.D0*t28*t16*s2l1+64.D0/
     #3.D0*t28*t22*s2l2

CCC Coefficients

      tcoff67 = taux34
      tcoff68 = taux34
      tcoff69 = taux35
      tcoff70 = taux35

CCC Jet functions

      tFjet67 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff67)
      tFjet68 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff68)
      tFjet69 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff69)
      tFjet70 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff70)

CCC Sum

      qgscnum18e1 = tFjet67*tcoff67+tFjet68*tcoff68+tFjet69*tcoff69+tFje
     #t70*tcoff70
      RETURN
      END

      double precision function qgscnum18e2(var,ndim,z,ux,l1,
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

      t3 = z**2
      t4 = l1**2
      t13 = dlog(1-z)
      taux34 = 16.D0/3.D0*z*(2*z*l1+t3*t4-2*z*t4+1-2*l1+t4)*(-33+2*NF-36
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
      taux35 = (-64.D0/3.D0*t9*s2l2-128.D0/3.D0*t13*t14)*s1l1+(64.D0/3.D
     #0*t9*s2l1+128.D0/3.D0*t13*t21)*s1l2-64.D0/3.D0*t29*t14*s2l1+64.D0/
     #3.D0*t29*t21*s2l2

CCC Coefficients

      tcoff67 = taux34
      tcoff68 = taux34
      tcoff69 = taux35
      tcoff70 = taux35

CCC Jet functions

      tFjet67 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff67)
      tFjet68 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff68)
      tFjet69 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff69)
      tFjet70 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff70)

CCC Sum

      qgscnum18e2 = tFjet67*tcoff67+tFjet68*tcoff68+tFjet69*tcoff69+tFje
     #t70*tcoff70
      RETURN
      END

      double precision function qgscnum18e3(var,ndim,z,ux,l1,
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

      taux34 = 0
      t5 = dlog(1-z)
      t8 = (1-l1+z*l1)*(-33+2*NF-36*t5)
      taux35 = -128.D0/3.D0*t8*s2l2*s1l1+128.D0/3.D0*t8*s2l1*s1l2-128.D0
     #/3.D0*t8*(-s1l2-s2l2-zp)*s2l1+128.D0/3.D0*t8*(-s1l1-s2l1-zp)*s2l2

CCC Coefficients

      tcoff67 = taux34
      tcoff68 = taux34
      tcoff69 = taux35
      tcoff70 = taux35

CCC Jet functions

      tFjet67 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tco
     #ff67)
      tFjet68 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*t
     #coff68)
      tFjet69 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tco
     #ff69)
      tFjet70 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*t
     #coff70)

CCC Sum

      qgscnum18e3 = tFjet67*tcoff67+tFjet68*tcoff68+tFjet69*tcoff69+tFje
     #t70*tcoff70
      RETURN
      END

      double precision function qgscnum19e0(var,ndim,z,ux,l1,
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
      taux36 = (3*yp*t1+3*t1/yp-3)*((128*t7*l1*s1l2+t13*s2l2+t18*t19)*s1
     #l1+(t13*s2l1+t18*t24)*s1l2+t29*t19*s2l1+t29*t24*s2l2+128*t24*t19)
      t1 = 1-yp
      t7 = zp*l1
      t8 = t7-l1-zp
      t11 = 1-l1-2*zp+t7
      t12 = -s1l2-s2l2-zp
      t18 = -s1l1-s2l1-zp
      t23 = -zp+1
      taux37 = (3*yp*t1+3*t1/yp-3)*((256*t8*s2l2+256*t11*t12)*s1l1+(-256
     #*t8*s2l1-256*t11*t18)*s1l2+256*t23*t12*s2l1-256*t23*t18*s2l2)

CCC Coefficients

      tcoff71 = taux36
      tcoff72 = taux36
      tcoff73 = taux37
      tcoff74 = taux37

CCC Jet functions

      tFjet71 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f71)
      tFjet72 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off72)
      tFjet73 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f73)
      tFjet74 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off74)

CCC Sum

      qgscnum19e0 = tFjet71*tcoff71+tFjet72*tcoff72+tFjet73*tcoff73+tFje
     #t74*tcoff74
      RETURN
      END

      double precision function qgscnum19e1(var,ndim,z,ux,l1,
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
      taux36 = (3*yp*t1+3*t1/yp-3)*(t9*s2l2*s1l1+t9*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t9*(-s1l2-s2l2-zp))*s2l1+t9*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(2*t26*l1-t26*t29-2*t7-2*t26+2*zp+2*zp*t29-t29-1))
      t1 = 1-yp
      t7 = 3*zp
      t8 = 3*l1
      t9 = zp*l1
      t10 = 3*t9
      t11 = t7+t8-t10+1
      t14 = -3+t8+6*zp-t10
      t15 = -s1l2-s2l2-zp
      t21 = -s1l1-s2l1-zp
      t26 = -l1+t7+t9-3
      taux37 = (3*yp*t1+3*t1/yp-3)*((256*t11*s2l2+256*t14*t15)*s1l1+(-25
     #6*t11*s2l1-256*t14*t21)*s1l2+256*t26*t15*s2l1-256*t26*t21*s2l2)

CCC Coefficients

      tcoff71 = taux36
      tcoff72 = taux36
      tcoff73 = taux37
      tcoff74 = taux37

CCC Jet functions

      tFjet71 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f71)
      tFjet72 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off72)
      tFjet73 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f73)
      tFjet74 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off74)

CCC Sum

      qgscnum19e1 = tFjet71*tcoff71+tFjet72*tcoff72+tFjet73*tcoff73+tFje
     #t74*tcoff74
      RETURN
      END

      double precision function qgscnum19e2(var,ndim,z,ux,l1,
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
      t10 = zp**2
      t11 = l1**2
      taux36 = -64*(3*yp*t1+3*t1/yp-3)*zp*(2*zp*l1+t10*t11-2*zp*t11+1-2*
     #l1+t11)
      t1 = 1-yp
      t7 = zp*l1
      t8 = 2*zp
      t9 = t7-3-l1-t8
      t12 = 1-l1-t8+t7
      t13 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t27 = -3*t7+1+3*l1-t8
      taux37 = (3*yp*t1+3*t1/yp-3)*((256*t9*s2l2+512*t12*t13)*s1l1+(-256
     #*t9*s2l1-512*t12*t20)*s1l2+256*t27*t13*s2l1-256*t27*t20*s2l2)

CCC Coefficients

      tcoff71 = taux36
      tcoff72 = taux36
      tcoff73 = taux37
      tcoff74 = taux37

CCC Jet functions

      tFjet71 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f71)
      tFjet72 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off72)
      tFjet73 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f73)
      tFjet74 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off74)

CCC Sum

      qgscnum19e2 = tFjet71*tcoff71+tFjet72*tcoff72+tFjet73*tcoff73+tFje
     #t74*tcoff74
      RETURN
      END

      double precision function qgscnum19e3(var,ndim,z,ux,l1,
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

      taux36 = 0
      t1 = 1-yp
      t8 = 1-l1+zp*l1
      taux37 = (3*yp*t1+3*t1/yp-3)*(512*t8*s2l2*s1l1-512*t8*s2l1*s1l2+51
     #2*t8*(-s1l2-s2l2-zp)*s2l1-512*t8*(-s1l1-s2l1-zp)*s2l2)

CCC Coefficients

      tcoff71 = taux36
      tcoff72 = taux36
      tcoff73 = taux37
      tcoff74 = taux37

CCC Jet functions

      tFjet71 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f71)
      tFjet72 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off72)
      tFjet73 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f73)
      tFjet74 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off74)

CCC Sum

      qgscnum19e3 = tFjet71*tcoff71+tFjet72*tcoff72+tFjet73*tcoff73+tFje
     #t74*tcoff74
      RETURN
      END

      double precision function qgscnum20e0(var,ndim,z,ux,l1,
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
      taux38 = 3*(128*t1*l1*s1l2+t7*s2l2+t12*t13)*s1l1+3*(t7*s2l1+t12*t1
     #9)*s1l2+3*t25*t13*s2l1+3*t25*t19*s2l2+384*t19*t13
      t1 = zp*l1
      t2 = t1-l1-zp
      t5 = 1-l1-2*zp+t1
      t6 = -s1l2-s2l2-zp
      t13 = -s1l1-s2l1-zp
      t19 = -zp+1
      taux39 = 3*(256*t2*s2l2+256*t5*t6)*s1l1+3*(-256*t2*s2l1-256*t5*t13
     #)*s1l2+768*t19*t6*s2l1-768*t19*t13*s2l2

CCC Coefficients

      tcoff75 = taux38
      tcoff76 = taux38
      tcoff77 = taux39
      tcoff78 = taux39

CCC Jet functions

      tFjet75 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f75)
      tFjet76 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off76)
      tFjet77 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f77)
      tFjet78 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off78)

CCC Sum

      qgscnum20e0 = tFjet75*tcoff75+tFjet76*tcoff76+tFjet77*tcoff77+tFje
     #t78*tcoff78
      RETURN
      END

      double precision function qgscnum20e1(var,ndim,z,ux,l1,
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
      taux38 = 3*t3*s2l2*s1l1+3*t3*s2l1*s1l2+3*(128*(zp-1)*(-1+l1)*s2l2+
     #t3*(-s1l2-s2l2-zp))*s2l1+3*t3*(-s1l1-s2l1-zp)*s2l2-192*zp*(2*t24*l
     #1-t24*t27-2*t1-2*t24+2*zp+2*zp*t27-t27-1)
      t1 = 3*zp
      t2 = 3*l1
      t3 = zp*l1
      t4 = 3*t3
      t5 = t1+t2-t4+1
      t8 = -3+t2+6*zp-t4
      t9 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t22 = -l1+t1+t3-3
      taux39 = 3*(256*t5*s2l2+256*t8*t9)*s1l1+3*(-256*t5*s2l1-256*t8*t16
     #)*s1l2+768*t22*t9*s2l1-768*t22*t16*s2l2

CCC Coefficients

      tcoff75 = taux38
      tcoff76 = taux38
      tcoff77 = taux39
      tcoff78 = taux39

CCC Jet functions

      tFjet75 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f75)
      tFjet76 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off76)
      tFjet77 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f77)
      tFjet78 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off78)

CCC Sum

      qgscnum20e1 = tFjet75*tcoff75+tFjet76*tcoff76+tFjet77*tcoff77+tFje
     #t78*tcoff78
      RETURN
      END

      double precision function qgscnum20e2(var,ndim,z,ux,l1,
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

      t3 = zp**2
      t4 = l1**2
      taux38 = -192*zp*(2*zp*l1+t3*t4-2*zp*t4+1-2*l1+t4)
      t1 = zp*l1
      t2 = 2*zp
      t3 = -3-l1+t1-t2
      t6 = 1-l1-t2+t1
      t7 = -s1l2-s2l2-zp
      t15 = -s1l1-s2l1-zp
      t23 = 1+3*l1-3*t1-t2
      taux39 = 3*(256*t3*s2l2+512*t6*t7)*s1l1+3*(-256*t3*s2l1-512*t6*t15
     #)*s1l2+768*t23*t7*s2l1-768*t23*t15*s2l2

CCC Coefficients

      tcoff75 = taux38
      tcoff76 = taux38
      tcoff77 = taux39
      tcoff78 = taux39

CCC Jet functions

      tFjet75 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f75)
      tFjet76 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off76)
      tFjet77 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f77)
      tFjet78 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off78)

CCC Sum

      qgscnum20e2 = tFjet75*tcoff75+tFjet76*tcoff76+tFjet77*tcoff77+tFje
     #t78*tcoff78
      RETURN
      END

      double precision function qgscnum20e3(var,ndim,z,ux,l1,
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

      taux38 = 0
      t2 = 1-l1+zp*l1
      taux39 = 1536*t2*s2l2*s1l1-1536*t2*s2l1*s1l2+1536*t2*(-s1l2-s2l2-z
     #p)*s2l1-1536*t2*(-s1l1-s2l1-zp)*s2l2

CCC Coefficients

      tcoff75 = taux38
      tcoff76 = taux38
      tcoff77 = taux39
      tcoff78 = taux39

CCC Jet functions

      tFjet75 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f75)
      tFjet76 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off76)
      tFjet77 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f77)
      tFjet78 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off78)

CCC Sum

      qgscnum20e3 = tFjet75*tcoff75+tFjet76*tcoff76+tFjet77*tcoff77+tFje
     #t78*tcoff78
      RETURN
      END

      double precision function qgscnum21e0(var,ndim,z,ux,l1,
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
      taux40 = (t1/4+t3/4)*((-128*t6*l1*s1l2+t11*s2l2+t15*t16)*s1l1+(t11
     #*s2l1+t15*t21)*s1l2+(128*t6*(-1+l1)*s2l2+t30*t16)*s2l1+t30*t21*s2l
     #2)
      t1 = yp**2
      t3 = (1-yp)**2
      t6 = 1+zp
      t8 = zp*l1
      t9 = -zp+t8-l1
      t10 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t21 = -1+l1-t8
      taux41 = (t1/4+t3/4)*((256*t6*s2l2-256*t9*t10)*s1l1+(-256*t6*s2l1+
     #256*t9*t16)*s1l2+256*t21*t10*s2l1-256*t21*t16*s2l2)

CCC Coefficients

      tcoff79 = taux40
      tcoff80 = taux40
      tcoff81 = taux41
      tcoff82 = taux41

CCC Jet functions

      tFjet79 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f79)
      tFjet80 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off80)
      tFjet81 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f81)
      tFjet82 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off82)

CCC Sum

      qgscnum21e0 = tFjet79*tcoff79+tFjet80*tcoff80+tFjet81*tcoff81+tFje
     #t82*tcoff82
      RETURN
      END

      double precision function qgscnum21e1(var,ndim,z,ux,l1,
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
      t20 = zp**2
      t23 = l1**2
      taux40 = (t1/4+t3/4)*(t9*s1l1+t12*s1l2+t9*s2l1+t12*s2l2+128*t11*t8
     #+64*zp*(-4*zp*l1-t20+2*t20*l1-2*t20*t23-1+4*zp*t23-2*t23+2*l1))
      t1 = yp**2
      t3 = (1-yp)**2
      t7 = -3-3*zp
      t9 = 4*l1
      t11 = 4*zp*l1
      t13 = t9-1-t11+4*zp
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = t11+3-t9
      taux41 = (t1/4+t3/4)*((256*t7*s2l2-256*t13*t14)*s1l1+(-256*t7*s2l1
     #+256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)

CCC Coefficients

      tcoff79 = taux40
      tcoff80 = taux40
      tcoff81 = taux41
      tcoff82 = taux41

CCC Jet functions

      tFjet79 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f79)
      tFjet80 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off80)
      tFjet81 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f81)
      tFjet82 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off82)

CCC Sum

      qgscnum21e1 = tFjet79*tcoff79+tFjet80*tcoff80+tFjet81*tcoff81+tFje
     #t82*tcoff82
      RETURN
      END

      double precision function qgscnum21e2(var,ndim,z,ux,l1,
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
      taux40 = 64*(t1/4+t3/4)*zp*(1+t7-2*zp)
      t1 = yp**2
      t3 = (1-yp)**2
      t6 = 1+zp
      t10 = 4*l1
      t12 = 4*zp*l1
      t13 = -5*zp+3-t10+t12
      t14 = -s1l2-s2l2-zp
      t21 = -s1l1-s2l1-zp
      t26 = -zp-1-t12+t10
      taux41 = (t1/4+t3/4)*((512*t6*s2l2-256*t13*t14)*s1l1+(-512*t6*s2l1
     #+256*t13*t21)*s1l2+256*t26*t14*s2l1-256*t26*t21*s2l2)

CCC Coefficients

      tcoff79 = taux40
      tcoff80 = taux40
      tcoff81 = taux41
      tcoff82 = taux41

CCC Jet functions

      tFjet79 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f79)
      tFjet80 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off80)
      tFjet81 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f81)
      tFjet82 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off82)

CCC Sum

      qgscnum21e2 = tFjet79*tcoff79+tFjet80*tcoff80+tFjet81*tcoff81+tFje
     #t82*tcoff82
      RETURN
      END

      double precision function qgscnum21e3(var,ndim,z,ux,l1,
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

      taux40 = 0
      t1 = yp**2
      t3 = (1-yp)**2
      t6 = zp-1
      t8 = t6*(-s1l2-s2l2-zp)
      t11 = t6*(-s1l1-s2l1-zp)
      taux41 = (t1/4+t3/4)*(-512*t8*s1l1+512*t11*s1l2+512*t8*s2l1-512*t1
     #1*s2l2)

CCC Coefficients

      tcoff79 = taux40
      tcoff80 = taux40
      tcoff81 = taux41
      tcoff82 = taux41

CCC Jet functions

      tFjet79 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tcof
     #f79)
      tFjet80 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,yp,Tcoeff*tc
     #off80)
      tFjet81 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tcof
     #f81)
      tFjet82 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0
     #,zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,yp,Tcoeff*tc
     #off82)

CCC Sum

      qgscnum21e3 = tFjet79*tcoff79+tFjet80*tcoff80+tFjet81*tcoff81+tFje
     #t82*tcoff82
      RETURN
      END

