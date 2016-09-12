CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCC Numerator functions CCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qNLORnum1e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v
      external fillspin
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13 = -(1.0d0-l1)*(1.0d0-z)
      s23 = -l1*(1.0d0-z)
      s1v = s23
      s2v = s13
      call fillspin(z,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -1+z
      t6 = 64+64*z
      t8 = z*l1
      t10 = 64*l1+64*z-64*t8
      t23 = -64*l1+64*t8+64
      taux1 = (-128*t1*l1*s1l2+t6*s2l2+t10*s3l2)*s1l1+(t6*s2l1+t10*s3l1)
     #*s1l2+(128*t1*(-1+l1)*s2l2+t23*s3l2)*s2l1+t23*s3l1*s2l2
      t1 = 1+z
      t3 = z*l1
      t4 = -z+t3-l1
      t14 = -1-t3+l1
      taux2 = (256*t1*s2l2-256*t4*s3l2)*s1l1+(-256*t1*s2l1+256*t4*s3l1)*
     #s1l2+256*t14*s3l2*s2l1-256*t14*s3l1*s2l2

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f1)
      tFjet2 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f2)

CCC Sum

      qNLORnum1e0 = tFjet1*tcoff1+tFjet2*tcoff2
      RETURN
      END

      double precision function qNLORnum1e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v
      external fillspin
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13 = -(1.0d0-l1)*(1.0d0-z)
      s23 = -l1*(1.0d0-z)
      s1v = s23
      s2v = s13
      call fillspin(z,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = 64-64*z
      t3 = t2*s3l2
      t5 = t2*s3l1
      t11 = z**2
      t14 = l1**2
      taux1 = t3*s1l1+t5*s1l2+t3*s2l1+t5*s2l2+128*s3l1*s3l2+64*z*(2*t11*
     #l1+4*z*t14-t11-2*t11*t14-4*z*l1-1-2*t14+2*l1)
      t2 = -3-3*z
      t4 = 4*l1
      t6 = 4*z*l1
      t8 = -1+t4-t6+4*z
      t18 = t6+3-t4
      taux2 = (256*t2*s2l2-256*t8*s3l2)*s1l1+(-256*t2*s2l1+256*t8*s3l1)*
     #s1l2+256*t18*s3l2*s2l1-256*t18*s3l1*s2l2

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f1)
      tFjet2 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f2)

CCC Sum

      qNLORnum1e1 = tFjet1*tcoff1+tFjet2*tcoff2
      RETURN
      END

      double precision function qNLORnum1e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v
      external fillspin
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13 = -(1.0d0-l1)*(1.0d0-z)
      s23 = -l1*(1.0d0-z)
      s1v = s23
      s2v = s13
      call fillspin(z,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      taux1 = 64*z*(1+t1-2*z)
      t1 = 1+z
      t5 = 4*l1
      t7 = 4*z*l1
      t8 = 3-5*z-t5+t7
      t19 = -1-z-t7+t5
      taux2 = (512*t1*s2l2-256*t8*s3l2)*s1l1+(-512*t1*s2l1+256*t8*s3l1)*
     #s1l2+256*t19*s3l2*s2l1-256*t19*s3l1*s2l2

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f1)
      tFjet2 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f2)

CCC Sum

      qNLORnum1e2 = tFjet1*tcoff1+tFjet2*tcoff2
      RETURN
      END

      double precision function qNLORnum2e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v
      external fillspin
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13 = -(1.0d0-l1)*(1.0d0-z)
      s23 = -l1*(1.0d0-z)
      s1v = s23
      s2v = s13
      call fillspin(z,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -1+z
      t5 = z*l1
      t7 = -64*l1-64*z+64*t5
      t12 = 64-64*l1+64*t5-128*z
      t21 = -64*t1
      taux3 = (128*t1*l1*s1l2+t7*s2l2+t12*s3l2)*s1l1+(t7*s2l1+t12*s3l1)*
     #s1l2+t21*s3l2*s2l1+t21*s3l1*s2l2+128*s3l1*s3l2
      t1 = z*l1
      t3 = -256*l1-256*z+256*t1
      t6 = 1-l1+t1-2*z
      t18 = 256-256*z
      taux4 = (t3*s2l2+256*t6*s3l2)*s1l1+(-t3*s2l1-256*t6*s3l1)*s1l2+t18
     #*s3l2*s2l1-t18*s3l1*s2l2

CCC Coefficients

      tcoff3 = taux3
      tcoff4 = taux3
      tcoff5 = taux4
      tcoff6 = taux4

CCC Jet functions

      tFjet3 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tcof
     #f3)
      tFjet4 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,
     #z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tc
     #off4)
      tFjet5 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tcof
     #f5)
      tFjet6 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,
     #z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tc
     #off6)

CCC Sum

      qNLORnum2e0 = tFjet3*tcoff3+tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tco
     #ff6
      RETURN
      END

      double precision function qNLORnum2e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v
      external fillspin
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13 = -(1.0d0-l1)*(1.0d0-z)
      s23 = -l1*(1.0d0-z)
      s1v = s23
      s2v = s13
      call fillspin(z,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z*l1
      t3 = -64*l1+64*t1+64
      t18 = z**2
      t21 = l1**2
      taux3 = t3*s2l2*s1l1+t3*s2l1*s1l2+(128*(-1+z)*(-1+l1)*s2l2+t3*s3l2
     #)*s2l1+t3*s3l1*s2l2-64*z*(-1+2*t18*l1-t18*t21-2*t18-2*t1+2*z+2*z*t
     #21-t21)
      t2 = 768*z
      t3 = z*l1
      t5 = 256+768*l1+t2-768*t3
      t10 = -3+3*l1-3*t3+6*z
      t23 = -256*l1+256*t3-768+t2
      taux4 = (t5*s2l2+256*t10*s3l2)*s1l1+(-t5*s2l1-256*t10*s3l1)*s1l2+t
     #23*s3l2*s2l1-t23*s3l1*s2l2

CCC Coefficients

      tcoff3 = taux3
      tcoff4 = taux3
      tcoff5 = taux4
      tcoff6 = taux4

CCC Jet functions

      tFjet3 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tcof
     #f3)
      tFjet4 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,
     #z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tc
     #off4)
      tFjet5 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tcof
     #f5)
      tFjet6 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,
     #z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tc
     #off6)

CCC Sum

      qNLORnum2e1 = tFjet3*tcoff3+tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tco
     #ff6
      RETURN
      END

      double precision function qNLORnum2e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v
      external fillspin
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13 = -(1.0d0-l1)*(1.0d0-z)
      s23 = -l1*(1.0d0-z)
      s1v = s23
      s2v = s13
      call fillspin(z,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = l1**2
      taux3 = -64*z*(t1*t2+2*z*l1-2*z*t2+1+t2-2*l1)
      t2 = 512*z
      t3 = z*l1
      t5 = -768-256*l1-t2+256*t3
      t8 = 1-l1+t3-2*z
      t21 = 768*l1-768*t3+256-t2
      taux4 = (t5*s2l2+512*t8*s3l2)*s1l1+(-t5*s2l1-512*t8*s3l1)*s1l2+t21
     #*s3l2*s2l1-t21*s3l1*s2l2

CCC Coefficients

      tcoff3 = taux3
      tcoff4 = taux3
      tcoff5 = taux4
      tcoff6 = taux4

CCC Jet functions

      tFjet3 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tcof
     #f3)
      tFjet4 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,
     #z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*tc
     #off4)
      tFjet5 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tcof
     #f5)
      tFjet6 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,
     #z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*tc
     #off6)

CCC Sum

      qNLORnum2e2 = tFjet3*tcoff3+tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tco
     #ff6
      RETURN
      END

