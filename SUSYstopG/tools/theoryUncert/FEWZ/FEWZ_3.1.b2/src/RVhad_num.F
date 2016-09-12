CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCC Numerator functions CCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RVnum1e0(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t1 = z*s2l1
      t4 = z**2
      t5 = s1l1*t4
      t15 = s1l1**2
      t20 = l1**2
      t25 = t4*z
      t30 = s2l1**2
      t38 = 384*t1*l1+384*t5*l1-384*t4*s2l1*l1-192*z-384*z*s1l1*l1-384*t
     #15*z-768*t4*l1+768*t4*t20-384*z*t20+384*t25*l1-384*t1-384*t5-384*t
     #30*z-192*t25-384*t25*t20+384*z*l1
      taux1 = t38/(1-l3*l1*(1-z))
      t1 = z**2
      t2 = t1*z
      t8 = z*s2l1
      t12 = s1l1*t1
      taux2 = -384*(4*t2*l1-4*z*s1l1*l1-4*t8*l1+4*t8+4*t12*l1+2*z-4*z*l1
     #-4*t12+4*t1*s2l1*l1-2*t2)/(1-l3*l1*(1-z))

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = -taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f1)
      tFjet2 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f2)
      tFjet3 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f3)
      tFjet4 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f4)

CCC Sum

      RVnum1e0 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff4
      RETURN
      END

      double precision function RVnum1e1(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t2 = z**2
      t4 = s1l1*s2l1
      t5 = t2*z
      t6 = t5*l1
      t9 = l1**2
      t10 = t5*t9
      t13 = z*t9
      t16 = t2*t9
      t19 = z*l1
      t22 = t2*l1
      t29 = z*s2l1
      t31 = 768*z-384*t2-768*t4*t6+768*t4*t10+1536*t4*t13-1920*t4*t16-15
     #36*t4*t19+1920*t4*t22+3456*t22-3456*t16+1728*t13-1728*t6+1152*t29
      t32 = s1l1*t2
      t34 = s2l1**2
      t35 = t34*z
      t41 = s2l1*t5
      t43 = s1l1**2
      t46 = t34*t2
      t50 = t34*t5
      t52 = z*s1l1
      t54 = s2l1*t2
      t58 = 1344*t32+1152*t35+1728*t10-1728*t19-192*t34*t9+576*t41-192*t
     #43*t9-576*t46+384*t43*l1+384*t50-192*t52-576*t54-192*t43+960*t5
      t62 = t43*t2
      t69 = t43*z
      t88 = 384*t4*t2+384*t62*l1+384*t52*t9+768*t35*t9+768*t69*t9+1536*t
     #46*l1-960*t46*t9-768*t35*l1+384*t4*l1+384*t41*t9+384*t29*t9-768*t5
     #4*t9-768*t50*l1
      t102 = s1l1*t5
      t118 = -384*t4*t9+384*t43*t5*t9-960*t62*t9+384*t50*t9-768*t69*l1-9
     #60*t41*l1+192*t102*l1-768*t32*t9+384*t102*t9-1728*t29*l1-1152*t32*
     #l1+2688*t54*l1+960*t52*l1+1152*t69
      taux1 = (t31+t58+t88+t118)/(1-l3*l1*(1-z))
      t1 = z*s2l1
      t2 = l1**2
      t3 = t2*l1
      t6 = z**2
      t7 = t6*z
      t8 = s1l1*t7
      t11 = z*s1l1
      t14 = s2l1*t7
      t17 = s2l1*t6
      t39 = 18*t1*t3+6*t8*t3+18*t11*t3+9*t14*l1+30*t17*t2-6*s1l1*t3-6*s2
     #l1*t3-8*s1l1*l1-3*t14-3*t11+5*t17-24*t1*t2-12*t14*t2-30*t11*t2+2*s
     #1l1+9*t7+t6+3*t8*l1
      t40 = t6*s1l1
      t74 = 24*t40*t2-6*t8*t2+31*t1*l1-32*t40*l1-38*t17*l1+37*t11*l1-2*t
     #6*l1-18*t7*l1-22*t1+21*t40+20*z*l1+6*s2l1*t2+12*s1l1*t2-2*s2l1*l1-
     #10*z-18*t40*t3-18*t17*t3+6*t14*t3
      taux2 = -384*(t39+t74)/(1-l3*l1*(1-z))

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = -taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f1)
      tFjet2 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f2)
      tFjet3 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f3)
      tFjet4 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f4)

CCC Sum

      RVnum1e1 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff4
      RETURN
      END

      double precision function RVnum1e2(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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
      t2 = t1*z
      t5 = s1l1*s2l1
      t6 = l1**2
      t7 = t2*t6
      t10 = t6*z
      t13 = t1*t6
      t16 = z*l1
      t19 = t1*l1
      t23 = t2*l1
      t25 = s1l1*t1
      t27 = s2l1**2
      t28 = t27*z
      t32 = -960*t2+192*t1-768*t5*t7+1536*t5*t10+384*t5*t13-1536*t5*t16-
     #384*t5*t19-960*z+1728*t23-384*t25-768*t28-1728*t7+1728*t16
      t35 = s2l1*t2
      t37 = s1l1**2
      t40 = t27*t1
      t44 = t27*t2
      t46 = z*s1l1
      t48 = s2l1*t1
      t58 = -576*t27*t6-384*t35-576*t37*t6+192*t40+1152*t37*l1-384*t44-1
     #536*t46-768*t48+768*t5*t23-1728*t10-576*t37-1152*t48*l1-1920*t5*z
      t62 = s1l1*t2
      t73 = t37*z
      t90 = 1536*t46*l1+384*t62*l1+1536*t25*t6-768*t62*t6-1920*t25*l1+11
     #52*t35*l1-1152*t73*l1-384*t44*t6+192*t37*t1*t6-384*t37*t2*t6-1152*
     #t5*t6+768*t44*l1+1536*t48*t6
      t91 = z*s2l1
      t114 = -768*t91*t6-768*t35*t6+1152*t5*l1-384*t28*l1+768*t73*t6+192
     #*t40*t6-384*t40*l1+768*t28*t6-768*t46*t6+3456*t13-3456*t19-384*t73
     #-768*t91
      taux1 = (t32+t58+t90+t114)/(1-l3*l1*(1-z))
      t3 = l1**2
      t10 = t3*l1
      t17 = z**2
      t18 = t17*z
      t19 = s2l1*t18
      t21 = z*s1l1
      t23 = s2l1*t17
      t25 = z*s2l1
      t28 = s1l1*t18
      t43 = -30*z*l1-22*s2l1*t3-44*s1l1*t3+11*s2l1*l1+22*s1l1*t10+22*s2l
     #1*t10+33*s1l1*l1+13*t19+18*t21-22*t23-66*t25*t10-22*t28*t10-66*t21
     #*t10-35*t19*l1-110*t23*t3+88*t25*t3+44*t19*t3+110*t21*t3
      t50 = s1l1*t17
      t75 = -2*t17-11*s1l1+15*z-13*t18-81*t25*l1+83*t50*l1+105*t23*l1-10
     #3*t21*l1+66*t50*t10+66*t23*t10-22*t19*t10-39*t50+4*t17*l1+26*t18*l
     #1+41*t25-13*t28*l1-88*t50*t3+22*t28*t3
      taux2 = -384*(t43+t75)/(1-l3*l1*(1-z))

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = -taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f1)
      tFjet2 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f2)
      tFjet3 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f3)
      tFjet4 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f4)

CCC Sum

      RVnum1e2 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff4
      RETURN
      END

      double precision function RVnum1e3(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t2 = z**2
      t4 = t2*z
      t10 = l1**2
      taux1 = (384*z-768*t2+384*t4-768*t2*l1+384*t4*l1-384*z*t10+768*t2*
     #t10-384*t4*t10+384*z*l1)/(1-l3*l1*(1-z))
      t3 = z**2
      t4 = s2l1*t3
      t6 = t3*z
      t9 = s1l1*t3
      t13 = z*s2l1
      t15 = l1**2
      t18 = z*s1l1
      t22 = t6*s1l1
      t35 = s2l1*t6
      t38 = t15*l1
      t43 = 8*z*l1+28*t4-16*t6*l1+32*t9+8*t3*l1-28*t13-64*t13*t15-80*t18
     #*t15-4*z-16*t22*t15+64*t9*t15-92*t4*l1+96*t18*l1+80*t13*l1-76*t9*l
     #1-32*t35*t15+16*t22*t38+48*t18*t38
      t75 = 20*s1l1-36*s1l1*l1-16*s2l1*t38-20*s2l1*l1-16*s1l1*t38+16*s2l
     #1*t15+32*s1l1*t15-16*t35-36*t18-48*t9*t38+16*t35*t38+48*t13*t38+80
     #*t4*t15+32*t35*l1-48*t4*t38-4*t3+8*t6+16*t22*l1
      taux2 = -384*(t43+t75)/(1-l3*l1*(1-z))

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = -taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f1)
      tFjet2 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f2)
      tFjet3 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f3)
      tFjet4 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f4)

CCC Sum

      RVnum1e3 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff4
      RETURN
      END

      double precision function RVnum1e4(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      taux1 = 0
      t1 = z**2
      t2 = s2l1*t1
      t4 = t1*z
      t7 = s1l1*t1
      t11 = z*s2l1
      t17 = l1**2
      t18 = t17*l1
      t29 = s2l1*t4
      t31 = z*s1l1
      t33 = s1l1*t4
      t42 = -4*t2+4*t4*l1-6*t7-8*t1*l1+2*t11+4*z*l1+2*s1l1*l1-4*s2l1*t18
     #+6*s2l1*l1-4*s1l1*t18+4*s2l1*t17+8*s1l1*t17+2*t29+12*t31-2*t33*l1-
     #4*t33*t17+16*t7*t17+2*t2*l1
      t75 = -6*t31*l1-10*t11*l1+6*t7*l1-8*t29*t17+4*t33*t18+12*t31*t18-1
     #2*t7*t18+4*t29*t18+12*t11*t18+20*t2*t17+2*t29*l1-12*t2*t18-16*t11*
     #t17-20*t31*t17+4*t1-2*t4-6*s1l1-2*z
      taux2 = -768*(t42+t75)/(1-l3*l1*(1-z))

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = -taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f1)
      tFjet2 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f2)
      tFjet3 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f3)
      tFjet4 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f4)

CCC Sum

      RVnum1e4 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff4
      RETURN
      END

      double precision function RVnum2e0(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t1 = l1**2
      t4 = s2l1**2
      t7 = z**2
      t10 = s1l1*t7
      t12 = z*s2l1
      t16 = t7*z
      t19 = s1l1**2
      taux3 = -384*z*t1-384*t4*z-768*t7*l1-384*t10-384*t12+768*t7*t1+384
     #*t16*l1-384*t19*z-384*t16*t1-192*z-192*t16+384*z*l1-384*z*s1l1*l1-
     #384*t7*s2l1*l1+384*t12*l1+384*t10*l1
      t4 = z**2
      t5 = t4*z
      t10 = z*s2l1
      t16 = s1l1*t4
      taux4 = 1536*z*s1l1*l1+768*t5-1536*t4*s2l1*l1-1536*t10-1536*t5*l1+
     #1536*z*l1+1536*t16-1536*t16*l1-768*z+1536*t10*l1

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = -taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f5)
      tFjet6 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f6)
      tFjet7 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f7)
      tFjet8 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f8)

CCC Sum

      RVnum2e0 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function RVnum2e1(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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
      t2 = t1*s2l1
      t5 = z*s2l1
      t8 = s1l1*t1
      t11 = z*t1
      t12 = s1l1*t11
      t15 = l1**2
      t20 = s2l1**2
      t21 = t20*t1
      t24 = t20*z
      t27 = s1l1**2
      t28 = t27*t1
      t38 = t27*z
      t41 = 2688*t2*l1-1728*t5*l1-1152*t8*l1+192*t12*l1-768*t8*t15+384*t
     #12*t15+1536*t21*l1+768*t24*t15+384*t28*l1-960*t28*t15-960*t21*t15+
     #384*t27*t11*t15+768*t38*t15
      t42 = s1l1*s2l1
      t45 = t20*t11
      t52 = t11*l1
      t55 = t1*l1
      t58 = z*l1
      t61 = z*t15
      t64 = t1*t15
      t67 = t11*t15
      t73 = 384*t42*l1+384*t45*t15-768*t38*l1-384*t1-192*t27-768*t42*t52
     #+1920*t42*t55-1536*t42*t58+1536*t42*t61-1920*t42*t64+768*t42*t67+7
     #68*z-1728*t52+1152*t38
      t83 = z*s1l1
      t88 = s2l1*t11
      t95 = 1728*t67+384*t27*l1-192*t27*t15-192*t20*t15-1728*t58-192*t83
     #-576*t2-576*t21+384*t45+576*t88+1152*t5-384*t42*t15-768*t45*l1
      t118 = -960*t88*l1+384*t88*t15-768*t2*t15+384*t5*t15-768*t24*l1+38
     #4*t42*t1+384*t83*t15+960*t83*l1+960*t11-3456*t64+1728*t61+1152*t24
     #+3456*t55+1344*t8
      taux3 = t41+t73+t95+t118
      t1 = z**2
      t5 = s1l1*t1
      t7 = z*s2l1
      t9 = t1*z
      t12 = l1**2
      t21 = t12*l1
      t28 = z*s1l1
      t30 = s2l1*t1
      t32 = s2l1*t9
      t37 = -384*t1+768*t1*l1-8064*t5+8448*t7+6912*t9*l1-2304*s2l1*t12-7
     #680*z*l1+3072*s1l1*l1-4608*s1l1*t12+2304*s1l1*t21+2304*s2l1*t21+76
     #8*s2l1*l1+1152*t28-1920*t30+1152*t32-768*s1l1+3840*z-3456*t9
      t56 = s1l1*t9
      t75 = 6912*t5*t21+6912*t30*t21-2304*t32*t21-6912*t7*t21+11520*t28*
     #t12-14208*t28*l1+14592*t30*l1-11904*t7*l1+12288*t5*l1-1152*t56*l1-
     #9216*t5*t12+2304*t56*t12-2304*t56*t21-6912*t28*t21-3456*t32*l1+460
     #8*t32*t12-11520*t30*t12+9216*t7*t12
      taux4 = t37+t75

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = -taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f5)
      tFjet6 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f6)
      tFjet7 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f7)
      tFjet8 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f8)

CCC Sum

      RVnum2e1 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function RVnum2e2(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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
      t2 = t1*z
      t4 = s2l1*t2
      t7 = l1**2
      t10 = s2l1*t1
      t13 = s1l1**2
      t16 = z*s2l1
      t19 = z*s1l1
      t24 = s1l1*s2l1
      t25 = t2*l1
      t28 = t1*l1
      t31 = z*l1
      t34 = t1*t7
      t36 = -960*t2+1152*t4*l1-768*t4*t7+1536*t10*t7-576*t13+192*t1-768*
     #t16*t7-768*t19*t7+1536*t19*l1+768*t24*t25-384*t24*t28-1536*t24*t31
     #+3456*t34
      t40 = t2*t7
      t42 = s2l1**2
      t43 = t2*t42
      t47 = s1l1*t1
      t56 = t42*t1
      t60 = 1728*t31-1536*t19-768*t10-1728*t40-384*t43-1152*t10*l1-1920*
     #t47*l1-576*t42*t7-384*t4+1152*t24*l1-960*z+192*t56-576*t13*t7
      t62 = t13*z
      t64 = t42*z
      t70 = z*t7
      t74 = s1l1*t2
      t85 = -384*t62-768*t64-3456*t28-384*t47-768*t16+1728*t25-1728*t70+
     #1152*t13*l1+384*t74*l1-384*t64*l1+768*t64*t7-384*t56*l1+192*t56*t7
      t114 = 1536*t47*t7-768*t74*t7-768*t24*t40+384*t24*t34+1536*t24*t70
     #-1152*t62*l1-1152*t24*t7-384*t43*t7-1920*t24*z+192*t13*t1*t7-384*t
     #13*t2*t7+768*t62*t7+768*t43*l1
      taux3 = t36+t60+t85+t114
      t1 = z**2
      t2 = t1*z
      t6 = s1l1*t1
      t8 = z*s2l1
      t12 = s2l1*t2
      t15 = l1**2
      t18 = s2l1*t1
      t23 = s1l1*t2
      t24 = t15*l1
      t27 = z*s1l1
      t42 = 4992*t2-1536*t1*l1+14976*t6-15744*t8-9984*t2*l1+13440*t12*l1
     #-16896*t12*t15+42240*t18*t15-33792*t8*t15+8448*t23*t24+25344*t27*t
     #24+768*t1-5760*z+8448*s2l1*t15+11520*z*l1-12672*s1l1*l1+16896*s1l1
     #*t15-8448*s1l1*t24
      t75 = -8448*s2l1*t24-4224*s2l1*l1-6912*t27+8448*t18-4992*t12+4992*
     #t23*l1+33792*t6*t15-8448*t23*t15-25344*t6*t24-25344*t18*t24+8448*t
     #12*t24+25344*t8*t24-42240*t27*t15+39552*t27*l1-40320*t18*l1+31104*
     #t8*l1-31872*t6*l1+4224*s1l1
      taux4 = t42+t75

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = -taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f5)
      tFjet6 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f6)
      tFjet7 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f7)
      tFjet8 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f8)

CCC Sum

      RVnum2e2 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function RVnum2e3(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t2 = z**2
      t6 = t2*z
      t9 = l1**2
      taux3 = 384*z-768*t2-768*t2*l1+384*t6*l1+768*t2*t9-384*z*t9-384*t6
     #*t9+384*z*l1+384*t6
      t1 = l1**2
      t2 = t1*l1
      t7 = z*s2l1
      t10 = z**2
      t11 = s2l1*t10
      t14 = t10*z
      t15 = s2l1*t14
      t31 = s1l1*t10
      t36 = z*s1l1
      t40 = s1l1*t14
      t43 = 6144*s2l1*t2+6144*s1l1*t2+24576*t7*t1-30720*t11*t1+12288*t15
     #*t1-12288*t15*l1-12288*s1l1*t1+13824*s1l1*l1-3072*z*l1-6144*s2l1*t
     #1+6144*t14*l1+10752*t7-12288*t31-3072*t10*l1-10752*t11+13824*t36-1
     #8432*t36*t2-6144*t40*t2
      t75 = 7680*s2l1*l1-18432*t7*t2-6144*t15*t2+18432*t11*t2+18432*t31*
     #t2-30720*t7*l1+35328*t11*l1-36864*t36*l1+30720*t36*t1+6144*t15+614
     #4*t40*t1-24576*t31*t1-6144*t40*l1+29184*t31*l1+1536*z-7680*s1l1-30
     #72*t14+1536*t10
      taux4 = t43+t75

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = -taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f5)
      tFjet6 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f6)
      tFjet7 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f7)
      tFjet8 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f8)

CCC Sum

      RVnum2e3 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function RVnum2e4(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      taux3 = 0
      t2 = z*s1l1
      t3 = l1**2
      t8 = z**2
      t9 = t8*s2l1
      t12 = z*s2l1
      t15 = s1l1*t8
      t18 = t8*z
      t19 = s1l1*t18
      t27 = s2l1*t18
      t42 = 4608*s1l1+15360*t2*t3+4608*t2*l1-1536*t9*l1+7680*t12*l1-4608
     #*t15*l1+1536*t19*l1-12288*t15*t3+3072*t19*t3+1536*t18-1536*t27*l1+
     #6144*t27*t3-15360*t9*t3+12288*t12*t3-3072*t8+1536*z-3072*s2l1*t3-3
     #072*z*l1
      t47 = t3*l1
      t75 = -1536*s1l1*l1-6144*s1l1*t3+3072*s1l1*t47+3072*s2l1*t47-4608*
     #s2l1*l1-9216*t2+3072*t9-1536*t27-3072*t19*t47-9216*t2*t47+9216*t15
     #*t47+9216*t9*t47-3072*t27*t47-9216*t12*t47+6144*t8*l1+4608*t15-153
     #6*t12-3072*t18*l1
      taux4 = t42+t75

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = -taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f5)
      tFjet6 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f6)
      tFjet7 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f7)
      tFjet8 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z,
     #s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tcof
     #f8)

CCC Sum

      RVnum2e4 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function RVnum3e0(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t1 = s1l1**2
      t4 = s2l1**2
      t7 = z**2
      t10 = s1l1*t7
      t12 = z*s2l1
      t14 = t7*z
      t17 = l1**2
      taux5 = -384*t1*z-384*t4*z-768*t7*l1-384*t10-384*t12+384*t14*l1-38
     #4*z*t17+768*t7*t17-192*z-192*t14+384*z*l1-384*z*s1l1*l1-384*t7*s2l
     #1*l1-384*t14*t17+384*t12*l1+384*t10*l1
      t1 = z**2
      t2 = t1*z
      t5 = z*s2l1
      t8 = s1l1*t1
      taux6 = -1536*t2*l1+1536*t5*l1-1536*t8*l1+1536*z*s1l1*l1+768*t2-15
     #36*t1*s2l1*l1+1536*t8-768*z-1536*t5+1536*z*l1

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f9)
      tFjet10 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff10)

CCC Sum

      RVnum3e0 = tFjet9*tcoff9+tFjet10*tcoff10
      RETURN
      END

      double precision function RVnum3e1(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t1 = s1l1**2
      t2 = z**2
      t3 = t1*t2
      t4 = l1**2
      t7 = t2*z
      t11 = s1l1*s2l1
      t16 = t2*s2l1
      t19 = z*s2l1
      t24 = t7*t4
      t26 = s2l1**2
      t27 = t26*t2
      t29 = t7*t26
      t33 = -960*t3*t4+384*t1*t7*t4+384*t11*t2+768*z-192*t1+2688*t16*l1-
     #1728*t19*l1-384*t2-576*t16+1728*t24-576*t27+384*t29+384*t1*l1
      t36 = s2l1*t7
      t44 = t1*z
      t52 = s1l1*t2
      t57 = s1l1*t7
      t61 = t26*z
      t63 = -192*t1*t4+576*t36-960*t27*t4-768*t29*l1+384*t29*t4-768*t44*
     #l1-384*t11*t4+768*t44*t4+960*t7-1152*t52*l1+384*t11*l1+192*t57*l1+
     #1152*t44+1152*t61
      t65 = t2*l1
      t69 = t7*l1
      t71 = z*t4
      t81 = z*l1
      t88 = z*s1l1
      t91 = 3456*t65+1344*t52+1152*t19-1728*t69+1728*t71-960*t36*l1+384*
     #t36*t4-768*t16*t4+1920*t11*t65-1536*t11*t81+384*t19*t4+384*t3*l1+3
     #84*t88*t4
      t94 = t2*t4
      t118 = 960*t88*l1-3456*t94-1728*t81-192*t26*t4-192*t88-768*t52*t4+
     #384*t57*t4-768*t61*l1+768*t61*t4+1536*t27*l1-768*t11*t69+768*t11*t
     #24-1920*t11*t94+1536*t11*t71
      taux5 = t33+t63+t91+t118
      t1 = z**2
      t2 = t1*z
      t3 = t2*s2l1
      t6 = l1**2
      t9 = s2l1*t1
      t12 = z*s2l1
      t15 = s1l1*t2
      t16 = t6*l1
      t19 = z*s1l1
      t22 = s1l1*t1
      t44 = -3456*t3*l1+4608*t3*t6-11520*t9*t6+9216*t12*t6-2304*t15*t16-
     #6912*t19*t16+6912*t22*t16+6912*t9*t16-2304*t3*t16+3840*z-768*s1l1-
     #3456*t2-6912*t12*t16+11520*t19*t6-14208*t19*l1+14592*t9*l1-11904*t
     #12*l1+12288*t22*l1
      t75 = -1152*t15*l1-9216*t22*t6-384*t1+2304*t15*t6+768*t1*l1-8064*t
     #22+8448*t12+6912*t2*l1-2304*s2l1*t6-7680*z*l1+3072*s1l1*l1-4608*s1
     #l1*t6+2304*s1l1*t16+2304*s2l1*t16+768*s2l1*l1+1152*t19-1920*t9+115
     #2*t3
      taux6 = t44+t75

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f9)
      tFjet10 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff10)

CCC Sum

      RVnum3e1 = tFjet9*tcoff9+tFjet10*tcoff10
      RETURN
      END

      double precision function RVnum3e2(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t1 = s1l1*s2l1
      t2 = z**2
      t3 = t2*l1
      t6 = z*l1
      t9 = t2*z
      t10 = t9*l1
      t13 = l1**2
      t14 = t9*t13
      t17 = t2*t13
      t20 = z*t13
      t24 = s1l1**2
      t28 = s2l1**2
      t31 = z*s1l1
      t33 = t2*s2l1
      t35 = -384*t1*t3-1536*t1*t6+768*t1*t10-768*t1*t14+384*t1*t17+1536*
     #t1*t20+192*t2-576*t24-960*z-960*t9-576*t28*t13-1536*t31-768*t33
      t37 = t28*t2
      t39 = t9*t28
      t45 = s2l1*t9
      t50 = t24*z
      t52 = t28*z
      t55 = s1l1*t2
      t57 = -1728*t14+192*t37-384*t39+1152*t24*l1-576*t24*t13-384*t45+17
     #28*t6+3456*t17-1728*t20-384*t50-768*t52-3456*t3-384*t55
      t59 = z*s2l1
      t80 = s1l1*t9
      t85 = -768*t59+1728*t10+1152*t45*l1-768*t45*t13+1536*t33*t13-768*t
     #59*t13-768*t31*t13+1536*t31*l1-1152*t33*l1-1920*t55*l1+1152*t1*l1+
     #384*t80*l1+1536*t55*t13
      t114 = -768*t80*t13-384*t52*l1+768*t52*t13-384*t37*l1+192*t37*t13+
     #768*t39*l1-384*t39*t13-1920*t1*z-1152*t50*l1-1152*t1*t13+768*t50*t
     #13+192*t24*t2*t13-384*t24*t9*t13
      taux5 = t35+t57+t85+t114
      t1 = l1**2
      t2 = t1*l1
      t7 = z*s2l1
      t10 = z**2
      t11 = s2l1*t10
      t14 = t10*z
      t15 = s2l1*t14
      t29 = s1l1*t10
      t33 = z*s1l1
      t41 = s1l1*t14
      t44 = -8448*s1l1*t2+16896*s1l1*t1-33792*t7*t1+42240*t11*t1-16896*t
     #15*t1+13440*t15*l1-12672*s1l1*l1+11520*z*l1+8448*s2l1*t1-9984*t14*
     #l1-15744*t7+14976*t29-1536*t10*l1-6912*t33-25344*t11*t2-25344*t29*
     #t2+25344*t33*t2+8448*t41*t2
      t75 = -4224*s2l1*l1-8448*s2l1*t2+25344*t7*t2+8448*t15*t2+33792*t29
     #*t1+4992*t41*l1-31872*t29*l1+31104*t7*l1-40320*t11*l1+39552*t33*l1
     #-4992*t15-42240*t33*t1+8448*t11-8448*t41*t1-5760*z+4224*s1l1+4992*
     #t14+768*t10
      taux6 = t44+t75

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f9)
      tFjet10 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff10)

CCC Sum

      RVnum3e2 = tFjet9*tcoff9+tFjet10*tcoff10
      RETURN
      END

      double precision function RVnum3e3(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t2 = z**2
      t6 = t2*z
      t9 = l1**2
      taux5 = 384*z-768*t2-768*t2*l1+384*t6*l1+768*t2*t9-384*z*t9-384*t6
     #*t9+384*z*l1+384*t6
      t2 = z**2
      t3 = t2*z
      t5 = s2l1*t3
      t8 = l1**2
      t11 = s2l1*t2
      t14 = z*s2l1
      t17 = s1l1*t3
      t18 = t8*l1
      t21 = z*s1l1
      t28 = s1l1*t2
      t41 = -7680*s1l1-3072*t3-12288*t5*l1+12288*t5*t8-30720*t11*t8+2457
     #6*t14*t8-6144*t17*t18-18432*t21*t18-3072*t2*l1+1536*t2+1536*z-1228
     #8*t28+10752*t14+6144*t3*l1-6144*s2l1*t8-3072*z*l1+13824*s1l1*l1-12
     #288*s1l1*t8
      t75 = 6144*s1l1*t18+6144*s2l1*t18+7680*s2l1*l1+13824*t21-10752*t11
     #+6144*t5+35328*t11*l1-30720*t14*l1+29184*t28*l1-6144*t17*l1-24576*
     #t28*t8+6144*t17*t8+18432*t28*t18+18432*t11*t18-6144*t5*t18-18432*t
     #14*t18+30720*t21*t8-36864*t21*l1
      taux6 = t41+t75

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f9)
      tFjet10 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff10)

CCC Sum

      RVnum3e3 = tFjet9*tcoff9+tFjet10*tcoff10
      RETURN
      END

      double precision function RVnum3e4(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      taux5 = 0
      t1 = z**2
      t2 = t1*z
      t3 = s2l1*t2
      t7 = s1l1*t1
      t9 = z*s2l1
      t13 = l1**2
      t22 = t13*l1
      t29 = z*s1l1
      t31 = s2l1*t1
      t37 = -1536*t3+6144*t1*l1+4608*t7-1536*t9-3072*t2*l1-3072*s2l1*t13
     #-3072*z*l1-1536*s1l1*l1-6144*s1l1*t13+3072*s1l1*t22+3072*s2l1*t22-
     #4608*s2l1*l1-9216*t29+3072*t31+1536*t2-3072*t1+1536*z+4608*s1l1
      t46 = s1l1*t2
      t75 = -1536*t3*l1+6144*t3*t13-15360*t31*t13+12288*t9*t13-3072*t46*
     #t22-9216*t29*t22+9216*t7*t22+9216*t31*t22-3072*t3*t22-9216*t9*t22+
     #15360*t29*t13+4608*t29*l1-1536*t31*l1+7680*t9*l1-4608*t7*l1+1536*t
     #46*l1-12288*t7*t13+3072*t46*t13
      taux6 = t37+t75

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z,
     #s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tcof
     #f9)
      tFjet10 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff10)

CCC Sum

      RVnum3e4 = tFjet9*tcoff9+tFjet10*tcoff10
      RETURN
      END

      double precision function RVnum4e0(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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
      t2 = t1**2
      t3 = t2*l1
      t4 = -1+z
      t5 = z*l1
      t7 = (1-l1+t5)**2
      t8 = 1/t7
      t9 = t4*t8
      t10 = l6**2
      t13 = 1024*t3*t9*t10
      t16 = 1024*t3*t9*l6
      t17 = t1*z
      t18 = t17*l1
      t20 = 2*z
      t21 = 3*l1
      t22 = l1**2
      t23 = t17*t22
      t25 = t2*t22
      t26 = t1*l1
      t27 = t1*t22
      t29 = z*t22
      t32 = t22*l1
      t33 = t17*t32
      t35 = t2*t32
      t36 = t1*t32
      t38 = z*t32
      t40 = 6*t18-1-t1-t20+t21+2*t23-t25-t26-4*t27+6*t29-3*t22-4*t33+t35
     #+6*t36-4*t38+t32
      t45 = -t13+t16-64*t40*t4*z*t8
      t46 = l5**2
      t48 = dsqrt(z)
      t51 = dsqrt(1-l5)
      t52 = t48*t1*t51
      t53 = dsqrt(l1)
      t55 = t26-2*t5-z+l1-1
      t58 = 1-l1
      t59 = dsqrt(t58)
      t60 = t59*t4
      t62 = t60*t8*l6
      t70 = -512*t52*t53*t55*t62+256*t52*t53*t55*t59*t9
      t71 = dsqrt(l5)
      t72 = t71*l5
      t82 = t13-t16-128*(-3*t26+t5+z-2*l1+1+t27-2*t29+t22)*t4*t1*t8
      t86 = t48*t17*t51*t53
      t92 = -512*t86*t62+256*t86*t60*t8
      t104 = t2*z
      t105 = t104*t22
      t107 = t104*t32
      t109 = -2304*t18+3264*t29-7680*t27+1856*t26+9024*t23-5760*t38+1267
     #2*t36+9216*t35-14976*t33-5184*t25+1152*t105-2304*t107
      t110 = t22**2
      t111 = t104*t110
      t112 = 1152*t111
      t113 = t2*t110
      t114 = 5184*t113
      t115 = z*t110
      t116 = 3456*t115
      t117 = t17*t110
      t118 = 9216*t117
      t119 = t1*t110
      t120 = 8064*t119
      t124 = 576*t110
      t125 = 1152*t32
      t127 = 640*t1
      t128 = t112-t114+t116+t118-t120+1152*t3-256*z-704*t5-t124+t125-576
     #*t22-t127
      t130 = s1l1**2
      t134 = 20736*t25
      t140 = 6912*t105
      t141 = t104*l1
      t143 = t127+2304*t111+19712*t23-t134-5504*t26-4864*t27-1024*t29-10
     #368*t113-6912*t107+t140-2304*t141+t125
      t157 = -39168*t33+26496*t35+27648*t36-9216*t38-1152*t110+3456*t3-1
     #6128*t119+18432*t117+6912*t115+1152*t2+1024*t5+3328*t18-2304*t17
      t166 = 39*t36
      t171 = 18*t113
      t172 = 72*t115
      t173 = 72*t117
      t174 = 108*t119
      t178 = 18*t110
      t182 = t171-t172-t173+t174+9*t3-65*t5+9*t17+t178+14*t32-36*t22-27*
     #t1
      t188 = 64*t1
      t197 = t188+t112+21568*t23-t134+1600*t26-8576*t27+832*t29-t114-460
     #8*t107+t140-4608*t141+1152*t104-24192*t33
      t207 = 17280*t35+14976*t36-3456*t38-t124+10368*t3-t120+t118+t116+6
     #4*z-1728*t2+64*t5-7424*t18+832*t17
      t209 = s2l1**2
      t221 = 255*t18-217*t29+509*t27-194*t26-447*t23+279*t38-441*t36-81*
     #t35+309*t33+135*t25+t171
      t222 = 99*t3
      t229 = -t172-t173+t174-t222+t20+27*t2+38*t5-45*t17+t178-66*t32+20*
     #t22+18*t1
      t235 = t110*t22
      t248 = t110*l1
      t249 = z*t248
      t253 = t2*t1
      t261 = t104*t248
      t266 = -40*t235-767*t249-425*t107+126*t105-9*t253*t22-103*t253*t11
     #0-18*t141+240*t235*z-667*t261+76*t32-40*t253*t235
      t279 = -2080*t33+1368*t35+1607*t36-591*t38-600*t235*t2-156*t110+t2
     #22-2367*t119+3028*t117+956*t115-9*t2
      t282 = t1*t248
      t284 = t17*t248
      t288 = t2*t248
      t298 = -20*t5-174*t18+1855*t282-2390*t284+107*t253*t248+1730*t288+
     #800*t17*t235+240*t235*t104+132*t248+45*t253*t32+17*t17
      t306 = (-l1-z+t5)**2
      t329 = 5*t5
      t331 = 2*t1
      t356 = -1-t331+3*t111+22*t23-9*t25+8*t26-25*t27+15*t29-15*t113-4*t
     #22-5*t249-3*t107+t105-t261+4*t32-37*t33
      t369 = 17*t35+t166-20*t38-3*t110+2*t3-30*t119+30*t117+15*t115+t21-
     #t329-8*t18+10*t282-10*t284+5*t288+t248+2*t17
      t377 = s13e0(-t58*t4)
      taux7 = (t45*t46+t70*t72+t82*l5+t92*t71+(t109+t128)*t130+((t143+t1
     #57)*s2l1+64*z*(-36*t18+104*t29-100*t27+88*t26+32*t23-15*t38-t166-2
     #7*t35+67*t33+4*l1+t182))*s1l1+(t197+t207)*t209+64*z*(t221+t229)*s2
     #l1-64*z*(-4*t1-600*t235*t1+744*t111+818*t23-501*t25+113*t26-604*t2
     #7+182*t29-2102*t113-12*t22+t266+t279+t298)*t8)/t306-(-t45*t46-t70*
     #t72-t82*l5-t92*t71+(192*t26+192*z-192*t5-t188)*t130+((128*t26+128*
     #z-128*t5)*s2l1-64*z*(-1-4*t29+2*t27-4*t26-l1-4*z+t329+2*t22+t331))
     #*s1l1+(64*t26+64*z-64*t5-64*t1)*t209+64*z*(2-l1-z+t5)*s2l1-64*z*(t
     #356+t369)*t8)*t377
      t1 = l1**2
      t2 = t1**2
      t3 = t2*l1
      t5 = t1*l1
      t7 = z**2
      t8 = t7*z
      t9 = t8*l1
      t11 = t1*z
      t13 = t7*t1
      t15 = t7*l1
      t17 = t8*t1
      t19 = z*t5
      t21 = t7*t5
      t23 = t7**2
      t24 = t23*t5
      t26 = t8*t5
      t28 = t23*t1
      t30 = t23*z
      t37 = 6912*t3+6912*t5+6784*t9+8064*t11-32128*t13-4352*t15+53888*t1
     #7-50560*t19+134400*t21+102400*t24-168960*t26-43648*t28+13824*t30*t
     #1-24192*t30*t5+20736*t30*t2
      t38 = t23*t2
      t40 = t2*z
      t42 = t8*t2
      t46 = t7*t2
      t58 = z*l1
      t64 = -96768*t38+76032*t40+179712*t42-3456*t30*l1-165888*t46+69120
     #*t7*t3-69120*t8*t3+34560*t23*t3-6912*t30*t3-34560*z*t3+1024*t58+34
     #56*t23-13824*t2-3456*t8+2048*t7
      t67 = -l1-z+t58
      t83 = t23*l1
      t87 = 27*t8
      t91 = -54*t38+216*t40+216*t42-324*t46+108*t83-27*t23-40*t58+t87-54
     #*t2+54*t5-16*t7
      t111 = -81*t9+180*t11-210*t13+150*t15+80*t17-78*t19+78*t21-26*t26+
     #24*l1+16*z-93*t58+t87+26*t5-50*t1-27*t7
      t115 = t67**2
      t118 = z*s1l1
      t121 = z*s2l1
      t124 = s1l1*t7
      t127 = t7*s2l1
      t132 = s1l1*t8
      t145 = -6*t118*l1-2*t121*l1+14*t124*l1+2*t127*l1+4*t118*t1-8*t132*
     #l1+4*t132*t1-8*t124*t1-6*t83+4*t132+2*t118-2*t127+4*t121-4*t124
      t160 = 10*t9-2*t11+10*t13-2*t15-14*t17+2*t19-6*t21-2*t24+6*t26+6*t
     #28-2*t58+2*t23+2*z-2*t8
      t165 = s13e0((1-l1)*(1-z))
      taux8 = ((t37+t64)*s1l1+128*t67*(-204*t9+123*t11-435*t13+136*t15+5
     #01*t17-324*t19+648*t21+162*t24-540*t26-189*t28+t91)*s2l1-128*z*t67
     #*t111)/t115-128*(t145+t160)*t165

CCC Coefficients

      tcoff11 = taux7
      tcoff12 = taux7
      tcoff13 = taux8
      tcoff14 = -taux8

CCC Jet functions

      tFjet11 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff11)
      tFjet12 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff12)
      tFjet13 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff13)
      tFjet14 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff14)

CCC Sum

      RVnum4e0 = tFjet11*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet14
     #*tcoff14
      RETURN
      END

      double precision function RVnum4e1(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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
      t2 = t1**2
      t3 = t2*l1
      t4 = -1+z
      t5 = z*l1
      t7 = (1-l1+t5)**2
      t8 = 1/t7
      t9 = t4*t8
      t10 = l6**2
      t12 = t3*t9*t10
      t13 = 1024*t12
      t15 = t3*t9*l6
      t16 = 1024*t15
      t17 = t1*z
      t18 = t17*l1
      t20 = 2*z
      t21 = 3*l1
      t22 = l1**2
      t23 = t17*t22
      t25 = t2*t22
      t26 = t1*l1
      t27 = t1*t22
      t29 = z*t22
      t32 = t22*l1
      t33 = t17*t32
      t35 = t2*t32
      t36 = t1*t32
      t38 = z*t32
      t40 = 6*t18-1-t1-t20+t21+2*t23-t25-t26-4*t27+6*t29-3*t22-4*t33+t35
     #+6*t36-4*t38+t32
      t46 = l5**2
      t47 = (t13-t16+64*t40*t4*z*t8)*t46
      t48 = dsqrt(z)
      t51 = dsqrt(1-l5)
      t52 = t48*t1*t51
      t53 = dsqrt(l1)
      t54 = 2*t5
      t55 = t26-t54-z+l1-1
      t58 = 1-l1
      t59 = dsqrt(t58)
      t60 = t59*t4
      t62 = t60*t8*l6
      t63 = t52*t53*t55*t62
      t68 = t52*t53*t55*t59*t9
      t71 = dsqrt(l5)
      t72 = t71*l5
      t73 = (512*t63-256*t68)*t72
      t76 = 2*t29
      t79 = t1*t8
      t83 = (-t13+t16+128*(-3*t26+t5+z-2*l1+1+t27-t76+t22)*t4*t79)*l5
      t86 = t48*t17*t51*t53
      t87 = t86*t62
      t90 = t86*t60*t8
      t93 = (512*t87-256*t90)*t71
      t94 = 256*l1
      t96 = t22**2
      t97 = 1728*t96
      t100 = 512*z
      t107 = t94+512*t17-t97+3712*t32-2240*t22+t100-1600*t1+6976*t23-165
     #76*t27-4096*t5-3456*t18+7296*t26
      t108 = t2*z
      t109 = t108*t32
      t111 = t108*t96
      t112 = 1152*t111
      t113 = t1*t96
      t114 = 5760*t113
      t116 = t108*t22
      t118 = t2*t96
      t119 = 2880*t118
      t120 = z*t96
      t121 = 5760*t120
      t127 = 2304*t109-t112-t114-13504*t38-1152*t116+t119+t121-4544*t35+
     #11328*t29+15936*t36+1664*t25-3904*t33
      t129 = s1l1**2
      t138 = t108*l1
      t142 = 512*t26
      t143 = -640*t17-3456*t96+3968*t32-512*t22-128*t1-10112*t27-1024*t5
     #+2304*t138-3328*t3+1536*t18+t142
      t148 = 6912*t116
      t156 = 6912*t109-2304*t111-11520*t113-17792*t38-t148+5760*t118+115
     #20*t120-14848*t35+5120*t29+22656*t36+12416*t25-896*t33
      t159 = 4*l1
      t161 = 36*t96
      t164 = 8*z
      t170 = t159+8*t17-t161+87*t32-55*t22+t164-31*t1+298*t23-435*t27-92
     #*t5+18*t3
      t171 = t17*t96
      t172 = 144*t171
      t175 = 216*t113
      t177 = 36*t118
      t178 = 144*t120
      t184 = t172-93*t18+163*t26-t175-351*t38-t177+t178+90*t35+264*t29+5
     #31*t36-72*t25-357*t33
      t197 = 64*t1
      t199 = t119+t121-1152*t108-t148-t114+1664*t2-64*z-64*t5+4160*t18-6
     #592*t23-576*t17-t197+13632*t25
      t211 = -832*t26-64*t27-64*t29+3008*t33-10304*t35+6720*t36-4288*t38
     #-t97+256*t32+4608*t138-7872*t3+4608*t109-t112
      t213 = s2l1**2
      t224 = 21*t17-t161+9*t32-12*t22-18*t2-t20-10*t1+363*t23-252*t27-22
     #*t5+90*t3
      t233 = t172-141*t18+73*t26-t175-153*t38-t177+t178+126*t35+63*t29+4
     #05*t36-162*t25-387*t33
      t242 = t96*l1
      t244 = t108*t242
      t261 = t17+20*t22+522*t25-83*t26+659*t27-221*t29+2659*t33-1764*t35
     #-2073*t36+795*t38+225*t96
      t263 = t2*t1
      t273 = t96*t22
      t293 = t1*t242
      t295 = t17*t242
      t297 = t2*t242
      t299 = z*t242
      t301 = 55*t263*t273-330*t273*z-1100*t17*t273-330*t273*t108-160*t26
     #3*t242-60*t263*t32+825*t273*t1-2650*t293+3450*t295-2525*t297+1085*
     #t299
      t309 = (-l1-z+t5)**2
      t315 = 5*t5
      t319 = 2*t1
      t323 = -1-15*t118+15*t120+t116-30*t113+t242-t244+t21-t315-8*t18+22
     #*t23+2*t17-t319-4*t22-9*t25+8*t26
      t340 = -25*t27+15*t29-37*t33+17*t35+39*t36-20*t38-3*t96+4*t32+2*t3
     #-3*t109+3*t111+30*t171+10*t293-10*t295+5*t297-5*t299
      t356 = 4*t29
      t357 = 2*t27
      t359 = 4*z
      t360 = 2*t22
      t372 = -t58*t4
      t373 = s13e1(t372)
      t375 = 2048*t12
      t376 = 2048*t15
      t390 = -12*t18+2+t319+t359-6*l1-4*t23+2*t25+2*t26+8*t27-12*t29+6*t
     #22+8*t33-2*t35-12*t36+8*t38-2*t32
      t411 = 128*t22
      t412 = 256*t1
      t413 = 128*t27
      t415 = 256*t29
      t419 = 256*z
      t428 = 6*t1
      t462 = 4+54*t118-66*t120-4*t116+124*t113-4*t242+4*t244-t359-14*l1+
     #31*t5+28*t18-84*t23-6*t17+t428+22*t22+32*t25
      t480 = -37*t26+110*t27-76*t29+139*t33-59*t35-163*t36+95*t38+14*t96
     #-22*t32-8*t3+10*t109-10*t111-116*t171-40*t293+40*t295-20*t297+20*t
     #299
      t486 = s13e0(t372)
      s1 = (t47+t73+t83+t93+(t107+t127)*t129+((t143+t156)*s2l1+64*z*(t17
     #0+t184))*s1l1+(t199+t211)*t213+64*z*(t224+t233)*s2l1-64*z*(2997*t1
     #18-1353*t120-130*t116+3322*t113-185*t242+985*t244+5*t2+24*t5+119*t
     #18-855*t23+t261+160*t263*t96-115*t32+10*t138-70*t3+558*t109-1093*t
     #111+5*t263*t22+55*t273+825*t273*t2-4258*t171+t301)*t8)/t309
      s2 = (t47+t93-64*z*(t323+t340)*t8+t73+(-t197+192*t26-192*t5+192*z)
     #*t129+t83+(64*t26-64*t1-64*t5+64*z)*t213+((128*t26+128*z-128*t5)*s
     #2l1-64*z*(-1-t356+t357-4*t26-l1-t359+t319+t315+t360))*s1l1+64*z*(2
     #-l1-z+t5)*s2l1)*t373-((-t375+t376+64*t390*z*t9)*t46+(-1024*t63+512
     #*t68)*t72+(t375-t376+128*(6*t26-t54-t20+t159-2-t357+t356-t360)*t4*
     #t79)*l5+(-1024*t87+512*t90)*t71+(128-t94+t411-t100+t412+t413+768*t
     #5-t142-t415)*t129+((-t94+256*t22-t419+t412+256*t27+896*t5-640*t26-
     #512*t29)*s2l1-64*z*(5*l1-5*t22+t164-t428-5*t27-15*t5+10*t26+10*t29
     #))*s1l1+(t411-t419+384*t1+t413+512*t5-t142-t415)*t213+64*z*(-6+8*l
     #1-t22+6*z+t319-t27-7*t5-t26+t76)*s2l1-64*z*(t462+t480)*t8)*t486
      taux7 = s1+s2
      t1 = z**2
      t2 = t1*z
      t3 = l1**2
      t4 = t3*l1
      t5 = t2*t4
      t7 = t1**2
      t8 = t7*t4
      t10 = t1*t4
      t12 = z*t4
      t14 = t2*t3
      t16 = t1*l1
      t18 = t1*t3
      t20 = z*t3
      t22 = t2*l1
      t24 = t7*z
      t25 = t3**2
      t32 = t7*t3
      t34 = z*l1
      t36 = t25*l1
      t41 = 661376*t5-392320*t8-540288*t10+210688*t12-234496*t14+14336*t
     #16+147584*t18-41216*t20-17920*t22-76032*t24*t25+91008*t24*t4-55296
     #*t24*t3+181504*t32-3328*t34+126720*z*t36+25344*t24*t36
      t48 = t7*l1
      t50 = t1*t25
      t54 = t2*t25
      t56 = z*t25
      t58 = t7*t25
      t67 = -126720*t7*t36+253440*t2*t36-253440*t1*t36-8064*t48+627456*t
     #50+14976*t24*l1-671744*t54-291584*t56+358016*t58+15872*t2+1920*t3-
     #8320*t1-30464*t4-12416*t7+53888*t25-25344*t36
      t70 = -l1-z+t34
      t82 = -144*t2+198*t25-223*t4+15*t3+117*t7+65*t1-1916*t14+1714*t18+
     #152*t34-432*t48-792*t54
      t94 = 847*t22-567*t16+1188*t50+1263*t12+198*t58-792*t56-594*t8-524
     #*t20-2451*t10+711*t32+2005*t5
      t115 = -89*l1-107*t2-78*t4+177*t3-65*z+134*t1-286*t14+749*t18+372*
     #t34+315*t22-598*t16+234*t12-640*t20-234*t10+78*t5
      t119 = t70**2
      t131 = z*s2l1
      t132 = t131*l1
      t134 = t1*s2l1
      t135 = t134*l1
      t137 = s1l1*t1
      t138 = t137*l1
      t140 = z*s1l1
      t141 = t140*l1
      t144 = 6*t5-2*t8-6*t10+2*t12-14*t14-2*t16+10*t18-2*t20+10*t22-2*t1
     #32+2*t135+14*t138-6*t141+6*t32
      t147 = t137*t3
      t149 = t2*s1l1
      t150 = t149*t3
      t152 = t149*l1
      t154 = t140*t3
      t164 = -2*t34-6*t48-8*t147+4*t150-8*t152+4*t154+2*z+2*t140+4*t149-
     #4*t137+4*t131-2*t134-2*t2+2*t7
      t168 = (1-l1)*(1-z)
      t169 = s13e1(t168)
      t183 = -30*t5+10*t8+30*t10-10*t12+74*t14+12*t16-52*t18+10*t20-56*t
     #22+22*t132-30*t135
      t196 = -90*t138+48*t141-32*t32+10*t34+34*t48+58*t147-26*t150+48*t1
     #52-38*t154-10*z-16*t140-24*t149
      t200 = s2l1*t2
      t219 = 28*t137-22*t131+10*t200*l1-8*t200*t3+20*t134*t3+6*t131*t4+2
     #*t200*t4-6*t137*t4+6*t140*t4-6*t134*t4+2*t149*t4
      t239 = -16*t131*t3+16*t134+12*t2-4*t200+4*t3*s2l1-2*t4*s1l1+6*t3*s
     #1l1-6*l1*s1l1-2*t4*s2l1-2*s2l1*l1-12*t7+2*s1l1
      t242 = s13e0(t168)
      taux8 = ((t41+t67)*s1l1+128*t70*(t82+t94)*s2l1-128*z*t70*t115)/t11
     #9+128*(t144+t164)*t169-128*(t183+t196+t219+t239)*t242

CCC Coefficients

      tcoff11 = taux7
      tcoff12 = taux7
      tcoff13 = taux8
      tcoff14 = -taux8

CCC Jet functions

      tFjet11 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff11)
      tFjet12 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff12)
      tFjet13 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff13)
      tFjet14 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff14)

CCC Sum

      RVnum4e1 = tFjet11*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet14
     #*tcoff14
      RETURN
      END

      double precision function RVnum4e2(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t1 = l1**2
      t2 = t1*l1
      t3 = 256*t2
      t4 = 512*t1
      t5 = 256*l1
      t6 = z**2
      t7 = 64*t6
      t8 = z*l1
      t10 = 256*z
      t11 = t6*z
      t13 = t1*z
      t15 = t6*t1
      t17 = t6*l1
      t19 = t11*l1
      t21 = t6**2
      t22 = t21*t1
      t24 = z*t2
      t25 = 704*t24
      t26 = t21*t2
      t27 = 64*t26
      t28 = t11*t2
      t29 = 64*t28
      t30 = t6*t2
      t31 = 576*t30
      t32 = -t3+t4-t5-t7+704*t8-t10+640*t11-1152*t13+704*t15+320*t17-768
     #*t19-64*t22+t25+t27+t29-t31
      t33 = s1l1**2
      t38 = t21*l1
      t44 = t11*t1
      t45 = 256*t44
      t51 = -512*t2+t4-512*t6+1152*t11+128*t38-1280*t13+768*t15+1408*t17
     #-1536*t19+t45-256*t22+1408*t24+128*t26+128*t28-1152*t30
      t53 = 5*t2
      t55 = 8*l1
      t58 = 8*z
      t61 = 19*t15
      t63 = 12*t19
      t65 = 15*t24
      t66 = 5*t28
      t67 = 15*t30
      t68 = -t53+13*t1-t55+8*t6+29*t8-t58+10*t11-29*t13+t61-9*t17-t63-3*
     #t44+t65+t66-t67
      t78 = 512*t17
      t81 = -t3-64*t21+384*t11+192*t38-128*t13+64*t15+t78-704*t19+t45-19
     #2*t22+t25+t27+t29-t31
      t82 = s2l1**2
      t84 = 2*t1
      t85 = 2*t6
      t88 = 6*t13
      t92 = 2*t44
      t93 = -t53+t84+t85+4*t8+8*t11-t88+6*t15+7*t17-11*t19-t92+t65+t66-t
     #67
      t98 = t1**2
      t99 = t98*l1
      t100 = t6*t99
      t104 = 5*t1
      t105 = t98*t1
      t109 = t21*z
      t110 = t109*t1
      t113 = t109*t2
      t115 = t109*t98
      t117 = t21*t98
      t119 = z*t98
      t121 = t11*t98
      t125 = t6*t98
      t127 = t21*t6
      t136 = -49*t113-139*t115+494*t117-287*t119-806*t121-20*t109*l1+679
     #*t125+11*t127*t98-10*t127*t1+315*t105*t21+21*t127*t105
      t148 = t21*t99
      t152 = 6*t11
      t155 = -126*t105*z-420*t11*t105-126*t105*t109-42*t127*t99+20*t127*
     #t2-710*t148+48*t98-58*t99+t152-10*t21-45*t19
      t156 = t109*t99
      t158 = t11*t99
      t168 = z*t99
      t170 = 268*t156+1000*t158-114*t22+76*t24-31*t26+193*t28-203*t30+52
     #*t38+71*t44+21*t105+332*t168
      t175 = (1-l1+t8)**2
      t176 = 1/t175
      t181 = (-l1-z+t8)**2
      t184 = -1+z
      t185 = t184*t176
      t186 = l6**2
      t188 = t38*t185*t186
      t189 = 1024*t188
      t191 = t38*t185*l6
      t192 = 1024*t191
      t194 = 2*z
      t195 = 3*l1
      t201 = 6*t19-1-t6-t194+t195+t92-t22-t17-4*t15+t88-3*t1-4*t28+t26+6
     #*t30-4*t24+t2
      t206 = l5**2
      t207 = (t189-t192+64*t201*z*t185)*t206
      t208 = dsqrt(z)
      t211 = dsqrt(1-l5)
      t213 = dsqrt(l1)
      t214 = t208*t11*t211*t213
      t215 = 1-l1
      t216 = dsqrt(t215)
      t217 = t216*t184
      t219 = t217*t176*l6
      t220 = t214*t219
      t223 = t214*t217*t176
      t226 = dsqrt(l5)
      t227 = (512*t220-256*t223)*t226
      t229 = 5*t8
      t240 = -1-25*t15-t229+10*t100+4*t2+15*t13+t195-4*t1+8*t17-t85+t110
     #-3*t113+3*t115-15*t117+15*t119+30*t121
      t255 = -30*t125+5*t148-3*t98+t99+2*t11-8*t19-t156-10*t158-9*t22-20
     #*t24+17*t26-37*t28+39*t30+2*t38+22*t44-5*t168
      t261 = t208*t6*t211
      t262 = 2*t8
      t263 = t17-t262-z+l1-1
      t266 = t261*t213*t263*t219
      t271 = t261*t213*t263*t216*t185
      t274 = t226*l5
      t275 = (512*t266-256*t271)*t274
      t276 = 192*t17
      t277 = 192*t8
      t278 = 192*z
      t282 = 2*l1
      t283 = 2*t13
      t286 = t6*t176
      t290 = (-t189+t192+128*(-3*t17+t8+z-t282+1+t15-t283+t1)*t184*t286)
     #*l5
      t297 = 4*t13
      t298 = 2*t15
      t299 = 4*t17
      t300 = 4*z
      t312 = -t215*t184
      t313 = s13e2(t312)
      t315 = 2048*t188
      t316 = 2048*t191
      t320 = 2*t17
      t329 = -t63+2+t85+t300-6*l1-4*t44+2*t22+t320+8*t15-12*t13+6*t1+8*t
     #28-2*t26-12*t30+8*t24-2*t2
      t335 = 128*t1
      t337 = 128*t15
      t339 = 256*t6
      t340 = 256*t13
      t361 = 6*t6
      t373 = 6*z
      t380 = 4*l1
      t400 = 4+110*t15+31*t8-40*t100-22*t2-76*t13-14*l1-t300+22*t1-37*t1
     #7+t361-4*t110+10*t113-10*t115+54*t117-66*t119
      t409 = 32*t22
      t417 = -116*t121+124*t125-20*t148+14*t98-4*t99-t152+28*t19+4*t156+
     #40*t158+t409+95*t24-59*t26+139*t28-163*t30-8*t38-84*t44+20*t168
      t423 = s13e1(t312)
      t425 = 64*t1
      t426 = 128*l1
      t429 = 192*t13
      t430 = 192*t15
      t431 = 64*t44
      t438 = 128*t19
      t470 = -5-119*t15-40*t8+50*t100+29*t2+92*t13+18*l1+t373-29*t1+43*t
     #17-7*t6+5*t110-8*t113+8*t115-50*t117+80*t119
      t486 = 120*t121-140*t125+25*t148-18*t98+5*t99+t152-31*t19-5*t156-5
     #0*t158-t409-113*t24+50*t26-131*t28+173*t30+10*t38+83*t44-25*t168
      t492 = s13e0(t312)
      s1 = (t32*t33+(t51*s2l1+64*z*t68)*s1l1+t81*t82+64*z*t93*s2l1-64*z*
     #(-t61-6*t8-790*t100-6*t2+11*t13-t104+315*t105*t6+19*t17-t6+66*t110
     #+t136+t155+t170)*t176)/t181+(t207+t227-64*z*(t240+t255)*t176+t275+
     #(-t7+t276-t277+t278)*t33+t290+(64*t17-64*t6-64*t8+64*z)*t82+((128*
     #t17+128*z-128*t8)*s2l1-64*z*(-1-t297+t298-t299-l1-t300+t85+t229+t8
     #4))*s1l1+64*z*(2-l1-z+t8)*s2l1)*t313
      s2 = s1
      s4 = ((-t315+t316+64*t329*z*t185)*t206+(128+t335-t5+768*t8+t337-51
     #2*z+t339-t78-t340)*t33+(-1024*t220+512*t223)*t226+(-1024*t266+512*
     #t271)*t274+((-t5-640*t17-t10-512*t13+256*t15+256*t1+t339+896*t8)*s
     #2l1-64*z*(5*l1-15*t8-5*t15+t58-t361+10*t13-t104+10*t17))*s1l1+(t33
     #5-t10+384*t6+t337-t78+512*t8-t340)*t82+64*z*(-6+t283-t1+t55-t17-t1
     #5+t373-7*t8+t85)*s2l1+(t315-t316+128*(6*t17-t262-t194+t380-2-t298+
     #t297-t84)*t184*t286)*l5-64*z*(t400+t417)*t176)*t423
      s5 = -(t207+t275+t290+t227+(t425-t426-192*t6-64*t8+64+t10-t429+t43
     #0+t276-t431)*t33+((t335-t426-640*t6+640*z-384*t13+384*t15+t438-128
     #*t44)*s2l1-64*z*(t84+t380+4*t6-t262-6-t297+t298-t320))*s1l1+(t425-
     #t339-t277+t278-64*t11-t429+t430+64*t17+t438-t431)*t82+64*z*(t84-t2
     #82-t361-t262+4-t297+t298+t299)*s2l1-64*z*(t470+t486)*t176)*t492
      s3 = s4+s5
      taux7 = s2+s3
      t1 = l1**2
      t2 = t1*z
      t4 = z**2
      t5 = t4*z
      t6 = t5*l1
      t8 = t4*l1
      t10 = t4*t1
      t12 = t4**2
      t13 = t12*z
      t14 = t1**2
      t15 = t14*l1
      t20 = z*l1
      t25 = t5*t1
      t27 = t1*l1
      t29 = z*t27
      t31 = t4*t27
      t33 = t12*t27
      t35 = t5*t27
      t37 = 62976*t2-12416*t6-4224*t8-194432*t10-18432*t13*t15-92160*z*t
     #15+2304*t20-17408*t5-6784*t1+10112*t4+282368*t25+34944*t27-217728*
     #t29+517248*t31+334080*t33-594816*t35
      t38 = t12*t1
      t46 = t12*t14
      t48 = z*t14
      t50 = t5*t14
      t54 = t4*t14
      t56 = t12*l1
      t67 = -199424*t38+55296*t13*t1-73728*t13*t27+55296*t13*t14-267776*
     #t46+241664*t48+518144*t50-18432*t13*l1-500736*t54+32768*t56+92160*
     #t12*t15+8192*t12-46592*t14-184320*t5*t15+18432*t15+184320*t4*t15
      t70 = -l1-z+t20
      t82 = -168*t20-1593*t10+1634*t25-930*t6+666*t8+588*t2-99*t4+236*t5
     #+220*t27-53*t1+576*t48
      t94 = 576*t50-144*t46-864*t54-576*t38-1516*t35+1956*t31-144*t12+43
     #2*t33-1092*t29+432*t56-144*t14
      t115 = -399*t20-512*t10+190*t25-278*t6+584*t8+454*t2-186*t4+89*z+1
     #04*t5+16*t27-132*t1+93*l1-16*t35+48*t31-48*t29
      t119 = t70**2
      t127 = s1l1*t4
      t129 = s1l1*t5
      t131 = t129*t1
      t134 = z*s2l1
      t137 = z*s1l1
      t138 = t137*t1
      t140 = s2l1*t4
      t143 = -2*t2+10*t6-2*t8+10*t10-2*t20-4*t127+4*t129+4*t131-2*t5+4*t
     #134+2*z+4*t138-2*t140-14*t25
      t152 = t140*l1
      t154 = t137*l1
      t156 = t134*l1
      t158 = t127*l1
      t160 = t129*l1
      t162 = t127*t1
      t164 = 2*t29-6*t31-2*t33+6*t35+6*t38-6*t56+2*t12+2*t137+2*t152-6*t
     #154-2*t156+14*t158-8*t160-8*t162
      t168 = (1-l1)*(1-z)
      t169 = s13e2(t168)
      t183 = 10*t2-56*t6+12*t8-52*t10+10*t20+28*t127-24*t129-26*t131+12*
     #t5-22*t134-10*z
      t188 = t140*t27
      t190 = t127*t27
      t192 = t129*t27
      t194 = t140*t1
      t200 = -38*t138+16*t140+74*t25+2*s1l1-6*t188-6*t190+2*t192+20*t194
     #-10*t29+30*t31+10*t33-30*t35
      t206 = t134*t27
      t208 = t137*t27
      t210 = s2l1*t5
      t211 = t210*t1
      t213 = t210*l1
      t215 = t134*t1
      t219 = -32*t38+34*t56-12*t12-16*t137+6*t206+6*t208-8*t211+10*t213-
     #16*t215-30*t152+48*t154
      t224 = s2l1*t1
      t226 = s2l1*t27
      t229 = t210*t27
      t231 = s1l1*t1
      t233 = s1l1*t27
      t235 = s2l1*l1
      t237 = s1l1*l1
      t239 = 22*t156-90*t158+48*t160+58*t162+4*t224-2*t226-4*t210+2*t229
     #+6*t231-2*t233-2*t235-6*t237
      t242 = s13e1(t168)
      t256 = -12*t2+110*t6-30*t8+78*t10-14*t20-75*t127+52*t129+54*t131-2
     #9*t5+41*t134+15*z
      t269 = 94*t138+4*t4-42*t140-120*t25-11*s1l1+18*t188+18*t190-6*t192
     #-78*t194+14*t29-42*t31-14*t33
      t283 = 42*t35+54*t38-66*t56+26*t12+50*t137-18*t206-18*t208+32*t211
     #-43*t213+60*t215+105*t152-123*t154
      t296 = -73*t156+195*t158-97*t160-128*t162-14*t224+6*t226+17*t210-6
     #*t229-20*t231+6*t233+11*t235+25*t237
      t299 = s13e0(t168)
      taux8 = ((t37+t67)*s1l1+128*t70*(t82+t94)*s2l1-128*z*t70*t115)/t11
     #9+128*(t143+t164)*t169+128*(t183+t200+t219+t239)*t242-128*(t256+t2
     #69+t283+t296)*t299

CCC Coefficients

      tcoff11 = taux7
      tcoff12 = taux7
      tcoff13 = taux8
      tcoff14 = -taux8

CCC Jet functions

      tFjet11 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff11)
      tFjet12 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff12)
      tFjet13 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff13)
      tFjet14 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff14)

CCC Sum

      RVnum4e2 = tFjet11*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet14
     #*tcoff14
      RETURN
      END

      double precision function RVnum4e3(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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
      t2 = t1**2
      t3 = t2*z
      t4 = l1**2
      t5 = t3*t4
      t7 = t2*t4
      t9 = t1*z
      t10 = t4*l1
      t11 = t9*t10
      t12 = 8*t11
      t13 = t2*t10
      t15 = t1*t10
      t17 = z*t10
      t19 = t9*t4
      t21 = t1*l1
      t23 = t1*t4
      t25 = z*t4
      t26 = 6*t25
      t27 = t9*l1
      t29 = t4**2
      t30 = z*t29
      t32 = t2*t29
      t34 = t3*t29
      t36 = t3*t10
      t38 = t29*l1
      t39 = z*t38
      t41 = t3*t38
      t43 = t2*t38
      t45 = -34*t5+75*t7-t12-33*t13+47*t15-36*t17-68*t19-21*t21+19*t23+t
     #26+46*t27+46*t30-55*t32+14*t34+28*t36-18*t39-18*t41+45*t43
      t46 = t9*t38
      t48 = t1*t38
      t50 = t2*t1
      t58 = t2*l1
      t60 = t1*t29
      t64 = t9*t29
      t66 = z*l1
      t67 = 2*t66
      t69 = 3*t4
      t75 = -60*t46+45*t48-7*t50*t10+3*t50*t38+5*t50*t4-t50*t29-37*t58-9
     #5*t60+10*t3*l1+100*t64+t67-10*t9-t69+5*t1+9*t10+5*t2-9*t29+3*t38
      t79 = (1-l1+t66)**2
      t80 = 1/t79
      t82 = (-l1-z+t66)**2
      t87 = -1+z
      t88 = t87*t80
      t89 = l6**2
      t91 = t58*t88*t89
      t92 = 1024*t91
      t94 = t58*t88*l6
      t95 = 1024*t94
      t97 = 2*z
      t98 = 3*l1
      t104 = 6*t27-1-t1-t97+t98+2*t19-t7-t21-4*t23+t26-t69-4*t11+t13+6*t
     #15-4*t17+t10
      t109 = l5**2
      t110 = (t92-t95+64*t104*z*t88)*t109
      t111 = dsqrt(z)
      t114 = dsqrt(1-l5)
      t116 = dsqrt(l1)
      t117 = t111*t9*t114*t116
      t118 = 1-l1
      t119 = dsqrt(t118)
      t120 = t119*t87
      t122 = t120*t80*l6
      t123 = t117*t122
      t126 = t117*t120*t80
      t129 = dsqrt(l5)
      t130 = (512*t123-256*t126)*t129
      t145 = -1+t5-9*t7-37*t11+17*t13+39*t15-20*t17+22*t19+8*t21-25*t23+
     #15*t25-8*t27+15*t30-15*t32+3*t34-3*t36
      t153 = 5*t66
      t154 = 2*t9
      t156 = 2*t1
      t159 = -5*t39-t41+5*t43-10*t46+10*t48+2*t58-30*t60+30*t64-t153+t98
     #+t154-4*t4-t156+4*t10-3*t29+t38
      t165 = t111*t1*t114
      t166 = t21-t67-z+l1-1
      t169 = t165*t116*t166*t122
      t174 = t165*t116*t166*t119*t88
      t177 = t129*l5
      t178 = (512*t169-256*t174)*t177
      t179 = 64*t1
      t180 = 192*t21
      t181 = 192*t66
      t182 = 192*z
      t184 = s1l1**2
      t187 = 2*l1
      t188 = 2*t25
      t191 = t1*t80
      t195 = (-t92+t95+128*(-3*t21+t66+z-t187+1+t23-t188+t4)*t87*t191)*l
     #5
      t198 = s2l1**2
      t203 = 4*t25
      t204 = 2*t23
      t205 = 4*t21
      t206 = 4*z
      t207 = 2*t4
      t219 = -t118*t87
      t220 = s13e3(t219)
      t222 = 2048*t91
      t223 = 2048*t94
      t228 = 2*t21
      t236 = -12*t27+2+t156+t206-6*l1-4*t19+2*t7+t228+8*t23-12*t25+6*t4+
     #t12-2*t13-12*t15+8*t17-2*t10
      t242 = 128*t4
      t243 = 256*l1
      t245 = 128*t23
      t247 = 256*t1
      t248 = 512*t21
      t249 = 256*t25
      t261 = 256*z
      t268 = 5*l1
      t272 = 6*t1
      t274 = 5*t4
      t286 = 6*z
      t293 = 4*l1
      t301 = 32*t7
      t314 = 4-4*t5+t301+139*t11-59*t13-163*t15+95*t17-84*t19-37*t21+110
     #*t23-76*t25+28*t27-t206-66*t30+54*t32-10*t34
      t326 = 6*t9
      t331 = 10*t36+20*t39+4*t41-20*t43+40*t46-40*t48-8*t58+124*t60-116*
     #t64+31*t66-14*l1-t326+22*t4+t272-22*t10+14*t29-4*t38
      t337 = s13e2(t219)
      t339 = 128*l1
      t342 = 64*t19
      t344 = 192*t25
      t345 = 192*t23
      t346 = 64*t4
      t351 = 128*t19
      t353 = 128*t27
      t354 = 384*t25
      t355 = 384*t23
      t359 = 4*t66
      t367 = 64*t9
      t388 = 2-2*t5+13*t7+20*t11-7*t13-28*t15+19*t17-32*t19-31*t21+38*t2
     #3-22*t25+19*t27-7*z-21*t30+9*t32-t34
      t398 = 7*t1
      t402 = t36+10*t39+2*t41-10*t43+20*t46-20*t48-4*t58+34*t60-26*t64+2
     #1*t66-t268-t154+t274+t398-5*t10+5*t29-2*t38
      t408 = s13e0(t219)
      t444 = -5+5*t5-t301-131*t11+50*t13+173*t15-113*t17+83*t19+43*t21-1
     #19*t23+92*t25-31*t27+t286+80*t30-50*t32+8*t34
      t460 = -8*t36-25*t39-5*t41+25*t43-50*t46+50*t48+10*t58-140*t60+120
     #*t64-40*t66+18*l1+t326-29*t4-t398+29*t10-18*t29+5*t38
      t466 = s13e1(t219)
      s1 = -64*z*(t45+t75)*t80/t82+(t110+t130-64*z*(t145+t159)*t80+t178+
     #(-t179+t180-t181+t182)*t184+t195+(64*t21-64*t1-64*t66+64*z)*t198+(
     #(128*t21+128*z-128*t66)*s2l1-64*z*(-1-t203+t204-t205-l1-t206+t156+
     #t153+t207))*s1l1+64*z*(2-l1-z+t66)*s2l1)*t220
      s2 = s1+((-t222+t223+64*t236*z*t88)*t109+(128+t242-t243+768*t66+t2
     #45-512*z+t247-t248-t249)*t184+(-1024*t123+512*t126)*t129+(-1024*t1
     #69+512*t174)*t177+((-t243-640*t21-t261-512*t25+256*t23+256*t4+t247
     #+896*t66)*s2l1-64*z*(t268-15*t66-5*t23+8*z-t272+10*t25-t274+10*t21
     #))*s1l1+(t242-t261+384*t1+t245-t248+512*t66-t249)*t198+64*z*(-6+t1
     #88-t4+8*l1-t21-t23+t286-7*t66+t156)*s2l1+(t222-t223+128*(6*t21-t67
     #-t97+t293-2-t204+t203-t207)*t87*t191)*l5-64*z*(t314+t331)*t80)*t33
     #7
      taux7 = s2-((-t339-64*z+256*t66-t342-128*t21-t344+t345+64+t346)*t1
     #84+((-t339+128*z+384*t66-t351-384*t21+t353-t354+t355-128*t1+t242)*
     #s2l1-64*z*(t187+t97-t359+t228-2))*s1l1+(128*t66-t342-256*t21+t353-
     #t344+t345+t179-t367+t346)*t198+64*z*(t187+t97-t359+t228-t156)*s2l1
     #-64*z*(t388+t402)*t80)*t408+(t110+t178+t195+t130+(t346-t339-192*t1
     #-64*t66+64+t261-t344+t345+t180-t342)*t184+((t242-t339-640*t1+640*z
     #-t354+t355+t353-t351)*s2l1-64*z*(t207+t293+4*t1-t67-6-t203+t204-t2
     #28))*s1l1+(t346-t247-t181+t182-t367-t344+t345+64*t21+t353-t342)*t1
     #98+64*z*(t207-t187-t272-t67+4-t203+t204+t205)*s2l1-64*z*(t444+t460
     #)*t80)*t466
      t1 = z**2
      t2 = t1*z
      t3 = l1**2
      t4 = t2*t3
      t6 = t3**2
      t7 = t6*l1
      t9 = t3*l1
      t11 = t1**2
      t13 = t1*t9
      t15 = t11*t9
      t17 = t2*t9
      t19 = t11*t3
      t21 = t11*z
      t26 = t11*t6
      t28 = z*t6
      t30 = t2*t6
      t34 = t1*t6
      t36 = -61952*t4+9216*t7+1280*t9+8192*t11+80640*t13+87040*t15-12416
     #0*t17+23936*t19-23040*t21*t9+27648*t21*t6-125824*t26+88576*t28+226
     #816*t30+4608*t21*l1-201984*t34
      t37 = t11*l1
      t49 = z*l1
      t55 = z*t9
      t57 = t2*l1
      t59 = z*t3
      t61 = t1*t3
      t63 = t1*l1
      t65 = -39424*t37+92160*t1*t7-92160*t2*t7+46080*t11*t7-9216*t21*t7-
     #46080*z*t7+1280*t49-6144*t2+4736*t3-1408*t1-15232*t6-21760*t55+563
     #20*t57-23552*t59+56832*t61-22784*t63
      t68 = -l1-z+t49
      t91 = -72*t26+288*t30-432*t34+288*t28-180*t19-695*t17+789*t13+37*t
     #3+47*t9+36*t11-72*t6
      t112 = -25*l1-41*z+28*t49+132*t4+77*t63-80*t57+238*t59-317*t61+32*
     #t1+14*t2-198*t55-66*t17+198*t13-53*t3+66*t9
      t116 = t68**2
      t121 = s1l1*t1
      t123 = s1l1*t2
      t125 = z*s1l1
      t127 = s2l1*t1
      t129 = z*s2l1
      t132 = 6*t13
      t133 = 2*t15
      t134 = 6*t17
      t138 = -14*t4+2*z-4*t121+4*t123+2*t125-2*t127+4*t129+2*t11-t132-t1
     #33+t134+6*t19-6*t37-2*t49
      t140 = t125*t3
      t142 = t127*l1
      t144 = t125*l1
      t146 = t129*l1
      t148 = t121*l1
      t150 = t123*l1
      t152 = 2*t55
      t157 = t121*t3
      t159 = t123*t3
      t161 = -2*t2+4*t140+2*t142-6*t144-2*t146+14*t148-8*t150+t152+10*t5
     #7-2*t59+10*t61-2*t63-8*t157+4*t159
      t165 = (1-l1)*(1-z)
      t166 = s13e3(t165)
      t177 = s2l1*t2
      t179 = s2l1*t9
      t180 = 2*t179
      t181 = s1l1*l1
      t183 = 74*t4-10*z+2*s1l1+28*t121-24*t123-16*t125+16*t127-22*t129-4
     #*t177-t180-6*t181
      t184 = s1l1*t3
      t186 = s1l1*t9
      t187 = 2*t186
      t188 = s2l1*l1
      t190 = s2l1*t3
      t200 = 6*t184-t187-2*t188+4*t190-12*t11+30*t13+10*t15-30*t17-32*t1
     #9+34*t37+10*t49+12*t2
      t202 = t129*t9
      t203 = 6*t202
      t204 = t129*t3
      t206 = t177*t3
      t208 = t127*t9
      t209 = 6*t208
      t210 = t177*l1
      t212 = t123*t9
      t213 = 2*t212
      t214 = t121*t9
      t215 = 6*t214
      t216 = t125*t9
      t217 = 6*t216
      t218 = t177*t9
      t219 = 2*t218
      t220 = t127*t3
      t223 = t203-16*t204-8*t206-t209+10*t210+t213-t215+t217+t219+20*t22
     #0-38*t140
      t236 = -30*t142+48*t144+22*t146-90*t148+48*t150-10*t55-56*t57+10*t
     #59-52*t61+12*t63+58*t157-26*t159
      t239 = s13e2(t165)
      t253 = -120*t4+15*z-11*s1l1-75*t121+52*t123+50*t125-42*t127+41*t12
     #9+17*t177+6*t179+25*t181
      t266 = -20*t184+6*t186+11*t188-14*t190+26*t11-42*t13-14*t15+42*t17
     #+54*t19-66*t37-14*t49-29*t2
      t280 = 4*t1-18*t202+60*t204+32*t206+18*t208-43*t210-6*t212+18*t214
     #-18*t216-6*t218-78*t220+94*t140
      t293 = 105*t142-123*t144-73*t146+195*t148-97*t150+14*t55+110*t57-1
     #2*t59+78*t61-30*t63-128*t157+54*t159
      t296 = s13e1(t165)
      t309 = 52*t4-4*z+20*s1l1+92*t121-48*t123-72*t125+40*t127-28*t129-2
     #0*t177-t180-31*t181
      t317 = 13*t184-t187-20*t188+11*t190-24*t11+t132+t133-t134-28*t19+5
     #0*t37-t49+36*t2
      t325 = -16*t1+t203-57*t204-35*t206-t209+53*t210+t213-t215+t217+t21
     #9+81*t220-63*t140
      t337 = -126*t142+118*t144+93*t146-159*t148+72*t150-t152-93*t57-4*t
     #59-20*t61+44*t63+87*t157-37*t159
      t340 = s13e0(t165)
      taux8 = ((t36+t65)*s1l1+128*t68*(60*t49+416*t4-163*t63+103*t57-18*
     #t59-255*t61+71*t1-112*t2-357*t55+216*t15+t91)*s2l1-128*z*t68*t112)
     #/t116+128*(t138+t161)*t166+128*(t183+t200+t223+t236)*t239+128*(t25
     #3+t266+t280+t293)*t296-128*(t309+t317+t325+t337)*t340

CCC Coefficients

      tcoff11 = taux7
      tcoff12 = taux7
      tcoff13 = taux8
      tcoff14 = -taux8

CCC Jet functions

      tFjet11 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff11)
      tFjet12 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff12)
      tFjet13 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff13)
      tFjet14 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff14)

CCC Sum

      RVnum4e3 = tFjet11*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet14
     #*tcoff14
      RETURN
      END

      double precision function RVnum4e4(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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
      t2 = t1**2
      t3 = t2*l1
      t4 = -1+z
      t5 = z*l1
      t7 = (1-l1+t5)**2
      t8 = 1/t7
      t9 = t4*t8
      t10 = l6**2
      t12 = t3*t9*t10
      t13 = 2048*t12
      t15 = t3*t9*l6
      t16 = 2048*t15
      t17 = t1*z
      t18 = t17*l1
      t20 = 2*t1
      t21 = 4*z
      t23 = l1**2
      t24 = t17*t23
      t26 = t2*t23
      t28 = t1*l1
      t29 = 2*t28
      t30 = t1*t23
      t32 = z*t23
      t34 = 6*t23
      t35 = t23*l1
      t36 = t17*t35
      t38 = t2*t35
      t40 = t1*t35
      t42 = z*t35
      t45 = -12*t18+2+t20+t21-6*l1-4*t24+2*t26+t29+8*t30-12*t32+t34+8*t3
     #6-2*t38-12*t40+8*t42-2*t35
      t50 = l5**2
      t52 = 256*l1
      t53 = 256*t1
      t54 = 512*t28
      t56 = 128*t30
      t57 = 256*t32
      t58 = 128*t23
      t61 = s1l1**2
      t63 = dsqrt(z)
      t66 = dsqrt(1-l5)
      t68 = dsqrt(l1)
      t69 = t63*t17*t66*t68
      t70 = 1-l1
      t71 = dsqrt(t70)
      t72 = t71*t4
      t74 = t72*t8*l6
      t75 = t69*t74
      t78 = t69*t72*t8
      t81 = dsqrt(l5)
      t84 = t63*t1*t66
      t85 = 2*t5
      t86 = t28-t85-z+l1-1
      t89 = t84*t68*t86*t74
      t94 = t84*t68*t86*t71*t9
      t97 = t81*l5
      t103 = 256*z
      t110 = 6*t1
      t112 = 5*t23
      t113 = 5*l1
      t123 = s2l1**2
      t125 = 2*t32
      t127 = 6*z
      t133 = 6*t28
      t134 = 2*z
      t135 = 4*l1
      t136 = 2*t30
      t137 = 4*t32
      t138 = 2*t23
      t141 = t1*t8
      t146 = t23**2
      t154 = t146*l1
      t155 = t2*t154
      t158 = t17*t154
      t160 = t2*z
      t161 = t160*t154
      t166 = 32*t26
      t167 = 4+14*t146-76*t32+139*t36-59*t38-163*t40+95*t42-22*t35-20*t1
     #55-4*t154+40*t158+4*t161+31*t5+28*t18-84*t24+t166
      t168 = z*t154
      t170 = t160*t23
      t173 = t160*t35
      t175 = t160*t146
      t177 = t2*t146
      t179 = z*t146
      t181 = t17*t146
      t183 = t1*t146
      t186 = t1*t154
      t188 = 6*t17
      t192 = 20*t168-4*t170-14*l1+10*t173-10*t175+54*t177-66*t179-116*t1
     #81+124*t183-8*t3-40*t186-t188+t110-t21+22*t23+110*t30-37*t28
      t199 = -t70*t4
      t200 = s13e3(t199)
      t203 = 192*t32
      t206 = 64*t24
      t207 = 192*t30
      t208 = 64*t23
      t209 = 128*l1
      t212 = 1024*t12
      t222 = 6*t18-1-t1-t134+3*l1+2*t24-t26-t28-4*t30+6*t32-3*t23-4*t36+
     #t38+6*t40-4*t42+t35
      t226 = 1024*t15
      t229 = 128*t24
      t230 = 384*t32
      t231 = 128*t18
      t232 = 384*t30
      t243 = 64*t17
      t258 = 2*l1
      t284 = -5-18*t146+92*t32-131*t36+50*t38+173*t40-113*t42+29*t35+25*
     #t155+5*t154-50*t158-5*t161-40*t5-31*t18+83*t24-t166
      t296 = 7*t1
      t300 = -25*t168+5*t170+18*l1-8*t173+8*t175-50*t177+80*t179+120*t18
     #1-140*t183+10*t3+50*t186+t188-t296+t127-29*t23-119*t30+43*t28
      t306 = s13e2(t199)
      t314 = 4*t5
      t335 = 22*t32
      t349 = 2+5*t146-t335+20*t36-7*t38-28*t40+19*t42-5*t35-10*t155-2*t1
     #54+20*t158+2*t161+21*t5+19*t18-32*t24+13*t26
      t362 = 10*t168-2*t170-t113+t173-t175+9*t177-21*t179-26*t181+34*t18
     #3-4*t3-20*t186-2*t17+t296-7*z+t112+38*t30-31*t28
      t368 = s13e1(t199)
      t388 = -4*t26-36*t36+18*t24+14*t38-2*t173-26*t42+t335+t258-2*t18+t
     #133-6*t5
      t391 = s13e0(t199)
      s2 = ((-t13+t16+64*t45*z*t9)*t50+(128-t52+t53-t54-512*z+t56-t57+t5
     #8+768*t5)*t61+(-1024*t75+512*t78)*t81+(-1024*t89+512*t94)*t97+((-6
     #40*t28-512*t32+256*t30+t53+896*t5-t103-t52+256*t23)*s2l1-64*z*(10*
     #t28+8*z+10*t32-t110-5*t30-t112+t113-15*t5))*s1l1+(-t103+t58-t54+t5
     #6+384*t1-t57+512*t5)*t123+64*z*(-6-t30+t125+8*l1-t28-t23+t127+t20-
     #7*t5)*s2l1+(t13-t16+128*(t133-t85-t134+t135-2-t136+t137-t138)*t4*t
     #141)*l5-64*z*(t167+t192)*t8)*t200
      s3 = ((-192*t1+64-t203+192*t28-64*t5-t206+t207+t103+t208-t209)*t61
     #+(t212+64*t222*z*t9-t226)*t50+((-t229-t209-t230+t231+t58+t232+640*
     #z-640*t1)*s2l1-64*z*(t136-t85-6+t138-t137-t29+4*t1+t135))*s1l1+(-t
     #243-t53-t206-192*t5-t203+192*z+t208+t207+64*t28+t231)*t123+(512*t7
     #5-256*t78)*t81+(512*t89-256*t94)*t97+64*z*(-t110+4*t28-t85+t136-t1
     #37+4+t138-t258)*s2l1+(-t212+128*(-3*t28+t5+z-t258+1+t30-t125+t23)*
     #t4*t141+t226)*l5-64*z*(t284+t300)*t8)*t306
      s1 = s2+s3
      taux7 = s1+(((-128*t1+t231-384*t28+t232-t230+128*z-t229+t58+384*t5
     #-t209)*s2l1-64*z*(t258-2+t29-t314+t134))*s1l1+(t208+t231+64*t1+t20
     #7-256*t28+128*t5-t206-t243-t203)*t123+64*z*(t29-t20+t258-t314+t134
     #)*s2l1+(t208-128*t28+t207-t203-t206+64-t209+256*t5-64*z)*t61-64*z*
     #(t349+t362)*t8)*t368+64*z*(-10*t177+2*t175-t34+6*t35-2*t146+20*t18
     #1-20*t183+10*t179+44*t40-30*t30+t388)*t8*t391
      t1 = z**2
      t2 = t1**2
      t3 = t2*l1
      t5 = l1**2
      t7 = t5**2
      t9 = t5*l1
      t12 = t1*z
      t23 = 12288*t3+1664*t5+4224*t7-5888*t9-4480*t1+10240*t12-5120*t2+4
     #224*t2*t7-16896*z*t7+25344*t1*t7-16896*t12*t7
      t24 = t2*t5
      t26 = t1*t5
      t28 = t1*t9
      t30 = t2*t9
      t32 = t12*t5
      t34 = t12*t9
      t36 = z*t5
      t38 = z*t9
      t40 = t12*l1
      t42 = z*l1
      t44 = t1*l1
      t46 = -5504*t24-7936*t26-35328*t28-5888*t30+13056*t32+23552*t34-12
     #80*t36+23552*t38-28416*t40-2816*t42+18944*t44
      t49 = -l1-z+t42
      t62 = 13*t5-33*t9-5*t1+57*t26-99*t28-22*t32+33*t34-48*t36+99*t38-1
     #1*t40-24*t42+35*t44
      t70 = 20*t12
      t71 = 15*z
      t74 = 16*t32
      t81 = 7*t5+12*t9+l1+40*t1-t70-t71+39*t26+36*t28-t74-12*t34-30*t36-
     #36*t38+48*t40+54*t42-103*t44
      t85 = t49**2
      t88 = l1*s1l1
      t90 = t5*s2l1
      t91 = 4*t90
      t92 = s1l1*t12
      t94 = z*s1l1
      t96 = s2l1*t1
      t98 = t38*s2l1
      t99 = 6*t98
      t100 = t34*s2l1
      t101 = 2*t100
      t102 = t28*s1l1
      t103 = 6*t102
      t107 = -6*t88+t91-24*t92-16*t94+16*t96+t99+t101-t103+34*t3-12*t2+1
     #0*t30
      t118 = s2l1*l1
      t119 = t118*t12
      t121 = t28*s2l1
      t122 = 6*t121
      t123 = 74*t32-30*t34+10*t36-10*t38-56*t40+10*t42-32*t24-52*t26+30*
     #t28+12*t44+10*t119-t122
      t125 = t34*s1l1
      t126 = 2*t125
      t127 = t36*s2l1
      t130 = t94*t5
      t132 = s2l1*t12
      t133 = 4*t132
      t136 = s1l1*t1
      t138 = z*s2l1
      t140 = t94*l1
      t142 = t96*l1
      t144 = t126-16*t127-10*z-38*t130-t133-2*t118+2*s1l1+28*t136-22*t13
     #8+48*t140-30*t142
      t145 = t138*l1
      t147 = t136*l1
      t149 = t92*t5
      t151 = t136*t5
      t153 = t92*l1
      t156 = t90*t12
      t158 = t26*s2l1
      t160 = t5*s1l1
      t162 = t9*s1l1
      t163 = 2*t162
      t164 = t9*s2l1
      t165 = 2*t164
      t166 = t38*s1l1
      t167 = 6*t166
      t168 = 22*t145-90*t147-26*t149+58*t151+48*t153+12*t12-8*t156+20*t1
     #58+6*t160-t163-t165+t167
      t173 = (1-l1)*(1-z)
      t174 = s13e3(t173)
      t188 = 25*t88-14*t90+52*t92+50*t94-42*t96-18*t98-6*t100+18*t102-66
     #*t3+26*t2-14*t30
      t201 = -120*t32+42*t34-12*t36+14*t38+110*t40-14*t42+54*t24+78*t26-
     #42*t28-30*t44-43*t119+18*t121
      t214 = -6*t125+60*t127+t71+94*t130+17*t132+11*t118-11*s1l1-75*t136
     #+41*t138-123*t140+105*t142-73*t145
      t227 = 195*t147+54*t149-128*t151-97*t153+4*t1-29*t12+32*t156-78*t1
     #58-20*t160+6*t162+6*t164-18*t166
      t230 = s13e2(t173)
      t246 = 4*z
      t248 = 16*t88-t91+16*t92+40*t94-8*t96-16*t3+8*t2-t74+44*t40+12*t42
     #+8*t24+8*t26-40*t44-16*t119+20*t127-t246+20*t130
      t260 = 16*t1
      t264 = t133+12*t118-12*s1l1-44*t136+4*t138-52*t140+44*t142-40*t145
     #+56*t147+12*t149-28*t151-20*t153+t260-t70+12*t156-28*t158-4*t160
      t266 = s13e0(t173)
      t277 = -31*t88+11*t90-48*t92-72*t94+40*t96+t99+t101-t103+50*t3-24*
     #t2+2*t30
      t288 = 52*t32-6*t34-4*t36-2*t38-93*t40-t42-28*t24-20*t26+6*t28+44*
     #t44+53*t119-t122
      t300 = t126-57*t127-t246-63*t130-20*t132-20*t118+20*s1l1+92*t136-2
     #8*t138+118*t140-126*t142+93*t145
      t309 = -159*t147-37*t149+87*t151+72*t153-t260+36*t12-35*t156+81*t1
     #58+13*t160-t163-t165+t167
      t312 = s13e1(t173)
      taux8 = ((t23+t46)*s1l1+128*t49*t62*s2l1-128*z*t49*t81)/t85+128*(t
     #107+t123+t144+t168)*t174+128*(t188+t201+t214+t227)*t230-128*(t248+
     #t264)*t266+128*(t277+t288+t300+t309)*t312

CCC Coefficients

      tcoff11 = taux7
      tcoff12 = taux7
      tcoff13 = taux8
      tcoff14 = -taux8

CCC Jet functions

      tFjet11 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff11)
      tFjet12 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff12)
      tFjet13 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff13)
      tFjet14 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff14)

CCC Sum

      RVnum4e4 = tFjet11*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet14
     #*tcoff14
      RETURN
      END

      double precision function RVnum5e0(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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
      t2 = t1**2
      t3 = t2*t1
      t4 = l1**2
      t5 = t4**2
      t6 = t5**2
      t7 = t3*t6
      t9 = t1*z
      t10 = t9*t6
      t12 = t4*l1
      t13 = t5*t12
      t14 = t1*t13
      t16 = t2*t6
      t18 = t2*z
      t19 = t18*t13
      t21 = t18*t6
      t23 = 1088*z
      t24 = 3584*l1
      t25 = t5*l1
      t26 = t9*t25
      t28 = t2*t25
      t30 = t18*t25
      t32 = z*t25
      t34 = t9*t13
      t36 = z*t13
      t38 = 9216*t7-102912*t10-342272.D0/3.D0*t14+99328*t16+270848.D0/3.
     #D0*t19-47872*t21+t23-t24+493696.D0/3.D0*t26-391168.D0/3.D0*t28+121
     #792.D0/3.D0*t30+40384.D0/3.D0*t32+648448.D0/3.D0*t34+70912.D0/3.D0
     #*t36
      t39 = z*t6
      t41 = t2*t13
      t43 = t1*t6
      t45 = t3*t13
      t47 = z*l1
      t48 = 2816*t47
      t49 = t6*l1
      t50 = t49*t9
      t51 = 20480*t50
      t52 = t49*z
      t53 = 2048*t52
      t54 = t49*t18
      t55 = 10240*t54
      t56 = t49*t2
      t57 = 20480*t56
      t58 = t49*t3
      t59 = 2048*t58
      t60 = t49*t1
      t61 = 10240*t60
      t62 = t18*t4
      t63 = 1088.D0/3.D0*t62
      t64 = 1088.D0/3.D0*t9
      t65 = 2560*t4
      t66 = -11008*t39-600832.D0/3.D0*t41+53248*t43-47104.D0/3.D0*t45-t4
     #8+t51+t53+t55-t57-t59-t61+t63+t64+t65
      t68 = 2368.D0/3.D0*t1
      t69 = 256*t12
      t70 = t1*t25
      t72 = t5*t4
      t73 = t72*t1
      t75 = t3*t25
      t77 = t72*t18
      t79 = t9*t72
      t81 = t72*z
      t83 = t3*t72
      t85 = t72*t2
      t87 = t2*l1
      t88 = 2176.D0/3.D0*t87
      t89 = t1*t5
      t91 = t9*t5
      t93 = 4160.D0/3.D0+t68-t69-254464.D0/3.D0*t70+389632.D0/3.D0*t73-1
     #0240.D0/3.D0*t75-257536.D0/3.D0*t77-737792.D0/3.D0*t79-75776.D0/3.
     #D0*t81+35840.D0/3.D0*t83+645632.D0/3.D0*t85+t88+95168.D0/3.D0*t89-
     #196288.D0/3.D0*t91
      t94 = z*t5
      t96 = t2*t5
      t98 = t18*t5
      t100 = t18*t12
      t101 = 5824.D0/3.D0*t100
      t102 = t2*t4
      t103 = 11456.D0/3.D0*t102
      t104 = t1*l1
      t105 = 2304*t104
      t106 = z*t4
      t107 = 1920*t106
      t108 = t1*t4
      t109 = 2240*t108
      t110 = t9*l1
      t111 = 1088*t110
      t112 = t9*t4
      t114 = t1*t12
      t116 = t2*t12
      t118 = t9*t12
      t120 = z*t12
      t121 = 64*t120
      t122 = -3072*t94+119360.D0/3.D0*t96-5824*t98-t101-t103-t105+t107+t
     #109-t111-3136.D0/3.D0*t112-6400*t114-1280.D0/3.D0*t116+14848*t118+
     #t121
      t126 = t6*t4
      t135 = t2*t9
      t140 = 10752*t10
      t142 = 7168*t16
      t144 = 1792*t21
      t145 = 576*z
      t146 = 5504.D0/3.D0*l1
      t150 = 2176.D0/3.D0-t140-32000*t14+t142+7808*t19-t144+t145-t146+81
     #152*t26-47232*t28+24064.D0/3.D0*t30
      t153 = 8064*t36
      t154 = 1792*t39
      t156 = 7168*t43
      t157 = 4352.D0/3.D0*t47
      t158 = 448.D0/3.D0*t62
      t159 = 448.D0/3.D0*t9
      t160 = 1408*t4
      t161 = 1280.D0/3.D0*t1
      t162 = 10240*t32+47616*t34+t153-t154-31488*t41+t156-t157+t158+t159
     #+t160+t161
      t170 = 896.D0/3.D0*t87
      t174 = -t69-156544.D0/3.D0*t70+57216*t73-12480*t77-86016*t79-13632
     #*t81+54912*t85+t170+76928.D0/3.D0*t89-130496.D0/3.D0*t91-3008*t94
      t177 = 2944.D0/3.D0*t100
      t178 = 5504.D0/3.D0*t102
      t179 = 3968.D0/3.D0*t104
      t180 = 1024*t106
      t181 = 1536*t108
      t182 = 1280.D0/3.D0*t110
      t187 = 57856.D0/3.D0*t96-2176.D0/3.D0*t98-t177-t178-t179+t180+t181
     #-t182-6400.D0/3.D0*t112-19456.D0/3.D0*t114-3328.D0/3.D0*t116+13568
     #*t118
      t190 = s2l1**2
      t197 = (t38+t66+t93+t122)*s2l1-30720*t126*t9+157696*t7-397312*t10-
     #202752*t14+1593344.D0/3.D0*t16+520704*t19-397312*t21+11264*t49*t13
     #5-2048*t126*t135+(t150+t162+t174+t187)*t190+630272.D0/3.D0*t26-864
     #256.D0/3.D0*t28+630272.D0/3.D0*t30+28160.D0/3.D0*t32+520704*t34
      t213 = 1536*t10
      t215 = 1024*t16
      t217 = 256*t21
      t218 = 1984.D0/3.D0*z
      t219 = 7936.D0/3.D0*l1
      t225 = -960-t213-12160.D0/3.D0*t14+t215+2944.D0/3.D0*t19-t217-t218
     #+t219+7040*t26-6912*t28+10688.D0/3.D0*t30+832*t32+17408.D0/3.D0*t3
     #4
      t227 = 256*t39
      t229 = 1024*t43
      t231 = 5632.D0/3.D0*t47
      t252 = 3712*t30+8320.D0/3.D0*t32+13824*t34+2304*t36-512*t39-9216*t
     #41+2048*t43+5888.D0/3.D0*t47-1792*t4-t161-13824*t70
      t276 = 896.D0/3.D0*t100+512.D0/3.D0*t102+4352.D0/3.D0*t104-1152*t1
     #06-5888.D0/3.D0*t108-1664.D0/3.D0*t110+5888.D0/3.D0*t112+2560.D0/3
     #.D0*t114-5632.D0/3.D0*t116-512.D0/3.D0*t118-384*t120
      t278 = -1024-3072*t10-9216*t14+2048*t16+2304*t19-512*t21-2176.D0/3
     #.D0*z+8192.D0/3.D0*l1+23040*t26-47104.D0/3.D0*t28+t252+48896.D0/3.
     #D0*t73-12416.D0/3.D0*t77-25088*t79-11648.D0/3.D0*t81+50432.D0/3.D0
     #*t85+14336.D0/3.D0*t89-29824.D0/3.D0*t91-1280.D0/3.D0*t94+23296.D0
     #/3.D0*t96-1664*t98+t276
      t280 = 1280*t62
      t281 = 1536*t4
      t282 = 1088.D0/3.D0*t1
      t283 = 128*t2
      t285 = t3*t12
      t286 = 512*t285
      t287 = 3200.D0/3.D0*t36-t227-11264.D0/3.D0*t41+t229-128.D0/3.D0*t4
     #5+t231+t278*s2l1-t280-t281-t282-t69+t283-11648.D0/3.D0*t70-t286
      t297 = 128*t3*t4
      t298 = t3*t5
      t300 = 768*t87
      t302 = 2112*t91
      t303 = 448.D0/3.D0*t94
      t304 = 18176.D0/3.D0*t73-640*t75-6272.D0/3.D0*t77-8960*t79-4736.D0
     #/3.D0*t81+256*t83+18944.D0/3.D0*t85+t297+2432.D0/3.D0*t298-t300+44
     #8.D0/3.D0*t89-t302+t303
      t309 = 256*t18*l1
      t311 = 3712.D0/3.D0*t104
      t312 = 960*t106
      t313 = 3904.D0/3.D0*t108
      t314 = 832.D0/3.D0*t110
      t315 = 1472*t112
      t316 = 3328.D0/3.D0*t114
      t319 = 1600.D0/3.D0*t120
      t320 = 18112.D0/3.D0*t96-12800.D0/3.D0*t98+9280.D0/3.D0*t100+t309+
     #2496*t102+t311-t312-t313-t314+t315+t316-13696.D0/3.D0*t116-1408*t1
     #18-t319
      t325 = 2*l6-1
      t326 = t325**2
      t329 = l5-1
      t331 = -1+l1
      t333 = (1-l1+t47)**2
      t334 = 1/t333
      t339 = dsqrt(z)
      t341 = dsqrt(l5)
      t343 = dsqrt(-t329)
      t346 = dsqrt(l1)
      t350 = l5*z*l1-l1+1-l5+l5*l1
      t353 = dsqrt(-t331)
      t358 = t350**2
      t361 = dsqrt(4*l1*t326*z*l5*t329*t331*t334-4*t325*t339*t341*t343*t
     #346*t350*t334*t353+t358*t334)
      t365 = 12*t12-78*t4+104*l1-41
      t372 = 32256*t36-77824.D0/3.D0*t39-700416*t41+157696*t43-202752*t4
     #5+168960*t50+11264*t52+168960*t54-225280*t56-67584*t58-67584*t60-5
     #120.D0/3.D0*t135*t5+28160.D0/3.D0*t135*t25+(t225+t287+t304+t320)*s
     #1l1-64.D0/3.D0*t361*t365*s2l1-5120.D0/3.D0*t62-69632.D0/3.D0*t135*
     #t72
      t379 = -2176.D0/3.D0-t213-4864*t14+t215+1408*t19-t217-t145+t146+19
     #200*t26-44672.D0/3.D0*t28+11776.D0/3.D0*t30
      t382 = 1152*t36
      t384 = 5120.D0/3.D0*t32+7680*t34+t382-t227-5376*t41+t229+t157-t158
     #-t159-t160-t161
      t395 = t69-29824.D0/3.D0*t70+29056.D0/3.D0*t73-9664.D0/3.D0*t77-16
     #384*t79-6208.D0/3.D0*t81+35968.D0/3.D0*t85-t170+15232.D0/3.D0*t89-
     #39488.D0/3.D0*t91-1216.D0/3.D0*t94
      t402 = 34304.D0/3.D0*t96-2688*t98+t177+t178+t179-t180-t181+t182-12
     #80*t112-1024.D0/3.D0*t114-17152.D0/3.D0*t116+15616.D0/3.D0*t118
      t412 = 448.D0/3.D0*z
      t413 = 896*l1
      t417 = 896.D0/3.D0-t140-32512*t14+t142+8320*t19-t144+t412-t413+924
     #16*t26-175744.D0/3.D0*t28+34432.D0/3.D0*t30
      t421 = 512*t47
      t422 = 128*t62
      t423 = 128*t9
      t424 = 384*t4
      t425 = 30592.D0/3.D0*t32+49152*t34+t153-t154-33024*t41+t156-t421-t
     #422-t423+t424+t69
      t433 = 256*t87
      t438 = -166528.D0/3.D0*t70+178304.D0/3.D0*t73-43712.D0/3.D0*t77-92
     #416*t79-41024.D0/3.D0*t81+183680.D0/3.D0*t85-t433+82816.D0/3.D0*t8
     #9-163648.D0/3.D0*t91-7936.D0/3.D0*t94+30976*t96
      t440 = 384*t100
      t441 = 1024*t102
      t442 = 256.D0/3.D0*t104
      t443 = 448.D0/3.D0*t106
      t444 = 640.D0/3.D0*t108
      t445 = 640*t110
      t450 = 256.D0/3.D0*t120
      t451 = -11072.D0/3.D0*t98+t440+t441-t442+t443-t444+t445-15616.D0/3
     #.D0*t112-5888*t114-25600.D0/3.D0*t116+20480*t118+t450
      t454 = s1l1**2
      t466 = -4160.D0/3.D0-t213-5376*t14+t215+1536*t19-t217-t23+t24+2316
     #8*t26-57344.D0/3.D0*t28+16448.D0/3.D0*t30
      t471 = 5696.D0/3.D0*t32+8448*t34+1280*t36-t227-5888*t41+t229+t48-t
     #63-t64-t65-t68+t69
      t483 = -34304.D0/3.D0*t70+33280.D0/3.D0*t73-11776.D0/3.D0*t77-1894
     #4*t79-7168.D0/3.D0*t81+42496.D0/3.D0*t85-t88+5824*t89-16832*t91-10
     #24.D0/3.D0*t94+49600.D0/3.D0*t96
      t489 = -4416*t98+t101+t103+t105-t107-t109+t111-2368*t112-1280.D0/3
     #.D0*t114-29440.D0/3.D0*t116+22016.D0/3.D0*t118-t121
      t502 = (t379+t384+t395+t402)*t190+40960*t126*t2-2048*z*t126+(t417+
     #t425+t438+t451)*t454-226304.D0/3.D0*t70+64.D0/3.D0*t361*t365*s2l1+
     #(t466+t471+t483+t489)*s2l1-10240.D0/3.D0*t285+471040.D0/3.D0*t73-2
     #26304.D0/3.D0*t75-1251328.D0/3.D0*t77-1251328.D0/3.D0*t79-69632.D0
     #/3.D0*t81+471040.D0/3.D0*t83+1699840.D0/3.D0*t85+66560.D0/3.D0*t29
     #8
      t522 = 960-11264*t7+101888*t10+95872*t14-105472*t16-120448*t19+545
     #28*t21-t278*s2l1+t218-t219-156544*t26+583936.D0/3.D0*t28-299456.D0
     #/3.D0*t30-12736.D0/3.D0*t32-211968*t34
      t528 = -16768*t36+9472*t39+228352*t41-49152*t43+24960*t45-t231-t51
     #-t53-t55+t57+t59+t61+t280+t281+t282
      t540 = t69-t283+149888.D0/3.D0*t70+t286-287488.D0/3.D0*t73+16000*t
     #75+429184.D0/3.D0*t77+238336*t79+40576.D0/3.D0*t81-27904*t83-81356
     #8.D0/3.D0*t85-t297-4224*t298+t300-31168.D0/3.D0*t89
      t548 = 58432*t91-t303-264896.D0/3.D0*t96+125440.D0/3.D0*t98-9920*t
     #100-t309-17728.D0/3.D0*t102-t311+t312+t313+t314-t315-t316+26752*t1
     #16-8832*t118+t319
      t569 = -896.D0/3.D0-t213-4352*t14+t215+896*t19-t217-t412+t413+7936
     #*t26-10624.D0/3.D0*t28+1408.D0/3.D0*t30
      t573 = 5248.D0/3.D0*t32+6144*t34+t382-t227-3840*t41+t229+t421+t422
     #+t423-t424-t69
      t584 = -19840.D0/3.D0*t70+22400.D0/3.D0*t73-3392.D0/3.D0*t77-9984*
     #t79-6080.D0/3.D0*t81+17024.D0/3.D0*t85+t433+9344.D0/3.D0*t89-t302-
     #768*t94-256*t96
      t590 = 832.D0/3.D0*t98-t440-t441+t442-t443+t444-t445+1792*t112-281
     #6.D0/3.D0*t114+5120.D0/3.D0*t116-5120.D0/3.D0*t118-t450
      t596 = 66560.D0/3.D0*t89-197632.D0/3.D0*t91-5120.D0/3.D0*t94+26624
     #0.D0/3.D0*t96-197632.D0/3.D0*t98+12800*t100+(t522+t528+t540+t548)*
     #s1l1-77824.D0/3.D0*t135*t6+32256*t135*t13-5120.D0/3.D0*t112-10240.
     #D0/3.D0*t114-40960.D0/3.D0*t116+12800*t118+12288*t126*t3+12288*t12
     #6*t1+(t569+t573+t584+t590)*t454-30720*t126*t18
      t600 = (l1+z-t47)**2
      taux9 = (t197+t372+t502+t596)/t600*t334
      t1 = l1**2
      t2 = t1**2
      t3 = t2**2
      t4 = t3*l1
      t5 = z**2
      t6 = t5**2
      t8 = 81920*t4*t6
      t9 = t6*z
      t10 = t9*t3
      t12 = t6*t3
      t14 = 256*z
      t15 = t5*z
      t16 = 256*t15
      t17 = 2048*t1
      t18 = t1*l1
      t19 = t5*t18
      t21 = t15*t1
      t23 = 256*t5
      t24 = 1024*t18
      t25 = 2560.D0/3.D0*l1
      t26 = t15*t18
      t28 = t9*t18
      t29 = 4864.D0/3.D0*t28
      t30 = 256-t8-192512*t10+401408*t12+t14-t16-t17-92672.D0/3.D0*t19-6
     #4768.D0/3.D0*t21+t23+t24+t25+306176.D0/3.D0*t26+t29
      t31 = t9*t2
      t33 = z*t2
      t35 = t15*t2
      t37 = t5*t2
      t39 = t6*t5
      t40 = t2*t1
      t41 = t39*t40
      t43 = t40*z
      t45 = t15*t40
      t47 = t40*t9
      t49 = t2*l1
      t50 = t39*t49
      t52 = t40*t5
      t54 = t5*t49
      t56 = t15*t49
      t58 = t9*t49
      t60 = z*t49
      t62 = -44800*t31-41984.D0/3.D0*t33-1028864.D0/3.D0*t35+457984.D0/3
     #.D0*t37+143360.D0/3.D0*t41-326656.D0/3.D0*t43-3184640.D0/3.D0*t45-
     #361472*t47-40960.D0/3.D0*t50+561152*t52-383488*t54+2277376.D0/3.D0
     #*t56+188672*t58+179456.D0/3.D0*t60
      t64 = t2*t18
      t65 = t15*t64
      t67 = z*t64
      t69 = z*t3
      t71 = t5*t3
      t73 = t39*t64
      t75 = t4*t15
      t76 = 81920*t75
      t78 = 8192*t4*z
      t80 = 40960*t4*t9
      t81 = t4*t39
      t82 = 8192*t81
      t83 = t4*t5
      t84 = 40960*t83
      t85 = t15*l1
      t86 = 7936.D0/3.D0*t85
      t87 = z*t1
      t88 = 2560*t87
      t89 = t5*t1
      t90 = 1280*t89
      t91 = t5*l1
      t92 = 512.D0/3.D0*t91
      t93 = 2695168.D0/3.D0*t65+297472.D0/3.D0*t67-45056*t69+217088*t71-
     #188416.D0/3.D0*t73+t76+t78+t80-t82-t84+t86-t88+t90-t92
      t94 = z*t18
      t95 = 2304*t94
      t96 = t6*t18
      t98 = t6*t1
      t99 = 12032.D0/3.D0*t98
      t100 = t6*t2
      t102 = t6*l1
      t103 = 512*t102
      t104 = t40*t6
      t106 = t6*t49
      t108 = t6*t64
      t110 = z*l1
      t111 = 2560.D0/3.D0*t110
      t112 = t5*t64
      t114 = t9*t64
      t116 = t9*t1
      t117 = 256*t116
      t118 = t39*t3
      t120 = t15*t3
      t122 = t95-156160.D0/3.D0*t96+t99+240896*t100-t103+2768896.D0/3.D0
     #*t104-1831424.D0/3.D0*t106-2477056.D0/3.D0*t108+t111-1430528.D0/3.
     #D0*t112+1103360.D0/3.D0*t114-t117+36864*t118-417792*t120
      t126 = t6*t15
      t141 = (t30+t62+t93+t122)*s2l1-112640.D0/3.D0*t126*t49+20480.D0/3.
     #D0*t126*t2+204800*t10-40960.D0/3.D0*t19-20480.D0/3.D0*t21-243712.D
     #0/3.D0*t126*t64+51200*t26-8192*t4*t126-51200*t28+544768.D0/3.D0*t3
     #1
      t145 = t39*t2
      t152 = t39*t18
      t156 = -20480.D0/3.D0*t33-544768.D0/3.D0*t35+266240.D0/3.D0*t37-26
     #6240.D0/3.D0*t145-1146880.D0/3.D0*t41-237568.D0/3.D0*t43-1581056.D
     #0/3.D0*t45+1581056.D0/3.D0*t47+741376.D0/3.D0*t50+40960.D0/3.D0*t1
     #52+1146880.D0/3.D0*t52-741376.D0/3.D0*t54
      t169 = 1169408.D0/3.D0*t56-1169408.D0/3.D0*t58+112640.D0/3.D0*t60+
     #1316864.D0/3.D0*t65+243712.D0/3.D0*t67-40960*t69+237568.D0/3.D0*t1
     #26*t40+163840*t71+1024000.D0/3.D0*t73+40960*t75+t78
      t172 = 512*t6
      t173 = 1536*t19
      t174 = 27392.D0/3.D0*t21
      t178 = -256-t172-t14+t17+t173+t174-t23-t24-t25-55808.D0/3.D0*t26-3
     #7120.D0/3.D0*t28+52736.D0/3.D0*t31
      t179 = 1792.D0/3.D0*t33
      t188 = 2048*t152
      t191 = -t179+70912.D0/3.D0*t35-2304*t37-9728.D0/3.D0*t145-1024*t41
     #-4096.D0/3.D0*t43+8192*t45+20480.D0/3.D0*t47+2560*t50+t188+1024.D0
     #/3.D0*t52+2048.D0/3.D0*t54
      t199 = 2304*t85
      t200 = 2304*t87
      t201 = 1280.D0/3.D0*t89
      t202 = 1024.D0/3.D0*t91
      t203 = 4352.D0/3.D0*t94
      t205 = -18432*t56-44288.D0/3.D0*t58+1792*t60-5120.D0/3.D0*t65+1024
     #.D0/3.D0*t67+512.D0/3.D0*t73-t199+t200+t201-t202-t203+78848.D0/3.D
     #0*t96
      t208 = 1024*t9*l1
      t210 = 3072*t102
      t212 = 512*t39*t1
      t218 = 5120*t116
      t219 = -35072.D0/3.D0*t98-t208-35072*t100+t210-t212-38912.D0/3.D0*
     #t104+28160*t106+8192.D0/3.D0*t108-t111-512.D0/3.D0*t112-4096.D0/3.
     #D0*t114+t218
      t238 = -256-t14+t16+t17+3584*t19+7936*t21-t23-t24-t25-13312*t26-t2
     #9+3840*t31+1024.D0/3.D0*t33+41728.D0/3.D0*t35-7424.D0/3.D0*t37-512
     #0.D0/3.D0*t43+2048*t45+7168.D0/3.D0*t47+8192.D0/3.D0*t52-1536*t54
      t250 = -8192*t56-4352*t58+4864.D0/3.D0*t60+512*t67-t86+t88-t90+t92
     #-t95+33280.D0/3.D0*t96-t99-15616*t100+t103-16384.D0/3.D0*t104+3737
     #6.D0/3.D0*t106+1024*t108-t111-1024*t112-512*t114+t117
      t262 = 256-t8-217088*t10+417792*t12+t172+t14-t17-t173-t174+t23+t24
     #+t25+178688.D0/3.D0*t26+39680*t28-503296.D0/3.D0*t31
      t276 = t179-746752.D0/3.D0*t35+43264*t37+16896*t145+111616*t41-139
     #264.D0/3.D0*t43-925696*t45-1712128.D0/3.D0*t47-64000*t50-t188+1076
     #224.D0/3.D0*t52-555008.D0/3.D0*t54+616448*t56+1199360.D0/3.D0*t58+
     #35584.D0/3.D0*t60
      t283 = 826368*t65+62464*t67-36864*t69+192512*t71-99840*t73+t76+t78
     #+t80-t82-t84+t199-t200-t201+t202+t203
      t294 = -345088.D0/3.D0*t96+25344*t98+t208+1092352.D0/3.D0*t100-t21
     #0+t212+3217408.D0/3.D0*t104-2337280.D0/3.D0*t106-901120*t108+t111-
     #367104*t112+479232*t114-t218+45056*t118-401408*t120
      t298 = -t80+32768*t81-32768*t83+(t178+t191+t205+t219)*s1l1-1024000
     #.D0/3.D0*t112-1316864.D0/3.D0*t114+20480.D0/3.D0*t116+(t238+t250)*
     #s2l1+40960*t126*t3-163840*t118-204800*t120+(t262+t276+t283+t294)*s
     #1l1
      t302 = (l1+z-t110)**2
      t306 = (1-l1+t110)**2
      taux10 = (t141+t156+t169+t298)/t302/t306

CCC Coefficients

      tcoff15 = taux9
      tcoff16 = taux9
      tcoff17 = taux10
      tcoff18 = -taux10

CCC Jet functions

      tFjet15 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff15)
      tFjet16 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff16)
      tFjet17 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff17)
      tFjet18 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff18)

CCC Sum

      RVnum5e0 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet18
     #*tcoff18
      RETURN
      END

      double precision function RVnum5e1(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t2 = 2*l6-1
      t3 = t2**2
      t6 = l5-1
      t8 = -1+l1
      t9 = z*l1
      t11 = (1-l1+t9)**2
      t12 = 1/t11
      t17 = dsqrt(z)
      t19 = dsqrt(l5)
      t21 = dsqrt(-t6)
      t24 = dsqrt(l1)
      t25 = z*l5
      t26 = t25*l1
      t27 = l5*l1
      t28 = t26-l1+1-l5+t27
      t31 = dsqrt(-t8)
      t36 = t28**2
      t39 = dsqrt(4*l1*t3*z*l5*t6*t8*t12-4*t2*t17*t19*t21*t24*t28*t12*t3
     #1+t36*t12)
      t40 = l1**2
      t41 = t40**2
      t42 = t41*l1
      t44 = 243*t39*t42
      t45 = 30*l5
      t46 = t40*l1
      t47 = 35*t46
      t48 = t17*t19
      t49 = t48*t21
      t51 = t24*t31*l6
      t52 = t49*t51
      t53 = 120*t52
      t54 = l5**2
      t55 = t54*t40
      t56 = 115*t55
      t57 = t54*l1
      t58 = 75*t57
      t59 = t24*l1
      t61 = t21*t59*t31
      t62 = t48*t61
      t63 = 240*t62
      t64 = t39*z
      t65 = t64*t41
      t67 = t64*t40
      t69 = t25*t41
      t70 = 120*t69
      t71 = t25*t46
      t72 = 330*t71
      t73 = t54*z
      t74 = t73*t41
      t75 = 120*t74
      t76 = 150*t27
      t77 = 15-t44-t45-t47-t53+t56-t58-t63-150*t65+600*t67+t70+t72-t75+t
     #76
      t78 = t19*l5
      t79 = t17*t78
      t80 = t79*t21
      t82 = t59*t31*l6
      t83 = t80*t82
      t85 = t24*t40
      t87 = t21*t85*t31
      t88 = t79*t87
      t90 = t48*t87
      t92 = t41*l6
      t93 = t25*t92
      t95 = l5*t40
      t98 = t17*z*t78
      t99 = t98*t21
      t100 = t24*t46
      t102 = t100*t31*l6
      t103 = t99*t102
      t106 = t21*t100*t31
      t107 = t79*t106
      t109 = t24*t41
      t111 = t109*t31*l6
      t114 = t79*t61
      t116 = t24*t42
      t118 = t116*t31*l6
      t121 = t46*l6
      t122 = t73*t121
      t124 = t73*t92
      t126 = l6**2
      t127 = t41*t126
      t128 = t73*t127
      t130 = t41*t40
      t131 = t130*l6
      t134 = -480*t83-220*t88+220*t90-320*t93-230*t95-80*t103-80*t107+33
     #6*t49*t111+240*t114-96*t49*t118+880*t122+320*t124-320*t128-192*t25
     #*t131
      t136 = t130*t126
      t138 = 192*t73*t136
      t140 = 192*t73*t131
      t142 = 192*t25*t136
      t143 = t42*l6
      t145 = 672*t25*t143
      t146 = t39*l1
      t148 = t39*t40
      t151 = t21*t116*t31
      t153 = 48*t79*t151
      t155 = 96*t80*t118
      t156 = t54*t41
      t157 = 62*t156
      t158 = l5*t41
      t159 = 124*t158
      t160 = t54*t46
      t161 = 35*t160
      t163 = t85*t31*l6
      t164 = t49*t163
      t165 = 440*t164
      t167 = 54*t39*t130
      t168 = t98*t87
      t169 = 180*t168
      t170 = -t138+t140+t142+t145+375*t146-585*t148-t153+t155-t157+t159-
     #t161-t165+t167-t169
      t171 = t98*t61
      t172 = 60*t171
      t173 = t99*t82
      t174 = 120*t173
      t175 = z**2
      t176 = t54*t175
      t177 = t176*t40
      t178 = 15*t177
      t179 = t73*t40
      t180 = 360*t179
      t181 = t25*t40
      t182 = 360*t181
      t183 = 90*t26
      t184 = t99*t163
      t185 = 360*t184
      t186 = t39*t41
      t188 = t39*t46
      t190 = t98*t106
      t191 = 40*t190
      t193 = 48*t48*t151
      t194 = 75*l1
      t195 = 54*t42
      t197 = t21*t109*t31
      t199 = 168*t48*t197
      t201 = 336*t80*t111
      t202 = t172-t174+t178+t180-t182+t183+t185+264*t186+210*t188+t191+t
     #193-t194+t195-t199-t201
      t207 = 240*t99*t111
      t208 = t48*t106
      t209 = 80*t208
      t210 = t40*t126
      t211 = t25*t210
      t212 = 960*t211
      t213 = l1*t126
      t214 = t25*t213
      t215 = 240*t214
      t216 = t73*t213
      t217 = 240*t216
      t218 = t40*l6
      t219 = t73*t218
      t220 = 960*t219
      t221 = l1*l6
      t222 = t73*t221
      t223 = 240*t222
      t224 = t73*t210
      t225 = 960*t224
      t226 = t25*t221
      t227 = 240*t226
      t228 = t25*t218
      t229 = 960*t228
      t230 = t46*t126
      t231 = t73*t230
      t232 = 880*t231
      t233 = t25*t230
      t234 = 880*t233
      t235 = t25*t127
      t236 = 320*t235
      t237 = -75*t39-t207+t209-t212+t215-t217-t220+t223+t225-t227+t229-t
     #232+t234+t236
      t238 = t25*t121
      t241 = t21*t24*t31
      t242 = t48*t241
      t244 = t79*t241
      t246 = t80*t51
      t252 = t176*t46
      t254 = t73*t46
      t258 = t49*t82
      t261 = t49*t102
      t263 = t80*t102
      t266 = -880*t238+60*t242-60*t244+120*t246+168*t79*t197+120*t98*t19
     #7-45*t252-330*t254+96*t99*t118+480*t258-62*t41-160*t261+160*t263+1
     #5*t54
      t269 = 54*t54*t42
      t271 = 108*l5*t42
      t273 = 12*t54*t130
      t275 = 24*l5*t130
      t276 = t80*t163
      t277 = 440*t276
      t278 = 115*t40
      t279 = 12*t130
      t280 = t39*t175
      t281 = t280*t40
      t284 = 48*t98*t151
      t285 = t64*t46
      t287 = t280*t46
      t289 = t64*l1
      t291 = t176*t41
      t292 = 10*t291
      t293 = t280*t41
      t295 = t269-t271-t273+t275+t277+t278-t279-75*t281-t284-570*t285+22
     #5*t287-150*t289+t292-60*t293
      t296 = t73*l1
      t306 = l5*t46
      t322 = t42*t126
      t327 = -90*t296+252*t73*t42-252*t25*t42-72*t73*t130+72*t25*t130+70
     #*t306-12*t176*t130-135*t280*t42+378*t64*t42+54*t280*t130-108*t64*t
     #130+30*t176*t42-672*t73*t143+672*t73*t322-672*t25*t322
      t334 = 2944*t46
      t335 = t175**2
      t336 = t335*z
      t337 = t336*t42
      t339 = z*t42
      t341 = t175*z
      t342 = t41*t46
      t343 = t341*t342
      t345 = z*t342
      t347 = t41**2
      t348 = z*t347
      t350 = t335*t342
      t352 = t175*t347
      t355 = t347*l1
      t356 = t355*t341
      t358 = t347*t40
      t359 = z*t358
      t361 = t358*t341
      t363 = t358*t335
      t365 = 1920+640*t175+t334-4736.D0/3.D0*t337+57088.D0/3.D0*t339+137
     #344.D0/3.D0*t343+46208.D0/3.D0*t345-24064.D0/3.D0*t348-73856.D0/3.
     #D0*t350+81664.D0/3.D0*t352-3840*t9+16896*t356-512*t359-3072*t361+2
     #048*t363
      t366 = t358*t175
      t368 = t355*z
      t370 = t355*t336
      t372 = t355*t335
      t374 = t335*t46
      t376 = t341*t46
      t378 = t335*t40
      t380 = t336*t46
      t382 = t336*t41
      t384 = t335*t41
      t386 = z*t41
      t388 = t341*t41
      t390 = t175*t41
      t392 = t335*l1
      t393 = 512*t392
      t394 = t130*t335
      t396 = 2048*t366+2816*t368+2816*t370-11264*t372-15872.D0/3.D0*t374
     #+52736.D0/3.D0*t376+2560*t378+2176.D0/3.D0*t380+128.D0/3.D0*t382+1
     #3952.D0/3.D0*t384-37760.D0/3.D0*t386-60160.D0/3.D0*t388+77696.D0/3
     #.D0*t390-t393+34048.D0/3.D0*t394
      t398 = t130*z
      t400 = t341*t130
      t402 = t130*t336
      t404 = t130*t175
      t406 = t175*t42
      t408 = t341*t42
      t410 = t335*t42
      t412 = z*t46
      t414 = t341*t40
      t417 = t358*t336
      t419 = t341*t347
      t421 = t175*t342
      t423 = t335*t347
      t425 = t336*t342
      t427 = -20608*t398-101888.D0/3.D0*t400-640.D0/3.D0*t402+118784.D0/
     #3.D0*t404-34560*t406+66688.D0/3.D0*t408-3072*t410+6400*t412-28544.
     #D0/3.D0*t414+512*t347-512*t417-38144*t419-116608.D0/3.D0*t421+7244
     #8.D0/3.D0*t423+14080.D0/3.D0*t425
      t428 = t336*t347
      t430 = t336*t40
      t431 = 256*t430
      t433 = t175*l1
      t435 = t355*t175
      t438 = t175*t40
      t441 = t341*l1
      t443 = t175*t46
      t445 = 1024.D0/3.D0*t41
      t450 = z*t40
      t452 = -17152.D0/3.D0*t428-t431+1280*z-2688*t433-11264*t435-5760*l
     #1+16640.D0/3.D0*t438-3840*t42+2432*t441-40832.D0/3.D0*t443+t445-71
     #68.D0/3.D0*t342-256*t341+5120.D0/3.D0*t40+4352*t130+512*t450
      t457 = 11584.D0/3.D0*t46
      t468 = 1024*t361
      t469 = 1536*t363
      t470 = 256*t366
      t471 = 512*t368
      t472 = (t365+t396+t427+t452)*s2l1+1664.D0/3.D0*t175+t457-384*t335+
     #10304.D0/3.D0*t337-2944*t339+23872.D0/3.D0*t343+7616.D0/3.D0*t345-
     #2048*t348-11840.D0/3.D0*t350+15808.D0/3.D0*t352-10880.D0/3.D0*t9-1
     #792*t356+t468-t469-t470+t471
      t485 = t335*t175
      t486 = t485*t41
      t489 = t485*t130
      t493 = -3840*t370+4864*t372+25408.D0/3.D0*t374-2048*t376-5696*t378
     #-19904.D0/3.D0*t380+14080.D0/3.D0*t382-15680.D0/3.D0*t384+5632.D0/
     #3.D0*t386+7552*t388-5824*t390+1920*t392-1472*t486+26240.D0/3.D0*t3
     #94+6272.D0/3.D0*t489+256.D0/3.D0*t398-14336.D0/3.D0*t400
      t496 = t485*t42
      t498 = t485*t46
      t499 = 1408*t498
      t506 = 1024*t417
      t512 = 3456*t430
      t513 = 5504.D0/3.D0-22528.D0/3.D0*t402-1856.D0/3.D0*t496+t499+6400
     #.D0/3.D0*t404+14336.D0/3.D0*t406-13888.D0/3.D0*t408-10048.D0/3.D0*
     #t410+2048*t412-3584.D0/3.D0*t414+t506-7168.D0/3.D0*t419-23296.D0/3
     #.D0*t421-3968*t423+5632.D0/3.D0*t425+12800.D0/3.D0*t428+t512
      t516 = t355*t485
      t518 = t485*t347
      t520 = t358*t485
      t521 = 256*t520
      t525 = 256*t42
      t528 = t485*t342
      t530 = 768*t41
      t533 = t336*l1
      t534 = 768*t533
      t535 = t485*t40
      t536 = 384*t535
      t538 = 3584.D0/3.D0*z-5120.D0/3.D0*t433+1024*t516-3392.D0/3.D0*t51
     #8-t521-768*t435-16640.D0/3.D0*l1-2944.D0/3.D0*t438-t525+1088.D0/3.
     #D0*t441+13376.D0/3.D0*t443-1984.D0/3.D0*t528-t530-256.D0/3.D0*t341
     #+1792.D0/3.D0*t40-t534-t536+512.D0/3.D0*t450
      t543 = 12*t42
      t544 = 78*t41
      t545 = 104*t46
      t553 = 4736.D0/3.D0*t46
      t563 = 9472*t359
      t564 = 59136*t361
      t565 = 64384*t363
      t566 = 32128*t366
      t571 = -2240.D0/3.D0*t175-t553-1628672.D0/3.D0*t337-124288.D0/3.D0
     #*t339-3627392.D0/3.D0*t343-360896.D0/3.D0*t345+301184.D0/3.D0*t348
     #+4421696.D0/3.D0*t350-1217408.D0/3.D0*t352+9856.D0/3.D0*t9-336000*
     #t356+t563+t564-t565-t566-47232*t368-257152*t370+382656*t372+65216*
     #t374
      t574 = t335*t341
      t575 = t355*t574
      t594 = -62720*t376-17600.D0/3.D0*t378-14976*t575-107840.D0/3.D0*t3
     #80+188608*t382-908800.D0/3.D0*t384+11136*t386+758336.D0/3.D0*t388-
     #93568*t390+128*t392-55680*t486-4098560.D0/3.D0*t394-339456*t489+89
     #792*t398+3301760.D0/3.D0*t400+943040*t402+182208*t496+7680*t498-14
     #19776.D0/3.D0*t404+794048.D0/3.D0*t406
      t598 = 7040.D0/3.D0*t412
      t600 = 1152*t358
      t603 = 41728*t417
      t614 = -4160.D0/3.D0-653184*t408+2434304.D0/3.D0*t410-t598+5184*t3
     #55-t600+28928.D0/3.D0*t414-9088*t347+t603+840064*t419+1650176.D0/3
     #.D0*t421-2981888.D0/3.D0*t423-3074240.D0/3.D0*t425+2052608.D0/3.D0
     #*t428+9472.D0/3.D0*t430-3200.D0/3.D0*z+9280.D0/3.D0*t433+95040*t51
     #6-256512*t518
      t615 = t358*t574
      t616 = 2304*t615
      t617 = 14976*t520
      t618 = t574*t42
      t621 = t574*t41
      t625 = 1280.D0/3.D0*t42
      t627 = t574*t342
      t629 = t574*t347
      t632 = t574*t130
      t635 = 1280*t41
      t637 = 128.D0/3.D0*t341
      t641 = t616-t617-21120*t618+172480*t435+3840*t621+12736.D0/3.D0*l1
     #-14848.D0/3.D0*t438+t625-896.D0/3.D0*t441-59136*t627+40704*t629+60
     #032.D0/3.D0*t443+48384*t632+381696*t528+t635+22976.D0/3.D0*t342-t6
     #37-2304*t40-9472.D0/3.D0*t130-5120.D0/3.D0*t450
      t644 = s2l1**2
      t647 = 2048.D0/3.D0*t46
      t648 = 64.D0/3.D0*t335
      t659 = -576-64*t175-t647-t648-108864*t337-112320*t339-2969216.D0/3
     #.D0*t343-221120*t345+447104.D0/3.D0*t348+876032*t350-1344128.D0/3.
     #D0*t352+2048.D0/3.D0*t9-314496*t356+t563+t564-t565-t566-56960*t368
      t666 = 1152*t380
      t679 = -201856*t370+325568*t372+24768*t374-195008.D0/3.D0*t376-928
     #0.D0/3.D0*t378-10368*t575-t666+20096*t382-364160.D0/3.D0*t384+3552
     #0*t386+218560*t388-441472.D0/3.D0*t390+2240.D0/3.D0*t392-2062336.D
     #0/3.D0*t394-57600*t489+604736.D0/3.D0*t398+868480*t400+288704*t402
     #+11520*t496
      t694 = 384*t430
      t698 = -1770752.D0/3.D0*t404+1081088.D0/3.D0*t406-521088*t408+1071
     #424.D0/3.D0*t410-7296*t412+7488*t355-t600+42112.D0/3.D0*t414-20608
     #*t347+t603+732544*t419+1927616.D0/3.D0*t421-2125568.D0/3.D0*t423-1
     #339328.D0/3.D0*t425+1223168.D0/3.D0*t428+t694-256*z+1088*t433+6969
     #6*t516
      t712 = 1024.D0/3.D0*t341
      t716 = -129792*t518+t616-t617+180928*t435+5696.D0/3.D0*l1-17920.D0
     #/3.D0*t438+14912*t42-5888.D0/3.D0*t441-13440*t627+17664*t629+38784
     #*t443+3840*t632+121152*t528-14336.D0/3.D0*t41+93056.D0/3.D0*t342+t
     #712+640.D0/3.D0*t40-82688.D0/3.D0*t130+1792*t450
      t719 = s1l1**2
      t723 = 15-t44-t45-t47-t53+t56-t58-t63-138*t65+620*t67+t70+t72-t75+
     #t76
      t727 = -t138+t140+t142+t145+386*t146-610*t148-t153+t155-t157+t159-
     #t161-t165+t167-t169
      t730 = t172-t174+t178+t180-t182+t183+t185+258*t186+231*t188+t191+t
     #193-t194+t195-t199-t201
      t734 = -76*t39-t207+t209-t212+t215-t217-t220+t223+t225-t227+t229-t
     #232+t234+t236
      t740 = 66*t293
      t741 = t269-t271-t273+t275+t277+t278-t279-76*t281-t284-600*t285+23
     #4*t287-152*t289+t292-t740
      t750 = 64.D0/3.D0*(t77+t134+t170+t202+t237+t266+t295+t327)*t12+(t4
     #72+t493+t513+t538)*s1l1-64.D0/3.D0*t175-64.D0/3.D0*t39*(t543-t544+
     #t545+112*t40-235*l1+79)*s2l1+(t571+t594+t614+t641)*t644+(t659+t679
     #+t698+t716)*t719-t648-64.D0/3.D0*(t723+t134+t727+t730+t734+t266+t7
     #41+t327)*t12-3204544.D0/3.D0*t337-40320*t339-2510784*t343
      t762 = -151488*t345+372032.D0/3.D0*t348+3401856*t350-756992*t352+1
     #92*t9-813120*t356+9856*t359+147840*t361-197120*t363-59136*t366-542
     #08*t368
      t766 = dlog(z)
      t767 = 1088.D0/3.D0*t175
      t773 = 256*t46
      t782 = 20480*t356
      t783 = 2048*t368
      t784 = 10240*t370
      t785 = 20480*t372
      t786 = 2368.D0/3.D0*t175-t773+121792.D0/3.D0*t337+40384.D0/3.D0*t3
     #39+648448.D0/3.D0*t343+70912.D0/3.D0*t345-11008*t348-600832.D0/3.D
     #0*t350+53248*t352-2816*t9+t782+t783+t784-t785
      t801 = -1280.D0/3.D0*t374+14848*t376-11456.D0/3.D0*t378-5824.D0/3.
     #D0*t380-5824*t382+119360.D0/3.D0*t384-3072*t386-196288.D0/3.D0*t38
     #8+95168.D0/3.D0*t390+2176.D0/3.D0*t392+645632.D0/3.D0*t394+35840.D
     #0/3.D0*t489-75776.D0/3.D0*t398-737792.D0/3.D0*t400
      t817 = -257536.D0/3.D0*t402-10240.D0/3.D0*t496+389632.D0/3.D0*t404
     #-254464.D0/3.D0*t406+493696.D0/3.D0*t408-391168.D0/3.D0*t410+64*t4
     #12-3136.D0/3.D0*t414-102912*t419-342272.D0/3.D0*t421+99328*t423+27
     #0848.D0/3.D0*t425-47872*t428+1088.D0/3.D0*t430
      t820 = 2048*t516
      t822 = 10240*t435
      t828 = 1088.D0/3.D0*t341
      t831 = 4160.D0/3.D0+1088*z-2304*t433-t820+9216*t518-t822-3584*l1+2
     #240*t438-1088*t441-6400*t443-47104.D0/3.D0*t528+t828+2560*t40+1920
     #*t450
      t835 = 1280.D0/3.D0*t175
      t839 = 8064*t345
      t840 = 1792*t348
      t842 = 7168*t352
      t844 = 2176.D0/3.D0+t835-t773+24064.D0/3.D0*t337+10240*t339+47616*
     #t343+t839-t840-31488*t350+t842-4352.D0/3.D0*t9
      t856 = -3328.D0/3.D0*t374+13568*t376-5504.D0/3.D0*t378-2944.D0/3.D
     #0*t380-2176.D0/3.D0*t382+57856.D0/3.D0*t384-3008*t386-130496.D0/3.
     #D0*t388+76928.D0/3.D0*t390+896.D0/3.D0*t392+54912*t394
      t866 = 10752*t419
      t868 = 7168*t423
      t869 = -13632*t398-86016*t400-12480*t402+57216*t404-156544.D0/3.D0
     #*t406+81152*t408-47232*t410-6400.D0/3.D0*t414-t866-32000*t421+t868
      t871 = 1792*t428
      t879 = 448.D0/3.D0*t341
      t882 = 7808*t425-t871+448.D0/3.D0*t430+576*z-3968.D0/3.D0*t433-550
     #4.D0/3.D0*l1+1536*t438-1280.D0/3.D0*t441-19456.D0/3.D0*t443+t879+1
     #408*t40+1024*t450
      t901 = t767-64*t574*t40+320*t574*t46+(t786+t801+t817+t831)*s2l1+(t
     #844+t856+t869+t882)*t644+216064*t337+8192*t339+1546816.D0/3.D0*t34
     #3+93632.D0/3.D0*t345-77120.D0/3.D0*t348-698368*t350+156800*t352-40
     #96.D0/3.D0*t9+168960*t356-2048*t359-30720*t361+40960*t363+12288*t3
     #66+11264*t368+168960*t370
      t905 = 128*t335
      t926 = 5632.D0/3.D0*t374+512.D0/3.D0*t376-512.D0/3.D0*t378-896.D0/
     #3.D0*t380+1664*t382-23296.D0/3.D0*t384+1280.D0/3.D0*t386+29824.D0/
     #3.D0*t388-14336.D0/3.D0*t390-50432.D0/3.D0*t394+11648.D0/3.D0*t398
      t950 = -2048*t423-2304*t425+512*t428+2176.D0/3.D0*z-4352.D0/3.D0*t
     #433-8192.D0/3.D0*l1+5888.D0/3.D0*t438+1664.D0/3.D0*t441-2560.D0/3.
     #D0*t443+1792*t40+1152*t450
      t962 = 960+t767+t773-t905+(1024+t835-3712*t337-8320.D0/3.D0*t339-1
     #3824*t343-2304*t345+512*t348+9216*t350-2048*t352-5888.D0/3.D0*t9+t
     #926+25088*t400+12416.D0/3.D0*t402-48896.D0/3.D0*t404+13824*t406-23
     #040*t408+47104.D0/3.D0*t410+384*t412-5888.D0/3.D0*t414+3072*t419+9
     #216*t421+t950)*s2l1-299456.D0/3.D0*t337-12736.D0/3.D0*t339-211968*
     #t343-16768*t345+9472*t348+228352*t350-49152*t352-5632.D0/3.D0*t9-t
     #782-t783
      t972 = 768*t392
      t976 = -t784+t785+26752*t374-8832*t376-17728.D0/3.D0*t378-9920*t38
     #0+125440.D0/3.D0*t382-264896.D0/3.D0*t384-448.D0/3.D0*t386+58432*t
     #388-31168.D0/3.D0*t390+t972-4224*t486-813568.D0/3.D0*t394-27904*t4
     #89
      t993 = 40576.D0/3.D0*t398+238336*t400+429184.D0/3.D0*t402+16000*t4
     #96+512*t498-287488.D0/3.D0*t404+149888.D0/3.D0*t406-156544*t408+58
     #3936.D0/3.D0*t410+1600.D0/3.D0*t412-1472*t414+101888*t419+95872*t4
     #21-105472*t423-120448*t425
      t1005 = 256*t533
      t1006 = 128*t535
      t1008 = 54528*t428+1280*t430+1984.D0/3.D0*z-3712.D0/3.D0*t433+t820
     #-11264*t518+t822-7936.D0/3.D0*l1+3904.D0/3.D0*t438+832.D0/3.D0*t44
     #1-3328.D0/3.D0*t443+24960*t528+1536*t40-t1005-t1006+960*t450
      t1028 = 1984.D0/3.D0-225280*t372-128*t485*l1+(t962+t976+t993+t1008
     #)*s1l1-14080*t374+34432.D0/3.D0*t376-64*t378+11264*t575+16896*t380
     #-216320.D0/3.D0*t382+268352.D0/3.D0*t384-4672.D0/3.D0*t386-184064.
     #D0/3.D0*t388+19648*t390+128.D0/3.D0*t392+25984*t486+1692544.D0/3.D
     #0*t394+159488*t489-64448.D0/3.D0*t398-408256*t400
      t1044 = 896.D0/3.D0+t773+34432.D0/3.D0*t337+30592.D0/3.D0*t339+491
     #52*t343+t839-t840-33024*t350+t842-512*t9-25600.D0/3.D0*t374
      t1056 = 20480*t376+1024*t378+384*t380-11072.D0/3.D0*t382+30976*t38
     #4-7936.D0/3.D0*t386-163648.D0/3.D0*t388+82816.D0/3.D0*t390-256*t39
     #2+183680.D0/3.D0*t394-41024.D0/3.D0*t398
      t1067 = -92416*t400-43712.D0/3.D0*t402+178304.D0/3.D0*t404-166528.
     #D0/3.D0*t406+92416*t408-175744.D0/3.D0*t410+256.D0/3.D0*t412-15616
     #.D0/3.D0*t414-t866-32512*t421+t868
      t1079 = 8320*t425-t871-128*t430+448.D0/3.D0*z-256.D0/3.D0*t433-896
     #*l1-640.D0/3.D0*t438+640*t441-5888*t443-128*t341+384*t40+448.D0/3.
     #D0*t450
      t1094 = -420416*t402-238208.D0/3.D0*t496-5632*t498+150144*t404-693
     #76*t406+605312.D0/3.D0*t408-287360*t410+192*t412+(t1044+t1056+t106
     #7+t1079)*t719-3904.D0/3.D0*t414-30720*t417-396096*t419-596480.D0/3
     #.D0*t421+1591552.D0/3.D0*t423+1564736.D0/3.D0*t425-397376*t428-103
     #04.D0/3.D0*t430+512*z-3328.D0/3.D0*t433-67584*t516
      t1114 = 157824*t518-2048*t615+12288*t520+10112*t618-67584*t435-236
     #8*t621-5248.D0/3.D0*l1+1088*t438-1216.D0/3.D0*t441+97216.D0/3.D0*t
     #627-77888.D0/3.D0*t629-10240.D0/3.D0*t443-70976.D0/3.D0*t632-61081
     #6.D0/3.D0*t528-64*t336+t879+1152*t40+448*t533+768*t535+832*t450
      t1134 = 256*t359
      t1135 = 1536*t361
      t1136 = 1024*t363
      t1137 = 1024*t366
      t1138 = 128.D0/3.D0*t175+t647-3968.D0/3.D0*t337+42304.D0/3.D0*t339
     #+65024.D0/3.D0*t343+9152*t345-14720.D0/3.D0*t348-8960*t350+47744.D
     #0/3.D0*t352-1472.D0/3.D0*t9+9216*t356-t1134-t1135+t1136+t1137
      t1140 = 1408*t370
      t1147 = 11712*t386
      t1152 = 1664*t368+t1140-5888*t372-14848.D0/3.D0*t374+25024.D0/3.D0
     #*t376+3072*t378+t666-2048.D0/3.D0*t382+4736.D0/3.D0*t384-t1147-569
     #6*t388+15488*t390-t972-1152*t394-38464.D0/3.D0*t398
      t1162 = 256*t347
      t1163 = 256*t417
      t1168 = 7424.D0/3.D0*t428
      t1169 = -10496*t400+4672.D0/3.D0*t402+59392.D0/3.D0*t404-17536*t40
     #6+3264*t408+12544.D0/3.D0*t410+23168.D0/3.D0*t412-17216.D0/3.D0*t4
     #14+t1162-t1163-20608*t419-64256.D0/3.D0*t421+35456.D0/3.D0*t423+30
     #08.D0/3.D0*t425-t1168
      t1183 = 1664.D0/3.D0-t694+704.D0/3.D0*z-896*t433-6400*t435-5120.D0
     #/3.D0*l1+9856.D0/3.D0*t438-4352*t42+5824.D0/3.D0*t441-9344*t443+85
     #76.D0/3.D0*t41-4352.D0/3.D0*t342-t828-1024.D0/3.D0*t40+10112.D0/3.
     #D0*t130-8320.D0/3.D0*t450
      t1192 = t17*t175
      t1193 = t1192*t19
      t1197 = t1192*t78
      t1198 = t1197*t21
      t1203 = t341*t126
      t1206 = t341*l6
      t1218 = l5*t175
      t1222 = t17*t341*t78
      t1226 = -16*t27*t1206-16*t57*t1203-8*t1197*t61+4*t1193*t241-2*t121
     #8+t176+2*t73-t160+4*t1222*t61-t433+4*t171
      t1249 = t54*t335
      t1251 = t54*t341
      t1262 = -8*t1222*t21*t82-4*t1197*t241+t1249*t40-2*t1251*t40+t176*l
     #1+4*t1218*t40-6*t252+4*t254-t1249*t46+4*t1251*t46+t54
      t1264 = t175+3*t55-3*t57-2*t25-4*t98*t241-8*t1193*t21*t51+16*t1198
     #*t82+8*t99*t51+16*t27*t1203+16*t57*t1206+t1226-8*t173+8*t1198*t51+
     #4*t177-6*t179-6*t57*t341+6*t27*t341-2*t181-2*t95*t341+4*t26-2*t121
     #8*l1+t1262
      t1272 = -64.D0/3.D0-813120*t370+1084160*t372-t766*(t901+t1028+t109
     #4+t1114)-64.D0/3.D0*t39*(t544-t545-118*t40-t543+244*l1-80)*s2l1+(t
     #1138+t1152+t1169+t1183)*t719+64.D0/3.D0*t1264*t335*t12+80000*t374-
     #164480.D0/3.D0*t376-14464.D0/3.D0*t378-54208*t575
      t1284 = -66880*t380+1035200.D0/3.D0*t382-1378688.D0/3.D0*t384+2048
     #0.D0/3.D0*t386+934976.D0/3.D0*t388-283520.D0/3.D0*t390+192*t392-34
     #1120.D0/3.D0*t486-8352512.D0/3.D0*t394-780544*t489+315968.D0/3.D0*
     #t398
      t1305 = 2176.D0/3.D0*t175+t553-6336*t337+7808.D0/3.D0*t339+10112.D
     #0/3.D0*t343+8384.D0/3.D0*t345-7040.D0/3.D0*t348-5504.D0/3.D0*t350+
     #24704.D0/3.D0*t352-9280.D0/3.D0*t9+7680*t356-t1134-t1135+t1136+t11
     #37
      t1319 = 4096.D0/3.D0+1152*t368+t1140-5376*t372+16256.D0/3.D0*t374-
     #1856.D0/3.D0*t376-3200.D0/3.D0*t378-2560.D0/3.D0*t380+10816.D0/3.D
     #0*t382-43904.D0/3.D0*t384-3456*t386+39232.D0/3.D0*t388-256.D0/3.D0
     #*t390+256.D0/3.D0*t392-43136.D0/3.D0*t394
      t1333 = -7232.D0/3.D0*t398+50944.D0/3.D0*t400+13888.D0/3.D0*t402-1
     #7408.D0/3.D0*t404+7040*t406-24512*t408+21504*t410+t598-6208.D0/3.D
     #0*t414+t1162-t1163-12928*t419-10880.D0/3.D0*t421+27776.D0/3.D0*t42
     #3+704.D0/3.D0*t425
      t1346 = -t1168+128.D0/3.D0*t430+3136.D0/3.D0*z-8704.D0/3.D0*t433-4
     #864*t435-12160.D0/3.D0*l1+14464.D0/3.D0*t438-t625+512*t441-13952.D
     #0/3.D0*t443-t635-2816.D0/3.D0*t342+t637+2176*t40+3712.D0/3.D0*t130
     #+4736.D0/3.D0*t450
      t1351 = 11456.D0/3.D0*t46
      t1361 = 2304*t359
      t1362 = 33536*t361
      t1363 = 44544*t363
      t1364 = 13568*t366
      t1367 = -2688-1408*t175-t1351+t905-2319296.D0/3.D0*t337-45248*t339
     #-4498688.D0/3.D0*t343-321856.D0/3.D0*t345+65216*t348+5667520.D0/3.
     #D0*t350-1119040.D0/3.D0*t352+19072.D0/3.D0*t9-275968*t356+t1361+t1
     #362-t1363-t1364-19968*t368-256256*t370
      t1386 = 354048*t372+231680.D0/3.D0*t374-70656*t376-4352.D0/3.D0*t3
     #78-16128*t575-52864*t380+273536*t382-1250944.D0/3.D0*t384+t1147+32
     #5824*t388-111360*t390-3136.D0/3.D0*t392-245696.D0/3.D0*t486-565312
     #0.D0/3.D0*t394-1373824.D0/3.D0*t489+287104.D0/3.D0*t398+4456768.D0
     #/3.D0*t400+1296128*t402+769792.D0/3.D0*t496
      t1393 = 11456.D0/3.D0*t412
      t1395 = 33536*t417
      t1406 = 36160.D0/3.D0*t498-599808*t404+339072*t406-2682304.D0/3.D0
     #*t408+3445760.D0/3.D0*t410-t1393+25280.D0/3.D0*t414+t1395+2658496.
     #D0/3.D0*t419+1877120.D0/3.D0*t421-1119616*t423-3959168.D0/3.D0*t42
     #5+2381248.D0/3.D0*t428+4736*t430-2048*z+15616.D0/3.D0*t433+99328*t
     #516-898624.D0/3.D0*t518+t616
      t1407 = 13568*t520
      t1419 = 2304*t41
      t1423 = -t1407-29760*t618+114944*t435+5760*t621+24832.D0/3.D0*l1-1
     #7152.D0/3.D0*t438-t525+2944.D0/3.D0*t441-73152*t627+47040*t629+200
     #96*t443+63936*t632+1454528.D0/3.D0*t528+t1419-384*t341-10496.D0/3.
     #D0*t40+t1005+t1006-2688*t450
      t1428 = 2014400*t400+6205760.D0/3.D0*t402+383040*t496+51200.D0/3.D
     #0*t498-2216192.D0/3.D0*t404+1032640.D0/3.D0*t406-1012928*t408+4325
     #888.D0/3.D0*t410+(t1305+t1319+t1333+t1346)*t644+(t1367+t1386+t1406
     #+t1423)*s2l1+4992*t414
      t1429 = 106*t46
      t1431 = 6*t41
      t1449 = -64.D0/3.D0*t39*(-t1429+134*l1+t1431+6*t40+t543-42)*s1l1+1
     #47840*t417+1913152*t419+967104*t421-7688704.D0/3.D0*t423-2539584*t
     #425+1919552*t428+22400.D0/3.D0*t430-64.D0/3.D0*z+192*t433-64.D0/3.
     #D0*t1264*t341*t12
      t1468 = 12*t46
      t1490 = -7+14*l5+42*t46+56*t52-55*t55+32*t57+100*t62-132*t65-276*t
     #67+72*t69-180*t71-72*t74-64*t27+200*t83+120*t88-120*t90-192*t93+11
     #0*t95+96*t103
      t1510 = -48*t107-100*t114-480*t122+192*t124-192*t128-176*t146+304*
     #t148-12*t156+24*t158+42*t160+240*t164+72*t168-28*t171+56*t173-7*t1
     #77-150*t179+150*t181-42*t26-144*t184
      t1531 = 66*t186-232*t188-48*t190+32*l1+38*t39+48*t208+400*t211-112
     #*t214+112*t216+400*t219-112*t222-400*t224+112*t226-400*t228+480*t2
     #31-480*t233+192*t235+480*t238-28*t242
      t1550 = 28*t244-56*t246+18*t252+180*t254-200*t258-12*t41-96*t261+9
     #6*t263-7*t54-240*t276-55*t40+38*t281+332*t285-100*t287+76*t289-12*
     #t291+t740+42*t296-84*t306
      t1562 = 64.D0/3.D0*t39*(t543+t1431-t1429+12*t40+125*l1-41)*s1l1+32
     #5248*t516-762112*t518+9856*t615-59136*t520-143360.D0/3.D0*t618+325
     #248*t435+8960*t621+192*l1-128*t438-t766*(-64.D0/3.D0*t39*(30*t40-6
     #0*l1+t1468+21)*s1l1-64.D0/3.D0*(t1490+t1510+t1531+t1550)*t12-64.D0
     #/3.D0*t39*(78*t40-104*l1-t1468+41)*s2l1)
      t1574 = 4096.D0/3.D0*t175+t1351-t905-28416*t337-10816.D0/3.D0*t339
     #-226880.D0/3.D0*t343-7232*t345+3072*t348+75072*t350-49472.D0/3.D0*
     #t352-17920.D0/3.D0*t9-9472*t356+t468-t469-t470-t471-7424*t370
      t1592 = 12032*t372+57728.D0/3.D0*t374-8064*t376-20608.D0/3.D0*t378
     #-6016*t380+45632.D0/3.D0*t382-126400.D0/3.D0*t384-6976.D0/3.D0*t38
     #6+112768.D0/3.D0*t388-28160.D0/3.D0*t390+1472*t392-3584.D0/3.D0*t4
     #86-270080.D0/3.D0*t394-16960.D0/3.D0*t489+24512.D0/3.D0*t398+29017
     #6.D0/3.D0*t400+116672.D0/3.D0*t402
      t1608 = 7936.D0/3.D0+8384.D0/3.D0*t496+1600.D0/3.D0*t498-142016.D0
     #/3.D0*t404+97280.D0/3.D0*t406-233792.D0/3.D0*t408+218048.D0/3.D0*t
     #410+t1393-448.D0/3.D0*t414+t506+108032.D0/3.D0*t419+37568*t421-398
     #08*t423-36864*t425+66560.D0/3.D0*t428+1664*t430+6016.D0/3.D0*z
      t1620 = -14464.D0/3.D0*t433+1792*t516-14912.D0/3.D0*t518-t521+3584
     #*t435-23680.D0/3.D0*l1+16384.D0/3.D0*t438+t525-1792.D0/3.D0*t441-2
     #176*t443+21248.D0/3.D0*t528-t1419+t712+9728.D0/3.D0*t40-t1005-t100
     #6+2432*t450
      t1646 = -1792.D0/3.D0*t175-t457+1024.D0/3.D0*t335+310592*t337+9516
     #8.D0/3.D0*t339+410048*t343+116480.D0/3.D0*t345-7744*t348-1187648.D
     #0/3.D0*t350+90560.D0/3.D0*t352+12032.D0/3.D0*t9-92928*t356+t1361+t
     #1362-t1363-t1364-5376*t368-112640*t370+135936*t372
      t1665 = -86592*t374+154624.D0/3.D0*t376+51520.D0/3.D0*t378-9216*t5
     #75+31168*t380-408640.D0/3.D0*t382+867328.D0/3.D0*t384-9344*t386-22
     #9696*t388+72064*t390-5888.D0/3.D0*t392+16832*t486+714752*t394+2392
     #96.D0/3.D0*t489-51264*t398-1992256.D0/3.D0*t400-1119808.D0/3.D0*t4
     #02-59328*t496-t499
      t1682 = -5632.D0/3.D0+888640.D0/3.D0*t404-618304.D0/3.D0*t406+5244
     #16*t408-1799680.D0/3.D0*t410-3584.D0/3.D0*t412-8704.D0/3.D0*t414+t
     #1395-87104.D0/3.D0*t419-608320.D0/3.D0*t421-29056*t423+523456.D0/3
     #.D0*t425+227008.D0/3.D0*t428-t512-3712.D0/3.D0*z+6272.D0/3.D0*t433
     #+49920*t516-157504.D0/3.D0*t518+t616
      t1710 = -5888.D0/3.D0-2048.D0/3.D0*t175-t334-872320.D0/3.D0*t337-1
     #27744.D0/3.D0*t339-5176192.D0/3.D0*t343-792704.D0/3.D0*t345+696832
     #.D0/3.D0*t348+5765888.D0/3.D0*t350-2372992.D0/3.D0*t352+4224*t9-65
     #0496*t356+18944*t359+118272*t361-128768*t363-64256*t366-104192*t36
     #8-459008*t370+708224*t372
      t1729 = -37888.D0/3.D0*t374-9856.D0/3.D0*t376+10240.D0/3.D0*t378-2
     #5344*t575+2688*t380+37504*t382-143360.D0/3.D0*t384+12160.D0/3.D0*t
     #386+37120*t388-41984.D0/3.D0*t390+t393-11520*t486-3822976.D0/3.D0*
     #t394-277632*t489+156800*t398+3309824.D0/3.D0*t400+2496640.D0/3.D0*
     #t402+90240*t496-1661312.D0/3.D0*t404
      t1749 = 159744*t406-1056640.D0/3.D0*t408+437376*t410-8960.D0/3.D0*
     #t412+12672*t355-2304*t358+21376.D0/3.D0*t414-28544*t347+83456*t417
     #+1467392*t419+2729344.D0/3.D0*t421-4822144.D0/3.D0*t423-3746560.D0
     #/3.D0*t425+3127552.D0/3.D0*t428+t431-3968.D0/3.D0*z+3072*t433+1647
     #36*t516-371328*t518
      t1768 = 4608*t615-29952*t520-7680*t618+353408*t435+6144*l1-4096*t4
     #38+7680*t42-2048*t441-62208*t627+56064*t629+6272.D0/3.D0*t443+3456
     #0*t632+435456*t528-t445+100480.D0/3.D0*t342+640.D0/3.D0*t341-5888.
     #D0/3.D0*t40-21632*t130-768*t450
      t1774 = -t1407+1920*t618+34304*t435+17792.D0/3.D0*l1+2176.D0/3.D0*
     #t438+t525+64.D0/3.D0*t441-4608*t627+12480*t629-12352*t443-2880*t63
     #2-20096*t528+(t1710+t1729+t1749+t1768)*s2l1+t530+t637-2560.D0/3.D0
     #*t40+t534+t536-1280.D0/3.D0*t450
      t1779 = (t1574+t1592+t1608+t1620)*s2l1+192*t441-157248*t627+375872
     #.D0/3.D0*t629+12800*t443+115136*t632+990144*t528-64.D0/3.D0*t341-1
     #28*t40+(t1646+t1665+t1682+t1774)*s1l1-128*t450
      t1784 = (l1+z-t9)**2
      taux9 = (t750+t762+t1272+t1284+t1428+t1449+t1562+t1779)/t1784*t12
      t1 = z**2
      t2 = t1**2
      t3 = t2*z
      t4 = l1**2
      t5 = t4*l1
      t6 = t4**2
      t7 = t6*t5
      t8 = t3*t7
      t10 = t6*t4
      t12 = t6**2
      t13 = t12*t4
      t15 = 1280*z
      t16 = 12800.D0/3.D0*l1
      t17 = t12*t5
      t18 = t1*z
      t19 = t2*t18
      t21 = 13824*t17*t19
      t22 = t6*l1
      t23 = t19*t22
      t25 = t3*t12
      t27 = t3*t4
      t29 = t18*t10
      t31 = -1280-28082432.D0/3.D0*t8-144640.D0/3.D0*t10-76032*t13-t15-t
     #16-t21-74880*t23+8571264*t25+11008*t27+34099328.D0/3.D0*t29
      t32 = t10*z
      t34 = t19*t12
      t36 = t19*t7
      t38 = t18*t4
      t40 = t1*t5
      t42 = t2*t5
      t44 = t18*t5
      t46 = t2*t4
      t48 = t12*l1
      t49 = t48*t19
      t51 = t3*t5
      t53 = t3*t6
      t55 = 1159296*t32+392832*t34-374400*t36+115712*t38+543232.D0/3.D0*
     #t40+1331584.D0/3.D0*t42-590336*t44-161792.D0/3.D0*t46-251136*t49-3
     #49312.D0/3.D0*t51+2536064.D0/3.D0*t53
      t57 = t2*t6
      t59 = z*t6
      t61 = t18*t6
      t63 = t1*t6
      t68 = 483840*t17*t18
      t70 = 483840*t17*t2
      t72 = 290304*t17*t3
      t73 = z*t13
      t75 = t13*t18
      t77 = -6307328.D0/3.D0*t57+94208*t59+6868864.D0/3.D0*t61-2840320.D
     #0/3.D0*t63-615424.D0/3.D0*t12+406784.D0/3.D0*t7-t68+t70-t72+546048
     #*t73+2868480*t75
      t78 = t13*t2
      t80 = t13*t1
      t82 = t19*t10
      t84 = t48*z
      t86 = t2*t1
      t87 = t86*t10
      t89 = t10*t3
      t91 = t86*t22
      t93 = t18*l1
      t94 = 49664.D0/3.D0*t93
      t95 = z*t4
      t96 = 38144.D0/3.D0*t95
      t97 = t1*t4
      t98 = 40192.D0/3.D0*t97
      t99 = t1*l1
      t100 = 4864.D0/3.D0*t99
      t101 = -2937600*t78-1679616*t80+220032*t82-1343744*t84-5764864.D0/
     #3.D0*t87+20225920.D0/3.D0*t89+2290816.D0/3.D0*t91-t94+t96-t98+t100
      t104 = z*t5
      t105 = 38656.D0/3.D0*t104
      t106 = t3*l1
      t107 = 1024*t106
      t108 = t2*l1
      t109 = 6400*t108
      t110 = t86*t4
      t111 = 512*t110
      t112 = 1536*t18
      t113 = 28160.D0/3.D0*t4
      t114 = 1280*t1
      t115 = 5632.D0/3.D0*t5
      t116 = 512*t2
      t118 = t86*t5
      t120 = -t105-t107+t109-t111+t112+t113-t114-t115-t116+9472*t22+6323
     #2.D0/3.D0*t118
      t121 = 4096*t6
      t123 = t10*t2
      t125 = t18*t12
      t127 = t19*t6
      t129 = t86*t6
      t131 = t10*t1
      t133 = t1*t22
      t135 = t18*t22
      t137 = t2*t22
      t139 = t3*t22
      t141 = -t121+172288*t48-11967488*t123+38845568.D0/3.D0*t125+11520*
     #t127-172800*t129-16651264.D0/3.D0*t131+8533120.D0/3.D0*t133-184928
     #00.D0/3.D0*t135+6200832*t137-3108864*t139
      t143 = z*t22
      t145 = t18*t7
      t147 = z*t7
      t149 = z*t12
      t151 = t2*t7
      t153 = t1*t12
      t155 = t86*t7
      t157 = z*l1
      t158 = 12800.D0/3.D0*t157
      t159 = t48*t18
      t161 = t48*t3
      t163 = t48*t86
      t165 = -465664*t143-43607168.D0/3.D0*t145-1814784*t147+1926784*t14
     #9+15592704*t151-6928896*t153+2939776*t155-t158-7800576*t159-509619
     #2*t161+1738496*t163
      t166 = t48*t1
      t168 = t2*t12
      t170 = t13*t3
      t172 = t13*t19
      t174 = t13*t86
      t177 = 96768*t17*z
      t179 = 290304*t17*t1
      t180 = t1*t7
      t182 = t48*t2
      t184 = 13824*t17
      t186 = 96768*t17*t86
      t187 = t86*t12
      t189 = 4385024*t166-13855744*t168+1804032*t170+89856*t172-615168*t
     #174-t177+t179+22252928.D0/3.D0*t180+8195840*t182+t184+t186-8548864
     #.D0/3.D0*t187
      t196 = dlog(z)
      t200 = 256*z
      t201 = 2560.D0/3.D0*l1
      t216 = 256-1280*t19*t5-1338112.D0/3.D0*t8+t200+t201+256*t19*t4-404
     #48*t23+512*t86*l1+206080*t25+14336*t27-521984*t29-235264.D0/3.D0*t
     #32+123136.D0/3.D0*t34-245504.D0/3.D0*t36-34304.D0/3.D0*t38-47104.D
     #0/3.D0*t40+8192*t42
      t227 = 8192*t84
      t232 = 2304*t95
      t234 = 1024.D0/3.D0*t99
      t235 = 57344*t44-11008.D0/3.D0*t46-8192*t49-209408.D0/3.D0*t51+210
     #176*t53-42752.D0/3.D0*t57-7424*t59-557312.D0/3.D0*t61+275200.D0/3.
     #D0*t63+242944.D0/3.D0*t82+t227-392192*t87+1636096.D0/3.D0*t89+7889
     #92.D0/3.D0*t91+4864.D0/3.D0*t93-t232-256*t97+t234
      t240 = 3072*t110
      t241 = 2048*t4
      t242 = 256*t1
      t243 = 1024*t5
      t254 = 1792*t104-1792*t106+1024*t108-t240-t241+t242+t243+22528*t11
     #8+256*t3-14336*t123-203520*t125+9472*t127-103936*t129+1141760.D0/3
     #.D0*t131-247808*t133+1166336.D0/3.D0*t135+52736.D0/3.D0*t137
      t263 = 2560.D0/3.D0*t157
      t265 = 40960*t161
      t281 = 256+479232*t8+t200+t201-217088*t25-5120*t27-925696*t29-1392
     #64.D0/3.D0*t32-27392.D0/3.D0*t38-1536*t40-345088.D0/3.D0*t42+17868
     #8.D0/3.D0*t44+25344*t46+39680*t51-503296.D0/3.D0*t53
      t293 = 1092352.D0/3.D0*t57+1792.D0/3.D0*t59-746752.D0/3.D0*t61+432
     #64*t63+t227+111616*t87-1712128.D0/3.D0*t89-64000*t91+2304*t93-t232
     #-1280.D0/3.D0*t97+t234+4352.D0/3.D0*t104+t107-3072*t108
      t305 = t111-t241+t242+t243+t116-2048*t118+3217408.D0/3.D0*t123-401
     #408*t125+16896*t129+1076224.D0/3.D0*t131-555008.D0/3.D0*t133+61644
     #8*t135-2337280.D0/3.D0*t137+1199360.D0/3.D0*t139+35584.D0/3.D0*t14
     #3
      t312 = 81920*t159
      t313 = 8192*t163
      t314 = 40960*t166
      t317 = 81920*t182
      t319 = 826368*t145+62464*t147-36864*t149-901120*t151+192512*t153-9
     #9840*t155+t263+t312+t265-t313-t314+417792*t168-367104*t180-t317+45
     #056*t187
      t335 = 256+1103360.D0/3.D0*t8+t200+t201-192512*t25-256*t27-3184640
     #.D0/3.D0*t29-326656.D0/3.D0*t32-64768.D0/3.D0*t38-92672.D0/3.D0*t4
     #0-156160.D0/3.D0*t42+306176.D0/3.D0*t44+12032.D0/3.D0*t46+4864.D0/
     #3.D0*t51
      t349 = -44800*t53+240896*t57-41984.D0/3.D0*t59-1028864.D0/3.D0*t61
     #+457984.D0/3.D0*t63+t227+143360.D0/3.D0*t87-361472*t89-40960.D0/3.
     #D0*t91+7936.D0/3.D0*t93-2560*t95+1280*t97-512.D0/3.D0*t99+2304*t10
     #4
      t362 = -512*t108-256*t18-t241+t242+t243+2768896.D0/3.D0*t123-41779
     #2*t125+561152*t131-383488*t133+2277376.D0/3.D0*t135-1831424.D0/3.D
     #0*t137+188672*t139+179456.D0/3.D0*t143+2695168.D0/3.D0*t145
      t371 = 297472.D0/3.D0*t147-45056*t149-2477056.D0/3.D0*t151+217088*
     #t153-188416.D0/3.D0*t155+t263+t312+t265-t313-t314+401408*t168-1430
     #528.D0/3.D0*t180-t317+36864*t187
      t376 = -418304*t139+111616.D0/3.D0*t143+1304320.D0/3.D0*t145+24243
     #2.D0/3.D0*t147-122624.D0/3.D0*t149+7168*t151+163328*t153+1034240.D
     #0/3.D0*t155+t263+40960*t159-t265+32768*t163-32768*t166-5120.D0/3.D
     #0*t168+(t281+t293+t305+t319)*s1l1-1018880.D0/3.D0*t180+(t335+t349+
     #t362+t371)*s2l1-164352*t187
      t390 = (t31+t55+t77+t101+t120+t141+t165+t189)*s2l1+5105920.D0/3.D0
     #*t8+499840.D0/3.D0*t23-t196*(t216+t235+t254+t376)-793600*t25-26880
     #*t27+8585728.D0/3.D0*t29+1187840.D0/3.D0*t32-181760*t34+1081472.D0
     #/3.D0*t36+30720*t38+194560.D0/3.D0*t40+36480*t42-734720.D0/3.D0*t4
     #4
      t417 = 1280+22272*t8+1280*t10+t15+t16-6656*t25-7168*t27-40960*t29-
     #256*t32-54272*t38-123392.D0/3.D0*t40-96768*t42+310528.D0/3.D0*t44+
     #109312.D0/3.D0*t46+77312.D0/3.D0*t51-157952.D0/3.D0*t53
      t422 = 256.D0/3.D0*t12
      t427 = 428800.D0/3.D0*t57-26624.D0/3.D0*t59-385792.D0/3.D0*t61+156
     #160.D0/3.D0*t63+t422-512*t7+9216*t87-141568.D0/3.D0*t89-24064.D0/3
     #.D0*t91+t94-t96+t98-t100+t105+t107-t109
      t438 = t111-t112-t113+t114+t115+t116-7936.D0/3.D0*t22-2304*t118+t1
     #21+69888*t123+7168.D0/3.D0*t125+5120*t129+24064.D0/3.D0*t131-94976
     #.D0/3.D0*t133+98560*t135-128000*t137
      t443 = 512*t149
      t447 = 1024*t159
      t448 = 1024*t161
      t449 = 512*t163
      t450 = 512*t166
      t454 = 190720.D0/3.D0*t139+3584*t143+5376*t145+2048.D0/3.D0*t147-t
     #443-64256.D0/3.D0*t151-1024*t153-6912*t155+t158-t447+t448-t449+t45
     #0+2816*t168+512*t180+8704.D0/3.D0*t187
      t458 = -3840*t46+36352*t49+602240.D0/3.D0*t51-2127104.D0/3.D0*t53-
     #155904*t57+102400.D0/3.D0*t59+2756096.D0/3.D0*t61-1296640.D0/3.D0*
     #t63-1054208.D0/3.D0*t82-40960*t84+4906496.D0/3.D0*t87-6138880.D0/3
     #.D0*t89-1051264*t91+(t417+t427+t438+t454)*s2l1
      t460 = 12*t22
      t461 = 50*l1
      t462 = 110*t4
      t463 = 25*t5
      t464 = 36*t6
      t467 = 2*l6-1
      t468 = t467**2
      t471 = l5-1
      t473 = -1+l1
      t475 = (1-l1+t157)**2
      t476 = 1/t475
      t481 = dsqrt(z)
      t483 = dsqrt(l5)
      t485 = dsqrt(-t471)
      t488 = dsqrt(l1)
      t492 = l5*z*l1-l1+1-l5+l5*l1
      t495 = dsqrt(-t473)
      t500 = t492**2
      t503 = dsqrt(4*l1*t468*z*l5*t471*t473*t476-4*t467*t481*t483*t485*t
     #488*t492*t476*t495+t500*t476)
      t509 = 31232*t27
      t519 = 1280+50944.D0/3.D0*t8+1792*t10+t15+t16-4608*t25-t509-55808*
     #t29-4096*t32-197888.D0/3.D0*t38-44032*t40-507392.D0/3.D0*t42+13798
     #4*t44+75008*t46+76544*t51-332288.D0/3.D0*t53
      t528 = 16640*t93
      t529 = 42752.D0/3.D0*t95
      t530 = 35584.D0/3.D0*t97
      t531 = 3328.D0/3.D0*t99
      t533 = 6144*t106
      t534 = 19456*t108
      t535 = 225280*t57-37120.D0/3.D0*t59-535808.D0/3.D0*t61+193792.D0/3
     #.D0*t63+t422-1792.D0/3.D0*t7+9984*t87-153088.D0/3.D0*t89-50944.D0/
     #3.D0*t91+t528-t529+t530-t531+45568.D0/3.D0*t104+t533-t534
      t537 = 512*t18
      t538 = 27392.D0/3.D0*t4
      t539 = 1536*t1
      t540 = 2560.D0/3.D0*t5
      t541 = 3072*t2
      t543 = 12288*t118
      t552 = t240-t537-t538+t539+t540+t541-11776.D0/3.D0*t22-t543+17152.
     #D0/3.D0*t6+80128*t123+13312.D0/3.D0*t125+19712*t129+57088.D0/3.D0*
     #t131-48640*t133+138496*t135-178944*t137
      t560 = 14336.D0/3.D0*t157
      t564 = 290816.D0/3.D0*t139+8448*t143+14080.D0/3.D0*t145+5632.D0/3.
     #D0*t147-t443-47104.D0/3.D0*t151-1536*t153-5376*t155+t560-t447+t448
     #-t449+t450-256*t168-5632.D0/3.D0*t180+7168.D0/3.D0*t187
      t582 = -1280-7591808*t8-183808*t10-89856*t13-t15-t16-t21-11520*t23
     #+7036800*t25+t509+34363520.D0/3.D0*t29
      t594 = 5068928.D0/3.D0*t32+217728*t34-156672*t36+75776*t38+92672*t
     #40+625024*t42-1270912.D0/3.D0*t44-132608*t46-175104*t49-213504*t51
     #+3000448.D0/3.D0*t53
      t604 = -2280960*t57+371200.D0/3.D0*t59+5756032.D0/3.D0*t61-680704*
     #t63-1132288.D0/3.D0*t12+340736*t7-t68+t70-t72+615168*t73+2937600*t
     #75
      t612 = -2868480*t78-1804032*t80+63360*t82-1723904*t84-3766016.D0/3
     #.D0*t87+17101952.D0/3.D0*t89+513408*t91-t528+t529-t530+t531
      t617 = -76288.D0/3.D0*t104-t533+t534-t240+t537+t538-t539-t540-t541
     #+173056.D0/3.D0*t22+t543
      t629 = -12544*t6+248320*t48-33782528.D0/3.D0*t123+41217152.D0/3.D0
     #*t125-110592*t129-18624512.D0/3.D0*t131+8102272.D0/3.D0*t133-17331
     #968.D0/3.D0*t135+18196480.D0/3.D0*t137-2935168*t139-1826176.D0/3.D
     #0*t143
      t641 = -15250432*t145-8250752.D0/3.D0*t147+2768000*t149+43560704.D
     #0/3.D0*t151-8411904*t153+1924480*t155-t560-8180736*t159-4411904*t1
     #61+1358336*t163+5069312*t166
      t649 = -12982016*t168+1679616*t170+76032*t172-546048*t174-t177+t17
     #9+8963712*t180+7815680*t182+t184+t186-5970688.D0/3.D0*t187
      t664 = 256.D0/3.D0*(-15+t460-t461+t462-t463-t464)*t503+(t519+t535+
     #t552+t564)*s1l1-256.D0/3.D0*(-15+t460-t461+t462-t463-t464-t196*(3-
     #24*t4+10*l1+12*t5))*t503+(t582+t594+t604+t612+t617+t629+t641+t649)
     #*s1l1-171520.D0/3.D0*t118-556032*t123+1139200*t125-90880.D0/3.D0*t
     #127+1123840.D0/3.D0*t129-5818880.D0/3.D0*t131+1232512*t133-2056576
     #*t135+379392*t137+1517056*t139
      t680 = -563200.D0/3.D0*t143-7272832.D0/3.D0*t145-1218560.D0/3.D0*t
     #147+204800*t149+484608*t151-842240*t153-4394624.D0/3.D0*t155-22784
     #0*t159+158720*t161-140800*t163+168448*t166-230400*t168+5244800.D0/
     #3.D0*t180+46080*t182+704000*t187
      t684 = (l1+z-t157)**2
      taux10 = (t390+t458+t664+t680)/t684*t476

CCC Coefficients

      tcoff15 = taux9
      tcoff16 = taux9
      tcoff17 = taux10
      tcoff18 = -taux10

CCC Jet functions

      tFjet15 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff15)
      tFjet16 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff16)
      tFjet17 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff17)
      tFjet18 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff18)

CCC Sum

      RVnum5e1 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet18
     #*tcoff18
      RETURN
      END

      double precision function RVnum5e2(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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
      t2 = t1**2
      t3 = t2*t1
      t4 = l1**2
      t5 = t4*l1
      t6 = t3*t5
      t9 = t4**2
      t10 = t9**2
      t11 = t10*t4
      t12 = t11*t2
      t13 = 1280*t12
      t15 = t10*l1
      t16 = t2*z
      t17 = t15*t16
      t18 = 1408*t17
      t20 = t9*l1
      t21 = 1024.D0/3.D0*t20
      t22 = t15*z
      t24 = t11*t1
      t25 = 2560*t24
      t26 = t15*t1
      t28 = t3*t10
      t29 = 128*t28
      t30 = t11*t16
      t31 = 256*t30
      t32 = t1*z
      t33 = t32*t10
      t35 = t9*t5
      t36 = t1*t35
      t38 = t2*t10
      t40 = t16*t35
      t42 = t16*t10
      t43 = 4352.D0/3.D0*t42
      t44 = t16*t4
      t46 = -640-1408.D0/3.D0*t6-t13+64*t9-t18-896.D0/3.D0*t1-t21-6016*t
     #22-t25+12544*t26+t29+t31+60416.D0/3.D0*t33+4992*t36-28288.D0/3.D0*
     #t38+5056*t40+t43-960*t44
      t47 = t16*l1
      t48 = 128*t47
      t50 = 7936.D0/3.D0*l1
      t52 = t32*l1
      t54 = z*t4
      t56 = t1*t4
      t58 = t1*l1
      t60 = t32*t4
      t62 = z*t5
      t64 = t1*t5
      t66 = t2*t5
      t68 = t32*t5
      t70 = t2*t4
      t72 = 768*t5
      t73 = 64*t2
      t76 = t48+1152*t15+t50-1408.D0/3.D0*z-128*t52-2048*t54-3584*t56+60
     #16.D0/3.D0*t58-5312.D0/3.D0*t60+256*t62-1792.D0/3.D0*t64-14848*t66
     #+13120*t68+9664.D0/3.D0*t70+t72+t73-5504.D0/3.D0*t10+1024*t35
      t78 = 320.D0/3.D0*t32
      t80 = t3*t4
      t81 = 64*t80
      t82 = 256*t11
      t83 = t9*t4
      t85 = t16*t5
      t87 = t16*t9
      t89 = t2*t9
      t91 = z*t9
      t93 = t32*t9
      t95 = t1*t9
      t97 = t2*l1
      t98 = 384*t97
      t99 = t3*t9
      t101 = t83*t2
      t103 = t3*t83
      t105 = t83*z
      t107 = t32*t83
      t109 = t83*t16
      t111 = t78-2688*t4+t81-t82+192*t83+4352*t85-12096*t87+37184*t89-53
     #12.D0/3.D0*t91-118912.D0/3.D0*t93+46208.D0/3.D0*t95-t98+4096.D0/3.
     #D0*t99+123712.D0/3.D0*t101+5056.D0/3.D0*t103-14080.D0/3.D0*t105-13
     #6960.D0/3.D0*t107-15872*t109
      t112 = t3*t20
      t114 = t83*t1
      t116 = t1*t20
      t118 = t32*t20
      t120 = t2*t20
      t122 = t16*t20
      t124 = z*t20
      t126 = t32*t35
      t128 = z*t35
      t130 = z*t10
      t132 = t2*t35
      t134 = t1*t10
      t136 = t3*t35
      t138 = z*l1
      t140 = t15*t32
      t142 = z*t11
      t143 = 1280*t142
      t144 = t11*t32
      t145 = 2560*t144
      t146 = t15*t2
      t148 = -2048*t112+23104*t114-30720*t116+183872.D0/3.D0*t118-161024
     #.D0/3.D0*t120+57280.D0/3.D0*t122+6144*t124+9088.D0/3.D0*t126-13504
     #.D0/3.D0*t128+29440.D0/3.D0*t130-26624.D0/3.D0*t132-60800.D0/3.D0*
     #t134-2176.D0/3.D0*t136+6208.D0/3.D0*t138-13056*t140+t143+t145+6784
     #*t146
      t151 = s2l1**2
      t153 = t2*t32
      t154 = t153*t20
      t156 = t153*t9
      t159 = 2*l6-1
      t160 = t159**2
      t163 = l5-1
      t165 = -1+l1
      t167 = (1-l1+t138)**2
      t168 = 1/t167
      t173 = dsqrt(z)
      t175 = dsqrt(l5)
      t177 = dsqrt(-t163)
      t180 = dsqrt(l1)
      t181 = z*l5
      t182 = t181*l1
      t183 = l5*l1
      t184 = t182-l1+1-l5+t183
      t187 = dsqrt(-t165)
      t192 = t184**2
      t195 = dsqrt(4*l1*t160*z*l5*t163*t165*t168-4*t159*t173*t175*t177*t
     #180*t184*t168*t187+t192*t168)
      t196 = 12*t9
      t205 = t175*l5
      t206 = t173*z*t205
      t207 = t206*t177
      t208 = t180*t5
      t210 = t208*t187*l6
      t211 = t207*t210
      t212 = 144*t211
      t213 = t195*z
      t214 = t213*t9
      t216 = t180*t4
      t218 = t216*t187*l6
      t219 = t207*t218
      t220 = 216*t219
      t221 = t173*t175
      t222 = t221*t177
      t223 = t222*t218
      t224 = 360*t223
      t225 = 48*t182
      t226 = t222*t210
      t227 = 144*t226
      t228 = t173*t205
      t229 = t228*t177
      t230 = t229*t210
      t231 = 144*t230
      t232 = l5**2
      t233 = z*t232
      t234 = l6**2
      t235 = t5*t234
      t236 = t233*t235
      t237 = 720*t236
      t238 = t9*t234
      t239 = t233*t238
      t240 = 288*t239
      t241 = t9*l6
      t242 = t233*t241
      t243 = 288*t242
      t244 = t5*l6
      t245 = t233*t244
      t246 = 720*t245
      t247 = t181*t241
      t248 = 288*t247
      t249 = 16*l5
      t250 = t195*t9
      t252 = t181*t244
      t253 = 720*t252
      t255 = t177*t216*t187
      t256 = t221*t255
      t257 = 180*t256
      t258 = t4*l6
      t259 = t181*t258
      t260 = 560*t259
      t261 = l1*l6
      t262 = t181*t261
      t263 = 128*t262
      t264 = -8+t212-180*t214-t220+t224-t225-t227+t231+t237-t240+t243-t2
     #46-t248+t249+90*t250+t253-t257-t260+t263
      t265 = t4*t234
      t266 = t233*t265
      t267 = 560*t266
      t268 = t233*t261
      t269 = 128*t268
      t270 = t233*t258
      t271 = 560*t270
      t272 = l1*t234
      t273 = t233*t272
      t274 = 128*t273
      t275 = t181*t272
      t276 = 128*t275
      t277 = t181*t238
      t278 = 288*t277
      t279 = t233*l1
      t280 = 48*t279
      t281 = t213*t4
      t283 = t181*t265
      t284 = 560*t283
      t285 = t181*t235
      t286 = 720*t285
      t287 = t213*l1
      t289 = t181*t4
      t290 = 210*t289
      t291 = t195*t1
      t292 = t291*t4
      t294 = 86*t183
      t295 = t213*t5
      t297 = t233*t9
      t298 = 108*t297
      t299 = t181*t5
      t300 = 270*t299
      t301 = t181*t9
      t302 = 108*t301
      t303 = t232*t1
      t304 = t303*t4
      t305 = 8*t304
      t306 = -t267-t269+t271+t274-t276+t278+t280-350*t281+t284-t286+80*t
     #287+t290+40*t292-t294+450*t295-t298-t300+t302-t305
      t308 = t233*t5
      t309 = 270*t308
      t310 = t303*t5
      t311 = 27*t310
      t312 = t233*t4
      t313 = 210*t312
      t314 = t303*t9
      t315 = 18*t314
      t316 = t291*t5
      t318 = l5*t5
      t319 = 126*t318
      t320 = l5*t9
      t321 = 36*t320
      t322 = t232*t9
      t323 = 18*t322
      t324 = l5*t4
      t325 = 160*t324
      t326 = t195*t5
      t329 = t180*t187*l6
      t330 = t229*t329
      t331 = 64*t330
      t332 = t291*t9
      t335 = t177*t180*t187
      t336 = t228*t335
      t337 = 32*t336
      t338 = 18*t9
      t340 = t177*t208*t187
      t341 = t228*t340
      t342 = 72*t341
      t343 = t180*l1
      t345 = t343*t187*l6
      t346 = t222*t345
      t347 = 280*t346
      t348 = t195*t4
      t351 = t177*t343*t187
      t352 = t221*t351
      t353 = 140*t352
      t354 = t206*t351
      t355 = 32*t354
      t356 = t309+t311-t313-t315-135*t316-t319+t321-t323+t325-315*t326-t
     #331+90*t332+t337-t338-t342-t347+400*t348+t353-t355
      t357 = t229*t218
      t358 = 360*t357
      t359 = t222*t329
      t360 = 64*t359
      t361 = t232*t5
      t362 = 63*t361
      t363 = 8*t232
      t364 = t228*t255
      t365 = 180*t364
      t366 = t207*t345
      t367 = 64*t366
      t368 = t232*l1
      t369 = 43*t368
      t370 = t232*t4
      t371 = 80*t370
      t372 = t195*l1
      t374 = t221*t335
      t375 = 32*t374
      t376 = 43*l1
      t377 = t206*t255
      t378 = 108*t377
      t379 = 63*t5
      t380 = t229*t345
      t381 = 280*t380
      t383 = t221*t340
      t384 = 72*t383
      t385 = t206*t340
      t386 = 72*t385
      t387 = 80*t4
      t388 = t228*t351
      t389 = 140*t388
      t390 = -t358+t360+t362-t363+t365+t367+t369-t371-215*t372-t375+t376
     #+t378+t379+t381+40*t195+t384-t386-t387-t389
      t404 = dlog(z)
      t405 = t404**2
      t408 = 12*t5
      t430 = -7+96*t211-132*t214-144*t219+240*t223-42*t182-96*t226+96*t2
     #30+480*t236-192*t239+192*t242-480*t245-192*t247+14*l5+66*t250+480*
     #t252-120*t256-400*t259+112*t262
      t450 = -400*t266-112*t268+400*t270+112*t273-112*t275+192*t277+42*t
     #279-276*t281+400*t283-480*t285+76*t287+150*t289+38*t292-64*t183+33
     #2*t295-72*t297-180*t299+72*t301-7*t304
      t463 = 66*t332
      t470 = 180*t308+18*t310-150*t312-12*t314-100*t316-84*t318+24*t320-
     #12*t322+110*t324-232*t326-56*t330+t463+28*t336-t196-48*t341-200*t3
     #46+304*t348+100*t352-28*t354
      t490 = -240*t357+56*t359+42*t361-7*t232+120*t364+56*t366+32*t368-5
     #5*t370-176*t372-28*t374+32*l1+72*t377+42*t5+200*t380+38*t195+48*t3
     #83-48*t385-55*t4-100*t388
      t503 = 64.D0/3.D0-10880*t6+(t46+t76+t111+t148)*t151+42240*t154-217
     #60.D0/3.D0*t156-64.D0/3.D0*t195*(t196+170*t4-160*l1-66*t5+39)*s2l1
     #+64.D0/3.D0*(t264+t306+t356+t390)*t168+204800*t12+256*t9-64.D0/3.D
     #0*t195*(26*t4-74*l1+6*t5+t196+27)*s1l1+t405*(-64.D0/3.D0*t195*(30*
     #t4-60*l1+t408+21)*s1l1-64.D0/3.D0*(t430+t450+t470+t490)*t168-64.D0
     #/3.D0*t195*(78*t4-104*l1-t408+41)*s2l1)/2
      t506 = 2560*t12
      t508 = t15*t3
      t509 = 2048*t508
      t512 = 1280*t24
      t515 = 1280*t30
      t516 = t11*t3
      t517 = 256*t516
      t522 = -3968.D0/3.D0-5632.D0/3.D0*t6+t506+9984*t17-t509-640*t1+179
     #2*t22+t512-9216*t26+18368.D0/3.D0*t28-t515+t517-174464.D0/3.D0*t33
     #-149824.D0/3.D0*t36+179072.D0/3.D0*t38+153472.D0/3.D0*t40
      t525 = 768*t47
      t527 = 2944.D0/3.D0*z
      t530 = 4096*t56
      t538 = 384*t2
      t539 = -91072.D0/3.D0*t42-13696.D0/3.D0*t44+t525+14848.D0/3.D0*l1-
     #t527+3328.D0/3.D0*t52-3200*t54-t530+9472.D0/3.D0*t58-5760*t60+128*
     #t62-4736.D0/3.D0*t64-106304.D0/3.D0*t66+70784.D0/3.D0*t68+12160*t7
     #0+t72+t538
      t541 = 128.D0/3.D0*t32
      t543 = 384*t80
      t558 = -t541-4224*t4+t543+42496.D0/3.D0*t85-90368.D0/3.D0*t87+2177
     #92.D0/3.D0*t89-2944.D0/3.D0*t91-186688.D0/3.D0*t93+18432*t95-8192.
     #D0/3.D0*t97+4096*t99+126592*t101+26560.D0/3.D0*t103-9664*t105-3792
     #64.D0/3.D0*t107-171584.D0/3.D0*t109-19136.D0/3.D0*t112
      t573 = 256*t142
      t575 = 173504.D0/3.D0*t114-43200*t116+106880*t118-332032.D0/3.D0*t
     #120+142336.D0/3.D0*t122+17024.D0/3.D0*t124+104704*t126+27520.D0/3.
     #D0*t128-16064.D0/3.D0*t130-317056.D0/3.D0*t132+84160.D0/3.D0*t134-
     #9408*t136+11392.D0/3.D0*t138+18944*t140-t573-t145-19456*t146
      t580 = 64.D0/3.D0*t1
      t588 = 256*t44
      t589 = -384-t13+2816.D0/3.D0*t9-t18+t580+896.D0/3.D0*t20-8064*t22-
     #t25+15616*t26+t29+t31+91136.D0/3.D0*t33+34688*t36-35968.D0/3.D0*t3
     #8+3264*t40+t43+t588
      t607 = 1664*t15+2432.D0/3.D0*l1-640.D0/3.D0*z-3136.D0/3.D0*t52-192
     #*t54+11840.D0/3.D0*t56-2176.D0/3.D0*t58-704.D0/3.D0*t60+8192.D0/3.
     #D0*t62-12416*t64-4096.D0/3.D0*t66+34496.D0/3.D0*t68-1664*t70-3712.
     #D0/3.D0*t5-13184.D0/3.D0*t10+17536.D0/3.D0*t35+704.D0/3.D0*t32
      t611 = 512*t85
      t617 = 512*t97
      t618 = 64*t99
      t625 = 448*t4-t82-11072.D0/3.D0*t83-t611-5312.D0/3.D0*t87+46208.D0
     #/3.D0*t89-19648.D0/3.D0*t91-102272.D0/3.D0*t93+74944.D0/3.D0*t95+t
     #617+t618+74432.D0/3.D0*t101+1856.D0/3.D0*t103+35584.D0/3.D0*t105-2
     #1760*t107-25600.D0/3.D0*t109-1024.D0/3.D0*t112
      t641 = -9920.D0/3.D0*t114-73856.D0/3.D0*t116+137536.D0/3.D0*t118-3
     #0720*t120+20992.D0/3.D0*t122+8384.D0/3.D0*t124-18304*t126-72128.D0
     #/3.D0*t128+60160.D0/3.D0*t130-2944.D0/3.D0*t132-106880.D0/3.D0*t13
     #4-1408.D0/3.D0*t136+1088.D0/3.D0*t138-15104*t140+t143+t145+7296*t1
     #46
      t644 = s1l1**2
      t650 = 12*t20
      t674 = 15-80*t211-138*t214+360*t219-440*t223+90*t182-160*t226+160*
     #t230-880*t236-320*t239+320*t242+880*t245-320*t247-30*l5
      t689 = 258*t250-880*t252+220*t256+960*t259-240*t262+960*t266+240*t
     #268-960*t270-240*t273+240*t275+320*t277-90*t279+620*t281-960*t283
      t705 = 880*t285-152*t287-360*t289-76*t292+150*t183-600*t295-120*t2
     #97+330*t299+120*t301+15*t304-330*t308-45*t310+360*t312+10*t314
      t707 = t180*t20
      t709 = t177*t707*t187
      t712 = t213*t20
      t715 = t83*t234
      t720 = t83*l6
      t730 = t20*l6
      t733 = 234*t316-48*t228*t709+378*t712+70*t318+192*t181*t715+124*t3
     #20-62*t322+192*t233*t720-230*t324+231*t326+120*t330-t463-192*t233*
     #t715-192*t181*t720+672*t181*t730
      t736 = t20*t234
      t741 = t180*t9
      t743 = t741*t187*l6
      t749 = t707*t187*l6
      t766 = -672*t181*t736+672*t233*t736-336*t229*t743-672*t233*t730+96
     #*t229*t749-60*t336-62*t9+24*l5*t83-80*t341+480*t346+30*t303*t20+25
     #2*t233*t20-252*t181*t20-12*t303*t83
      t781 = t177*t741*t187
      t787 = t291*t20
      t789 = t195*t83
      t793 = -72*t233*t83-240*t207*t743+96*t207*t749-610*t348+72*t181*t8
     #3+54*t20-240*t352-48*t206*t709-168*t221*t781+60*t354+54*t232*t20-1
     #35*t787+54*t789+168*t228*t781
      t795 = t195*t20
      t806 = t291*t83
      t809 = t213*t83
      t814 = -243*t795+440*t357-120*t359-35*t361+48*t221*t709+15*t232+33
     #6*t222*t743-220*t364+54*t806-120*t366-108*t809-75*t368+115*t370+38
     #6*t372
      t834 = 60*t374-12*t232*t83-75*l1-108*l5*t20-180*t377-35*t5-480*t38
     #0-76*t195+80*t383+40*t385+115*t4+120*t206*t781+240*t388-96*t222*t7
     #49-12*t83
      t844 = 21120*t154
      t845 = 3840*t156
      t846 = 5760*t12
      t852 = 38016*t24
      t855 = 13824*t30
      t856 = t11*t153
      t857 = 2304*t856
      t858 = 10368*t516
      t862 = 1984.D0/3.D0-3584*t6+t844-t845-t846+192*t9+96768*t17-67392*
     #t508+320*t1+t21-85248*t22-t852+177984*t26+185472*t28-t855-t857+t85
     #8+271616*t33+308672*t36+74240*t38
      t869 = t153*t10
      t870 = 40704*t869
      t871 = t153*t35
      t872 = 59136*t871
      t884 = 457984*t40-285312*t42+4160.D0/3.D0*t44-t48+15552*t15-8512.D
     #0/3.D0*l1+1472.D0/3.D0*z-t870+t872-256.D0/3.D0*t52+7808.D0/3.D0*t5
     #4+12416.D0/3.D0*t56-6592.D0/3.D0*t58-11008.D0/3.D0*t60-896*t62-359
     #68.D0/3.D0*t64-8960.D0/3.D0*t66+17792*t68-2048*t70-1408*t5
      t889 = 3456*t11
      t891 = t15*t153
      t892 = 14976*t891
      t904 = -t73-26624*t10+61888.D0/3.D0*t35-t78+9728.D0/3.D0*t4-t81-t8
     #89-19136.D0/3.D0*t83+t892+4608*t85-194944.D0/3.D0*t87+54528*t89-11
     #840.D0/3.D0*t91-140096.D0/3.D0*t93+28800*t95+512.D0/3.D0*t97+10246
     #4.D0/3.D0*t99+821056.D0/3.D0*t101+245632*t103+152768.D0/3.D0*t105
      t919 = t153*t83
      t920 = 48384*t919
      t925 = 18432*t142
      t926 = 34560*t144
      t928 = -136320*t128+155776*t130-218304*t132-t920-334464*t134-83628
     #8.D0/3.D0*t136-6784.D0/3.D0*t138-155520*t140+t925+t926+2880*t146
      t935 = -8+t212-292*t214-t220+t224-t225-t227+t231+t237-t240+t243-t2
     #46-t248+t249+188*t250+t253-t257-t260+t263-t267
      t941 = t290+41*t292-t294+520*t295-t298-t300+t302-t305+t309+t311-t3
     #13
      t952 = -t315-144*t316+84*t712-t319+t321-t323+t325-406*t326-t331+11
     #6*t332+t337-t338-t342-t347+445*t348+t353-t355-30*t787+12*t789-54*t
     #795
      t958 = -226*t372-t375+t376+t378+t379+t381+41*t195+t384-t386-t387-t
     #389
      t964 = 640*t6
      t965 = 128*t6
      t980 = -1024-t965-t506+2432.D0/3.D0*t9-2816*t17-1024.D0/3.D0*t1-32
     #00.D0/3.D0*t20-14080*t22-5120*t24+28160*t26+256*t28+512*t30+171520
     #.D0/3.D0*t33+55808*t36-78848.D0/3.D0*t38+2560*t40+12544.D0/3.D0*t4
     #2
      t981 = 384*t44
      t986 = 1152*t54
      t988 = 2304*t58
      t990 = 384*t62
      t998 = t981+2816*t15+3456*l1-2048.D0/3.D0*z-12416.D0/3.D0*t52-t986
     #-17408.D0/3.D0*t56+t988+36224.D0/3.D0*t60-t990+20096.D0/3.D0*t64+1
     #7536.D0/3.D0*t66-14080*t68-15104.D0/3.D0*t70-128*t5-17920.D0/3.D0*
     #t10+17152.D0/3.D0*t35
      t1001 = 2304*t4
      t1004 = 384*t85
      t1010 = 1152*t97
      t1017 = 640*t32-t1001-512*t11-1664*t83-t1004-13312.D0/3.D0*t87+857
     #6*t89-8192.D0/3.D0*t91-15104.D0/3.D0*t93+9088.D0/3.D0*t95+t1010+24
     #32.D0/3.D0*t99+63104.D0/3.D0*t101+2176*t103+30592.D0/3.D0*t105-435
     #2*t107-39040.D0/3.D0*t109
      t1036 = -1920*t112-43136.D0/3.D0*t114-32512.D0/3.D0*t116+78080.D0/
     #3.D0*t118-86528.D0/3.D0*t120+39040.D0/3.D0*t122+10880.D0/3.D0*t124
     #-135680.D0/3.D0*t126-29696*t128+90368.D0/3.D0*t130+12032*t132-1784
     #32.D0/3.D0*t134-3584.D0/3.D0*t136+2304*t138-28160*t140+2560*t142+5
     #120*t144+14080*t146
      t1042 = 384*t1
      t1049 = -3200.D0/3.D0-t964+t506+(t980+t998+t1017+t1036)*s2l1+4096*
     #t17-768*t508-t1042+1024*t22+t512-4864*t26-832.D0/3.D0*t28-t515+t51
     #7-28544.D0/3.D0*t33+2944.D0/3.D0*t36+17792.D0/3.D0*t38
      t1052 = 1408*t44
      t1053 = 512*t47
      t1055 = 2176.D0/3.D0*z
      t1066 = -48064.D0/3.D0*t40-2752.D0/3.D0*t42-t1052+t1053+9344.D0/3.
     #D0*l1-t1055-3072*t52-12544.D0/3.D0*t56+4736.D0/3.D0*t58+38656.D0/3
     #.D0*t60-7936.D0/3.D0*t62+12544*t64+64640.D0/3.D0*t66-32576*t68-384
     #0*t70-t72+256*t2
      t1068 = 1024.D0/3.D0*t32
      t1069 = 1152*t4
      t1070 = 256*t80
      t1084 = t1068-t1069+t1070-7040.D0/3.D0*t85+17152*t87-157888.D0/3.D
     #0*t89+11008.D0/3.D0*t91+56768*t93-71168.D0/3.D0*t95-t98-2944.D0/3.
     #D0*t99-181120.D0/3.D0*t101-21952.D0/3.D0*t103+4672.D0/3.D0*t105+13
     #9648.D0/3.D0*t107+105280.D0/3.D0*t109+5248*t112
      t1097 = 5888.D0/3.D0*t138
      t1100 = -46528.D0/3.D0*t114+78080.D0/3.D0*t116-199936.D0/3.D0*t118
     #+74112*t120-34880*t122-10816.D0/3.D0*t124-11392*t126+1856.D0/3.D0*
     #t128-4544.D0/3.D0*t130+64768.D0/3.D0*t132+18880.D0/3.D0*t134+4224*
     #t136+t1097+9216*t140-t573-t145-8704*t146
      t1104 = 36*t5
      t1115 = 1280.D0/3.D0*t1
      t1116 = 10752*t33
      t1118 = 7168*t38
      t1120 = 1792*t42
      t1125 = 2176.D0/3.D0+t1115-t1116-32000*t36+t1118+7808*t40-t1120+44
     #8.D0/3.D0*t44-5504.D0/3.D0*l1+576*z-1280.D0/3.D0*t52
      t1133 = 5504.D0/3.D0*t70
      t1134 = 256*t5
      t1135 = 448.D0/3.D0*t32
      t1137 = 1024*t54+1536*t56-3968.D0/3.D0*t58-6400.D0/3.D0*t60-19456.
     #D0/3.D0*t64-3328.D0/3.D0*t66+13568*t68-t1133-t1134+t1135+1408*t4
      t1150 = -2944.D0/3.D0*t85-2176.D0/3.D0*t87+57856.D0/3.D0*t89-3008*
     #t91-130496.D0/3.D0*t93+76928.D0/3.D0*t95+896.D0/3.D0*t97+54912*t10
     #1-13632*t105-86016*t107-12480*t109
      t1158 = 8064*t128
      t1159 = 1792*t130
      t1161 = 7168*t134
      t1163 = 57216*t114-156544.D0/3.D0*t116+81152*t118-47232*t120+24064
     #.D0/3.D0*t122+10240*t124+47616*t126+t1158-t1159-31488*t132+t1161-4
     #352.D0/3.D0*t138
      t1167 = t153*t5
      t1169 = t153*t4
      t1174 = 10240*t17
      t1175 = 1088.D0/3.D0*t1
      t1176 = 2048*t22
      t1177 = 10240*t26
      t1185 = 256*t47
      t1186 = 960+512*t6-t1174+t509+t1175-t1176+t1177-11264*t28+101888*t
     #33+95872*t36-105472*t38-120448*t40+54528*t42+1280*t44-t1185
      t1198 = 128*t2
      t1200 = -t50+1984.D0/3.D0*z+832.D0/3.D0*t52+960*t54+3904.D0/3.D0*t
     #56-3712.D0/3.D0*t58-1472*t60+1600.D0/3.D0*t62-3328.D0/3.D0*t64+267
     #52*t66-8832*t68-17728.D0/3.D0*t70+t1134-t1198+1536*t4
      t1202 = 128*t80
      t1220 = t986+5888.D0/3.D0*t56-4352.D0/3.D0*t58-5888.D0/3.D0*t60+t9
     #90-2560.D0/3.D0*t64+5632.D0/3.D0*t66+512.D0/3.D0*t68-512.D0/3.D0*t
     #70+1792*t4-896.D0/3.D0*t85
      t1243 = 13824*t116-23040*t118+47104.D0/3.D0*t120-3712*t122-8320.D0
     #/3.D0*t124-13824*t126-2304*t128+512*t130+9216*t132-2048*t134-t1097
      t1260 = -t1202+(1024+t1115+3072*t33+9216*t36-2048*t38-2304*t40+512
     #*t42-8192.D0/3.D0*l1+t1055+1664.D0/3.D0*t52+t1220+1664*t87-23296.D
     #0/3.D0*t89+1280.D0/3.D0*t91+29824.D0/3.D0*t93-14336.D0/3.D0*t95-50
     #432.D0/3.D0*t101+11648.D0/3.D0*t105+25088*t107+12416.D0/3.D0*t109-
     #48896.D0/3.D0*t114+t1243)*s2l1-9920*t85+125440.D0/3.D0*t87-264896.
     #D0/3.D0*t89-448.D0/3.D0*t91+58432*t93-31168.D0/3.D0*t95+768*t97-42
     #24*t99-813568.D0/3.D0*t101-27904*t103+40576.D0/3.D0*t105+238336*t1
     #07+429184.D0/3.D0*t109
      t1275 = 20480*t140
      t1276 = 20480*t146
      t1277 = 16000*t112-287488.D0/3.D0*t114+149888.D0/3.D0*t116-156544*
     #t118+583936.D0/3.D0*t120-299456.D0/3.D0*t122-12736.D0/3.D0*t124-21
     #1968*t126-16768*t128+9472*t130+228352*t132-49152*t134+24960*t136-5
     #632.D0/3.D0*t138-t1275+t1276
      t1289 = 896.D0/3.D0-t1116-32512*t36+t1118+8320*t40-t1120-128*t44-8
     #96*l1+448.D0/3.D0*z+640*t52+448.D0/3.D0*t54
      t1300 = -640.D0/3.D0*t56-256.D0/3.D0*t58-15616.D0/3.D0*t60+256.D0/
     #3.D0*t62-5888*t64-25600.D0/3.D0*t66+20480*t68+1024*t70+t1134-128*t
     #32+384*t4
      t1312 = t1004-11072.D0/3.D0*t87+30976*t89-7936.D0/3.D0*t91-163648.
     #D0/3.D0*t93+82816.D0/3.D0*t95-256*t97+183680.D0/3.D0*t101-41024.D0
     #/3.D0*t105-92416*t107-43712.D0/3.D0*t109
      t1322 = 178304.D0/3.D0*t114-166528.D0/3.D0*t116+92416*t118-175744.
     #D0/3.D0*t120+34432.D0/3.D0*t122+30592.D0/3.D0*t124+49152*t126+t115
     #8-t1159-33024*t132+t1161-512*t138
      t1328 = t3*l1
      t1333 = 1984.D0/3.D0-5632*t6+10112*t154-2368*t156+40960*t12+(t1125
     #+t1137+t1150+t1163)*t151+320*t1167-64*t1169-64*t16+168960*t17+(t11
     #86+t1200+t1260+t1277)*s1l1-67584*t508+t1175+(t1289+t1300+t1312+t13
     #22)*t644+11264*t22+12288*t24-128*t1328-67584*t26+157824*t28-30720*
     #t30
      t1355 = t1174-t509+2368.D0/3.D0*t1+t1176-t1177+9216*t28-102912*t33
     #-342272.D0/3.D0*t36+99328*t38+270848.D0/3.D0*t40-47872*t42+1088.D0
     #/3.D0*t44-3584*l1+1088*z
      t1367 = 4160.D0/3.D0-1088*t52+1920*t54+2240*t56-t988-3136.D0/3.D0*
     #t60+64*t62-6400*t64-1280.D0/3.D0*t66+14848*t68-11456.D0/3.D0*t70-t
     #1134+1088.D0/3.D0*t32+2560*t4
      t1383 = -5824.D0/3.D0*t85-5824*t87+119360.D0/3.D0*t89-3072*t91-196
     #288.D0/3.D0*t93+95168.D0/3.D0*t95+2176.D0/3.D0*t97+645632.D0/3.D0*
     #t101+35840.D0/3.D0*t103-75776.D0/3.D0*t105-737792.D0/3.D0*t107-257
     #536.D0/3.D0*t109-10240.D0/3.D0*t112+389632.D0/3.D0*t114
      t1396 = -254464.D0/3.D0*t116+493696.D0/3.D0*t118-391168.D0/3.D0*t1
     #20+121792.D0/3.D0*t122+40384.D0/3.D0*t124+648448.D0/3.D0*t126+7091
     #2.D0/3.D0*t128-11008*t130-600832.D0/3.D0*t132+53248*t134-47104.D0/
     #3.D0*t136-2816*t138+t1275-t1276
      t1408 = -2048*t856+12288*t516-396096*t33-596480.D0/3.D0*t36+159155
     #2.D0/3.D0*t38+1564736.D0/3.D0*t40-397376*t42-10304.D0/3.D0*t44+448
     #*t47-5248.D0/3.D0*l1+512*z+(t1355+t1367+t1383+t1396)*s2l1-77888.D0
     #/3.D0*t869+97216.D0/3.D0*t871-1216.D0/3.D0*t52+832*t54+1088*t56-33
     #28.D0/3.D0*t58-3904.D0/3.D0*t60+192*t62
      t1428 = -10240.D0/3.D0*t64-14080*t66+34432.D0/3.D0*t68-64*t70+t113
     #5+t1069+768*t80+11264*t891+16896*t85-216320.D0/3.D0*t87+268352.D0/
     #3.D0*t89-4672.D0/3.D0*t91-184064.D0/3.D0*t93+19648*t95+128.D0/3.D0
     #*t97+25984*t99+1692544.D0/3.D0*t101+159488*t103-64448.D0/3.D0*t105
     #-408256*t107
      t1449 = -420416*t109-238208.D0/3.D0*t112+150144*t114-69376*t116+60
     #5312.D0/3.D0*t118-287360*t120+216064*t122+8192*t124+1546816.D0/3.D
     #0*t126+93632.D0/3.D0*t128-77120.D0/3.D0*t130-698368*t132-70976.D0/
     #3.D0*t919+156800*t134-610816.D0/3.D0*t136-4096.D0/3.D0*t138+168960
     #*t140-2048*t142-30720*t144-225280*t146
      t1455 = 841280*t17+(t522+t539+t558+t575)*s2l1+(t589+t607+t625+t641
     #)*t644-t404*(-64.D0/3.D0*t195*(-106*t5+134*l1+6*t9+6*t4+t650-42)*s
     #1l1-64.D0/3.D0*t195*(78*t9-104*t5-118*t4-t650+244*l1-80)*s2l1-64.D
     #0/3.D0*(t674+t689+t705+t733+t766+t793+t814+t834)*t168)-335104*t508
     #+(t862+t884+t904+119936.D0/3.D0*t107-1290304.D0/3.D0*t109-125824*t
     #112-376064.D0/3.D0*t114-8448*t116+165248.D0/3.D0*t118-529024.D0/3.
     #D0*t120+233600*t122+1408.D0/3.D0*t124-639104.D0/3.D0*t126+t928)*t1
     #51+t580-64.D0/3.D0*(t935-t269+t271+t274-t276+t278+t280-370*t281+t2
     #84-t286+82*t287+t941+t952-t358+t360+t362-t363+t365+12*t806+t367-24
     #*t809+t369-t371+t958)*t168+(t1049+t1066+t1084+t1100)*s1l1-64.D0/3.
     #D0*t195*(t1104-144*t4+151*l1-38)*s2l1+t405*(t1333+t1408+t1428+t144
     #9)/2+57024*t22
      t1468 = 61952*t24-340736*t26+2327552.D0/3.D0*t28-152960*t30-10112*
     #t856+60928*t516-5943104.D0/3.D0*t33-1006848*t36+7899904.D0/3.D0*t3
     #8+2515392*t40-5883328.D0/3.D0*t42
      t1472 = 92160*t12
      t1473 = 512*t9
      t1476 = 2048.D0/3.D0*t1
      t1478 = 27648*t24
      t1481 = 69120*t30
      t1482 = 4608*t856
      t1483 = 27648*t516
      t1486 = 4096.D0/3.D0+512.D0/3.D0*t6+24960*t154-t845+t1472+t1473+37
     #7600*t17-158464*t508+t1476+20224*t22+t1478-128768*t26+378496*t28-t
     #1481-t1482+t1483-701056*t33-791936.D0/3.D0*t36
      t1505 = 1041152*t38+3165056.D0/3.D0*t40-862720*t42+21376.D0/3.D0*t
     #44-t525-16000.D0/3.D0*l1+1024*z-68736*t869+90240*t871-4480.D0/3.D0
     #*t52+12928.D0/3.D0*t54+15616.D0/3.D0*t56-10624.D0/3.D0*t58+896*t60
     #-1408*t62-37504.D0/3.D0*t64+40832*t66+4096.D0/3.D0*t68-42880.D0/3.
     #D0*t70
      t1508 = 256.D0/3.D0*t32
      t1524 = -2048*t5-t538+t1508+16000.D0/3.D0*t4-t543+27648*t891-68096
     #.D0/3.D0*t85-32384.D0/3.D0*t87-56704.D0/3.D0*t89-16000.D0/3.D0*t91
     #-76928.D0/3.D0*t93+110720.D0/3.D0*t95+2304*t97+72832.D0/3.D0*t99+7
     #79648*t101+1047424.D0/3.D0*t103-30592*t105-1450240.D0/3.D0*t107-72
     #5248*t109
      t1541 = 4608*t142
      t1542 = 69120*t144
      t1544 = -138112*t112+176128*t114-261760.D0/3.D0*t116+164608*t118-2
     #34240*t120+754816.D0/3.D0*t122+18432*t124+2316928.D0/3.D0*t126+115
     #840.D0/3.D0*t128-36608*t130-3628288.D0/3.D0*t132-65664*t919+249472
     #*t134-1448320.D0/3.D0*t136-12544.D0/3.D0*t138+340480*t140-t1541-t1
     #542-478720*t146
      t1551 = t32*t234
      t1554 = t173*t1
      t1555 = t1554*t175
      t1558 = t32*l6
      t1569 = -4*t182+16*t368*t1551-4*t1555*t335+16*t183*t1558-16*t368*t
     #1558-16*t183*t1551+2*t289-4*t304-4*t308+6*t310
      t1570 = 6*t312
      t1571 = 2*t233
      t1572 = t1554*t205
      t1573 = t1572*t177
      t1575 = 16*t1573*t345
      t1576 = l5*t1
      t1577 = 2*t1576
      t1579 = 4*t206*t335
      t1581 = 8*t1573*t329
      t1582 = 2*t181
      t1583 = 4*t354
      t1586 = 8*t1555*t177*t329
      t1587 = t1570-t303-t1571-t1575+t1577+t1579-t1581-t1+t1582-t1583+t1
     #586
      t1590 = 8*t207*t329
      t1592 = t173*t32*t205
      t1594 = 4*t1592*t351
      t1595 = 8*t366
      t1596 = 3*t368
      t1597 = 3*t370
      t1599 = 8*t1572*t351
      t1601 = 2*t1576*l1
      t1611 = t232*t32
      t1614 = t232*t2
      t1625 = 2*t324*t32-6*t183*t32+6*t368*t32+4*t1572*t335+2*t1611*t4-t
     #1614*t4-t303*l1-4*t1576*t4-4*t1611*t5+t1614*t5+8*t1592*t177*t345
      t1627 = t1569+t1587+t361-t1590-t232-t1594+t1595+t1596-t1597+t1599+
     #t58+t1601+t1625
      t1629 = t1570-t303-t1571-t1575+t1577+t404*t1627+t1579-t1581-t1+t15
     #82-t1583
      t1631 = t1586+t361-t1590-t232-t1594+t1595+t1596-t1597+t1599+t58+t1
     #601
      t1648 = -8960.D0/3.D0*t44+(t1486+t1505+t1524+t1544)*s2l1-192*l1+64
     #.D0/3.D0*z+64.D0/3.D0*(t1569+t1629+t1631+t1625)*t2*t168+64.D0/3.D0
     #*t195*(t1104-65*l1+26)*s1l1-382144.D0/3.D0*t869+156096*t871-192*t5
     #2+1664.D0/3.D0*t54+1664.D0/3.D0*t56-192*t58
      t1668 = 1216.D0/3.D0-t846-20608.D0/3.D0*t9+55296*t17-46656*t508+33
     #728*t20-117504*t22-t852+240192*t26+81792*t28-t855-t857+t858+617216
     #*t33+965504*t36-213760*t38+35968*t40-77952*t42
      t1672 = 17664*t869
      t1673 = 13440*t871
      t1685 = 64.D0/3.D0*t2
      t1687 = -t588+22464*t15-3008.D0/3.D0*l1+704.D0/3.D0*z-t1672+t1673+
     #3200.D0/3.D0*t52+1600*t54-26176.D0/3.D0*t56+1600.D0/3.D0*t58-1792.
     #D0/3.D0*t60-54272.D0/3.D0*t62+51136*t64+11840.D0/3.D0*t66-92032.D0
     #/3.D0*t68+5056.D0/3.D0*t70+1792.D0/3.D0*t5+t1685-61184*t10
      t1690 = 640.D0/3.D0*t32
      t1693 = 10368*t891
      t1705 = 268288.D0/3.D0*t35-t1690+320.D0/3.D0*t4-t889-222784.D0/3.D
     #0*t83+t1693+t611+14272.D0/3.D0*t87-165184.D0/3.D0*t89+255872.D0/3.
     #D0*t91+514240.D0/3.D0*t93-603776.D0/3.D0*t95-1472.D0/3.D0*t97-t618
     #-1072960.D0/3.D0*t101+80384.D0/3.D0*t103+1266880.D0/3.D0*t105+2579
     #072.D0/3.D0*t107+13632*t109
      t1717 = 3840*t919
      t1723 = -3712*t112-886656*t114+1568384.D0/3.D0*t116-492032*t118+58
     #5152.D0/3.D0*t120-54272.D0/3.D0*t122-713216.D0/3.D0*t124-2809088.D
     #0/3.D0*t126-471040*t128+317056*t130+371584*t132-t1717-645504*t134-
     #205568.D0/3.D0*t136-1664.D0/3.D0*t138-224640*t140+t925+t926+60480*
     #t146
      t1739 = 44544*t12
      t1744 = 256*t20
      t1745 = 7680*t154
      t1761 = -5888.D0/3.D0-t1745-128768*t12-1024.D0/3.D0*t9-459008*t17+
     #164736*t508-t1476+7680*t20-104192*t22-64256*t24+353408*t26-371328*
     #t28+83456*t30+t1482-29952*t516+1467392*t33+2729344.D0/3.D0*t36-482
     #2144.D0/3.D0*t38-3746560.D0/3.D0*t40
      t1766 = 56064*t869
      t1767 = 62208*t871
      t1779 = 3127552.D0/3.D0*t42+t588+12672*t15+6144*l1-3968.D0/3.D0*z+
     #t1766-t1767-2048*t52-768*t54-t530+3072*t58+21376.D0/3.D0*t60-8960.
     #D0/3.D0*t62+6272.D0/3.D0*t64-37888.D0/3.D0*t66-9856.D0/3.D0*t68+10
     #240.D0/3.D0*t70-2944*t5-28544*t10
      t1785 = 25344*t891
      t1798 = 100480.D0/3.D0*t35+t1690-5888.D0/3.D0*t4-2304*t11-21632*t8
     #3-t1785+2688*t85+37504*t87-143360.D0/3.D0*t89+12160.D0/3.D0*t91+37
     #120*t93-41984.D0/3.D0*t95+t617-11520*t99-3822976.D0/3.D0*t101-2776
     #32*t103+156800*t105+3309824.D0/3.D0*t107+2496640.D0/3.D0*t109
      t1810 = 34560*t919
      t1818 = 90240*t112-1661312.D0/3.D0*t114+159744*t116-1056640.D0/3.D
     #0*t118+437376*t120-872320.D0/3.D0*t122-127744.D0/3.D0*t124-5176192
     #.D0/3.D0*t126-792704.D0/3.D0*t128+696832.D0/3.D0*t130+5765888.D0/3
     #.D0*t132+t1810-2372992.D0/3.D0*t134+435456*t136+4224*t138-650496*t
     #140+18944*t142+118272*t144+708224*t146
      t1823 = 13568*t24
      t1826 = 33536*t30
      t1827 = 13568*t516
      t1830 = -5632.D0/3.D0-1408*t6+1920*t154-t1739+768*t9-112640*t17+49
     #920*t508-1792.D0/3.D0*t1+t1744+(t1761+t1779+t1798+t1818)*s2l1-5376
     #*t22-t1823+34304*t26-157504.D0/3.D0*t28+t1826+t857-t1827-87104.D0/
     #3.D0*t33-608320.D0/3.D0*t36
      t1849 = -29056*t38+523456.D0/3.D0*t40+227008.D0/3.D0*t42-3456*t44+
     #t525+17792.D0/3.D0*l1-3712.D0/3.D0*z+12480*t869-4608*t871+64.D0/3.
     #D0*t52-1280.D0/3.D0*t54+2176.D0/3.D0*t56+6272.D0/3.D0*t58-8704.D0/
     #3.D0*t60-3584.D0/3.D0*t62-12352*t64-86592*t66+154624.D0/3.D0*t68+5
     #1520.D0/3.D0*t70
      t1868 = -11584.D0/3.D0*t5+1024.D0/3.D0*t2+t541-2560.D0/3.D0*t4+t54
     #3-9216*t891+31168*t85-408640.D0/3.D0*t87+867328.D0/3.D0*t89-9344*t
     #91-229696*t93+72064*t95-5888.D0/3.D0*t97+16832*t99+714752*t101+239
     #296.D0/3.D0*t103-51264*t105-1992256.D0/3.D0*t107-1119808.D0/3.D0*t
     #109
      t1885 = 2304*t142
      t1886 = 33536*t144
      t1888 = -59328*t112+888640.D0/3.D0*t114-618304.D0/3.D0*t116+524416
     #*t118-1799680.D0/3.D0*t120+310592*t122+95168.D0/3.D0*t124+410048*t
     #126+116480.D0/3.D0*t128-7744*t130-1187648.D0/3.D0*t132-2880*t919+9
     #0560.D0/3.D0*t134-20096*t136+12032.D0/3.D0*t138-92928*t140+t1885+t
     #1886+135936*t146
      t1904 = -2688+36160.D0/3.D0*t6-29760*t154+5760*t156-t1739+2304*t9-
     #256256*t17+99328*t508-1408*t1-t1744-19968*t22-t1823+114944*t26-898
     #624.D0/3.D0*t28+t1826+t857-t1827+2658496.D0/3.D0*t33+1877120.D0/3.
     #D0*t36
      t1923 = -1119616*t38-3959168.D0/3.D0*t40+2381248.D0/3.D0*t42+4736*
     #t44+t1185+24832.D0/3.D0*l1-2048*z+47040*t869-73152*t871+2944.D0/3.
     #D0*t52-2688*t54-17152.D0/3.D0*t56+15616.D0/3.D0*t58+25280.D0/3.D0*
     #t60-11456.D0/3.D0*t62+20096*t64+231680.D0/3.D0*t66-70656*t68-4352.
     #D0/3.D0*t70
      t1942 = -11456.D0/3.D0*t5+t1198-384*t32-10496.D0/3.D0*t4+t1202-161
     #28*t891-52864*t85+273536*t87-1250944.D0/3.D0*t89+11712*t91+325824*
     #t93-111360*t95-3136.D0/3.D0*t97-245696.D0/3.D0*t99-5653120.D0/3.D0
     #*t101-1373824.D0/3.D0*t103+287104.D0/3.D0*t105+4456768.D0/3.D0*t10
     #7+1296128*t109
      t1960 = 769792.D0/3.D0*t112-599808*t114+339072*t116-2682304.D0/3.D
     #0*t118+3445760.D0/3.D0*t120-2319296.D0/3.D0*t122-45248*t124-449868
     #8.D0/3.D0*t126-321856.D0/3.D0*t128+65216*t130+5667520.D0/3.D0*t132
     #+63936*t919-1119040.D0/3.D0*t134+1454528.D0/3.D0*t136+19072.D0/3.D
     #0*t138-275968*t140+t1885+t1886+354048*t146
      t1969 = 64384*t12
      t1976 = 32128*t24
      t1979 = 41728*t30
      t1980 = 14976*t516
      t1985 = -576-t1969-14336.D0/3.D0*t9-201856*t17+69696*t508-64*t1+14
     #912*t20-56960*t22-t1976+180928*t26-129792*t28+t1979+t857-t1980+732
     #544*t33+1927616.D0/3.D0*t36-2125568.D0/3.D0*t38-1339328.D0/3.D0*t4
     #0
      t2001 = 1223168.D0/3.D0*t42+t981+7488*t15+5696.D0/3.D0*l1-256*z+t1
     #672-t1673-5888.D0/3.D0*t52+1792*t54-17920.D0/3.D0*t56+1088*t58+421
     #12.D0/3.D0*t60-7296*t62+38784*t64+24768*t66-195008.D0/3.D0*t68-928
     #0.D0/3.D0*t70-2048.D0/3.D0*t5-t1685
      t2006 = 1152*t11
      t2020 = -20608*t10+93056.D0/3.D0*t35+t1068+640.D0/3.D0*t4-t2006-82
     #688.D0/3.D0*t83-t1693-1152*t85+20096*t87-364160.D0/3.D0*t89+35520*
     #t91+218560*t93-441472.D0/3.D0*t95+2240.D0/3.D0*t97-2062336.D0/3.D0
     #*t101-57600*t103+604736.D0/3.D0*t105+868480*t107+288704*t109
      t2036 = 9472*t142
      t2037 = 59136*t144
      t2039 = 11520*t112-1770752.D0/3.D0*t114+1081088.D0/3.D0*t116-52108
     #8*t118+1071424.D0/3.D0*t120-108864*t122-112320*t124-2969216.D0/3.D
     #0*t126-221120*t128+447104.D0/3.D0*t130+876032*t132+t1717-1344128.D
     #0/3.D0*t134+121152*t136+2048.D0/3.D0*t138-314496*t140+t2036+t2037+
     #325568*t146
      t2045 = -3904.D0/3.D0+71744.D0/3.D0*t6-145024.D0/3.D0*t154+32384.D
     #0/3.D0*t156-187904*t12-1152*t1167+256*t1169+256*t16+896*t9-779328*
     #t17+(t1830+t1849+t1868+t1888)*s1l1+(t1904+t1923+t1942+t1960)*s2l1+
     #313984*t508-1984.D0/3.D0*t1-51392*t22-56064*t24+512*t1328+(t1985+t
     #2001+t2020+t2039)*t644+308352*t26-2207360.D0/3.D0*t28
      t2069 = -4160.D0/3.D0+7680*t6-t844+t845-t1969+1280*t9-257152*t17+9
     #5040*t508-2240.D0/3.D0*t1+1280.D0/3.D0*t20-47232*t22-t1976+172480*
     #t26-256512*t28+t1979+t857-t1980+840064*t33+1650176.D0/3.D0*t36
      t2088 = -2981888.D0/3.D0*t38-3074240.D0/3.D0*t40+2052608.D0/3.D0*t
     #42+9472.D0/3.D0*t44+5184*t15+12736.D0/3.D0*l1-3200.D0/3.D0*z+t870-
     #t872-896.D0/3.D0*t52-5120.D0/3.D0*t54-14848.D0/3.D0*t56+9280.D0/3.
     #D0*t58+28928.D0/3.D0*t60-7040.D0/3.D0*t62+60032.D0/3.D0*t64+65216*
     #t66-62720*t68-17600.D0/3.D0*t70-4736.D0/3.D0*t5
      t2105 = -9088*t10+22976.D0/3.D0*t35-t541-t1001-t2006-9472.D0/3.D0*
     #t83-t892-107840.D0/3.D0*t85+188608*t87-908800.D0/3.D0*t89+11136*t9
     #1+758336.D0/3.D0*t93-93568*t95+128*t97-55680*t99-4098560.D0/3.D0*t
     #101-339456*t103+89792*t105+3301760.D0/3.D0*t107
      t2123 = 943040*t109+182208*t112-1419776.D0/3.D0*t114+794048.D0/3.D
     #0*t116-653184*t118+2434304.D0/3.D0*t120-1628672.D0/3.D0*t122-12428
     #8.D0/3.D0*t124-3627392.D0/3.D0*t126-360896.D0/3.D0*t128+301184.D0/
     #3.D0*t130+4421696.D0/3.D0*t132+t920-1217408.D0/3.D0*t134+381696*t1
     #36+9856.D0/3.D0*t138-336000*t140+t2036+t2037+382656*t146
      t2135 = 4032*l1+(t2069+t2088+t2105+t2123)*t151-t527+366016.D0/3.D0
     #*t869-459584.D0/3.D0*t871-256*t52-2176.D0/3.D0*t54-3072*t56+2624*t
     #58+25600.D0/3.D0*t60-2496*t62
      t2158 = -78912*t85+1054336.D0/3.D0*t87-1348096.D0/3.D0*t89+7680*t9
     #1+886144.D0/3.D0*t93-267584.D0/3.D0*t95+2944.D0/3.D0*t97-121280*t9
     #9-7934464.D0/3.D0*t101-2276672.D0/3.D0*t103+290240.D0/3.D0*t105
      t2182 = -424768.D0/3.D0*t128+351296.D0/3.D0*t130+3234112*t132+1128
     #32*t919-715520*t134+2871296.D0/3.D0*t136+3072*t138-772288*t140+934
     #4*t142+140416*t144+1033472*t146
      t2205 = 3200.D0/3.D0+t965+t1745-11520*t12-896.D0/3.D0*t9+152064*t1
     #7-114048*t508+t1042+13440*t20-202752*t22-76032*t24+418176*t26+2568
     #96*t28-27648*t30-t1482+20736*t516+848128*t33+1084672*t36-129664*t3
     #8
      t2221 = 436352*t40-350464*t42-t981+38016*t15-3840*l1+t1055-t1766+t
     #1767+11264.D0/3.D0*t52+6784.D0/3.D0*t54+6912*t56-2688*t58-19072*t6
     #0-7808.D0/3.D0*t62-67456.D0/3.D0*t64-37632*t66+62080*t68+9728*t70-
     #1152*t5
      t2240 = -84352*t10+283520.D0/3.D0*t35-1792.D0/3.D0*t32+10240.D0/3.
     #D0*t4-6912*t11-54016*t83+t1785+4224*t85-96128.D0/3.D0*t87+245120.D
     #0/3.D0*t89+41728.D0/3.D0*t91-88704*t93+51328.D0/3.D0*t95-t1010+729
     #6*t99+17920.D0/3.D0*t101+571904.D0/3.D0*t103+947840.D0/3.D0*t105+1
     #663744.D0/3.D0*t107
      t2259 = -949376.D0/3.D0*t109-183040.D0/3.D0*t112-1986304.D0/3.D0*t
     #114+521600.D0/3.D0*t116-251264.D0/3.D0*t118-267776.D0/3.D0*t120+40
     #3072.D0/3.D0*t122-285952.D0/3.D0*t124-2898688.D0/3.D0*t126-519808*
     #t128+453376*t130+108928*t132-t1810-937856*t134-901888.D0/3.D0*t136
     #-2688*t138-380160*t140+36864*t142+t1542+63360*t146
      t2268 = t964+3840*t154+t1472+t1473+382720*t17-145664*t508+t1115+30
     #464*t22+t1478+(t2205+t2221+t2240+t2259)*s2l1-175360*t26+314496*t28
     #-t1481-t1482+t1483-1097856*t33-2245888.D0/3.D0*t36+1322752*t38
      t2286 = 3366784.D0/3.D0*t40-888320*t42+t1052-t1053-10496.D0/3.D0*l
     #1+768*z-45696*t869+44544*t871+2688*t52+3328.D0/3.D0*t54+6144*t56-5
     #888.D0/3.D0*t58-62848.D0/3.D0*t60+4096.D0/3.D0*t62-41984*t64-18483
     #2.D0/3.D0*t66+103296*t68+15488.D0/3.D0*t70-512*t5
      t2305 = 3328.D0/3.D0-640.D0/3.D0*t2-896.D0/3.D0*t32+6784.D0/3.D0*t
     #4-t1070+23040*t891+27520.D0/3.D0*t85-286336.D0/3.D0*t87+890240.D0/
     #3.D0*t89-53632.D0/3.D0*t91-1109120.D0/3.D0*t93+519040.D0/3.D0*t95+
     #1280.D0/3.D0*t97+29824.D0/3.D0*t99+1552512*t101+663040.D0/3.D0*t10
     #3-133376*t105-4486144.D0/3.D0*t107-840832*t109
      t2323 = -73088*t112+717184*t114-1316864.D0/3.D0*t116+906496*t118-2
     #545024.D0/3.D0*t120+1136768.D0/3.D0*t122+213376.D0/3.D0*t124+48611
     #84.D0/3.D0*t126+423296.D0/3.D0*t128-87808*t130-5474816.D0/3.D0*t13
     #2-21120*t919+482432*t134-1064192.D0/3.D0*t136-7040.D0/3.D0*t138+41
     #9840*t140-t1541-t1542-535040*t146
      t2327 = -13696.D0/3.D0*t60-640*t62-45440.D0/3.D0*t64-50560*t66+141
     #440.D0/3.D0*t68+t1133+(t1668+t1687+t1705+t1723)*t644-640*t5+t1685-
     #t404*(141696*t30+5440448.D0/3.D0*t33+2723392.D0/3.D0*t36-2440576*t
     #38-7302016.D0/3.D0*t40+1839296*t42+13568*t44+1891008*t107+1988928*
     #t109+t2135-1664*t47+380352*t112-2058176.D0/3.D0*t114+44800.D0/3.D0
     #*t64+258112.D0/3.D0*t66+950912.D0/3.D0*t116+t2045-174656.D0/3.D0*t
     #68-2841536.D0/3.D0*t118-57088*t516-8512*t70+4128704.D0/3.D0*t120+t
     #2158-3121600.D0/3.D0*t122-36608*t124-7110656.D0/3.D0*t126-52800*t8
     #91-t1508+t2182-t1685-2816*t80-3584.D0/3.D0*t4-2240*t5+9600*t856)+(
     #t2268+t2286+t2305+t2323)*s1l1
      t2342 = 64.D0/3.D0*t32+1664.D0/3.D0*t4-64.D0/3.D0*t1627*t32*t168+5
     #5616*t891+37760*t85-755456.D0/3.D0*t87+1070912.D0/3.D0*t89-24832.D
     #0/3.D0*t91-831488.D0/3.D0*t93+299008.D0/3.D0*t95-192*t97+84800*t99
      t2356 = 7950464.D0/3.D0*t101+728448*t103-351040.D0/3.D0*t105-19816
     #32*t107-5803712.D0/3.D0*t109-977984.D0/3.D0*t112+2291456.D0/3.D0*t
     #114-1064192.D0/3.D0*t116+2851712.D0/3.D0*t118-3803840.D0/3.D0*t120
     #+2726144.D0/3.D0*t122+141440.D0/3.D0*t124
      t2369 = 2552256*t126+163008*t128-131264*t130-3401856*t132-327616.D
     #0/3.D0*t919+2374912.D0/3.D0*t134-978048*t136-192*t138+848320*t140-
     #10368*t142-154240*t144-1126400*t146
      t2374 = (l1+z-t138)**2
      taux9 = (t503+t1455+t1468+t1648+t2327+t2342+t2356+t2369)/t2374*t16
     #8
      t1 = l1**2
      t2 = t1*l1
      t3 = t1**2
      t4 = t3*t2
      t5 = z*t4
      t7 = z**2
      t8 = t7**2
      t9 = t8*t2
      t11 = t8*t4
      t13 = t7*z
      t14 = t8*t13
      t15 = t3*t1
      t16 = t14*t15
      t18 = t3**2
      t19 = t7*t18
      t21 = t8*t7
      t22 = t21*t4
      t24 = t18*l1
      t25 = t24*t13
      t27 = t24*z
      t29 = t8*z
      t30 = t24*t29
      t32 = t24*t8
      t34 = t24*t21
      t36 = t24*t7
      t38 = t21*t18
      t40 = t13*t18
      t42 = 1858304.D0/3.D0*t5-415360.D0/3.D0*t9-1718016*t11+1321472.D0/
     #3.D0*t16+1328640*t19+5180800.D0/3.D0*t22+392704*t25+62464*t27-1474
     #56*t30-162816*t32+166400*t34-265728*t36-832000*t38-1963520*t40
      t43 = t7*t4
      t45 = t8*t18
      t47 = t29*t4
      t49 = t29*t18
      t51 = t3*l1
      t52 = t14*t51
      t54 = t29*t1
      t56 = t13*t4
      t58 = dlog(z)
      t59 = t58**2
      t66 = z*l1
      t67 = 2560.D0/3.D0*t66
      t69 = 8192*t27
      t70 = 40960*t30
      t77 = 256+242432.D0/3.D0*t5+8192*t9+7168*t11+242944.D0/3.D0*t16+16
     #3328*t19+1034240.D0/3.D0*t22+t67+40960*t25+t69-t70+32768*t34-32768
     #*t36-164352*t38-203520*t40-1018880.D0/3.D0*t43-5120.D0/3.D0*t45
      t83 = z*t51
      t85 = t13*t2
      t87 = t8*t1
      t89 = t24*t14
      t91 = t29*t2
      t93 = t29*t3
      t95 = t8*t3
      t97 = z*t3
      t99 = t13*t3
      t101 = t7*t3
      t103 = t13*l1
      t105 = z*t1
      t106 = 2304*t105
      t107 = t7*t1
      t109 = -1338112.D0/3.D0*t47+206080*t49-40448*t52+14336*t54+1304320
     #.D0/3.D0*t56+111616.D0/3.D0*t83+57344*t85-11008.D0/3.D0*t87-8192*t
     #89-209408.D0/3.D0*t91+210176*t93-42752.D0/3.D0*t95-7424*t97-557312
     #.D0/3.D0*t99+275200.D0/3.D0*t101+4864.D0/3.D0*t103-t106-256*t107
      t111 = t7*l1
      t112 = 1024.D0/3.D0*t111
      t113 = z*t2
      t115 = t29*l1
      t117 = t8*l1
      t119 = t21*t3
      t121 = t21*t1
      t122 = 3072*t121
      t123 = t15*t8
      t125 = t21*t15
      t127 = t15*z
      t129 = t13*t15
      t131 = t15*t29
      t133 = t21*t51
      t135 = t21*t2
      t137 = t15*t7
      t139 = t7*t51
      t141 = t13*t51
      t143 = t8*t51
      t145 = t112+1792*t113-1792*t115+1024*t117-103936*t119-t122-14336*t
     #123-392192*t125-235264.D0/3.D0*t127-521984*t129+1636096.D0/3.D0*t1
     #31+788992.D0/3.D0*t133+22528*t135+1141760.D0/3.D0*t137-247808*t139
     #+1166336.D0/3.D0*t141+52736.D0/3.D0*t143
      t146 = t29*t51
      t153 = 81920*t25
      t154 = 81920*t32
      t155 = 8192*t34
      t156 = 40960*t36
      t159 = 256+62464*t5-345088.D0/3.D0*t9-901120*t11+192512*t19-99840*
     #t22+t67+t153+t69+t70-t154-t155-t156+45056*t38-401408*t40
      t175 = -367104*t43+417792*t45+479232*t47-217088*t49-5120*t54+82636
     #8*t56+35584.D0/3.D0*t83+178688.D0/3.D0*t85+25344*t87+39680*t91-503
     #296.D0/3.D0*t93+1092352.D0/3.D0*t95+1792.D0/3.D0*t97-746752.D0/3.D
     #0*t99+43264*t101
      t178 = 1280.D0/3.D0*t107
      t180 = 1024*t115
      t183 = 512*t121
      t190 = 2304*t103-t106-t178+t112+4352.D0/3.D0*t113+t180-3072*t117+1
     #6896*t119+t183+3217408.D0/3.D0*t123+111616*t125-139264.D0/3.D0*t12
     #7-925696*t129-1712128.D0/3.D0*t131-64000*t133
      t197 = t7*t2
      t199 = 2048*t1
      t200 = 256*z
      t201 = 2560.D0/3.D0*l1
      t202 = t13*t1
      t204 = z*t18
      t206 = 256*t7
      t207 = 1024*t2
      t208 = 512*t8
      t209 = -2048*t135+1076224.D0/3.D0*t137-555008.D0/3.D0*t139+616448*
     #t141-2337280.D0/3.D0*t143+1199360.D0/3.D0*t146-1536*t197-t199+t200
     #+t201-27392.D0/3.D0*t202-36864*t204+t206+t207+t208
      t215 = t14*t4
      t217 = t14*t3
      t219 = t21*l1
      t228 = 256+297472.D0/3.D0*t5-156160.D0/3.D0*t9-2477056.D0/3.D0*t11
     #+217088*t19-188416.D0/3.D0*t22+t67+t153+t69+t70-t154-t155-t156+368
     #64*t38
      t243 = -417792*t40-1430528.D0/3.D0*t43+401408*t45+1103360.D0/3.D0*
     #t47-192512*t49-256*t54+2695168.D0/3.D0*t56+179456.D0/3.D0*t83+3061
     #76.D0/3.D0*t85+12032.D0/3.D0*t87+4864.D0/3.D0*t91-44800*t93+240896
     #*t95-41984.D0/3.D0*t97
      t259 = -1028864.D0/3.D0*t99+457984.D0/3.D0*t101+7936.D0/3.D0*t103-
     #2560*t105+1280*t107-512.D0/3.D0*t111+2304*t113-512*t117+2768896.D0
     #/3.D0*t123+143360.D0/3.D0*t125-326656.D0/3.D0*t127-3184640.D0/3.D0
     #*t129-361472*t131-40960.D0/3.D0*t133
      t269 = 561152*t137-383488*t139+2277376.D0/3.D0*t141-1831424.D0/3.D
     #0*t143+188672*t146-92672.D0/3.D0*t197-256*t13-t199+t200+t201-64768
     #.D0/3.D0*t202-45056*t204+t206+t207
      t274 = t14*t18
      t276 = t14*t2
      t278 = t14*t1
      t280 = -418304*t146+(t159+t175+t190+t209)*s1l1-47104.D0/3.D0*t197-
     #t199+256*t29+t200+t201-245504.D0/3.D0*t215+9472*t217+512*t219-3430
     #4.D0/3.D0*t202+(t228+t243+t259+t269)*s2l1-122624.D0/3.D0*t204+t206
     #+t207+123136.D0/3.D0*t274-1280*t276+256*t278
      t290 = 28928.D0/3.D0*t66
      t291 = 3072*t25
      t292 = 2560*t27
      t293 = 512*t30
      t294 = 2816*t32
      t295 = 768*t34
      t296 = 1280*t36
      t300 = -2048-143488.D0/3.D0*t5+383104*t9+52864.D0/3.D0*t11-14208*t
     #19+8448*t22-t290+t291-t292-t293-t294+t295+t296-3712*t38-7168*t40+1
     #90976.D0/3.D0*t43
      t304 = 67712*t54
      t315 = 45824*t103
      t316 = 111232.D0/3.D0*t105
      t318 = 10112*t45-64640.D0/3.D0*t47+3072*t49+t304-97024.D0/3.D0*t56
     #-343168.D0/3.D0*t83-1155200.D0/3.D0*t85-518528.D0/3.D0*t87-489472.
     #D0/3.D0*t91+226944*t93-1498880.D0/3.D0*t95+104320*t97+1524736.D0/3
     #.D0*t99-924032.D0/3.D0*t101-t315+t316-216320.D0/3.D0*t107
      t320 = 39808.D0/3.D0*t111
      t322 = 13312*t115
      t323 = 44928*t117
      t325 = 6656*t121
      t332 = 25728*t135
      t337 = t320-209408.D0/3.D0*t113-t322+t323-118016.D0/3.D0*t119-t325
     #-470912.D0/3.D0*t123-50048.D0/3.D0*t125+263936.D0/3.D0*t127+555008
     #.D0/3.D0*t129+86784*t131+94336.D0/3.D0*t133+t332-166272*t137+84390
     #4.D0/3.D0*t139-1228288.D0/3.D0*t141+382976*t143
      t342 = 2944*t13
      t343 = 39808.D0/3.D0*t1
      t344 = 768*t24
      t346 = 2048*z
      t347 = 20480.D0/3.D0*l1
      t352 = 3456*t7
      t353 = 15232.D0/3.D0*t2
      t354 = 6656*t8
      t355 = -185984*t146+204032*t197+24320*t51+35968.D0/3.D0*t4+t342+t3
     #43+t344-68224.D0/3.D0*t3-t346-t347-19840*t15+540416.D0/3.D0*t202-4
     #480*t18+16384*t204-t352+t353-t354
      t365 = -2747776*t43+814080*t45-1582080*t47+737280*t49-626560.D0/3.
     #D0*t52+78080.D0/3.D0*t54+12460160.D0/3.D0*t56+t59*(t77+t109+t145+t
     #280)/2+(t300+t318+t337+t355)*s1l1+858880.D0/3.D0*t83+362240*t85+44
     #800.D0/3.D0*t87-45568*t89-572800.D0/3.D0*t91+665856*t93
      t377 = t18*t2
      t378 = t377*t21
      t379 = 354816*t378
      t380 = t377*t14
      t381 = 50688*t380
      t383 = 2048+29630464.D0/3.D0*t5-1154304*t9-133551232.D0/3.D0*t11-2
     #59200*t16+29224832*t19-18558848.D0/3.D0*t22-t379+t381+t290+2922547
     #2*t25
      t384 = t18*t1
      t385 = z*t384
      t387 = t384*t13
      t389 = t384*t8
      t391 = t384*t7
      t393 = t377*z
      t394 = 354816*t393
      t395 = t377*t7
      t396 = 1064448*t395
      t402 = -2255616*t385-10771200*t387+10517760*t389+6614784*t391+t394
     #-t396+6272512*t27+15792640*t30-27869440*t32-4915456*t34-18234112*t
     #36
      t404 = t377*t13
      t405 = 1774080*t404
      t406 = t377*t8
      t407 = 1774080*t406
      t408 = t377*t29
      t409 = 1064448*t408
      t411 = t384*t29
      t413 = t384*t14
      t415 = t384*t21
      t421 = t405-t407+t409+6921344*t38-6158592*t411-278784*t413+2002176
     #*t415-46620800*t40-29895552*t43+43573120*t45+23233024*t47
      t432 = -23742080*t49+49920*t52-t304+144670720.D0/3.D0*t56+2474240*
     #t83+1038080*t85+246656*t87+646656*t89+391424*t91-5969152.D0/3.D0*t
     #93+14239744.D0/3.D0*t95
      t441 = -1847680.D0/3.D0*t97-14317952.D0/3.D0*t99+7153024.D0/3.D0*t
     #101+t315-t316+78080*t107-t320+368128.D0/3.D0*t113+t322-t323+243712
     #*t119
      t451 = t325+92506496.D0/3.D0*t123+10778240.D0/3.D0*t125-18719360.D
     #0/3.D0*t127-100259456.D0/3.D0*t129-15339648*t131-1284352*t133-t332
     #+20048512*t137-8710400*t139+46508288.D0/3.D0*t141
      t460 = -14572544*t143+20415232.D0/3.D0*t146-435712*t197-273920*t51
     #-4032896.D0/3.D0*t4-t342-t343-918272*t24+178304.D0/3.D0*t3+t346+t3
     #47
      t467 = 50688*t377
      t469 = 609792*t215+329472*t384+781696*t15-212992*t202+1426560*t18-
     #9963904*t204+t352-t353+t354-t467-819072*t274
      t483 = 577280*t95-156160.D0/3.D0*t97-4347904.D0/3.D0*t99+643840*t1
     #01+(t383+t402+t421+t432+t441+t451+t460+t469)*s1l1-1280000.D0/3.D0*
     #t119+1987584*t123-5747200.D0/3.D0*t125-1811456.D0/3.D0*t127-144972
     #80.D0/3.D0*t129+1904640*t131+3651200.D0/3.D0*t133+64000*t135+30192
     #64*t137
      t494 = 12800.D0/3.D0*t66
      t505 = -1280-1212160.D0/3.D0*t5-27776.D0/3.D0*t9+1441280.D0/3.D0*t
     #11-1082624.D0/3.D0*t16-839680*t19-4416640.D0/3.D0*t22-t494-228864*
     #t25-40960*t27+155648*t30+49152*t32-139776*t34+168448*t36+701440*t3
     #8+1137920*t40+1739648*t43-711680.D0/3.D0*t45
      t523 = 1705728*t47-784640*t49+182400*t52-69888*t54-7239296.D0/3.D0
     #*t56-555008.D0/3.D0*t83-269568*t85+52480.D0/3.D0*t87+36352*t89+303
     #744*t91-2570752.D0/3.D0*t93-253952.D0/3.D0*t95+34048*t97+2796032.D
     #0/3.D0*t99-1348352.D0/3.D0*t101-24320.D0/3.D0*t103+32000.D0/3.D0*t
     #105+t178
      t543 = -5120.D0/3.D0*t111-17920.D0/3.D0*t113+10496*t115-5888*t117+
     #1373440.D0/3.D0*t119+17920*t121-517632*t123+1677568*t125+1176320.D
     #0/3.D0*t127+2844416*t129-6314752.D0/3.D0*t131-3396224.D0/3.D0*t133
     #-322816.D0/3.D0*t135-5792000.D0/3.D0*t137+1236096*t139-6150784.D0/
     #3.D0*t141+307968*t143+4934912.D0/3.D0*t146
      t551 = 96768*t378
      t552 = 13824*t380
      t555 = -1280-8250752.D0/3.D0*t5+625024*t9+43560704.D0/3.D0*t11+633
     #60*t16-8411904*t19+1924480*t22+t551-t552-14336.D0/3.D0*t66-8180736
     #*t25
      t560 = 96768*t393
      t561 = 290304*t395
      t567 = 615168*t385+2937600*t387-2868480*t389-1804032*t391-t560+t56
     #1-1723904*t27-4411904*t30+7815680*t32+1358336*t34+5069312*t36
      t569 = 483840*t404
      t570 = 483840*t406
      t571 = 290304*t408
      t580 = -t569+t570-t571-5970688.D0/3.D0*t38+1679616*t411+76032*t413
     #-546048*t415+41217152.D0/3.D0*t40+8963712*t43-12982016*t45-7591808
     #*t47
      t592 = 7036800*t49-11520*t52+31232*t54-15250432*t56-1826176.D0/3.D
     #0*t83-1270912.D0/3.D0*t85-132608*t87-175104*t89-213504*t91+3000448
     #.D0/3.D0*t93-2280960*t95
      t606 = 371200.D0/3.D0*t97+5756032.D0/3.D0*t99-680704*t101-16640*t1
     #03+42752.D0/3.D0*t105-35584.D0/3.D0*t107+3328.D0/3.D0*t111-76288.D
     #0/3.D0*t113-6144*t115+19456*t117-110592*t119
      t617 = -t122-33782528.D0/3.D0*t123-3766016.D0/3.D0*t125+5068928.D0
     #/3.D0*t127+34363520.D0/3.D0*t129+17101952.D0/3.D0*t131+513408*t133
     #+12288*t135-18624512.D0/3.D0*t137+8102272.D0/3.D0*t139-17331968.D0
     #/3.D0*t141
      t628 = 1280*z
      t629 = 12800.D0/3.D0*l1
      t630 = 18196480.D0/3.D0*t143-2935168*t146+92672*t197+173056.D0/3.D
     #0*t51+340736*t4+512*t13+27392.D0/3.D0*t1+248320*t24-12544*t3-t628-
     #t629
      t640 = 13824*t377
      t642 = -156672*t215-89856*t384-183808*t15+75776*t202-1132288.D0/3.
     #D0*t18+2768000*t204-1536*t7-2560.D0/3.D0*t2-3072*t8+t640+217728*t2
     #74
      t648 = 28160.D0/3.D0*t1
      t661 = -1280-1814784*t5+1331584.D0/3.D0*t9+15592704*t11+220032*t16
     #-6928896*t19+2939776*t22+t551-t552-t494-7800576*t25
      t671 = 546048*t385+2868480*t387-2937600*t389-1679616*t391-t560+t56
     #1-1343744*t27-5096192*t30+8195840*t32+1738496*t34+4385024*t36
      t681 = -t569+t570-t571-8548864.D0/3.D0*t38+1804032*t411+89856*t413
     #-615168*t415+38845568.D0/3.D0*t40+22252928.D0/3.D0*t43-13855744*t4
     #5-28082432.D0/3.D0*t47
      t693 = 8571264*t49-74880*t52+11008*t54-43607168.D0/3.D0*t56-465664
     #*t83-590336*t85-161792.D0/3.D0*t87-251136*t89-349312.D0/3.D0*t91+2
     #536064.D0/3.D0*t93-6307328.D0/3.D0*t95
      t706 = 94208*t97+6868864.D0/3.D0*t99-2840320.D0/3.D0*t101-49664.D0
     #/3.D0*t103+38144.D0/3.D0*t105-40192.D0/3.D0*t107+4864.D0/3.D0*t111
     #-38656.D0/3.D0*t113-t180+6400*t117-172800*t119
      t717 = -t183-11967488*t123-5764864.D0/3.D0*t125+1159296*t127+34099
     #328.D0/3.D0*t129+20225920.D0/3.D0*t131+2290816.D0/3.D0*t133+63232.
     #D0/3.D0*t135-16651264.D0/3.D0*t137+8533120.D0/3.D0*t139-18492800.D
     #0/3.D0*t141
      t727 = 6200832*t143-3108864*t146+543232.D0/3.D0*t197+9472*t51+4067
     #84.D0/3.D0*t4+1536*t13+t648+172288*t24-4096*t3-t628-t629
      t735 = 1280*t7
      t738 = -374400*t215+11520*t217-76032*t384-144640.D0/3.D0*t15+11571
     #2*t202-615424.D0/3.D0*t18+1926784*t204-t735-5632.D0/3.D0*t2-t208+t
     #640+392832*t274
      t749 = 77824*t197+(t555+t567+t580+t592+t606+t617+t630+t642)*s1l1+1
     #024*t51+t648-1536*t29-3072*t3-t628-t629+1090688.D0/3.D0*t215-13542
     #4.D0/3.D0*t217-3072*t219+(t661+t671+t681+t693+t706+t717+t727+t738)
     #*s2l1+156416.D0/3.D0*t202+613120.D0/3.D0*t204-t735-6400.D0/3.D0*t2
     #-546560.D0/3.D0*t274+7424*t276-1536*t278
      t754 = 26*t2
      t755 = 80*l1
      t756 = 72*t3
      t757 = 24*t51
      t758 = 172*t1
      t774 = 2*l6-1
      t775 = t774**2
      t778 = l5-1
      t780 = -1+l1
      t782 = (1-l1+t66)**2
      t783 = 1/t782
      t788 = dsqrt(z)
      t790 = dsqrt(l5)
      t792 = dsqrt(-t778)
      t795 = dsqrt(l1)
      t799 = l5*z*l1-l1+1-l5+l5*l1
      t802 = dsqrt(-t780)
      t807 = t799**2
      t810 = dsqrt(4*l1*t775*z*l5*t778*t780*t783-4*t774*t788*t790*t792*t
     #795*t799*t783*t802+t807*t783)
      t818 = 20480.D0/3.D0*t66
      t822 = -2048-81536.D0/3.D0*t5+221696*t9+17408.D0/3.D0*t11-6272*t19
     #+32896.D0/3.D0*t22-t818+t291-t292-t293-t294+t295+t296-4480*t38-133
     #12*t40+62080.D0/3.D0*t43
      t837 = 104192.D0/3.D0*t103
      t838 = 68224.D0/3.D0*t105
      t839 = 116864.D0/3.D0*t107
      t840 = 10368*t45-77440.D0/3.D0*t47+5120*t49+24704*t54+21760.D0/3.D
     #0*t56-156032.D0/3.D0*t83-721024.D0/3.D0*t85-273152.D0/3.D0*t87-212
     #992.D0/3.D0*t91+119936*t93-912256.D0/3.D0*t95+141952.D0/3.D0*t97+9
     #05728.D0/3.D0*t99-169472*t101-t837+t838-t839
      t842 = 15232.D0/3.D0*t111
      t843 = 101888.D0/3.D0*t113
      t844 = 4352*t115
      t845 = 19328*t117
      t847 = 2176*t121
      t859 = t842-t843-t844+t845-44416.D0/3.D0*t119-t847-97664*t123-1574
     #4*t125+41728*t127+223232.D0/3.D0*t129+72960*t131+16896*t133+8576*t
     #135-64128*t137+405248.D0/3.D0*t139-666112.D0/3.D0*t141+721792.D0/3
     #.D0*t143
      t864 = 3200*t13
      t865 = 44032.D0/3.D0*t1
      t866 = 12032*t3
      t871 = 2048*t7
      t872 = 2432.D0/3.D0*t2
      t873 = 2176*t8
      t874 = -363136.D0/3.D0*t146+354944.D0/3.D0*t197+38272.D0/3.D0*t51+
     #24832.D0/3.D0*t4+t864+t865+t344-t866-t346-t347-34688.D0/3.D0*t15+3
     #61216.D0/3.D0*t202-3712*t18+12288*t204-t871-t872-t873
      t892 = 2048+17800832.D0/3.D0*t5-3317888.D0/3.D0*t9-147740672.D0/3.
     #D0*t11-868992*t16+23489408*t19-10365952*t22-t379+t381+t818+2783155
     #2*t25
      t902 = -2002176*t385-10517760*t387+10771200*t389+6158592*t391+t394
     #-t396+4878592*t27+18301696*t30-29263360*t32-6309376*t34-15725056*t
     #36
      t912 = t405-t407+t409+10173568*t38-6614784*t411-329472*t413+225561
     #6*t415-43525760*t40-68690432.D0/3.D0*t43+46949760*t45+92898176.D0/
     #3.D0*t47
      t924 = -29643392*t49+309120*t52-46464*t54+134215040.D0/3.D0*t56+36
     #86272.D0/3.D0*t83+1145728*t85+461312.D0/3.D0*t87+925440*t89+449664
     #*t91-2771712*t93+16072832.D0/3.D0*t95
      t931 = -241024*t97-5009536*t99+6020096.D0/3.D0*t101+t837-t838+t839
     #-t842+t843+t844-t845+2086016.D0/3.D0*t119
      t942 = t847+35004032*t123+20481664.D0/3.D0*t125-3407232*t127-94402
     #688.D0/3.D0*t129-64479104.D0/3.D0*t131-8510080.D0/3.D0*t133-261248
     #.D0/3.D0*t135+45543296.D0/3.D0*t137-20461952.D0/3.D0*t139+45638528
     #.D0/3.D0*t141
      t950 = -50508928.D0/3.D0*t143+9783168*t146-364928*t197-148352.D0/3
     #.D0*t51-1671680.D0/3.D0*t4-t864-t865-639488*t24+t866+t346+t347
      t959 = 1428864*t215-49920*t217+278784*t384+673408.D0/3.D0*t15-6300
     #16.D0/3.D0*t202+787072*t18-6764928*t204+t871+t872+t873-t467-146572
     #8*t274
      t965 = -5653120.D0/3.D0*t139+3382912*t141-4128256.D0/3.D0*t143-141
     #5936*t146-t58*(t505+t523+t543+t749)-281600.D0/3.D0*t197-256.D0/3.D
     #0*(t754+t755+t756-t757-t758+24-t58*(110*t1-25*t2-15-50*l1+12*t51-3
     #6*t3)+t59*(3+12*t2+10*l1-24*t1)/2)*t810+(t822+t840+t859+t874)*s2l1
     #-1355648.D0/3.D0*t215+113920.D0/3.D0*t217-125440.D0/3.D0*t202-3123
     #20*t204+256.D0/3.D0*(t754+t755+t756-t757-t758+24)*t810+(t892+t902+
     #t912+t924+t931+t942+t950+t959)*s2l1+227840*t274
      t969 = (l1+z-t66)**2
      taux10 = (t42+t365+t483+t965)/t969*t783

CCC Coefficients

      tcoff15 = taux9
      tcoff16 = taux9
      tcoff17 = taux10
      tcoff18 = -taux10

CCC Jet functions

      tFjet15 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff15)
      tFjet16 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff16)
      tFjet17 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff17)
      tFjet18 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff18)

CCC Sum

      RVnum5e2 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet18
     #*tcoff18
      RETURN
      END

      double precision function RVnum5e3(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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
      t2 = t1*z
      t3 = l1**2
      t4 = t3**2
      t5 = t2*t4
      t8 = t3*l1
      t9 = t4*t8
      t10 = t1*t9
      t11 = 1280*t10
      t12 = t1**2
      t13 = t12*z
      t14 = t13*t9
      t15 = 1280*t14
      t16 = t13*t3
      t17 = 384*t16
      t18 = t12*t1
      t19 = t18*t4
      t21 = t12*l1
      t22 = 384*t21
      t24 = t4**2
      t25 = t18*t24
      t27 = t2*t24
      t30 = t12*t24
      t33 = t13*t24
      t35 = 128*t16
      t36 = 2432.D0/3.D0*t4
      t38 = 128*t21
      t39 = t4*t3
      t40 = 2176*t39
      t41 = -69376.D0/3.D0*t5+256*t25-5120*t27-17920*t10+3840*t30+7168*t
     #14-1536*t33-t35+t36+2432.D0/3.D0*t19+t38+t40
      t42 = t39*t12
      t44 = t18*t39
      t46 = t39*z
      t48 = t2*t39
      t50 = t39*t13
      t52 = t4*l1
      t53 = t18*t52
      t55 = t18*t8
      t56 = 128*t55
      t57 = t39*t1
      t59 = t1*t52
      t61 = t2*t52
      t63 = t12*t52
      t65 = t13*t52
      t67 = z*t52
      t69 = 33664*t42+2176*t44-13312*t46-45056*t48-13312*t50-1920*t53-t5
     #6+33664*t57-97664.D0/3.D0*t59+44032*t61-97664.D0/3.D0*t63+37376.D0
     #/3.D0*t65+37376.D0/3.D0*t67
      t71 = t2*t9
      t73 = t18*t9
      t75 = 1920*t52
      t76 = 256*t24
      t77 = z*t3
      t78 = 128*t77
      t79 = t1*l1
      t80 = 128*t79
      t81 = 3584.D0/3.D0*t9
      t82 = 128*t2
      t83 = t1*t24
      t85 = t2*t3
      t86 = 1792.D0/3.D0*t85
      t87 = t2*l1
      t88 = 2048.D0/3.D0*t87
      t89 = t1*t3
      t91 = z*t8
      t93 = 71680.D0/3.D0*t71-3584.D0/3.D0*t73-t75+t76-t78+t80-t81+t82+3
     #840*t83+t86-t88-128.D0/3.D0*t89+4352.D0/3.D0*t91
      t94 = t1*t8
      t96 = t12*t8
      t98 = t2*t8
      t100 = t12*t3
      t102 = t13*t8
      t104 = t13*t4
      t106 = t12*t4
      t108 = z*t4
      t110 = t1*t4
      t112 = t12*t9
      t114 = z*t9
      t116 = z*t24
      t118 = 128*t8
      t119 = -3968*t94-3968*t96+15872.D0/3.D0*t98-128.D0/3.D0*t100+4352.
     #D0/3.D0*t102-18304.D0/3.D0*t104+50560.D0/3.D0*t106-18304.D0/3.D0*t
     #108+50560.D0/3.D0*t110-17920*t112+7168*t114-1536*t116-t118
      t135 = 58112.D0/3.D0*t5+t11-t15+t17-2048.D0/3.D0*t19+t22+(t41+t69+
     #t93+t119)*s2l1-31232.D0/3.D0*t42-2816.D0/3.D0*t44+3584.D0/3.D0*t46
     #+32768.D0/3.D0*t48+14848.D0/3.D0*t50+3712.D0/3.D0*t53+t56-17152.D0
     #/3.D0*t57+31360.D0/3.D0*t59-58624.D0/3.D0*t61+17664*t63-22912.D0/3
     #.D0*t65-2176*t67
      t136 = 2560*t71
      t137 = 256*t73
      t138 = 256*t79
      t140 = 3200.D0/3.D0*t87
      t152 = 2560*t112
      t153 = 256*t114
      t154 = -t136+t137+t78+t138+t82+13312.D0/3.D0*t85-t140-5632.D0/3.D0
     #*t89-2432.D0/3.D0*t91+17536.D0/3.D0*t94+8960*t96-35072.D0/3.D0*t98
     #-2816*t100-2432*t102+6016*t104-16384*t106+1920*t108-10240*t110+t15
     #2-t153
      t157 = t24*l1
      t158 = t157*t1
      t160 = t157*t18
      t162 = dlog(z)
      t163 = dsqrt(z)
      t165 = dsqrt(l5)
      t166 = t165*l5
      t167 = t163*z*t166
      t168 = 1-l5
      t169 = dsqrt(t168)
      t170 = dsqrt(l1)
      t172 = 1-l1
      t173 = dsqrt(t172)
      t174 = t169*t170*t173
      t177 = l5**2
      t178 = t177*l1
      t179 = l6**2
      t180 = t2*t179
      t183 = l5*l1
      t184 = t2*l6
      t191 = t177*t12
      t195 = l5*t1
      t198 = l5*t3
      t204 = t177*t2
      t210 = t177*t1
      t215 = t163*t2*t166
      t216 = t170*l1
      t218 = t169*t216*t173
      t222 = t163*t1
      t223 = t222*t166
      t230 = t216*t173*l6
      t233 = -4*t204*t8+t191*t8-4*t195*t3-t210*l1+2*t204*t3-4*t215*t218-
     #t210+2*t195+4*t223*t174+8*t223*t218+8*t215*t169*t230
      t235 = t167*t169
      t236 = t235*t230
      t238 = t222*t165
      t241 = t223*t169
      t244 = z*l5
      t248 = t170*t173*l6
      t253 = t210*t3
      t255 = t177*z
      t256 = t255*t3
      t259 = t210*t8
      t261 = t255*t8
      t263 = t244*l1
      t265 = t244*t3
      t270 = t177*t8
      t271 = t167*t218
      t273 = t177*t3
      t276 = 6*t259-4*t261-4*t263+2*t265-t177-8*t235*t248-2*t255+t270-4*
     #t271-3*t273+3*t178
      t278 = 4*t167*t174+16*t178*t180+16*t183*t184-16*t178*t184-16*t183*
     #t180-t191*t3+6*t178*t2+2*t195*l1+2*t198*t2-6*t183*t2+t233+8*t236-4
     #*t238*t174-16*t241*t230+2*t244+8*t238*t169*t248+t79-t1-8*t241*t248
     #-4*t253+6*t256+t276
      t280 = t162**2
      t285 = z*l1
      t287 = (1-l1+t285)**2
      t288 = 1/t287
      t295 = t18*l1
      t297 = t12*t2
      t298 = t297*t39
      t301 = t157*t2
      t303 = t24*t3
      t304 = z*t303
      t306 = t303*t2
      t308 = t157*z
      t311 = t303*t13
      t313 = t303*t297
      t315 = t303*t18
      t317 = t297*t4
      t322 = 704-231360*t5-301312*t158-306944*t160-320*t13-640*t295-2942
     #08.D0/3.D0*t298-1920*t285+756800*t301-9088*t304-137600*t306+49984*
     #t308+2123008.D0/3.D0*t25-138880*t311-9344*t313+55808*t315-23168.D0
     #/3.D0*t317-5277184.D0/3.D0*t27-879680*t10+7075072.D0/3.D0*t30
      t326 = t157*t297
      t336 = 4608*t304
      t337 = 69120*t306
      t340 = 69120*t311
      t341 = 4608*t313
      t342 = 27648*t315
      t343 = 3840*t317
      t348 = -76928.D0/3.D0*t5-128768*t158-158464*t160-65664*t298-12544.
     #D0/3.D0*t285+340480*t301-t336-t337+20224*t308+378496*t25-t340-t341
     #+t342-t343-701056*t27-791936.D0/3.D0*t10+1041152*t30+3165056.D0/3.
     #D0*t14
      t352 = 512*t4
      t355 = t18*t3
      t356 = 384*t355
      t368 = 4096.D0/3.D0-862720*t33+21376.D0/3.D0*t16+27648*t326+t352+7
     #2832.D0/3.D0*t19+2304*t21-t356+779648*t42+1047424.D0/3.D0*t44-3059
     #2*t46-1450240.D0/3.D0*t48-725248*t50-138112*t53+512.D0/3.D0*t55+17
     #6128*t57-261760.D0/3.D0*t59+164608*t61-234240*t63
      t373 = t303*t1
      t374 = 27648*t373
      t376 = t297*t9
      t378 = t297*t24
      t384 = t297*t52
      t386 = t303*t12
      t387 = 92160*t386
      t389 = 256.D0/3.D0*t2
      t394 = 754816.D0/3.D0*t65+18432*t67+2316928.D0/3.D0*t71+t374-14483
     #20.D0/3.D0*t73+90240*t376-68736*t378-16000.D0/3.D0*l1+1024*z+12928
     #.D0/3.D0*t77-10624.D0/3.D0*t79+24960*t384+t387+16000.D0/3.D0*t3+t3
     #89+249472*t83+896*t85-4480.D0/3.D0*t87+15616.D0/3.D0*t89
      t400 = t13*l1
      t401 = 768*t400
      t407 = 2048.D0/3.D0*t1
      t411 = t157*t12
      t414 = t157*t13
      t417 = -1408*t91-37504.D0/3.D0*t94+40832*t96+4096.D0/3.D0*t98-4288
     #0.D0/3.D0*t100-t401-68096.D0/3.D0*t102-32384.D0/3.D0*t104-56704.D0
     #/3.D0*t106-16000.D0/3.D0*t108+110720.D0/3.D0*t110+t407-3628288.D0/
     #3.D0*t112+115840.D0/3.D0*t114-36608*t116-478720*t411-2048*t8+37760
     #0*t414-384*t12
      t422 = t297*t8
      t423 = 1152*t422
      t424 = t297*t3
      t440 = -1109120.D0/3.D0*t5-175360*t158-145664*t160-21120*t298-7040
     #.D0/3.D0*t285+419840*t301-t336-t337+30464*t308+314496*t25-t340-t34
     #1+t342-1097856*t27-2245888.D0/3.D0*t10+1322752*t30+3366784.D0/3.D0
     #*t14-888320*t33
      t447 = 34560*t298
      t460 = -88704*t5+418176*t158-114048*t160-6912*t303-t447-2688*t285-
     #380160*t301+36864*t304+t337-202752*t308+256896*t25-27648*t311-t341
     #+20736*t315+848128*t27+1084672*t10-129664*t30+436352*t14-350464*t3
     #3
      t461 = 25344*t326
      t478 = -t17+t461-896.D0/3.D0*t4+7296*t19-1152*t21-54016*t39+17920.
     #D0/3.D0*t42+571904.D0/3.D0*t44+947840.D0/3.D0*t46+1663744.D0/3.D0*
     #t48-949376.D0/3.D0*t50-183040.D0/3.D0*t53+t56-1986304.D0/3.D0*t57+
     #521600.D0/3.D0*t59-251264.D0/3.D0*t61-267776.D0/3.D0*t63+403072.D0
     #/3.D0*t65-285952.D0/3.D0*t67
      t483 = 62208*t376
      t485 = 56064*t378
      t487 = 2176.D0/3.D0*z
      t491 = 7680*t384
      t498 = 3200.D0/3.D0-2898688.D0/3.D0*t71-76032*t373-901888.D0/3.D0*
     #t73+t483+13440*t52-t485-3840*l1+t487-84352*t24+6784.D0/3.D0*t77-26
     #88*t79+t491-11520*t386+10240.D0/3.D0*t3+283520.D0/3.D0*t9-1792.D0/
     #3.D0*t2-937856*t83-19072*t85
      t519 = 11264.D0/3.D0*t87+6912*t89-7808.D0/3.D0*t91-67456.D0/3.D0*t
     #94-37632*t96+62080*t98+9728*t100+4224*t102-96128.D0/3.D0*t104+2451
     #20.D0/3.D0*t106+41728.D0/3.D0*t108+51328.D0/3.D0*t110+384*t1+10892
     #8*t112-519808*t114+453376*t116+63360*t411-1152*t8+152064*t414+3801
     #6*t157
      t538 = 1408*t16+23040*t326+t352+(t460+t478+t498+t519)*s2l1+29824.D
     #0/3.D0*t19+1280.D0/3.D0*t21-256*t355+1552512*t42+663040.D0/3.D0*t4
     #4-133376*t46-4486144.D0/3.D0*t48-840832*t50-73088*t53+640*t55+7171
     #84*t57-1316864.D0/3.D0*t59+906496*t61-2545024.D0/3.D0*t63+1136768.
     #D0/3.D0*t65
      t557 = 213376.D0/3.D0*t67+4861184.D0/3.D0*t71+t374-1064192.D0/3.D0
     #*t73+44544*t376-45696*t378-10496.D0/3.D0*l1+768*z+3328.D0/3.D0*t77
     #-5888.D0/3.D0*t79+3840*t384+t387+6784.D0/3.D0*t3-896.D0/3.D0*t2+48
     #2432*t83-62848.D0/3.D0*t85+2688*t87+6144*t89+4096.D0/3.D0*t91
      t568 = 1280.D0/3.D0*t1
      t576 = 3328.D0/3.D0-41984*t94-184832.D0/3.D0*t96+103296*t98+15488.
     #D0/3.D0*t100-512*t400+27520.D0/3.D0*t102-286336.D0/3.D0*t104+89024
     #0.D0/3.D0*t106-53632.D0/3.D0*t108+519040.D0/3.D0*t110+t568-5474816
     #.D0/3.D0*t112+423296.D0/3.D0*t114-87808*t116-535040*t411-512*t8+38
     #2720*t414-640.D0/3.D0*t12
      t590 = (t440+t538+t557+t576)*s1l1+3008*t355+6896000.D0/3.D0*t42+19
     #40608.D0/3.D0*t44-306112.D0/3.D0*t46-5124160.D0/3.D0*t48-5089216.D
     #0/3.D0*t50-854528.D0/3.D0*t53-44672.D0/3.D0*t55+1977088.D0/3.D0*t5
     #7-304192*t59
      t610 = 3456*t303
      t611 = 48384*t298
      t614 = 18432*t304
      t615 = 34560*t306
      t618 = 13824*t311
      t619 = 2304*t313
      t620 = 10368*t315
      t626 = -140096.D0/3.D0*t5+177984*t158-67392*t160-t610-t611-6784.D0
     #/3.D0*t285-155520*t301+t614+t615-85248*t308+185472*t25-t618-t619+t
     #620-t343+271616*t27+308672*t10+74240*t30+457984*t14-285312*t33
      t628 = 14976*t326
      t633 = 64*t355
      t646 = 1984.D0/3.D0+4160.D0/3.D0*t16+t628+192*t4+102464.D0/3.D0*t1
     #9+512.D0/3.D0*t21-19136.D0/3.D0*t39-t633+821056.D0/3.D0*t42+245632
     #*t44+152768.D0/3.D0*t46+119936.D0/3.D0*t48-1290304.D0/3.D0*t50-125
     #824*t53-3584*t55-376064.D0/3.D0*t57-8448*t59+165248.D0/3.D0*t61-52
     #9024.D0/3.D0*t63+233600*t65
      t650 = 38016*t373
      t652 = 59136*t376
      t653 = 1024.D0/3.D0*t52
      t654 = 40704*t378
      t660 = 21120*t384
      t661 = 5760*t386
      t666 = 11008.D0/3.D0*t85
      t668 = 1408.D0/3.D0*t67-639104.D0/3.D0*t71-t650-836288.D0/3.D0*t73
     #+t652+t653-t654-8512.D0/3.D0*l1+1472.D0/3.D0*z-26624*t24+7808.D0/3
     #.D0*t77-6592.D0/3.D0*t79+t660-t661+9728.D0/3.D0*t3+61888.D0/3.D0*t
     #9-320.D0/3.D0*t2-334464*t83-t666-256.D0/3.D0*t87
      t675 = 128*t400
      t690 = 64*t12
      t691 = -11840.D0/3.D0*t108+28800*t110+320*t1-218304*t112-136320*t1
     #14+155776*t116+2880*t411-1408*t8+96768*t414+15552*t157-t690
      t695 = s2l1**2
      t703 = 1600.D0/3.D0*z+4160.D0/3.D0*t77-1984*t79+(t626+t646+t668+12
     #416.D0/3.D0*t89-896*t91-35968.D0/3.D0*t94-8960.D0/3.D0*t96+17792*t
     #98-2048*t100-t675+4608*t102-194944.D0/3.D0*t104+54528*t106+t691)*t
     #695+113600.D0/3.D0*t384+184320*t386+6272.D0/3.D0*t3-t82+2092544.D0
     #/3.D0*t83-27328.D0/3.D0*t85+1152*t87
      t718 = 1088.D0/3.D0*t1
      t728 = 3840*t298
      t737 = 514240.D0/3.D0*t5+240192*t158-46656*t160-t610-t728-1664.D0/
     #3.D0*t285-224640*t301+t614+t615-117504*t308+81792*t25-t618-t619+t6
     #20+617216*t27+965504*t10-213760*t30+35968*t14
      t739 = 256*t16
      t740 = 10368*t326
      t742 = 64*t19
      t756 = 1216.D0/3.D0-77952*t33-t739+t740-20608.D0/3.D0*t4-t742-1472
     #.D0/3.D0*t21-222784.D0/3.D0*t39-1072960.D0/3.D0*t42+80384.D0/3.D0*
     #t44+1266880.D0/3.D0*t46+2579072.D0/3.D0*t48+13632*t50-3712*t53-886
     #656*t57+1568384.D0/3.D0*t59-492032*t61+585152.D0/3.D0*t63-54272.D0
     #/3.D0*t65
      t761 = 13440*t376
      t763 = 17664*t378
      t771 = 640.D0/3.D0*t2
      t773 = -713216.D0/3.D0*t67-2809088.D0/3.D0*t71-t650-205568.D0/3.D0
     #*t73+t761+33728*t52-t763-3008.D0/3.D0*l1+704.D0/3.D0*z-61184*t24+1
     #600*t77+1600.D0/3.D0*t79-t661+320.D0/3.D0*t3+268288.D0/3.D0*t9-t77
     #1-645504*t83-t86+t140
      t792 = 64.D0/3.D0*t12
      t793 = -26176.D0/3.D0*t89-54272.D0/3.D0*t91+51136*t94+11840.D0/3.D
     #0*t96-92032.D0/3.D0*t98+5056.D0/3.D0*t100+512*t102+14272.D0/3.D0*t
     #104-165184.D0/3.D0*t106+255872.D0/3.D0*t108-603776.D0/3.D0*t110+37
     #1584*t112-471040*t114+317056*t116+60480*t411+1792.D0/3.D0*t8+55296
     #*t414+22464*t157+t792
      t796 = s1l1**2
      t799 = -26240.D0/3.D0*t108+270784.D0/3.D0*t110+t718-9035264.D0/3.D
     #0*t112+426880.D0/3.D0*t114-344704.D0/3.D0*t116-1013760*t411-640*t8
     #+763840*t414+(t737+t756+t773+t793)*t796+256.D0/3.D0*t12
      t806 = 2304*t304
      t814 = 163264.D0/3.D0*t5+(t135+t154)*s1l1-76032*t158-76032*t160+64
     #.D0/3.D0*(-t162*t278-t280*t278/2)*t12*t288-t162*(1920*t400+54784*t
     #373-128192.D0/3.D0*t96+46144*t98+14848.D0/3.D0*t100+112832.D0/3.D0
     #*t102+t423+t322+231296.D0/3.D0*t19+(t348+t368+t394+t417)*s2l1+4268
     #80.D0/3.D0*t376-5316992.D0/3.D0*t33-320*t424-117248*t378+51392*t32
     #6-634688.D0/3.D0*t104+288576*t106-57472.D0/3.D0*t94-4288.D0/3.D0*t
     #21+2248448*t14+2241536*t71-883136*t73+2397440.D0/3.D0*t61-3195584.
     #D0/3.D0*t63+t703-23744.D0/3.D0*t16+12992.D0/3.D0*t89-2496*l1+256*t
     #4+775296*t65+41920*t67+t590+t799+512*t91)-8448*t298+190080*t301-t8
     #06-t615+12672*t308+165248*t25-34560*t311-t619+13824*t315+1920*t317
     #-413696*t27-157248*t10
      t821 = 3328*t25
      t822 = 7168*t27
      t824 = 15872*t30
      t826 = 12288*t33
      t829 = 3200.D0/3.D0*t21
      t830 = 128*t355
      t831 = -35072.D0/3.D0*t5+t821-t822+3328*t10+t824+52736*t14-t826+75
     #52.D0/3.D0*t16+4160.D0/3.D0*t19+t829-t830
      t843 = 354944.D0/3.D0*t42+22592*t44+20992.D0/3.D0*t46-50688*t48-26
     #9312.D0/3.D0*t50-42880.D0/3.D0*t53+4352.D0/3.D0*t55-22336.D0/3.D0*
     #t57+33344.D0/3.D0*t59+119488.D0/3.D0*t61-321728.D0/3.D0*t63
      t849 = 768*t83
      t850 = 384*t87
      t851 = 384*t89
      t852 = 128*t91
      t855 = 227008.D0/3.D0*t65-15232.D0/3.D0*t67+30464*t71-14336*t73-t8
     #49+t666-t850-t851-t852+12416.D0/3.D0*t94-14016*t96
      t858 = 256*t400
      t866 = 1024*t116
      t867 = 128*t12
      t868 = -3840*t98+896*t100-t858+6464.D0/3.D0*t102-92288.D0/3.D0*t10
     #4+53696*t106+1536*t108-30208.D0/3.D0*t110-67840*t112-4352*t114+t86
     #6-t867
      t873 = 55552.D0/3.D0*t5
      t874 = 128*t285
      t875 = 128*t25
      t876 = 2560*t27
      t878 = 1920*t30
      t880 = 768*t33
      t881 = 4096.D0/3.D0*t4
      t882 = 5056.D0/3.D0*t39
      t884 = -t873+t874+t875-t876-9600*t10+t878+3072*t14-t880+t881+t742+
     #t882+45376.D0/3.D0*t42
      t887 = 23808*t48
      t892 = 80128.D0/3.D0*t61
      t896 = 35840.D0/3.D0*t71
      t898 = 1856.D0/3.D0*t44-27520.D0/3.D0*t46-t887-14720.D0/3.D0*t50-1
     #024.D0/3.D0*t53+61376.D0/3.D0*t57-24576*t59+t892-45056.D0/3.D0*t63
     #+3968*t65+33920.D0/3.D0*t67+t896-1408.D0/3.D0*t73
      t900 = 2048*t52
      t901 = 128*t24
      t902 = 3200.D0/3.D0*t77
      t903 = 1792.D0/3.D0*t79
      t904 = 64*t3
      t905 = 2176.D0/3.D0*t9
      t906 = 1920*t83
      t907 = 6400.D0/3.D0*t85
      t908 = 256*t87
      t912 = -t900+t901-t902-t903+t904-t905+t906-t907+t908+2880*t89+1177
     #6.D0/3.D0*t91-27008.D0/3.D0*t94
      t914 = 24320.D0/3.D0*t98
      t915 = 384*t100
      t916 = 256*t102
      t921 = 64*t1
      t924 = 768*t116
      t925 = 1408.D0/3.D0*t8
      t926 = -2816*t96+t914+t915+t916-4864.D0/3.D0*t104+8704*t106-8448*t
     #108+55360.D0/3.D0*t110+t921-8320*t112+4096*t114-t924-t925
      t942 = 138368*t5+30208*t27+95104*t10-19712*t30-22144*t14+4864*t33+
     #t35-t36-384*t19-t38-t40-483712.D0/3.D0*t42
      t953 = 124288.D0/3.D0*t67
      t955 = -640.D0/3.D0*t44+143104.D0/3.D0*t46+762112.D0/3.D0*t48+1126
     #40.D0/3.D0*t50+1280.D0/3.D0*t53+t56-526976.D0/3.D0*t57+169344*t59-
     #745856.D0/3.D0*t61+437248.D0/3.D0*t63-27392*t65-t953-410368.D0/3.D
     #0*t71
      t962 = 128.D0/3.D0*t73+t75-t76+t78-t80+t81-t82-20992*t83+18688.D0/
     #3.D0*t85+t88-1664*t89-9472.D0/3.D0*t91
      t975 = 72064.D0/3.D0*t94+33664.D0/3.D0*t96-43264*t98+5248.D0/3.D0*
     #t100+2048.D0/3.D0*t102+18944.D0/3.D0*t104-66816*t106+52864.D0/3.D0
     #*t108-270592.D0/3.D0*t110+89344*t112-26752*t114+5888*t116+t118
      t979 = -118144.D0/3.D0*t5+t821-t822+3584*t10+t824+57856*t14-t826-t
     #17+20608.D0/3.D0*t19-t22+(t942+t955+t962+t975)*s2l1
      t991 = 147840*t42+83392.D0/3.D0*t44+8960*t46-197632.D0/3.D0*t48-33
     #2032.D0/3.D0*t50-66752.D0/3.D0*t53-t56-24128.D0/3.D0*t57+32128.D0/
     #3.D0*t59+201536.D0/3.D0*t61-472768.D0/3.D0*t63+329024.D0/3.D0*t65
      t999 = -23168.D0/3.D0*t67+34048*t71-15616*t73-t78-t138-t82-t849-29
     #44*t85+t140+1024*t89-128.D0/3.D0*t91
      t1011 = 3584.D0/3.D0*t94-39360*t96+39872.D0/3.D0*t98+27008.D0/3.D0
     #*t100+44416.D0/3.D0*t102-176768.D0/3.D0*t104+300608.D0/3.D0*t106+8
     #320.D0/3.D0*t108-22336.D0/3.D0*t110-75008*t112-4864*t114+t866
      t1024 = 557824*t30+456192*t14-418176*t33+5120*t16+12672*t326+(t831
     #+t843+t855+t868)*s2l1-67520.D0/3.D0*t19+(t884+t898+t912+t926)*t796
     #+(t979+t991+t999+t1011)*s1l1+937088.D0/3.D0*t42+230656.D0/3.D0*t44
     #-4736.D0/3.D0*t46-592000.D0/3.D0*t48-231424*t50+26752.D0/3.D0*t53+
     #7040*t55+49536*t57+29376*t59
      t1034 = t280*t162
      t1036 = 2*l6-1
      t1037 = t1036**2
      t1050 = t263-l1+1-l5+t183
      t1056 = t1050**2
      t1059 = dsqrt(4*l1*t1037*z*l5*t168*t172*t288-4*t1036*t163*t165*t16
     #9*t170*t1050*t288*t173+t1056*t288)
      t1062 = 12*t8
      t1066 = 12*t4
      t1068 = t163*t165
      t1069 = t1068*t169
      t1070 = t170*t8
      t1072 = t1070*t173*l6
      t1073 = t1069*t1072
      t1075 = t1068*t174
      t1078 = l5*t4
      t1080 = t1059*t8
      t1082 = t163*t166
      t1083 = t170*t3
      t1085 = t169*t1083*t173
      t1086 = t1082*t1085
      t1088 = t1069*t248
      t1090 = t1082*t174
      t1092 = t8*l6
      t1093 = t244*t1092
      t1096 = t1083*t173*l6
      t1097 = t1069*t1096
      t1099 = t235*t1096
      t1102 = t169*t1070*t173
      t1103 = t167*t1102
      t1105 = t4*l6
      t1106 = t244*t1105
      t1108 = t4*t179
      t1109 = t244*t1108
      t1111 = t255*t1092
      t1113 = t255*t1105
      t1115 = -7-t1066+56*t236-96*t1073-28*t1075+38*t1059+24*t1078-232*t
     #1080+120*t1086+56*t1088+28*t1090+480*t1093+240*t1097-144*t1099-48*
     #t1103-192*t1106+192*t1109-480*t1111+192*t1113
      t1116 = t255*t1108
      t1118 = t8*t179
      t1119 = t255*t1118
      t1121 = t244*t1118
      t1123 = l1*l6
      t1124 = t244*t1123
      t1126 = l1*t179
      t1127 = t244*t1126
      t1129 = t255*t1126
      t1131 = t3*l6
      t1132 = t255*t1131
      t1134 = t255*t1123
      t1136 = t3*t179
      t1137 = t255*t1136
      t1139 = t244*t1131
      t1141 = t1068*t1102
      t1145 = t1082*t169
      t1146 = t1145*t1096
      t1149 = t1059*t1
      t1150 = t1149*t8
      t1152 = t1059*z
      t1153 = t1152*l1
      t1155 = t1145*t230
      t1157 = t177*t4
      t1159 = -192*t1116+480*t1119-480*t1121+112*t1124-112*t1127+112*t11
     #29+400*t1132-112*t1134-400*t1137-400*t1139+48*t1141+32*l1+14*l5-24
     #0*t1146-55*t3-100*t1150+76*t1153+200*t1155-12*t1157
      t1165 = t167*t1085
      t1168 = t1068*t1085
      t1170 = t1069*t230
      t1172 = t210*t4
      t1174 = t1149*t4
      t1175 = 66*t1174
      t1176 = t255*l1
      t1178 = t244*t4
      t1180 = t244*t8
      t1182 = t255*t4
      t1186 = t1068*t218
      t1189 = t1082*t218
      t1191 = -7*t253-150*t256+18*t259+180*t261+72*t1165-42*t263-120*t11
     #68-200*t1170-12*t1172+t1175+42*t1176+72*t1178-180*t1180-72*t1182+1
     #50*t265+42*t8+100*t1186-7*t177-100*t1189
      t1193 = t1145*t1072
      t1195 = t235*t1072
      t1197 = t244*t1136
      t1199 = t1059*t4
      t1202 = t1152*t4
      t1204 = t1152*t3
      t1206 = t1149*t3
      t1208 = t1152*t8
      t1210 = t1082*t1102
      t1212 = t1145*t248
      t1215 = t1059*l1
      t1217 = t1059*t3
      t1222 = l5*t8
      t1224 = -64*t183+96*t1193+96*t1195+400*t1197+66*t1199+42*t270-132*
     #t1202-276*t1204+38*t1206+332*t1208-48*t1210-56*t1212-28*t271-176*t
     #1215+304*t1217+110*t198-55*t273+32*t178-84*t1222
      t1252 = -16384*t5-t11+t15-5632.D0/3.D0*t16+1920*t19-t829+t830+3276
     #8.D0/3.D0*t42+3584.D0/3.D0*t44-2816.D0/3.D0*t46-31232.D0/3.D0*t48-
     #17152.D0/3.D0*t50-2176*t53-2432.D0/3.D0*t55+14848.D0/3.D0*t57-2291
     #2.D0/3.D0*t59+17664*t61-58624.D0/3.D0*t63+31360.D0/3.D0*t65
      t1264 = 3712.D0/3.D0*t67+t136-t137-2816*t85+t850+t851+t852-2432*t9
     #4-35072.D0/3.D0*t96+8960*t98+13312.D0/3.D0*t100+t858+17536.D0/3.D0
     #*t102-10240*t104+58112.D0/3.D0*t106-2048.D0/3.D0*t108+6016*t110-t1
     #52+t153+t867
      t1269 = 15104*t27
      t1271 = 9856*t30
      t1273 = 2432*t33
      t1276 = 108224*t5-t874+t1269+52672*t10-t1271-11968*t14+t1273-t881-
     #t742-t882-97216*t42-256.D0/3.D0*t44
      t1288 = 64.D0/3.D0*t73
      t1289 = 31040*t46+465152.D0/3.D0*t48+22464*t50+128*t53-109568*t57+
     #363904.D0/3.D0*t59-171008*t61+99584*t63-57664.D0/3.D0*t65-32832*t6
     #7-225152.D0/3.D0*t71+t1288
      t1291 = 10496*t83
      t1296 = t900-t901+t902+t903-t904+t905-t1291+22400.D0/3.D0*t85-t908
     #-4800*t89-17536.D0/3.D0*t91+79808.D0/3.D0*t94
      t1305 = 2944*t116
      t1306 = 39808.D0/3.D0*t96-39360*t98-t915-t916+6528*t104-53824*t106
     #+56384.D0/3.D0*t108-228736.D0/3.D0*t110-t921+48640*t112-15040*t114
     #+t1305+t925
      t1318 = 64*t4
      t1320 = 1792.D0/3.D0*t21
      t1321 = 1856.D0/3.D0*t39
      t1322 = 45248*t5+t1269+42432*t10-t1271-10176*t14+t1273+2560.D0/3.D
     #0*t16-t1318-1216.D0/3.D0*t19+t1320-t1321-t633
      t1335 = -197056.D0/3.D0*t42-128*t44+52672.D0/3.D0*t46+305408.D0/3.
     #D0*t48+46400.D0/3.D0*t50+320*t53+256*t55-205184.D0/3.D0*t57+169280
     #.D0/3.D0*t59-263680.D0/3.D0*t61+156928.D0/3.D0*t63-28928.D0/3.D0*t
     #65
      t1339 = 1408.D0/3.D0*t9
      t1341 = 256*t91
      t1343 = -35584.D0/3.D0*t67-185216.D0/3.D0*t71+t1288+t653-t901+t133
     #9-t1291+3840*t85-t908-t851-t1341+18688.D0/3.D0*t94
      t1354 = 5056*t96-46720.D0/3.D0*t98-960*t100-t675-832*t102+6080.D0/
     #3.D0*t104-66304.D0/3.D0*t106+3328*t108-25984*t110+40704*t112-11712
     #*t114+t1305-t690
      t1359 = -49408.D0/3.D0*t61-58816.D0/3.D0*t63+8768*t65-24320.D0/3.D
     #0*t67+436032*t71+13824*t373-173376*t73+25920*t376-t1034*(-64.D0/3.
     #D0*t1059*(30*t3-60*l1+t1062+21)*s1l1-64.D0/3.D0*(t1115+t1159+t1191
     #+t1224)*t288-64.D0/3.D0*t1059*(78*t3-104*l1-t1062+41)*s2l1)/6+(t12
     #52+t1264)*s2l1-26880*t378+(t1276+t1289+t1296+t1306)*t796-2880*t384
     #+46080*t386+161664*t83+10240.D0/3.D0*t85+(t1322+t1335+t1343+t1354)
     #*t695+20480.D0/3.D0*t94
      t1388 = -3904.D0/3.D0+886144.D0/3.D0*t5+308352*t158+313984*t160+25
     #6*t13+512*t295+112832*t298+3072*t285-772288*t301+9344*t304+140416*
     #t306-51392*t308-2207360.D0/3.D0*t25+141696*t311+9600*t313-57088*t3
     #15+32384.D0/3.D0*t317+5440448.D0/3.D0*t27+2723392.D0/3.D0*t10-2440
     #576*t30
      t1397 = 1152*t303
      t1400 = 9472*t304
      t1401 = 59136*t306
      t1404 = 41728*t311
      t1405 = 14976*t315
      t1410 = 758336.D0/3.D0*t5+172480*t158+95040*t160-t1397+t611+9856.D
     #0/3.D0*t285-336000*t301+t1400+t1401-47232*t308-256512*t25+t1404+t6
     #19-t1405+t343+840064*t27+1650176.D0/3.D0*t10-2981888.D0/3.D0*t30-3
     #074240.D0/3.D0*t14
      t1428 = 2052608.D0/3.D0*t33+9472.D0/3.D0*t16-t628+1280*t4-55680*t1
     #9+t38-9472.D0/3.D0*t39-4098560.D0/3.D0*t42-339456*t44+89792*t46+33
     #01760.D0/3.D0*t48+943040*t50+182208*t53+7680*t55-1419776.D0/3.D0*t
     #57+794048.D0/3.D0*t59-653184*t61+2434304.D0/3.D0*t63-1628672.D0/3.
     #D0*t65-t953
      t1431 = 32128*t373
      t1439 = 64384*t386
      t1442 = 128.D0/3.D0*t2
      t1445 = -4160.D0/3.D0-3627392.D0/3.D0*t71-t1431+381696*t73-t652+12
     #80.D0/3.D0*t52+t654+12736.D0/3.D0*l1-3200.D0/3.D0*z-9088*t24-5120.
     #D0/3.D0*t77+9280.D0/3.D0*t79-t660-t1439-2304*t3+22976.D0/3.D0*t9-t
     #1442-1217408.D0/3.D0*t83+28928.D0/3.D0*t85
      t1466 = -896.D0/3.D0*t87-14848.D0/3.D0*t89-7040.D0/3.D0*t91+60032.
     #D0/3.D0*t94+65216*t96-62720*t98-17600.D0/3.D0*t100-107840.D0/3.D0*
     #t102+188608*t104-908800.D0/3.D0*t106+11136*t108-93568*t110-2240.D0
     #/3.D0*t1+4421696.D0/3.D0*t112-360896.D0/3.D0*t114+301184.D0/3.D0*t
     #116+382656*t411-4736.D0/3.D0*t8-257152*t414+5184*t157
      t1476 = 33536*t306
      t1479 = 33536*t311
      t1480 = 13568*t315
      t1487 = -229696*t5+34304*t158+49920*t160-2880*t298+12032.D0/3.D0*t
     #285-92928*t301+t806+t1476-5376*t308-157504.D0/3.D0*t25+t1479+t619-
     #t1480-87104.D0/3.D0*t27-608320.D0/3.D0*t10-29056*t30+523456.D0/3.D
     #0*t14+227008.D0/3.D0*t33-3456*t16
      t1509 = 37120*t5+353408*t158+164736*t160-2304*t303+t447+4224*t285-
     #650496*t301+18944*t304+118272*t306-104192*t308-371328*t25+83456*t3
     #11+t341-29952*t315+1467392*t27+2729344.D0/3.D0*t10-4822144.D0/3.D0
     #*t30-3746560.D0/3.D0*t14+3127552.D0/3.D0*t33
      t1526 = -5888.D0/3.D0+t739-t461-1024.D0/3.D0*t4-11520*t19+512*t21-
     #21632*t39-3822976.D0/3.D0*t42-277632*t44+156800*t46+3309824.D0/3.D
     #0*t48+2496640.D0/3.D0*t50+90240*t53-1661312.D0/3.D0*t57+159744*t59
     #-1056640.D0/3.D0*t61+437376*t63-872320.D0/3.D0*t65-127744.D0/3.D0*
     #t67
      t1543 = -5176192.D0/3.D0*t71-64256*t373+435456*t73-t483+7680*t52+t
     #485+6144*l1-3968.D0/3.D0*z-28544*t24-768*t77+3072*t79-t491-128768*
     #t386-5888.D0/3.D0*t3+100480.D0/3.D0*t9+t771-2372992.D0/3.D0*t83+21
     #376.D0/3.D0*t85-2048*t87
      t1562 = -4096*t89-8960.D0/3.D0*t91+6272.D0/3.D0*t94-37888.D0/3.D0*
     #t96-9856.D0/3.D0*t98+10240.D0/3.D0*t100+2688*t102+37504*t104-14336
     #0.D0/3.D0*t106+12160.D0/3.D0*t108-41984.D0/3.D0*t110-t407+5765888.
     #D0/3.D0*t112-792704.D0/3.D0*t114+696832.D0/3.D0*t116+708224*t411-2
     #944*t8-459008*t414+12672*t157
      t1579 = -9216*t326+768*t4+16832*t19-5888.D0/3.D0*t21+(t1509+t1526+
     #t1543+t1562)*s2l1+t356+714752*t42+239296.D0/3.D0*t44-51264*t46-199
     #2256.D0/3.D0*t48-1119808.D0/3.D0*t50-59328*t53-1408*t55+888640.D0/
     #3.D0*t57-618304.D0/3.D0*t59+524416*t61-1799680.D0/3.D0*t63+310592*
     #t65+95168.D0/3.D0*t67
      t1582 = 13568*t373
      t1585 = 256*t52
      t1592 = 44544*t386
      t1598 = -5632.D0/3.D0+410048*t71-t1582-20096*t73-4608*t376+t1585+1
     #2480*t378+17792.D0/3.D0*l1-3712.D0/3.D0*z-1280.D0/3.D0*t77+6272.D0
     #/3.D0*t79+1920*t384-t1592-2560.D0/3.D0*t3+t1442+90560.D0/3.D0*t83-
     #8704.D0/3.D0*t85+64.D0/3.D0*t87+2176.D0/3.D0*t89
      t1617 = -3584.D0/3.D0*t91-12352*t94-86592*t96+154624.D0/3.D0*t98+5
     #1520.D0/3.D0*t100+t401+31168*t102-408640.D0/3.D0*t104+867328.D0/3.
     #D0*t106-9344*t108+72064*t110-1792.D0/3.D0*t1-1187648.D0/3.D0*t112+
     #116480.D0/3.D0*t114-7744*t116+135936*t411-11584.D0/3.D0*t8-112640*
     #t414+1024.D0/3.D0*t12
      t1635 = -576+218560*t5+180928*t158+69696*t160-t1397+t728+2048.D0/3
     #.D0*t285-314496*t301+t1400+t1401-56960*t308-129792*t25+t1404+t619-
     #t1405+732544*t27+1927616.D0/3.D0*t10-2125568.D0/3.D0*t30
      t1653 = -1339328.D0/3.D0*t14+1223168.D0/3.D0*t33+t17-t740-14336.D0
     #/3.D0*t4+2240.D0/3.D0*t21-82688.D0/3.D0*t39-2062336.D0/3.D0*t42-57
     #600*t44+604736.D0/3.D0*t46+868480*t48+288704*t50+11520*t53-1770752
     #.D0/3.D0*t57+1081088.D0/3.D0*t59-521088*t61+1071424.D0/3.D0*t63-10
     #8864*t65-112320*t67
      t1670 = -2969216.D0/3.D0*t71-t1431+121152*t73-t761+14912*t52+t763+
     #5696.D0/3.D0*l1-256*z-20608*t24+1792*t77+1088*t79-t1439+640.D0/3.D
     #0*t3+93056.D0/3.D0*t9+1024.D0/3.D0*t2-1344128.D0/3.D0*t83+42112.D0
     #/3.D0*t85-5888.D0/3.D0*t87-17920.D0/3.D0*t89
      t1688 = -7296*t91+38784*t94+24768*t96-195008.D0/3.D0*t98-9280.D0/3
     #.D0*t100-1152*t102+20096*t104-364160.D0/3.D0*t106+35520*t108-44147
     #2.D0/3.D0*t110-t921+876032*t112-221120*t114+447104.D0/3.D0*t116+32
     #5568*t411-2048.D0/3.D0*t8-201856*t414+7488*t157-t792
      t1701 = 256*t424+(t1635+t1653+t1670+t1688)*t796-2816*t355-7934464.
     #D0/3.D0*t42-2276672.D0/3.D0*t44+290240.D0/3.D0*t46+1891008*t48+198
     #8928*t50+380352*t53+71744.D0/3.D0*t55-2058176.D0/3.D0*t57
      t1724 = 4032*l1-2944.D0/3.D0*z-2176.D0/3.D0*t77+2624*t79-145024.D0
     #/3.D0*t384-187904*t386-3584.D0/3.D0*t3-t389-715520*t83+25600.D0/3.
     #D0*t85-t908
      t1754 = -2688+325824*t5+114944*t158+99328*t160+63936*t298+19072.D0
     #/3.D0*t285-275968*t301+t806+t1476-19968*t308-898624.D0/3.D0*t25+t1
     #479+t619-t1480+5760*t317+2658496.D0/3.D0*t27+1877120.D0/3.D0*t10-1
     #119616*t30-3959168.D0/3.D0*t14
      t1773 = 2381248.D0/3.D0*t33+4736*t16-16128*t326+2304*t4-245696.D0/
     #3.D0*t19-3136.D0/3.D0*t21+t830-5653120.D0/3.D0*t42-1373824.D0/3.D0
     #*t44+287104.D0/3.D0*t46+4456768.D0/3.D0*t48+1296128*t50+769792.D0/
     #3.D0*t53+36160.D0/3.D0*t55-599808*t57+339072*t59-2682304.D0/3.D0*t
     #61+3445760.D0/3.D0*t63-2319296.D0/3.D0*t65
      t1791 = -45248*t67-4498688.D0/3.D0*t71-t1582+1454528.D0/3.D0*t73-7
     #3152*t376-t1585+47040*t378+24832.D0/3.D0*l1-2048*z-2688*t77+15616.
     #D0/3.D0*t79-29760*t384-t1592-10496.D0/3.D0*t3-384*t2-1119040.D0/3.
     #D0*t83+25280.D0/3.D0*t85+2944.D0/3.D0*t87-17152.D0/3.D0*t89
      t1809 = -11456.D0/3.D0*t91+20096*t94+231680.D0/3.D0*t96-70656*t98-
     #4352.D0/3.D0*t100+t858-52864*t102+273536*t104-1250944.D0/3.D0*t106
     #+11712*t108-111360*t110-1408*t1+5667520.D0/3.D0*t112-321856.D0/3.D
     #0*t114+65216*t116+354048*t411-11456.D0/3.D0*t8-256256*t414+t867
      t1818 = 7680*t108-267584.D0/3.D0*t110-1984.D0/3.D0*t1+3234112*t112
     #+(t1754+t1773+t1791+t1809)*s2l1-424768.D0/3.D0*t114+351296.D0/3.D0
     #*t116+1033472*t411-2240*t8-779328*t414-t792
      t1831 = -t873+t875-t876-8320*t10+t878+4096*t14-t880-3200.D0/3.D0*t
     #16+t1318+4096.D0/3.D0*t19-t1320+t1321
      t1842 = t633+61376.D0/3.D0*t42+5056.D0/3.D0*t44-14720.D0/3.D0*t46-
     #t887-27520.D0/3.D0*t50-2048*t53-1408.D0/3.D0*t55+45376.D0/3.D0*t57
     #-45056.D0/3.D0*t59+t892-24576*t63+33920.D0/3.D0*t65
      t1847 = 3968*t67+t896-2176.D0/3.D0*t73-t653+t901-t1339+t906-t907+t
     #908+t851+t1341-2816*t94
      t1857 = -27008.D0/3.D0*t96+t914+2880*t100+t675+11776.D0/3.D0*t102-
     #8448*t104+55360.D0/3.D0*t106-4864.D0/3.D0*t108+8704*t110-9600*t112
     #+3072*t114-t924+t690
      t1893 = -8-18*t4+64*t236-144*t1073-32*t1075+41*t1059+36*t1078-406*
     #t1080+180*t1086+64*t1088+32*t1090+720*t1093+360*t1097-216*t1099-72
     #*t1103-288*t1106+288*t1109-720*t1111+288*t1113-288*t1116
      t1911 = t1149*t39
      t1913 = t1152*t39
      t1915 = t1149*t52
      t1919 = 43*l1+16*l5-360*t1146-80*t3-144*t1150+82*t1153+12*t1911-24
     #*t1913-30*t1915+280*t1155-18*t1157
      t1926 = t1059*t39
      t1932 = t1152*t52
      t1944 = -8*t253-210*t256+27*t259+270*t261+12*t1926+108*t1165-48*t2
     #63-180*t1168-280*t1170+84*t1932-18*t1172+116*t1174+48*t1176+108*t1
     #178-270*t1180-108*t1182+210*t265+63*t8+140*t1186-8*t177
      t1963 = t1059*t52
      t1968 = 520*t1208-72*t1210-64*t1212-32*t271-226*t1215+445*t1217+16
     #0*t198-54*t1963-80*t273+43*t178-126*t1222
      t1996 = -184064.D0/3.D0*t5-67584*t158-67584*t160-64*t13-128*t295-7
     #0976.D0/3.D0*t298-4096.D0/3.D0*t285+168960*t301-2048*t304-30720*t3
     #06+11264*t308+157824*t25-30720*t311-2048*t313+12288*t315-2368*t317
     #-396096*t27-596480.D0/3.D0*t10+1591552.D0/3.D0*t30+1564736.D0/3.D0
     #*t14
      t2006 = 10752*t27
      t2008 = 7168*t30
      t2010 = 1792*t33
      t2014 = 2176.D0/3.D0-130496.D0/3.D0*t5-4352.D0/3.D0*t285-t2006-320
     #00*t10+t2008+7808*t14-t2010+448.D0/3.D0*t16+896.D0/3.D0*t21+54912*
     #t42
      t2026 = -13632*t46-86016*t48-12480*t50+57216*t57-156544.D0/3.D0*t5
     #9+81152*t61-47232*t63+24064.D0/3.D0*t65+10240*t67+47616*t71-5504.D
     #0/3.D0*l1
      t2032 = 448.D0/3.D0*t2
      t2033 = 7168*t83
      t2038 = 576*z+1024*t77-3968.D0/3.D0*t79+1408*t3+t2032+t2033-t907-1
     #280.D0/3.D0*t87+1536*t89-19456.D0/3.D0*t94-3328.D0/3.D0*t96
      t2047 = 8064*t114
      t2048 = 1792*t116
      t2049 = 256*t8
      t2050 = 13568*t98-5504.D0/3.D0*t100-2944.D0/3.D0*t102-2176.D0/3.D0
     #*t104+57856.D0/3.D0*t106-3008*t108+76928.D0/3.D0*t110+t568-31488*t
     #112+t2047-t2048-t2049
      t2065 = 1984.D0/3.D0-397376*t33-10304.D0/3.D0*t16+11264*t326+25984
     #*t19+128.D0/3.D0*t21+320*t422-64*t424+(t2014+t2026+t2038+t2050)*t6
     #95+768*t355+1692544.D0/3.D0*t42+159488*t44-64448.D0/3.D0*t46-40825
     #6*t48-420416*t50-238208.D0/3.D0*t53-5632*t55+150144*t57-69376*t59+
     #605312.D0/3.D0*t61
      t2072 = 10240*t158
      t2073 = 2048*t160
      t2075 = 20480*t301
      t2076 = 2048*t308
      t2085 = 960+58432*t5+t2072+t2073-5632.D0/3.D0*t285-t2075-t2076-112
     #64*t25+101888*t27+95872*t10-105472*t30-120448*t14+54528*t33+1280*t
     #16-4224*t19
      t2100 = 768*t21-t830-813568.D0/3.D0*t42-27904*t44+40576.D0/3.D0*t4
     #6+238336*t48+429184.D0/3.D0*t50+16000*t53+512*t55-287488.D0/3.D0*t
     #57+149888.D0/3.D0*t59-156544*t61+583936.D0/3.D0*t63-299456.D0/3.D0
     #*t65-12736.D0/3.D0*t67
      t2117 = -211968*t71+24960*t73-7936.D0/3.D0*l1+1984.D0/3.D0*z+960*t
     #77-3712.D0/3.D0*t79+1536*t3-49152*t83-1472*t85+832.D0/3.D0*t87+390
     #4.D0/3.D0*t89+1600.D0/3.D0*t91-3328.D0/3.D0*t94+26752*t96-8832*t98
      t2127 = 20480*t411
      t2128 = 10240*t414
      t2149 = 25088*t48+12416.D0/3.D0*t50-48896.D0/3.D0*t57+13824*t59-23
     #040*t61+47104.D0/3.D0*t63-3712*t65-8320.D0/3.D0*t67-13824*t71-8192
     #.D0/3.D0*l1+t487
      t2172 = 512.D0/3.D0*t98-512.D0/3.D0*t100-896.D0/3.D0*t102+1664*t10
     #4-23296.D0/3.D0*t106+1280.D0/3.D0*t108-14336.D0/3.D0*t110+t568+921
     #6*t112-2304*t114+512*t116
      t2176 = -17728.D0/3.D0*t100-t858-9920*t102+125440.D0/3.D0*t104-264
     #896.D0/3.D0*t106-448.D0/3.D0*t108-31168.D0/3.D0*t110+t718+228352*t
     #112-16768*t114+9472*t116+t2127+t2049-t2128+(1024+29824.D0/3.D0*t5-
     #5888.D0/3.D0*t285+3072*t27+9216*t10-2048*t30-2304*t14+512*t33-5043
     #2.D0/3.D0*t42+11648.D0/3.D0*t46+t2149+1152*t77-4352.D0/3.D0*t79+17
     #92*t3-2048*t83-5888.D0/3.D0*t85+1664.D0/3.D0*t87+5888.D0/3.D0*t89+
     #384*t91-2560.D0/3.D0*t94+5632.D0/3.D0*t96+t2172)*s2l1-t867
      t2194 = -287360*t63+216064*t65+8192*t67+1546816.D0/3.D0*t71+(t2085
     #+t2100+t2117+t2176)*s1l1+12288*t373-610816.D0/3.D0*t73+97216.D0/3.
     #D0*t376-77888.D0/3.D0*t378-5248.D0/3.D0*l1+512*z+832*t77-3328.D0/3
     #.D0*t79+10112*t384+40960*t386+1152*t3+t2032+156800*t83-3904.D0/3.D
     #0*t85-1216.D0/3.D0*t87
      t2218 = 896.D0/3.D0-163648.D0/3.D0*t5-512*t285-t2006-32512*t10+t20
     #08+8320*t14-t2010-t35-256*t21+183680.D0/3.D0*t42
      t2230 = -41024.D0/3.D0*t46-92416*t48-43712.D0/3.D0*t50+178304.D0/3
     #.D0*t57-166528.D0/3.D0*t59+92416*t61-175744.D0/3.D0*t63+34432.D0/3
     #.D0*t65+30592.D0/3.D0*t67+49152*t71-896*l1
      t2241 = 448.D0/3.D0*z+448.D0/3.D0*t77-256.D0/3.D0*t79+384*t3-t82+t
     #2033-15616.D0/3.D0*t85+640*t87-640.D0/3.D0*t89+256.D0/3.D0*t91-588
     #8*t94
      t2251 = -25600.D0/3.D0*t96+20480*t98+1024*t100+384*t102-11072.D0/3
     #.D0*t104+30976*t106-7936.D0/3.D0*t108+82816.D0/3.D0*t110-33024*t11
     #2+t2047-t2048+t2049
      t2265 = -196288.D0/3.D0*t5-t2072-t2073-2816*t285+t2075+t2076+9216*
     #t25-102912*t27-342272.D0/3.D0*t10+99328*t30+270848.D0/3.D0*t14-478
     #72*t33+1088.D0/3.D0*t16+2176.D0/3.D0*t21
      t2280 = 645632.D0/3.D0*t42+35840.D0/3.D0*t44-75776.D0/3.D0*t46-737
     #792.D0/3.D0*t48-257536.D0/3.D0*t50-10240.D0/3.D0*t53+389632.D0/3.D
     #0*t57-254464.D0/3.D0*t59+493696.D0/3.D0*t61-391168.D0/3.D0*t63+121
     #792.D0/3.D0*t65+40384.D0/3.D0*t67+648448.D0/3.D0*t71-47104.D0/3.D0
     #*t73
      t2296 = -3584*l1+1088*z+1920*t77-2304*t79+2560*t3+1088.D0/3.D0*t2+
     #53248*t83-3136.D0/3.D0*t85-1088*t87+2240*t89+64*t91-6400*t94-1280.
     #D0/3.D0*t96+14848*t98
      t2307 = 4160.D0/3.D0-11456.D0/3.D0*t100-5824.D0/3.D0*t102-5824*t10
     #4+119360.D0/3.D0*t106-3072*t108+95168.D0/3.D0*t110+2368.D0/3.D0*t1
     #-600832.D0/3.D0*t112+70912.D0/3.D0*t114-11008*t116-t2127-t2049+t21
     #28
      t2311 = 1088*t89+192*t91-10240.D0/3.D0*t94-14080*t96+34432.D0/3.D0
     #*t98-64*t100+448*t400+16896*t102-216320.D0/3.D0*t104+268352.D0/3.D
     #0*t106-4672.D0/3.D0*t108+19648*t110+t718-698368*t112+93632.D0/3.D0
     #*t114-77120.D0/3.D0*t116-225280*t411+168960*t414+(t2218+t2230+t224
     #1+t2251)*t796+(t2265+t2280+t2296+t2307)*s2l1
      t2320 = 12*t52
      t2332 = t170*t4
      t2334 = t169*t2332*t173
      t2338 = t2332*t173*l6
      t2341 = t170*t52
      t2343 = t169*t2341*t173
      t2346 = t52*l6
      t2349 = t52*t179
      t2360 = 15-62*t4-168*t1068*t2334+336*t1069*t2338-48*t1082*t2343-67
     #2*t255*t2346+672*t255*t2349-672*t244*t2349-120*t236-160*t1073+60*t
     #1075-76*t1059+124*t1078+231*t1080
      t2371 = t2341*t173*l6
      t2379 = -220*t1086-120*t1088-12*t39+120*t167*t2334-60*t1090-880*t1
     #093-440*t1097+360*t1099-96*t1069*t2371+40*t1103-320*t1106+320*t110
     #9+880*t1111+320*t1113
      t2395 = -320*t1116-880*t1119+880*t1121-240*t1124+240*t1127-240*t11
     #29-960*t1132+240*t1134+960*t1137+960*t1139+54*t52+80*t1141-75*l1-3
     #0*l5
      t2418 = 96*t1145*t2371+440*t1146+115*t3+234*t1150-152*t1153+54*t19
     #11-108*t1913-135*t1915+30*t210*t52+252*t255*t52-252*t244*t52-12*t2
     #10*t39-72*t255*t39+72*t244*t39-480*t1155
      t2435 = -240*t235*t2338-62*t1157+15*t253+360*t256-45*t259-330*t261
     #+54*t1926-180*t1165+90*t263+220*t1168+480*t1170+378*t1932+10*t1172
     #-t1175
      t2452 = -90*t1176+120*t1178+330*t1180-120*t1182-360*t265-35*t8-240
     #*t1186+15*t177+240*t1189+150*t183+160*t1193+168*t1082*t2334-80*t11
     #95-336*t1145*t2338
      t2465 = t39*l6
      t2468 = t39*t179
      t2481 = 24*l5*t39+54*t177*t52-108*l5*t52-12*t177*t39-960*t1197+672
     #*t244*t2346-192*t244*t2465-192*t255*t2468+192*t255*t2465+192*t244*
     #t2468-48*t167*t2343+96*t235*t2371+258*t1199-35*t270
      t2498 = -138*t1202+620*t1204-76*t1206-600*t1208-80*t1210+120*t1212
     #+60*t271+48*t1068*t2343+386*t1215-610*t1217-230*t198-243*t1963+115
     #*t273-75*t178+70*t1222
      s1 = 65920.D0/3.D0*t96-64640.D0/3.D0*t98-5120.D0/3.D0*t100-73600.D
     #0/3.D0*t102+147136.D0/3.D0*t104-162752.D0/3.D0*t106+10240.D0/3.D0*
     #t108-83840.D0/3.D0*t110-609408*t112
      s2 = s1+t280*((t1487+t1579+t1598+t1617)*s1l1-1664*t400+950912.D0/3
     #.D0*t59+(t1410+t1428+t1445+t1466)*t695-56064*t373+258112.D0/3.D0*t
     #96-174656.D0/3.D0*t98-8512*t100-78912*t102-t423-121280*t19-459584.
     #D0/3.D0*t376+1839296*t33+366016.D0/3.D0*t378-52800*t326+1054336.D0
     #/3.D0*t104-1348096.D0/3.D0*t106+44800.D0/3.D0*t94+2944.D0/3.D0*t21
     #-7302016.D0/3.D0*t14-7110656.D0/3.D0*t71+2871296.D0/3.D0*t73-28415
     #36.D0/3.D0*t61+4128704.D0/3.D0*t63+13568*t16+t1818-3072*t89+896*t4
     #+t1388-3121600.D0/3.D0*t65-36608*t67+t1701+t1724-2496*t91)/2+21888
     #*t114-25984*t116-253440*t411
      t2507 = s2+(t1831+t1842+t1847+t1857)*t695+190080*t414-t162*(-64.D0
     #/3.D0*t1059*(t1066+170*t3-160*l1-66*t8+39)*s2l1-64.D0/3.D0*t1059*(
     #26*t3-74*l1+6*t8+t1066+27)*s1l1-64.D0/3.D0*(t1893+720*t1119-720*t1
     #121+128*t1124-128*t1127+128*t1129+560*t1132-128*t1134-560*t1137-56
     #0*t1139+72*t1141+t1919+t1944-140*t1189-86*t183+144*t1193+144*t1195
     #+560*t1197+188*t1199+63*t270-292*t1202-370*t1204+41*t1206+t1968)*t
     #288)-t1034*(t1996+t2065+t2194+t2311)/6+t280*(-64.D0/3.D0*t1059*(-1
     #06*t8+134*l1+6*t4+6*t3+t2320-42)*s1l1-64.D0/3.D0*t1059*(78*t4-104*
     #t8-118*t3-t2320+244*l1-80)*s2l1-64.D0/3.D0*(t2360+t2379+t2395+t241
     #8+t2435+t2452+t2481+t2498)*t288)/2
      t2511 = (l1+z-t285)**2
      taux9 = (t814+t1024+t1359+t2507)/t2511*t288
      t1 = z**2
      t2 = t1**2
      t3 = l1**2
      t4 = t3**2
      t5 = t2*t4
      t7 = t4**2
      t8 = t7*l1
      t9 = t1*z
      t10 = t8*t9
      t12 = dlog(z)
      t13 = t12**2
      t16 = 1280*t1
      t17 = t3*l1
      t19 = t2*t1
      t20 = t19*l1
      t27 = t7*t3
      t28 = t2*z
      t29 = t27*t28
      t32 = t7*t17
      t33 = 13824*t32
      t34 = t4*l1
      t35 = t9*t34
      t37 = t4*t3
      t39 = -1280-2280960*t5-8180736*t10-3072*t2-1536*t1-2560.D0/3.D0*t1
     #7+1679616*t29+248320*t8+t33-17331968.D0/3.D0*t35-183808*t37
      t40 = z*l1
      t43 = t1*t34
      t45 = t2*t9
      t46 = t32*t45
      t47 = 13824*t46
      t48 = t9*t7
      t50 = t9*t4
      t52 = t9*l1
      t54 = z*t3
      t56 = t1*t3
      t58 = t1*l1
      t60 = z*t17
      t62 = -14336.D0/3.D0*t40-12544*t4+8102272.D0/3.D0*t43-t47+41217152
     #.D0/3.D0*t48+5756032.D0/3.D0*t50-16640*t52+42752.D0/3.D0*t54-35584
     #.D0/3.D0*t56+3328.D0/3.D0*t58-76288.D0/3.D0*t60
      t64 = t37*t1
      t66 = t27*t2
      t68 = t4*t17
      t69 = t1*t68
      t71 = t2*t7
      t73 = t28*t68
      t75 = t28*t7
      t77 = t45*t34
      t79 = t32*z
      t80 = 96768*t79
      t81 = t32*t1
      t82 = 290304*t81
      t83 = t28*t3
      t85 = t27*t1
      t87 = -18624512.D0/3.D0*t64-2868480*t66+8963712*t69-12982016*t71-7
     #591808*t73+7036800*t75-11520*t77-t80+t82+31232*t83-1804032*t85
      t88 = t8*z
      t90 = t8*t28
      t92 = t8*t2
      t94 = t8*t19
      t96 = t8*t1
      t98 = t19*t7
      t101 = t19*t17
      t103 = z*t4
      t105 = t45*t7
      t107 = t45*t68
      t109 = -1723904*t88-4411904*t90+7815680*t92+1358336*t94+5069312*t9
     #6-5970688.D0/3.D0*t98+173056.D0/3.D0*t34+12288*t101+371200.D0/3.D0
     #*t103+217728*t105-156672*t107
      t112 = t9*t3
      t114 = t1*t17
      t116 = t2*t17
      t118 = t9*t17
      t120 = t2*t3
      t122 = t8*t45
      t124 = t28*t17
      t126 = t28*t4
      t128 = t28*l1
      t130 = t28*t34
      t132 = z*t34
      t134 = 75776*t112+92672*t114+625024*t116-1270912.D0/3.D0*t118-1326
     #08*t120-175104*t122-213504*t124+3000448.D0/3.D0*t126-6144*t128-293
     #5168*t130-1826176.D0/3.D0*t132
      t135 = t9*t68
      t137 = z*t68
      t139 = z*t7
      t141 = t2*t68
      t143 = t45*t37
      t145 = t1*t7
      t147 = t32*t9
      t148 = 483840*t147
      t149 = t32*t2
      t150 = 483840*t149
      t151 = t32*t28
      t152 = 290304*t151
      t153 = t19*t68
      t156 = -15250432*t135-8250752.D0/3.D0*t137+2768000*t139+43560704.D
     #0/3.D0*t141+63360*t143-8411904*t145-t148+t150-t152+1924480*t153-89
     #856*t27
      t158 = t27*t9
      t160 = t32*t19
      t161 = 96768*t160
      t162 = t27*t19
      t164 = 1280*z
      t165 = t1*t4
      t167 = t2*l1
      t169 = t19*t4
      t171 = t19*t3
      t172 = 3072*t171
      t173 = t37*t2
      t175 = t19*t37
      t177 = t37*z
      t179 = 2937600*t158+t161-546048*t162-t164-680704*t165+19456*t167-1
     #10592*t169-t172-33782528.D0/3.D0*t173-3766016.D0/3.D0*t175+5068928
     #.D0/3.D0*t177
      t180 = t9*t37
      t182 = t37*t28
      t184 = t19*t34
      t186 = t2*t34
      t188 = 12800.D0/3.D0*l1
      t190 = z*t27
      t195 = t27*t45
      t197 = 34363520.D0/3.D0*t180+17101952.D0/3.D0*t182+513408*t184+181
     #96480.D0/3.D0*t186-t188-1132288.D0/3.D0*t7+615168*t190+340736*t68+
     #512*t9+27392.D0/3.D0*t3+76032*t195
      t204 = 12800.D0/3.D0*t40
      t208 = t45*t4
      t213 = 1280.D0/3.D0*t56
      t214 = -1280-253952.D0/3.D0*t5-228864*t10-t16-6400.D0/3.D0*t17-307
     #2*t20+(t39+t62+t87+t109+t134+t156+t179+t197)*s1l1-6150784.D0/3.D0*
     #t35-1536*t28-t204-3072*t4+1236096*t43+1137920*t48-135424.D0/3.D0*t
     #208+2796032.D0/3.D0*t50-24320.D0/3.D0*t52+32000.D0/3.D0*t54+t213
      t230 = t45*t17
      t232 = t45*t3
      t235 = -5120.D0/3.D0*t58-17920.D0/3.D0*t60-5792000.D0/3.D0*t64+173
     #9648*t69-711680.D0/3.D0*t71+1705728*t73-784640*t75+182400*t77-6988
     #8*t83-40960*t88+155648*t90+49152*t92-139776*t94+168448*t96+701440*
     #t98+7424*t230-1536*t232+1024*t34
      t252 = 512*t2
      t258 = -1280-6307328.D0/3.D0*t5-7800576*t10-t252-t16-5632.D0/3.D0*
     #t17+1804032*t29+172288*t8+t33-18492800.D0/3.D0*t35-144640.D0/3.D0*
     #t37
      t268 = -t204-4096*t4+8533120.D0/3.D0*t43-t47+38845568.D0/3.D0*t48+
     #11520*t208+6868864.D0/3.D0*t50-49664.D0/3.D0*t52+38144.D0/3.D0*t54
     #-40192.D0/3.D0*t56+4864.D0/3.D0*t58
      t279 = -38656.D0/3.D0*t60-16651264.D0/3.D0*t64-2937600*t66+2225292
     #8.D0/3.D0*t69-13855744*t71-28082432.D0/3.D0*t73+8571264*t75-74880*
     #t77-t80+t82+11008*t83
      t291 = -1679616*t85-1343744*t88-5096192*t90+8195840*t92+1738496*t9
     #4+4385024*t96-8548864.D0/3.D0*t98+9472*t34+63232.D0/3.D0*t101+9420
     #8*t103+392832*t105
      t303 = 1024*t128
      t305 = -374400*t107+115712*t112+543232.D0/3.D0*t114+1331584.D0/3.D
     #0*t116-590336*t118-161792.D0/3.D0*t120-251136*t122-349312.D0/3.D0*
     #t124+2536064.D0/3.D0*t126-t303-3108864*t130
      t314 = -465664*t132-43607168.D0/3.D0*t135-1814784*t137+1926784*t13
     #9+15592704*t141+220032*t143-6928896*t145-t148+t150-t152+2939776*t1
     #53
      t322 = 512*t171
      t325 = -76032*t27+2868480*t158+t161-615168*t162-t164-2840320.D0/3.
     #D0*t165+6400*t167-172800*t169-t322-11967488*t173-5764864.D0/3.D0*t
     #175
      t335 = 28160.D0/3.D0*t3
      t337 = 1159296*t177+34099328.D0/3.D0*t180+20225920.D0/3.D0*t182+22
     #90816.D0/3.D0*t184+6200832*t186-t188-615424.D0/3.D0*t7+546048*t190
     #+406784.D0/3.D0*t68+1536*t9+t335+89856*t195
      t346 = -322816.D0/3.D0*t101+34048*t103-546560.D0/3.D0*t105+1090688
     #.D0/3.D0*t107+156416.D0/3.D0*t112+77824*t114-27776.D0/3.D0*t116-26
     #9568*t118+52480.D0/3.D0*t120+36352*t122+303744*t124-2570752.D0/3.D
     #0*t126+10496*t128+(t258+t268+t279+t291+t305+t314+t325+t337)*s2l1+4
     #934912.D0/3.D0*t130-555008.D0/3.D0*t132-7239296.D0/3.D0*t135-12121
     #60.D0/3.D0*t137
      t363 = 613120.D0/3.D0*t139+1441280.D0/3.D0*t141-1082624.D0/3.D0*t1
     #43-839680*t145-4416640.D0/3.D0*t153-t164-1348352.D0/3.D0*t165-5888
     #*t167+1373440.D0/3.D0*t169+17920*t171-517632*t173+1677568*t175+117
     #6320.D0/3.D0*t177+2844416*t180-6314752.D0/3.D0*t182-3396224.D0/3.D
     #0*t184+307968*t186-t188+t335
      t369 = 1024*t10
      t370 = 2560*t2
      t371 = 1024*t1
      t372 = 13568.D0/3.D0*t17
      t373 = 2816*t8
      t376 = 10240.D0/3.D0*t40
      t377 = 24448.D0/3.D0*t4
      t381 = 25472*t52
      t382 = 53248.D0/3.D0*t54
      t383 = 1024+36352*t5-t369+t370+t371+t372-t373-19456.D0/3.D0*t35+88
     #576.D0/3.D0*t37+t376+t377-418432.D0/3.D0*t43+82432.D0/3.D0*t48-389
     #632.D0/3.D0*t50+t381-t382
      t384 = 118912.D0/3.D0*t56
      t385 = 18944.D0/3.D0*t58
      t386 = 46464*t60
      t393 = 10240*t88
      t394 = 9216*t90
      t395 = 12032*t92
      t396 = 2304*t94
      t397 = 11520*t96
      t401 = t384-t385+t386+280448.D0/3.D0*t64-66048*t69-76416*t71-12633
     #6*t73+52480*t75-66304.D0/3.D0*t83+t393-t394+t395+t396-t397-37760.D
     #0/3.D0*t98-21888*t34-21632.D0/3.D0*t101
      t411 = 5120*t128
      t420 = -257920.D0/3.D0*t103-259456.D0/3.D0*t112-328448.D0/3.D0*t11
     #4-100096*t116+459392.D0/3.D0*t118+60928*t120+105344.D0/3.D0*t124+3
     #456*t126+t411-295424.D0/3.D0*t130+112640*t132-279040.D0/3.D0*t135+
     #259712.D0/3.D0*t137-43776*t139+585088.D0/3.D0*t141+40064*t145+8704
     #0.D0/3.D0*t153
      t421 = 1024*z
      t423 = 52352.D0/3.D0*t167
      t425 = 2560*t171
      t433 = 10240.D0/3.D0*l1
      t436 = 2560*t9
      t437 = 27136.D0/3.D0*t3
      t438 = t421+476288.D0/3.D0*t165-t423+3584.D0/3.D0*t169+t425-748288
     #.D0/3.D0*t173-105088.D0/3.D0*t175-337280.D0/3.D0*t177+113920*t180+
     #479872.D0/3.D0*t182+59264.D0/3.D0*t184+409472.D0/3.D0*t186+t433+38
     #272.D0/3.D0*t7-25216*t68-t436-t437
      t445 = 2048*t1
      t454 = 50688*t32
      t457 = 2048+14239744.D0/3.D0*t5+29225472*t10+6656*t2+3456*t1-15232
     #.D0/3.D0*t17-6158592*t29-918272*t8-t454+46508288.D0/3.D0*t35+78169
     #6*t37
      t461 = 50688*t46
      t469 = 28928.D0/3.D0*t40+178304.D0/3.D0*t4-8710400*t43+t461-466208
     #00*t48-14317952.D0/3.D0*t50+45824*t52-111232.D0/3.D0*t54+78080*t56
     #-39808.D0/3.D0*t58+368128.D0/3.D0*t60
      t478 = 354816*t79
      t479 = 1064448*t81
      t482 = 20048512*t64+10517760*t66-29895552*t69+43573120*t71+2323302
     #4*t73-23742080*t75+49920*t77+t478-t479-67712*t83+6614784*t85
      t494 = 6272512*t88+15792640*t90-27869440*t92-4915456*t94-18234112*
     #t96+6921344*t98-273920*t34-25728*t101-1847680.D0/3.D0*t103-819072*
     #t105+609792*t107
      t508 = -212992*t112-435712*t114-1154304*t116+1038080*t118+246656*t
     #120+646656*t122+391424*t124-5969152.D0/3.D0*t126+13312*t128+204152
     #32.D0/3.D0*t130+2474240*t132
      t515 = 1774080*t147
      t516 = 1774080*t149
      t517 = 1064448*t151
      t520 = 144670720.D0/3.D0*t135+29630464.D0/3.D0*t137-9963904*t139-1
     #33551232.D0/3.D0*t141-259200*t143+29224832*t145+t515-t516+t517-185
     #58848.D0/3.D0*t153+329472*t27
      t523 = 354816*t160
      t525 = 2048*z
      t533 = -10771200*t158-t523+2002176*t162+t525+7153024.D0/3.D0*t165-
     #44928*t167+243712*t169+6656*t171+92506496.D0/3.D0*t173+10778240.D0
     #/3.D0*t175-18719360.D0/3.D0*t177
      t538 = 20480.D0/3.D0*l1
      t545 = -100259456.D0/3.D0*t180-15339648*t182-1284352*t184-14572544
     #*t186+t538+1426560*t7-2255616*t190-4032896.D0/3.D0*t68-2944*t9-398
     #08.D0/3.D0*t3-278784*t195
      t553 = 20480.D0/3.D0*t40
      t561 = 2048+966400.D0/3.D0*t5+383488*t10-384*t2+t445+6656.D0/3.D0*
     #t17+(t457+t469+t482+t494+t508+t520+t533+t545)*s1l1+6656*t20+954252
     #8.D0/3.D0*t35+3328*t28+t553+6144*t4-5390848.D0/3.D0*t43-1908480*t4
     #8+198656.D0/3.D0*t208-3993088.D0/3.D0*t50+30464.D0/3.D0*t52-49792.
     #D0/3.D0*t54
      t580 = 2048+16072832.D0/3.D0*t5+27831552*t10+2176*t2+t445+2432.D0/
     #3.D0*t17-6614784*t29-639488*t8-t454+45638528.D0/3.D0*t35+673408.D0
     #/3.D0*t37
      t590 = t553+12032*t4-20461952.D0/3.D0*t43+t461-43525760*t48-49920*
     #t208-5009536*t50+104192.D0/3.D0*t52-68224.D0/3.D0*t54+116864.D0/3.
     #D0*t56-15232.D0/3.D0*t58
      t601 = 101888.D0/3.D0*t60+45543296.D0/3.D0*t64+10771200*t66-686904
     #32.D0/3.D0*t69+46949760*t71+92898176.D0/3.D0*t73-29643392*t75+3091
     #20*t77+t478-t479-46464*t83
      t613 = 6158592*t85+4878592*t88+18301696*t90-29263360*t92-6309376*t
     #94-15725056*t96+10173568*t98-148352.D0/3.D0*t34-261248.D0/3.D0*t10
     #1-241024*t103-1465728*t105
      t627 = 1428864*t107-630016.D0/3.D0*t112-364928*t114-3317888.D0/3.D
     #0*t116+1145728*t118+461312.D0/3.D0*t120+925440*t122+449664*t124-27
     #71712*t126+4352*t128+9783168*t130
      t636 = 3686272.D0/3.D0*t132+134215040.D0/3.D0*t135+17800832.D0/3.D
     #0*t137-6764928*t139-147740672.D0/3.D0*t141-868992*t143+23489408*t1
     #45+t515-t516+t517-10365952*t153
      t647 = 278784*t27-10517760*t158-t523+2255616*t162+t525+6020096.D0/
     #3.D0*t165-19328*t167+2086016.D0/3.D0*t169+2176*t171+35004032*t173+
     #20481664.D0/3.D0*t175
      t657 = 44032.D0/3.D0*t3
      t659 = -3407232*t177-94402688.D0/3.D0*t180-64479104.D0/3.D0*t182-8
     #510080.D0/3.D0*t184-50508928.D0/3.D0*t186+t538+787072*t7-2002176*t
     #190-1671680.D0/3.D0*t68-3200*t9-t657-329472*t195
      t672 = -13696.D0/3.D0*t56+8960.D0/3.D0*t58+19712.D0/3.D0*t60+29030
     #40*t64-7980032.D0/3.D0*t69+2360320.D0/3.D0*t71-4843648.D0/3.D0*t73
     #+733440*t75-710528.D0/3.D0*t77+117376*t83+(t580+t590+t601+t613+t62
     #7+t636+t647+t659)*s2l1+60928*t88-144896*t90-159744*t92+164864*t94-
     #259072*t96-827904*t98-15104*t230+3328*t232
      t693 = -2048*t34+495104.D0/3.D0*t101-130432.D0/3.D0*t103+685312.D0
     #/3.D0*t105-456576*t107-56320*t112-283136.D0/3.D0*t114+41984.D0/3.D
     #0*t116+333952*t118-143872.D0/3.D0*t120-45568*t122-1228160.D0/3.D0*
     #t124+2947456.D0/3.D0*t126-22528*t128-1698176*t130+798848.D0/3.D0*t
     #132+12034688.D0/3.D0*t135+1799296.D0/3.D0*t137-912128.D0/3.D0*t139
      t710 = -4859392.D0/3.D0*t141+1368320.D0/3.D0*t143+1291776*t145+521
     #8816.D0/3.D0*t153+t525+1846912.D0/3.D0*t165+15360*t167-587392*t169
     #-37248*t171+5337088.D0/3.D0*t173-5972992.D0/3.D0*t175-577792*t177-
     #13859072.D0/3.D0*t180+2050304*t182+4102144.D0/3.D0*t184-3272192.D0
     #/3.D0*t186+t538+128*t9-t657
      t717 = 36*t4
      t719 = 12*t34
      t723 = 106*t3
      t724 = 40*l1
      t725 = 83*t17
      t726 = t13*t12
      t742 = 2*l6-1
      t743 = t742**2
      t746 = l5-1
      t748 = -1+l1
      t750 = (1-l1+t40)**2
      t751 = 1/t750
      t756 = dsqrt(z)
      t758 = dsqrt(l5)
      t760 = dsqrt(-t746)
      t763 = dsqrt(l1)
      t767 = l5*z*l1-l1+1-l5+l5*l1
      t770 = dsqrt(-t748)
      t775 = t767**2
      t778 = dsqrt(4*l1*t743*z*l5*t746*t748*t751-4*t742*t756*t758*t760*t
     #763*t767*t751*t770+t775*t751)
      t786 = 36864*t32
      t789 = -1024-12591488.D0/3.D0*t5-20465664*t10-t370-t371-t372+48107
     #52*t29+490496*t8+t786-35483392.D0/3.D0*t35-958336.D0/3.D0*t37
      t791 = 36864*t46
      t795 = -t376-t377+5645568*t43-t791+98058752.D0/3.D0*t48+61440*t208
     #+11176192.D0/3.D0*t50-t381+t382-t384+t385
      t804 = 258048*t79
      t805 = 774144*t81
      t807 = -t386-12566784*t64-7833600*t66+55110272.D0/3.D0*t69-3500313
     #6*t71-71995904.D0/3.D0*t73+22216448*t75-337920*t77-t804+t805+15078
     #4.D0/3.D0*t83
      t819 = -4478976*t85-3664384*t88-13430272*t90+21452800*t92+4652032*
     #t94+11653120*t96-23256832.D0/3.D0*t98+87168*t34+280192.D0/3.D0*t10
     #1+817280.D0/3.D0*t103+1148928*t105
      t832 = -1222656*t107+135552*t112+839168.D0/3.D0*t114+2512768.D0/3.
     #D0*t116-2339072.D0/3.D0*t118-115968*t120-688128*t122-439168*t124+7
     #445120.D0/3.D0*t126-t411-24872320.D0/3.D0*t130
      t840 = 1290240*t147
      t841 = 1290240*t149
      t842 = 774144*t151
      t844 = -3775360.D0/3.D0*t132-103308544.D0/3.D0*t135-15686912.D0/3.
     #D0*t137+5438208*t139+112702592.D0/3.D0*t141+835584*t143-18047744*t
     #145-t840+t841-t842+25040512.D0/3.D0*t153
      t848 = 258048*t160
      t854 = -202752*t27+7649280*t158+t848-1640448*t162-t421-4779008.D0/
     #3.D0*t165+t423-716288*t169-t425-81891968.D0/3.D0*t173-5928192*t175
      t864 = 9827456.D0/3.D0*t177+73864192.D0/3.D0*t180+17378944*t182+80
     #68352.D0/3.D0*t184+39864448.D0/3.D0*t186-t433-2060032.D0/3.D0*t7+1
     #456128*t190+601984*t68+t436+t437+239616*t195
      t872 = 6144*t2
      t873 = 3584*t1
      t874 = 20096.D0/3.D0*t17
      t877 = 25600.D0/3.D0*t40
      t882 = 173440.D0/3.D0*t52
      t883 = 48384*t54
      t884 = 1024+1318016.D0/3.D0*t5-t369+t872+t873-t874-t373+1767424.D0
     #/3.D0*t35+164224.D0/3.D0*t37+t877+93184.D0/3.D0*t4-1778048.D0/3.D0
     #*t43-61952.D0/3.D0*t48-1997824.D0/3.D0*t50+t882-t883
      t886 = 28800*t58
      t893 = 60928*t83
      t896 = 21504*t101
      t897 = 406016.D0/3.D0*t56-t886+133760*t60+1342976.D0/3.D0*t64-2366
     #72*t69-41344*t71-63616*t73+38656*t75-t893+t393-t394+t395+t396-t397
     #-30848.D0/3.D0*t98-154240.D0/3.D0*t34-t896
      t907 = 12288*t128
      t916 = -238976*t103-671360.D0/3.D0*t112-1052416.D0/3.D0*t114-10726
     #40.D0/3.D0*t116+1464704.D0/3.D0*t118+508288.D0/3.D0*t120+407168.D0
     #/3.D0*t124-468608.D0/3.D0*t126+t907+224000.D0/3.D0*t130+883456.D0/
     #3.D0*t132+375808.D0/3.D0*t135+475520.D0/3.D0*t137-59648*t139+10777
     #6.D0/3.D0*t141+77696*t145+18560*t153
      t918 = 45056*t167
      t920 = 6144*t171
      t930 = 5632*t9
      t931 = 19456.D0/3.D0*t3
      t932 = t421+1669888.D0/3.D0*t165-t918+80000.D0/3.D0*t169+t920+2543
     #36.D0/3.D0*t173-40192.D0/3.D0*t175-777344.D0/3.D0*t177-343296*t180
     #+85888.D0/3.D0*t182-25472.D0/3.D0*t184-909440.D0/3.D0*t186+t433+46
     #720.D0/3.D0*t7-113920.D0/3.D0*t68-t930-t931
      t943 = -1024-3702016*t5-21479424*t10-t872-t873+t874+4478976*t29+69
     #3248*t8+t786-40774144.D0/3.D0*t35-2764288.D0/3.D0*t37
      t950 = -t877-289024.D0/3.D0*t4+26797952.D0/3.D0*t43-t791+105348608
     #.D0/3.D0*t48+4546048*t50-t882+t883-463616.D0/3.D0*t56+t886-664960.
     #D0/3.D0*t60
      t960 = -54159488.D0/3.D0*t64-7649280*t66+73541504.D0/3.D0*t69-3253
     #2224*t71-53535488.D0/3.D0*t73+17807616*t75-61440*t77-t804+t805+t89
     #3-4810752*t85
      t971 = -4678144*t88-11605504*t90+20439040*t92+3638272*t94+13477888
     #*t96-15942400.D0/3.D0*t98+1219840.D0/3.D0*t34+t896+2903168.D0/3.D0
     #*t103+663552*t105-559104*t107
      t984 = 264320*t112+2141056.D0/3.D0*t114+2677760.D0/3.D0*t116-33751
     #04.D0/3.D0*t118-209536*t120-485376*t122-852608.D0/3.D0*t124+434316
     #8.D0/3.D0*t126-t907-5211008*t130-9183104.D0/3.D0*t132
      t993 = -37769472*t135-25980416.D0/3.D0*t137+7890176*t139+102158720
     #.D0/3.D0*t141+276480*t143-22454016*t145-t840+t841-t842+14935424.D0
     #/3.D0*t153-239616*t27
      t1002 = 7833600*t158+t848-1456128*t162-t421-2965376*t165+t918-1964
     #80*t169-t920-72364544.D0/3.D0*t173-8995712.D0/3.D0*t175+19024256.D
     #0/3.D0*t177
      t1011 = 27542528*t180+35802752.D0/3.D0*t182+3225472.D0/3.D0*t184+3
     #4517504.D0/3.D0*t186-t433-3531520.D0/3.D0*t7+1640448*t190+3865984.
     #D0/3.D0*t68+t930+t931+202752*t195
      t1016 = -709376.D0/3.D0*t5-124928*t10+t13*(t214+t235+t346+t363)/2+
     #(t383+t401+t420+t438)*s2l1-t12*(t561+t672+t693+t710)-3352832.D0/3.
     #D0*t35-256.D0/3.D0*(t13*(-15-50*l1+110*t3-t717-25*t17+t719)/2-t719
     #+t717+t723-t724-t725-12-t726*(3+10*l1+12*t17-24*t3)/6-t12*(24+26*t
     #17+80*l1+72*t4-24*t34-172*t3))*t778+1834112.D0/3.D0*t43+(t789+t795
     #+t807+t819+t832+t844+t854+t864)*s2l1+624640*t48-33280.D0/3.D0*t208
     #+(t884+t897+t916+t932)*s1l1+1475840.D0/3.D0*t50+(t943+t950+t960+t9
     #71+t984+t993+t1002+t1011)*s1l1
      t1032 = 256*t1
      t1033 = 1024*t17
      t1036 = 81920*t10
      t1038 = 2560.D0/3.D0*t40
      t1043 = 2304*t54
      t1044 = 1024.D0/3.D0*t58
      t1045 = 256+1092352.D0/3.D0*t5+t1036+t252+t1032+t1033+616448*t35+t
     #1038-555008.D0/3.D0*t43-401408*t48-746752.D0/3.D0*t50+2304*t52-t10
     #43-t213+t1044
      t1053 = 8192*t88
      t1054 = 40960*t90
      t1055 = 81920*t92
      t1056 = 8192*t94
      t1057 = 40960*t96
      t1061 = 4352.D0/3.D0*t60+1076224.D0/3.D0*t64-367104*t69+417792*t71
     #+479232*t73-217088*t75-5120*t83+t1053+t1054-t1055-t1056-t1057+4505
     #6*t98-2048*t101+1792.D0/3.D0*t103
      t1077 = -27392.D0/3.D0*t112-1536*t114-345088.D0/3.D0*t116+178688.D
     #0/3.D0*t118+25344*t120+39680*t124-503296.D0/3.D0*t126+t303+1199360
     #.D0/3.D0*t130+35584.D0/3.D0*t132+826368*t135+62464*t137-36864*t139
     #-901120*t141+192512*t145
      t1079 = 256*z
      t1090 = 2560.D0/3.D0*l1
      t1091 = 2048*t3
      t1092 = -99840*t153+t1079+43264*t165-3072*t167+16896*t169+t322+321
     #7408.D0/3.D0*t173+111616*t175-139264.D0/3.D0*t177-925696*t180-1712
     #128.D0/3.D0*t182-64000*t184-2337280.D0/3.D0*t186+t1090-t1091
      t1104 = 256-42752.D0/3.D0*t5+40960*t10+t1032+t1033+512*t20+(t1045+
     #t1061+t1077+t1092)*s1l1+1166336.D0/3.D0*t35+256*t28+t1038-247808*t
     #43-203520*t48+9472*t208-557312.D0/3.D0*t50+4864.D0/3.D0*t52-t1043-
     #256*t56
      t1120 = t1044+1792*t60+1141760.D0/3.D0*t64-1018880.D0/3.D0*t69-512
     #0.D0/3.D0*t71-1338112.D0/3.D0*t73+206080*t75-40448*t77+14336*t83+t
     #1053-t1054+32768*t94-32768*t96-164352*t98-1280*t230+256*t232+22528
     #*t101-7424*t103
      t1139 = 123136.D0/3.D0*t105-245504.D0/3.D0*t107-34304.D0/3.D0*t112
     #-47104.D0/3.D0*t114+8192*t116+57344*t118-11008.D0/3.D0*t120-8192*t
     #122-209408.D0/3.D0*t124+210176*t126-1792*t128-418304*t130+111616.D
     #0/3.D0*t132+1304320.D0/3.D0*t135+242432.D0/3.D0*t137-122624.D0/3.D
     #0*t139+7168*t141
      t1162 = 256+240896*t5+t1036+t1032+t1033+2277376.D0/3.D0*t35+t1038-
     #383488*t43-417792*t48-1028864.D0/3.D0*t50+7936.D0/3.D0*t52-2560*t5
     #4+1280*t56-512.D0/3.D0*t58
      t1172 = 2304*t60+561152*t64-1430528.D0/3.D0*t69+401408*t71+1103360
     #.D0/3.D0*t73-192512*t75-256*t83+t1053+t1054-t1055-t1056-t1057+3686
     #4*t98-41984.D0/3.D0*t103
      t1188 = -64768.D0/3.D0*t112-92672.D0/3.D0*t114-156160.D0/3.D0*t116
     #+306176.D0/3.D0*t118+12032.D0/3.D0*t120+4864.D0/3.D0*t124-44800*t1
     #26+188672*t130+179456.D0/3.D0*t132+2695168.D0/3.D0*t135+297472.D0/
     #3.D0*t137-45056*t139-2477056.D0/3.D0*t141+217088*t145
      t1200 = -188416.D0/3.D0*t153+t1079+457984.D0/3.D0*t165-512*t167+27
     #68896.D0/3.D0*t173+143360.D0/3.D0*t175-326656.D0/3.D0*t177-3184640
     #.D0/3.D0*t180-361472*t182-40960.D0/3.D0*t184-1831424.D0/3.D0*t186+
     #t1090-256*t9-t1091
      t1204 = 242944.D0/3.D0*t143+163328*t145+1034240.D0/3.D0*t153+t1079
     #+275200.D0/3.D0*t165+1024*t167-103936*t169-t172-14336*t173-392192*
     #t175-235264.D0/3.D0*t177-521984*t180+1636096.D0/3.D0*t182+788992.D
     #0/3.D0*t184+52736.D0/3.D0*t186+(t1162+t1172+t1188+t1200)*s2l1+t109
     #0-t1091
      t1210 = -968192*t64+875648*t69-271360*t71+438912*t73-207360*t75+18
     #3040.D0/3.D0*t77-17920.D0/3.D0*t83-19968*t88+41472*t90+54272*t92-4
     #8640*t94+84480*t96+243200*t98-t726*(t1104+t1120+t1139+t1204)/6-512
     #00.D0/3.D0*t101
      t1227 = 16640*t103-66560*t105+396032.D0/3.D0*t107+44800.D0/3.D0*t1
     #12+94720.D0/3.D0*t114+62080*t116-379520.D0/3.D0*t118-21760.D0/3.D0
     #*t120+13312*t122+131840.D0/3.D0*t124-159744*t126+360832*t130-91520
     #*t132-3983872.D0/3.D0*t135-198016*t137
      t1245 = 99840*t139+1746688.D0/3.D0*t141-386048.D0/3.D0*t143-422400
     #*t145-1508480.D0/3.D0*t153-212480*t165+352000.D0/3.D0*t169-2102272
     #.D0/3.D0*t173+1656320.D0/3.D0*t175+193024*t177+4691968.D0/3.D0*t18
     #0-511488*t182-1032320.D0/3.D0*t184+1553920.D0/3.D0*t186+256.D0/3.D
     #0*(-t719+t717+t723-t724-t725-12)*t778
      t1249 = (l1+z-t40)**2
      taux10 = (t1016+t1210+t1227+t1245)/t1249*t751

CCC Coefficients

      tcoff15 = taux9
      tcoff16 = taux9
      tcoff17 = taux10
      tcoff18 = -taux10

CCC Jet functions

      tFjet15 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff15)
      tFjet16 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff16)
      tFjet17 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff17)
      tFjet18 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff18)

CCC Sum

      RVnum5e3 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet18
     #*tcoff18
      RETURN
      END

      double precision function RVnum5e4(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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
      t2 = t1**2
      t3 = t2*z
      t4 = l1**2
      t5 = t3*t4
      t7 = dlog(z)
      t8 = t7**2
      t9 = t8*t7
      t11 = 2*l6-1
      t12 = t11**2
      t15 = l5-1
      t17 = -1+l1
      t18 = z*l1
      t20 = (1-l1+t18)**2
      t21 = 1/t20
      t26 = dsqrt(z)
      t28 = dsqrt(l5)
      t30 = dsqrt(-t15)
      t33 = dsqrt(l1)
      t34 = z*l5
      t35 = t34*l1
      t36 = l5*l1
      t37 = t35-l1+1-l5+t36
      t40 = dsqrt(-t17)
      t45 = t37**2
      t48 = dsqrt(4*l1*t12*z*l5*t15*t17*t21-4*t11*t26*t28*t30*t33*t37*t2
     #1*t40+t45*t21)
      t49 = t4*l1
      t52 = t4**2
      t55 = t52*l1
      t56 = 12*t55
      t67 = t26*t28
      t68 = t67*t30
      t69 = t33*t52
      t71 = t69*t40*l6
      t75 = t28*l5
      t76 = t26*z*t75
      t77 = t76*t30
      t80 = l5**2
      t81 = z*t80
      t82 = t55*l6
      t87 = l6**2
      t88 = t55*t87
      t93 = t33*t55
      t95 = t93*t40*l6
      t98 = t52*t4
      t99 = t98*t87
      t102 = t98*l6
      t109 = t33*l1
      t111 = t109*t40*l6
      t112 = t77*t111
      t116 = 15+336*t68*t71-240*t77*t71-672*t81*t82+672*t34*t82-672*t34*
     #t88+672*t81*t88-96*t68*t95-192*t81*t99-192*t34*t102+192*t34*t99+19
     #2*t81*t102-120*t112+96*t77*t95
      t118 = t68*t111
      t122 = t48*t98
      t124 = t81*t4
      t127 = t34*t4
      t129 = t80*t1
      t130 = t129*t4
      t133 = t26*t75
      t134 = t133*t30
      t137 = l5*t4
      t139 = t48*l1
      t141 = t80*l1
      t143 = -75*l1+480*t118-35*t49+115*t4+54*t122+360*t124+90*t35-360*t
     #127+15*t130-62*t52-336*t134*t71-230*t137+386*t139-75*t141
      t145 = t48*z
      t146 = t145*t49
      t148 = t48*t1
      t149 = t148*t49
      t151 = t145*l1
      t153 = t145*t98
      t155 = t129*t52
      t157 = t33*t49
      t159 = t157*t40*l6
      t160 = t68*t159
      t162 = l5*t49
      t164 = t33*t4
      t166 = t164*t40*l6
      t167 = t68*t166
      t170 = t30*t93*t40
      t183 = -600*t146+234*t149-152*t151-108*t153+10*t155-160*t160+70*t1
     #62-440*t167-48*t133*t170+30*t129*t55+252*t81*t55-252*t34*t55-12*t1
     #29*t98-72*t81*t98
      t186 = t81*t49
      t189 = t30*t69*t40
      t200 = t30*t164*t40
      t201 = t76*t200
      t204 = t30*t109*t40
      t205 = t76*t204
      t207 = t133*t200
      t212 = t30*t157*t40
      t213 = t133*t212
      t218 = 72*t34*t98-330*t186-168*t67*t189+150*t36-12*t98+96*t134*t95
     #-76*t48-12*t80*t98-180*t201+60*t205-220*t207+168*t133*t189-80*t213
     #-30*l5-48*t76*t170
      t221 = t49*t87
      t222 = t81*t221
      t224 = t34*t221
      t226 = t52*t87
      t227 = t34*t226
      t229 = t49*l6
      t230 = t34*t229
      t232 = t52*l6
      t233 = t34*t232
      t235 = t81*t229
      t237 = t81*t232
      t239 = t81*t226
      t241 = t4*t87
      t242 = t34*t241
      t244 = l1*t87
      t245 = t34*t244
      t247 = t81*t244
      t249 = t4*l6
      t250 = t81*t249
      t252 = l1*l6
      t253 = t81*t252
      t255 = t81*t241
      t257 = -880*t222+880*t224+320*t227-880*t230-320*t233+880*t235+320*
     #t237-320*t239-960*t242+240*t245-240*t247-960*t250+240*t253+960*t25
     #5
      t258 = t34*t252
      t260 = t34*t249
      t262 = t80*t4
      t265 = t134*t159
      t268 = t30*t33*t40
      t269 = t67*t268
      t271 = t148*t4
      t273 = t34*t52
      t275 = t34*t49
      t277 = t81*t52
      t279 = t48*t49
      t281 = t80*t49
      t283 = t76*t212
      t285 = t145*t55
      t287 = -240*t258+960*t260+115*t262+54*t55+160*t265+60*t269-76*t271
     #+120*t273+330*t275-120*t277+231*t279-35*t281+40*t283+378*t285
      t289 = t145*t52
      t291 = t148*t98
      t293 = t67*t204
      t298 = t33*t40*l6
      t299 = t68*t298
      t303 = t48*t52
      t305 = t67*t200
      t307 = t48*t4
      t310 = t134*t111
      t318 = -138*t289+54*t291-240*t293+48*t67*t170-120*t299+120*t76*t18
     #9+258*t303+220*t305-610*t307+15*t80-480*t310+54*t80*t55-108*l5*t55
     #+24*l5*t98
      t319 = t133*t204
      t321 = t134*t298
      t323 = t80*t52
      t325 = l5*t52
      t327 = t48*t55
      t329 = t133*t268
      t331 = t77*t166
      t333 = t77*t159
      t335 = t148*t55
      t337 = t129*t49
      t339 = t67*t212
      t341 = t134*t166
      t343 = t148*t52
      t344 = 66*t343
      t345 = t81*l1
      t347 = t145*t4
      t349 = 240*t319+120*t321-62*t323+124*t325-243*t327-60*t329+360*t33
     #1-80*t333-135*t335-45*t337+80*t339+440*t341-t344-90*t345+620*t347
      t358 = t52*t49
      t359 = t3*t358
      t361 = t52**2
      t362 = t3*t361
      t364 = t2*t55
      t366 = t2*t49
      t368 = t1*z
      t369 = t368*t49
      t371 = t2*t4
      t372 = 10240.D0/3.D0*t371
      t373 = t3*t49
      t375 = t3*t52
      t377 = t2*t52
      t379 = -5120.D0/3.D0*t5-t9*(-64.D0/3.D0*t48*(-106*t49+134*l1+6*t52
     #+6*t4+t56-42)*s1l1-64.D0/3.D0*t48*(78*t52-104*t49-118*t4-t56+244*l
     #1-80)*s2l1-64.D0/3.D0*(t116+t143+t183+t218+t257+t287+t318+t349)*t2
     #1)/6+92160*t359-20480*t362-487424.D0/3.D0*t364-25600*t366+12800*t3
     #69+t372+58880.D0/3.D0*t373-77824*t375+262144.D0/3.D0*t377
      t380 = t368*t52
      t382 = t1*t52
      t384 = t2*t1
      t385 = t384*t52
      t387 = t98*t2
      t389 = t384*t98
      t391 = t368*t98
      t393 = t98*t3
      t395 = t384*t55
      t397 = t384*t49
      t399 = t98*t1
      t401 = t1*t55
      t403 = t368*t55
      t405 = -131072.D0/3.D0*t380+7680*t382+89600.D0/3.D0*t385+507904.D0
     #/3.D0*t387+244736.D0/3.D0*t389-253952.D0/3.D0*t391-501760.D0/3.D0*
     #t393-211456.D0/3.D0*t395-5120*t397+50176.D0/3.D0*t399-46592.D0/3.D
     #0*t401+243712.D0/3.D0*t403
      t407 = t26*t1
      t408 = t407*t75
      t409 = t408*t30
      t412 = t368*t87
      t418 = t368*l6
      t424 = t407*t28
      t428 = t80*t368
      t432 = l5*t1
      t434 = t80*t2
      t450 = 2*t432+t434*t49-4*t432*t4+2*t432*l1+2*t137*t368-6*t36*t368+
     #2*t428*t4-t1-t434*t4+6*t141*t368-t129*l1
      t460 = t26*t368*t75
      t466 = t1*l1
      t483 = 4*t408*t268-4*t424*t268-4*t186+8*t460*t30*t111+4*t76*t268-4
     #*t205-16*t409*t111-3*t262+t281-t80+6*t337
      t485 = -8*t409*t298-16*t36*t412+2*t34+16*t141*t412+16*t36*t418-16*
     #t141*t418+8*t112+8*t424*t30*t298-4*t428*t49-t129+t450-2*t81-8*t77*
     #t298+6*t124-4*t35+2*t127-4*t130-4*t460*t204+3*t141+8*t408*t204+t46
     #6+t483
      t494 = t8**2
      t497 = 12*t49
      t510 = 12*t52
      t519 = -7+56*t112+32*l1-200*t118+42*t49-55*t4-150*t124-42*t35+150*
     #t127-7*t130-t510+110*t137-176*t139+32*t141+332*t146-100*t149+76*t1
     #51-12*t155-96*t160
      t539 = -84*t162+240*t167+180*t186-64*t36+38*t48+72*t201-28*t205+12
     #0*t207-48*t213+14*l5+480*t222-480*t224+192*t227+480*t230-192*t233-
     #480*t235+192*t237-192*t239+400*t242
      t560 = -112*t245+112*t247+400*t250-112*t253-400*t255+112*t258-400*
     #t260-55*t262+96*t265-28*t269+38*t271+72*t273-180*t275-72*t277-232*
     #t279+42*t281-48*t283-132*t289+100*t293
      t579 = 56*t299+66*t303-120*t305+304*t307-7*t80+200*t310-100*t319-5
     #6*t321-12*t323+24*t325+28*t329-144*t331+96*t333+18*t337+48*t339-24
     #0*t341+t344+42*t345-276*t347
      t592 = t368*t358
      t595 = t361*l1
      t596 = t595*z
      t598 = t595*t3
      t600 = t595*t2
      t602 = t595*t384
      t606 = t595*t368
      t608 = t361*t4
      t609 = z*t608
      t611 = t608*t368
      t613 = t608*t2
      t615 = t608*t1
      t619 = 2432*t362
      t621 = 64*t2
      t623 = z*t55
      t630 = 2560.D0/3.D0*t5-10176*t359+t619+1408.D0/3.D0*t358-t621+1569
     #28.D0/3.D0*t364-35584.D0/3.D0*t623+5056*t366-46720.D0/3.D0*t369-96
     #0*t371-832*t373+6080.D0/3.D0*t375
      t643 = -66304.D0/3.D0*t377+45248*t380-25984*t382-1216.D0/3.D0*t385
     #-197056.D0/3.D0*t387-128*t389+305408.D0/3.D0*t391+46400.D0/3.D0*t3
     #93+320*t395+256*t397-205184.D0/3.D0*t399+169280.D0/3.D0*t401
      t646 = t3*l1
      t647 = 128*t646
      t648 = z*t52
      t650 = t2*l1
      t652 = t384*t4
      t653 = 64*t652
      t654 = t98*z
      t657 = z*t358
      t659 = z*t361
      t660 = 2944*t659
      t662 = t2*t358
      t664 = t1*t361
      t665 = 10496*t664
      t666 = -263680.D0/3.D0*t403-t647+3328*t648+1792.D0/3.D0*t650-t653+
     #52672.D0/3.D0*t654-185216.D0/3.D0*t592-11712*t657+t660-64*t52+4070
     #4*t662-t665
      t667 = t384*t358
      t668 = 64.D0/3.D0*t667
      t669 = t368*t361
      t670 = 15104*t669
      t671 = t1*t358
      t673 = t2*t361
      t674 = 9856*t673
      t675 = t368*l1
      t676 = 256*t675
      t677 = t1*t4
      t678 = 384*t677
      t679 = z*t49
      t682 = 1024.D0/3.D0*t55
      t683 = t1*t49
      t685 = t368*t4
      t687 = t3*t55
      t689 = 128*t361
      t690 = t668+t670+42432*t671-t674-t676-t678-256*t679-1856.D0/3.D0*t
     #98+t682+18688.D0/3.D0*t683+3840*t685-28928.D0/3.D0*t687-t689
      t693 = s2l1**2
      t696 = 128*t368
      t700 = 25280.D0/3.D0*t5+15488*t596+232320*t598-309760*t600-92928*t
     #602+618624*t359-524672*t362+232320*t606-2816*t609-42240*t611+56320
     #*t613+16896*t615+(t630+t643+t666+t690)*t693-320*t2+t696-496064.D0/
     #3.D0*t364-3904*t623+49024.D0/3.D0*t366
      t720 = -18880*t369-6272*t371-65600.D0/3.D0*t373+12032*t375+2240.D0
     #/3.D0*t377+17408*t380-17728*t382-12288*t385+1603712.D0/3.D0*t387+4
     #08320.D0/3.D0*t389-1077824.D0/3.D0*t391-1180096.D0/3.D0*t393-25984
     #*t395+22528.D0/3.D0*t397+326272.D0/3.D0*t399-16576.D0/3.D0*t401+87
     #104*t403-1216*t646+128*t3
      t722 = 384*t5
      t724 = 12288*t362
      t732 = -t722+57856*t359-t724-t696-472768.D0/3.D0*t364-23168.D0/3.D
     #0*t623-39360*t366+39872.D0/3.D0*t369+27008.D0/3.D0*t371+44416.D0/3
     #.D0*t373-176768.D0/3.D0*t375
      t742 = 128*t397
      t745 = 300608.D0/3.D0*t377-118144.D0/3.D0*t380-22336.D0/3.D0*t382+
     #20608.D0/3.D0*t385+147840*t387+83392.D0/3.D0*t389-197632.D0/3.D0*t
     #391-332032.D0/3.D0*t393-66752.D0/3.D0*t395-t742-24128.D0/3.D0*t399
     #+32128.D0/3.D0*t401
      t753 = 1024*t659
      t755 = 768*t664
      t757 = t384*t361
      t758 = 3328*t757
      t759 = 201536.D0/3.D0*t403+8320.D0/3.D0*t648-384*t650+8960*t654+34
     #048*t592-4864*t657+t753-75008*t662-t755-15616*t667+t758
      t760 = 7168*t669
      t762 = 15872*t673
      t763 = 3200.D0/3.D0*t675
      t764 = z*t4
      t765 = 128*t764
      t767 = 256*t466
      t771 = 128*t5
      t777 = 124288.D0/3.D0*t623
      t782 = t771-22144*t359+4864*t362+3584.D0/3.D0*t358+128*t49-t696+43
     #7248.D0/3.D0*t364-t777+33664.D0/3.D0*t366-43264*t369+5248.D0/3.D0*
     #t371+2048.D0/3.D0*t373
      t795 = 18944.D0/3.D0*t375-66816*t377+138368*t380-270592.D0/3.D0*t3
     #82-384*t385-483712.D0/3.D0*t387-640.D0/3.D0*t389+762112.D0/3.D0*t3
     #91+112640.D0/3.D0*t393+1280.D0/3.D0*t395+t742-526976.D0/3.D0*t399+
     #169344*t401
      t799 = 128*t650
      t809 = -745856.D0/3.D0*t403+52864.D0/3.D0*t648-t799+143104.D0/3.D0
     #*t654-410368.D0/3.D0*t592-26752*t657+5888*t659-2432.D0/3.D0*t52+89
     #344*t662-20992*t664+128.D0/3.D0*t667+30208*t669
      t822 = 95104*t671-19712*t673+2048.D0/3.D0*t675+t765-1664*t677-128*
     #t466-9472.D0/3.D0*t679-2176*t98+1920*t55+72064.D0/3.D0*t683+18688.
     #D0/3.D0*t685-27392*t687-256*t361
      t827 = -t760+3584*t671+t762+t763-t765+1024*t677-t767-128.D0/3.D0*t
     #679+3584.D0/3.D0*t683-2944*t685+(t782+t795+t809+t822)*s2l1+329024.
     #D0/3.D0*t687
      t833 = t2*t368
      t834 = t833*t49
      t835 = 320*t834
      t836 = t833*t4
      t848 = 64*t1
      t854 = -128*t18-11968*t359+t619+2176.D0/3.D0*t358+1408.D0/3.D0*t49
     #-64*t4-t848+99584*t364-32832*t623+39808.D0/3.D0*t366-39360*t369-38
     #4*t371
      t860 = 64*t385
      t867 = -256*t373+6528*t375-53824*t377+108224*t380-228736.D0/3.D0*t
     #382-t860-97216*t387-256.D0/3.D0*t389+465152.D0/3.D0*t391+22464*t39
     #3+128*t395-109568*t399
      t877 = 363904.D0/3.D0*t401-171008*t403+56384.D0/3.D0*t648+31040*t6
     #54-225152.D0/3.D0*t592-15040*t657+t660-4096.D0/3.D0*t52+48640*t662
     #-t665+t668+t670
      t888 = 52672*t671-t674-t676+3200.D0/3.D0*t764-4800*t677+1792.D0/3.
     #D0*t466-17536.D0/3.D0*t679-5056.D0/3.D0*t98+2048*t55+79808.D0/3.D0
     #*t683+22400.D0/3.D0*t685-57664.D0/3.D0*t687-t689
      t891 = s1l1**2
      t894 = t833*t98
      t902 = (t732+t745+t759+t827)*s1l1+8128.D0/3.D0*t648+6400.D0/3.D0*t
     #650-t835+128*t836-1216*t652-29888.D0/3.D0*t654+598464*t592+31488*t
     #657-98432.D0/3.D0*t659+(t854+t867+t877+t888)*t891-829056*t662-1683
     #2*t894+203648*t664-235584*t667+207232*t757-520192*t669-219456*t671
     #+2101504.D0/3.D0*t673
      t906 = t833*t52
      t908 = t833*t55
      t910 = t595*t833
      t914 = 128*t2
      t922 = 7552.D0/3.D0*t5+52736*t359-t724-t914-321728.D0/3.D0*t364-15
     #232.D0/3.D0*t623-14016*t366-3840*t369+896*t371+6464.D0/3.D0*t373-9
     #2288.D0/3.D0*t375
      t934 = 53696*t377-35072.D0/3.D0*t380-30208.D0/3.D0*t382+4160.D0/3.
     #D0*t385+354944.D0/3.D0*t387+22592*t389-50688*t391-269312.D0/3.D0*t
     #393-42880.D0/3.D0*t395+4352.D0/3.D0*t397-22336.D0/3.D0*t399
      t938 = 256*t646
      t941 = 128*t652
      t946 = 33344.D0/3.D0*t401+119488.D0/3.D0*t403-t938+1536*t648+3200.
     #D0/3.D0*t650-t941+20992.D0/3.D0*t654+30464*t592-4352*t657+t753-678
     #40*t662
      t952 = 11008.D0/3.D0*t685
      t954 = -t755-14336*t667+t758-t760+3328*t671+t762-384*t675-t678-128
     #*t679+12416.D0/3.D0*t683+t952+227008.D0/3.D0*t687
      t958 = t384*l1
      t960 = t608*t384
      t964 = t833*t361
      t966 = t833*t358
      t969 = t595*t1
      t971 = t608*t3
      t973 = t608*t833
      t975 = -1216*t675+t765-1216*t677+t767-320*t679+1216*t906+3968.D0/3
     #.D0*t908+15488*t910+(t922+t934+t946+t954)*s2l1+256*t958+16896*t960
     #+7296*t683+6720*t685-101120.D0/3.D0*t964+35520*t966+337024.D0/3.D0
     #*t687-92928*t969-42240*t971-2816*t973
      t1010 = -8+64*t112+43*l1-280*t118+63*t49-80*t4+12*t122-210*t124-48
     #*t35+210*t127-8*t130-18*t52+160*t137-226*t139+43*t141+520*t146-144
     #*t149+82*t151-24*t153-18*t155
      t1033 = 16*l5+720*t222-720*t224+288*t227+720*t230-288*t233-720*t23
     #5+288*t237-288*t239+560*t242-128*t245
      t1056 = 128*t247+560*t250-128*t253-560*t255+128*t258-560*t260-80*t
     #262+144*t265-32*t269+41*t271+108*t273-270*t275-108*t277-406*t279+6
     #3*t281-72*t283+84*t285-292*t289+12*t291+140*t293
      t1079 = -54*t327+32*t329-216*t331+144*t333-30*t335+27*t337+72*t339
     #-360*t341+116*t343+48*t345-370*t347
      t1090 = 2048*t664
      t1094 = 64.D0/3.D0*(t8*t485/2+t9*t485/6)*t2*t21+t494*(-64.D0/3.D0*
     #t48*(30*t4-60*l1+t497+21)*s1l1-64.D0/3.D0*(t519+t539+t560+t579)*t2
     #1-64.D0/3.D0*t48*(78*t4-104*l1-t497+41)*s2l1)/24+46080*t592-t7*(t7
     #00+t720+t902+t975)+t8*(-64.D0/3.D0*t48*(t510+170*t4-160*l1-66*t49+
     #39)*s2l1-64.D0/3.D0*t48*(26*t4-74*l1+6*t49+t510+27)*s1l1-64.D0/3.D
     #0*(t1010-144*t160-126*t162+360*t167+270*t186-86*t36+41*t48+108*t20
     #1-32*t205+180*t207-72*t213+t1033+t1056+64*t299+188*t303-180*t305+4
     #45*t307-8*t80+280*t310-140*t319-64*t321-18*t323+36*t325+t1079)*t21
     #)/2-92160*t662-47104.D0/3.D0*t894+t1090-46080*t667+10240*t757-1024
     #0*t669
      t1095 = 9216*t671
      t1114 = 64.D0/3.D0*t2
      t1115 = 256.D0/3.D0*t368
      t1118 = -3904.D0/3.D0+13568*t5-51392*t596-779328*t598+1033472*t600
     #+313984*t602+3072*t18-7302016.D0/3.D0*t359+1839296*t362-772288*t60
     #6+9344*t609+140416*t611-187904*t613-56064*t615+4032*l1-2240*t49-t1
     #114-t1115-3584.D0/3.D0*t4-1984.D0/3.D0*t1
      t1141 = -121280*t385-7934464.D0/3.D0*t387-2276672.D0/3.D0*t389+189
     #1008*t391+1988928*t393+380352*t395+71744.D0/3.D0*t397-2058176.D0/3
     #.D0*t399+950912.D0/3.D0*t401-2841536.D0/3.D0*t403-1664*t646
      t1147 = 1152*t834
      t1165 = 896*t52+3234112*t662+112832*t894-715520*t664+2871296.D0/3.
     #D0*t667-2207360.D0/3.D0*t757+5440448.D0/3.D0*t669+2723392.D0/3.D0*
     #t671-2440576*t673-t676-2176.D0/3.D0*t764
      t1182 = 9472*t609
      t1183 = 59136*t611
      t1184 = 64384*t613
      t1185 = 32128*t615
      t1189 = 128.D0/3.D0*t368
      t1192 = 9472.D0/3.D0*t5-47232*t596-257152*t598+382656*t600+95040*t
     #602+9856.D0/3.D0*t18-3074240.D0/3.D0*t359+2052608.D0/3.D0*t362-336
     #000*t606+t1182+t1183-t1184-t1185+12736.D0/3.D0*l1+22976.D0/3.D0*t3
     #58-4736.D0/3.D0*t49-t1189-2304*t4-2240.D0/3.D0*t1
      t1212 = 2434304.D0/3.D0*t364-t777+65216*t366-62720*t369-17600.D0/3
     #.D0*t371-107840.D0/3.D0*t373+188608*t375-908800.D0/3.D0*t377+75833
     #6.D0/3.D0*t380-93568*t382-55680*t385-4098560.D0/3.D0*t387-339456*t
     #389+3301760.D0/3.D0*t391+943040*t393+182208*t395+7680*t397-1419776
     #.D0/3.D0*t399+794048.D0/3.D0*t401-653184*t403
      t1222 = 48384*t894
      t1232 = 5184*t595+11136*t648+t799+89792*t654-3627392.D0/3.D0*t592-
     #360896.D0/3.D0*t657+301184.D0/3.D0*t659+1280*t52+4421696.D0/3.D0*t
     #662+t1222-1217408.D0/3.D0*t664+381696*t667-256512*t757+840064*t669
     #+1650176.D0/3.D0*t671-2981888.D0/3.D0*t673-896.D0/3.D0*t675-5120.D
     #0/3.D0*t764-14848.D0/3.D0*t677
      t1235 = 3840*t906
      t1236 = 21120*t908
      t1237 = 14976*t910
      t1239 = 1152*t608
      t1241 = 14976*t960
      t1244 = 40704*t964
      t1245 = 59136*t966
      t1249 = 41728*t971
      t1250 = 2304*t973
      t1252 = -4160.D0/3.D0+9280.D0/3.D0*t466-7040.D0/3.D0*t679+t1235-t1
     #236-t1237-9472.D0/3.D0*t98-t1239+1280.D0/3.D0*t55-t1241+60032.D0/3
     #.D0*t683+28928.D0/3.D0*t685+t1244-t1245-1628672.D0/3.D0*t687-3200.
     #D0/3.D0*z+172480*t969+t1249+t1250-9088*t361
      t1273 = 2304*t609
      t1274 = 33536*t611
      t1275 = 44544*t613
      t1276 = 13568*t615
      t1281 = -2688+4736*t5-19968*t596-256256*t598+354048*t600+99328*t60
     #2+19072.D0/3.D0*t18-3959168.D0/3.D0*t359+2381248.D0/3.D0*t362-2759
     #68*t606+t1273+t1274-t1275-t1276+24832.D0/3.D0*l1-11456.D0/3.D0*t49
     #+t914-384*t368-10496.D0/3.D0*t4
      t1301 = -1408*t1+3445760.D0/3.D0*t364-45248*t623+231680.D0/3.D0*t3
     #66-70656*t369-4352.D0/3.D0*t371-52864*t373+273536*t375-1250944.D0/
     #3.D0*t377+325824*t380-111360*t382-245696.D0/3.D0*t385-5653120.D0/3
     #.D0*t387-1373824.D0/3.D0*t389+4456768.D0/3.D0*t391+1296128*t393+76
     #9792.D0/3.D0*t395+36160.D0/3.D0*t397-599808*t399
      t1320 = 339072*t401-2682304.D0/3.D0*t403+t938+11712*t648-3136.D0/3
     #.D0*t650+t941+287104.D0/3.D0*t654-4498688.D0/3.D0*t592-321856.D0/3
     #.D0*t657+65216*t659+2304*t52+5667520.D0/3.D0*t662+63936*t894-11190
     #40.D0/3.D0*t664+1454528.D0/3.D0*t667-898624.D0/3.D0*t757+2658496.D
     #0/3.D0*t669+1877120.D0/3.D0*t671-1119616*t673
      t1329 = 256*t55
      t1330 = 13568*t960
      t1338 = 33536*t971
      t1339 = 2944.D0/3.D0*t675-2688*t764-17152.D0/3.D0*t677+15616.D0/3.
     #D0*t466-11456.D0/3.D0*t679+5760*t906-29760*t908-16128*t910-t1329-t
     #1330+20096*t683+25280.D0/3.D0*t685+47040*t964-73152*t966-2319296.D
     #0/3.D0*t687-2048*z+114944*t969+t1338+t1250
      t1359 = -3456*t5-5376*t596-112640*t598+135936*t600+49920*t602+1203
     #2.D0/3.D0*t18+523456.D0/3.D0*t359+227008.D0/3.D0*t362-92928*t606+t
     #1273+t1274-t1275-t1276+17792.D0/3.D0*l1-11584.D0/3.D0*t49+1024.D0/
     #3.D0*t2+t1189-2560.D0/3.D0*t4-1792.D0/3.D0*t1
      t1379 = -1799680.D0/3.D0*t364+95168.D0/3.D0*t623-86592*t366+154624
     #.D0/3.D0*t369+51520.D0/3.D0*t371+31168*t373-408640.D0/3.D0*t375+86
     #7328.D0/3.D0*t377-229696*t380+72064*t382+16832*t385+714752*t387+23
     #9296.D0/3.D0*t389-1992256.D0/3.D0*t391-1119808.D0/3.D0*t393-59328*
     #t395-1408*t397+888640.D0/3.D0*t399-618304.D0/3.D0*t401
      t1382 = 768*t646
      t1385 = 384*t652
      t1400 = 524416*t403+t1382-9344*t648-5888.D0/3.D0*t650+t1385-51264*
     #t654+410048*t592+116480.D0/3.D0*t657-7744*t659+768*t52-1187648.D0/
     #3.D0*t662-2880*t894+90560.D0/3.D0*t664-20096*t667-157504.D0/3.D0*t
     #757-87104.D0/3.D0*t669-608320.D0/3.D0*t671-29056*t673+64.D0/3.D0*t
     #675
      t1414 = 256*t5
      t1430 = 640.D0/3.D0*t368
      t1432 = 2048.D0/3.D0*t1
      t1433 = t1414-104192*t596-459008*t598+708224*t600+164736*t602+4224
     #*t18-3746560.D0/3.D0*t359+3127552.D0/3.D0*t362-650496*t606+18944*t
     #609+118272*t611-128768*t613-64256*t615+6144*l1+100480.D0/3.D0*t358
     #-2944*t49+t1430-5888.D0/3.D0*t4-t1432
      t1452 = 437376*t364-127744.D0/3.D0*t623-37888.D0/3.D0*t366-9856.D0
     #/3.D0*t369+t372+2688*t373+37504*t375-143360.D0/3.D0*t377+37120*t38
     #0-41984.D0/3.D0*t382-11520*t385-3822976.D0/3.D0*t387-277632*t389+3
     #309824.D0/3.D0*t391+2496640.D0/3.D0*t393+90240*t395-1661312.D0/3.D
     #0*t399+159744*t401-1056640.D0/3.D0*t403
      t1463 = 34560*t894
      t1472 = -5888.D0/3.D0+12672*t595+12160.D0/3.D0*t648+512*t650+15680
     #0*t654-5176192.D0/3.D0*t592-792704.D0/3.D0*t657+696832.D0/3.D0*t65
     #9-1024.D0/3.D0*t52+5765888.D0/3.D0*t662+t1463-2372992.D0/3.D0*t664
     #+435456*t667-371328*t757+1467392*t669+2729344.D0/3.D0*t671-4822144
     #.D0/3.D0*t673-2048*t675-768*t764
      t1476 = 7680*t908
      t1477 = 25344*t910
      t1484 = 56064*t964
      t1485 = 62208*t966
      t1490 = 4608*t973
      t1492 = -4096*t677+3072*t466-8960.D0/3.D0*t679-t1476-t1477-21632*t
     #98-2304*t608+7680*t55-29952*t960+6272.D0/3.D0*t683+21376.D0/3.D0*t
     #685+t1484-t1485-872320.D0/3.D0*t687-3968.D0/3.D0*z+353408*t969+834
     #56*t971+t1490-28544*t361
      t1496 = -5632.D0/3.D0-1280.D0/3.D0*t764+2176.D0/3.D0*t677+6272.D0/
     #3.D0*t466-3584.D0/3.D0*t679+1920*t908-9216*t910+t1329-t1330-12352*
     #t683-8704.D0/3.D0*t685+12480*t964-4608*t966+310592*t687-3712.D0/3.
     #D0*z+34304*t969+(t1433+t1452+t1472+t1492)*s2l1+t1338+t1250
      t1513 = -576+t722-56960*t596-201856*t598+325568*t600+69696*t602+20
     #48.D0/3.D0*t18-1339328.D0/3.D0*t359+1223168.D0/3.D0*t362-314496*t6
     #06+t1182+t1183-t1184-t1185+5696.D0/3.D0*l1+93056.D0/3.D0*t358-2048
     #.D0/3.D0*t49-t1114
      t1532 = 1024.D0/3.D0*t368+640.D0/3.D0*t4-t848+1071424.D0/3.D0*t364
     #-112320*t623+24768*t366-195008.D0/3.D0*t369-9280.D0/3.D0*t371-1152
     #*t373+20096*t375-364160.D0/3.D0*t377+218560*t380-441472.D0/3.D0*t3
     #82-2062336.D0/3.D0*t387-57600*t389+868480*t391+288704*t393+11520*t
     #395-1770752.D0/3.D0*t399
      t1545 = 3840*t894
      t1553 = 1081088.D0/3.D0*t401-521088*t403+7488*t595+35520*t648+2240
     #.D0/3.D0*t650+604736.D0/3.D0*t654-2969216.D0/3.D0*t592-221120*t657
     #+447104.D0/3.D0*t659-14336.D0/3.D0*t52+876032*t662+t1545-1344128.D
     #0/3.D0*t664+121152*t667-129792*t757+732544*t669+1927616.D0/3.D0*t6
     #71-2125568.D0/3.D0*t673-5888.D0/3.D0*t675
      t1558 = 10368*t910
      t1563 = 17664*t964
      t1564 = 13440*t966
      t1569 = 1792*t764-17920.D0/3.D0*t677+1088*t466-7296*t679-t1558-826
     #88.D0/3.D0*t98-t1239+14912*t55-t1241+38784*t683+42112.D0/3.D0*t685
     #+t1563-t1564-108864*t687-256*z+180928*t969+t1249+t1250-20608*t361
      t1573 = 25600.D0/3.D0*t685+366016.D0/3.D0*t964-459584.D0/3.D0*t966
     #-3121600.D0/3.D0*t687+(t1281+t1301+t1320+t1339)*s2l1-2944.D0/3.D0*
     #z+308352*t969+(t1359+t1379+t1400+t1496)*s1l1+141696*t971+9600*t973
     #+(t1513+t1532+t1553+t1569)*t891
      t1596 = 1088.D0/3.D0*t1
      t1597 = 704-23744.D0/3.D0*t5+49984*t596+763840*t598-1013760*t600-3
     #06944*t602-1920*t18+2248448*t359-5316992.D0/3.D0*t362+756800*t606-
     #9088*t609-137600*t611+184320*t613+54784*t615-2496*l1-640*t49+256.D
     #0/3.D0*t2-t696+6272.D0/3.D0*t4+t1596
      t1627 = 18432*t609
      t1628 = 34560*t611
      t1629 = 5760*t613
      t1630 = 38016*t615
      t1634 = -t1414-117504*t596+55296*t598+60480*t600-46656*t602-1664.D
     #0/3.D0*t18+35968*t359-77952*t362-224640*t606+t1627+t1628-t1629-t16
     #30-3008.D0/3.D0*l1+268288.D0/3.D0*t358+1792.D0/3.D0*t49+t1114-t143
     #0
      t1653 = 320.D0/3.D0*t4+585152.D0/3.D0*t364-713216.D0/3.D0*t623+118
     #40.D0/3.D0*t366-92032.D0/3.D0*t369+5056.D0/3.D0*t371+512*t373+1427
     #2.D0/3.D0*t375-165184.D0/3.D0*t377+514240.D0/3.D0*t380-603776.D0/3
     #.D0*t382-t860-1072960.D0/3.D0*t387+80384.D0/3.D0*t389+2579072.D0/3
     #.D0*t391+13632*t393-3712*t395-886656*t399+1568384.D0/3.D0*t401
      t1671 = 1216.D0/3.D0-492032*t403+22464*t595+255872.D0/3.D0*t648-14
     #72.D0/3.D0*t650+1266880.D0/3.D0*t654-2809088.D0/3.D0*t592-471040*t
     #657+317056*t659-20608.D0/3.D0*t52+371584*t662-t1545-645504*t664-20
     #5568.D0/3.D0*t667+81792*t757+617216*t669+965504*t671-213760*t673+t
     #763
      t1677 = 3456*t608
      t1679 = 10368*t960
      t1685 = 13824*t971
      t1687 = 1600*t764-26176.D0/3.D0*t677+1600.D0/3.D0*t466-54272.D0/3.
     #D0*t679+t1558-222784.D0/3.D0*t98-t1677+33728*t55+t1679+51136*t683-
     #1792.D0/3.D0*t685-t1563+t1564-54272.D0/3.D0*t687+704.D0/3.D0*z+240
     #192*t969-t1685-t1250-61184*t361
      t1691 = 231296.D0/3.D0*t385+6896000.D0/3.D0*t387+1940608.D0/3.D0*t
     #389-5124160.D0/3.D0*t391-5089216.D0/3.D0*t393-854528.D0/3.D0*t395-
     #44672.D0/3.D0*t397+1977088.D0/3.D0*t399-304192*t401+2397440.D0/3.D
     #0*t403+(t1634+t1653+t1671+t1687)*t891
      t1715 = 4608*t609
      t1716 = 69120*t611
      t1717 = 92160*t613
      t1718 = 27648*t615
      t1724 = 1408*t5+30464*t596+382720*t598-535040*t600-145664*t602-704
     #0.D0/3.D0*t18+3366784.D0/3.D0*t359-888320*t362+419840*t606-t1715-t
     #1716+t1717+t1718-10496.D0/3.D0*l1-512*t49-640.D0/3.D0*t2-896.D0/3.
     #D0*t368+6784.D0/3.D0*t4
      t1725 = 1280.D0/3.D0*t1
      t1744 = t1725-2545024.D0/3.D0*t364+213376.D0/3.D0*t623-184832.D0/3
     #.D0*t366+103296*t369+15488.D0/3.D0*t371+27520.D0/3.D0*t373-286336.
     #D0/3.D0*t375+890240.D0/3.D0*t377-1109120.D0/3.D0*t380+519040.D0/3.
     #D0*t382+29824.D0/3.D0*t385+1552512*t387+663040.D0/3.D0*t389-448614
     #4.D0/3.D0*t391-840832*t393-73088*t395+640*t397+717184*t399
      t1756 = 512*t52
      t1765 = -1316864.D0/3.D0*t401+906496*t403-512*t646-53632.D0/3.D0*t
     #648+1280.D0/3.D0*t650-256*t652-133376*t654+4861184.D0/3.D0*t592+42
     #3296.D0/3.D0*t657-87808*t659+t1756-5474816.D0/3.D0*t662-21120*t894
     #+482432*t664-1064192.D0/3.D0*t667+314496*t757-1097856*t669-2245888
     #.D0/3.D0*t671+1322752*t673
      t1773 = 27648*t960
      t1792 = -t722-202752*t596+152064*t598+63360*t600-114048*t602-2688*
     #t18+436352*t359-350464*t362-380160*t606+36864*t609+t1716-11520*t61
     #3-76032*t615-3840*l1+283520.D0/3.D0*t358-1152*t49-1792.D0/3.D0*t36
     #8+10240.D0/3.D0*t4+384*t1
      t1811 = -267776.D0/3.D0*t364-285952.D0/3.D0*t623-37632*t366+62080*
     #t369+9728*t371+4224*t373-96128.D0/3.D0*t375+245120.D0/3.D0*t377-88
     #704*t380+51328.D0/3.D0*t382+7296*t385+17920.D0/3.D0*t387+571904.D0
     #/3.D0*t389+1663744.D0/3.D0*t391-949376.D0/3.D0*t393-183040.D0/3.D0
     #*t395+t742-1986304.D0/3.D0*t399+521600.D0/3.D0*t401
      t1831 = -251264.D0/3.D0*t403+38016*t595+41728.D0/3.D0*t648-1152*t6
     #50+947840.D0/3.D0*t654-2898688.D0/3.D0*t592-519808*t657+453376*t65
     #9-896.D0/3.D0*t52+108928*t662-t1463-937856*t664-901888.D0/3.D0*t66
     #7+256896*t757+848128*t669+1084672*t671-129664*t673+11264.D0/3.D0*t
     #675+6784.D0/3.D0*t764
      t1842 = 2176.D0/3.D0*z
      t1846 = 3200.D0/3.D0+6912*t677-2688*t466-7808.D0/3.D0*t679+t1476+t
     #1477-54016*t98-6912*t608+13440*t55+20736*t960-67456.D0/3.D0*t683-1
     #9072*t685-t1484+t1485+403072.D0/3.D0*t687+t1842+418176*t969-27648*
     #t971-t1490-84352*t361
      t1856 = 69120*t971
      t1857 = 3328.D0/3.D0+2688*t675+3328.D0/3.D0*t764+6144*t677-5888.D0
     #/3.D0*t466+4096.D0/3.D0*t679+3840*t908+23040*t910+t1773-41984*t683
     #+(t1792+t1811+t1831+t1846)*s2l1-62848.D0/3.D0*t685-45696*t964+4454
     #4*t966+1136768.D0/3.D0*t687+768*z-175360*t969-t1856-t1490
      t1869 = -344704.D0/3.D0*t659+256*t52+(t1724+t1744+t1765+t1857)*s1l
     #1-9035264.D0/3.D0*t662-294208.D0/3.D0*t894+2092544.D0/3.D0*t664-88
     #3136*t667+2123008.D0/3.D0*t757-5277184.D0/3.D0*t669-879680*t671+70
     #75072.D0/3.D0*t673
      t1904 = 4160.D0/3.D0*t5-85248*t596+96768*t598+2880*t600-67392*t602
     #-6784.D0/3.D0*t18+457984*t359-285312*t362-155520*t606+t1627+t1628-
     #t1629-t1630-8512.D0/3.D0*l1+61888.D0/3.D0*t358-1408*t49-t621-320.D
     #0/3.D0*t368+9728.D0/3.D0*t4+320*t1
      t1925 = -529024.D0/3.D0*t364+1408.D0/3.D0*t623-8960.D0/3.D0*t366+1
     #7792*t369-2048*t371+4608*t373-194944.D0/3.D0*t375+54528*t377-14009
     #6.D0/3.D0*t380+28800*t382+102464.D0/3.D0*t385+821056.D0/3.D0*t387+
     #245632*t389+119936.D0/3.D0*t391-1290304.D0/3.D0*t393-125824*t395-3
     #584*t397-376064.D0/3.D0*t399-8448*t401+165248.D0/3.D0*t403
      t1943 = 1984.D0/3.D0-t647+15552*t595-11840.D0/3.D0*t648+512.D0/3.D
     #0*t650-t653+152768.D0/3.D0*t654-639104.D0/3.D0*t592-136320*t657+15
     #5776*t659+192*t52-218304*t662-t1222-334464*t664-836288.D0/3.D0*t66
     #7+185472*t757+271616*t669+308672*t671+74240*t673-256.D0/3.D0*t675
      t1955 = t1679-35968.D0/3.D0*t683-t952-t1244+t1245+233600*t687+1472
     #.D0/3.D0*z+177984*t969-t1685-t1250-26624*t361
      t1973 = 21376.D0/3.D0*t5+20224*t596+377600*t598-478720*t600-158464
     #*t602-12544.D0/3.D0*t18+3165056.D0/3.D0*t359-862720*t362+340480*t6
     #06-t1715-t1716+t1717+t1718-16000.D0/3.D0*l1-2048*t49-384*t2+t1115+
     #16000.D0/3.D0*t4
      t1992 = t1432-234240*t364+18432*t623+40832*t366+4096.D0/3.D0*t369-
     #42880.D0/3.D0*t371-68096.D0/3.D0*t373-32384.D0/3.D0*t375-56704.D0/
     #3.D0*t377-76928.D0/3.D0*t380+110720.D0/3.D0*t382+72832.D0/3.D0*t38
     #5+779648*t387+1047424.D0/3.D0*t389-1450240.D0/3.D0*t391-725248*t39
     #3-138112*t395+512.D0/3.D0*t397+176128*t399
      t2009 = 4096.D0/3.D0-261760.D0/3.D0*t401+164608*t403-t1382-16000.D
     #0/3.D0*t648+2304*t650-t1385-30592*t654+2316928.D0/3.D0*t592+115840
     #.D0/3.D0*t657-36608*t659+t1756-3628288.D0/3.D0*t662-65664*t894+249
     #472*t664-1448320.D0/3.D0*t667+378496*t757-701056*t669-791936.D0/3.
     #D0*t671
      t2025 = 1041152*t673-4480.D0/3.D0*t675+12928.D0/3.D0*t764+15616.D0
     #/3.D0*t677-10624.D0/3.D0*t466-1408*t679-t1235+24960*t908+27648*t91
     #0+t1773-37504.D0/3.D0*t683+896*t685-68736*t964+90240*t966+754816.D
     #0/3.D0*t687+1024*z-128768*t969-t1856-t1490
      t2031 = -57472.D0/3.D0*t683-27328.D0/3.D0*t685-117248*t964+426880.
     #D0/3.D0*t966+775296*t687+1600.D0/3.D0*z-301312*t969+(t1904+t1925+t
     #1943+7808.D0/3.D0*t764+12416.D0/3.D0*t677-6592.D0/3.D0*t466-896*t6
     #79-t1235+t1236+t1237-19136.D0/3.D0*t98-t1677+t682+t1955)*t693+(t19
     #73+t1992+t2009+t2025)*s2l1-138880*t971-9344*t973
      t2057 = 1792*t362
      t2059 = 256*t49
      t2064 = -t771-512*t18+8320*t359-t2057-896*l1+t2059-t696+384*t4-175
     #744.D0/3.D0*t364+30592.D0/3.D0*t623-25600.D0/3.D0*t366
      t2076 = 20480*t369+1024*t371+384*t373-11072.D0/3.D0*t375+30976*t37
     #7-163648.D0/3.D0*t380+82816.D0/3.D0*t382+183680.D0/3.D0*t387-92416
     #*t391-43712.D0/3.D0*t393+178304.D0/3.D0*t399
      t2084 = 8064*t657
      t2085 = 1792*t659
      t2087 = 7168*t664
      t2088 = 896.D0/3.D0-166528.D0/3.D0*t401+92416*t403-7936.D0/3.D0*t6
     #48-256*t650-41024.D0/3.D0*t654+49152*t592+t2084-t2085-33024*t662+t
     #2087
      t2089 = 10752*t669
      t2091 = 7168*t673
      t2101 = -t2089-32512*t671+t2091+640*t675+448.D0/3.D0*t764-640.D0/3
     #.D0*t677-256.D0/3.D0*t466+256.D0/3.D0*t679-5888*t683-15616.D0/3.D0
     #*t685+34432.D0/3.D0*t687+448.D0/3.D0*z
      t2106 = 448.D0/3.D0*t368
      t2110 = -10304.D0/3.D0*t5+11264*t596+168960*t598-225280*t600-67584
     #*t602-4096.D0/3.D0*t18+1564736.D0/3.D0*t359-397376*t362+168960*t60
     #6-2048*t609-30720*t611+40960*t613+12288*t615+(t2064+t2076+t2088+t2
     #101)*t891-5248.D0/3.D0*l1+t2106+1152*t4+t1596-287360*t364+8192*t62
     #3
      t2131 = -14080*t366+34432.D0/3.D0*t369-64*t371+16896*t373-216320.D
     #0/3.D0*t375+268352.D0/3.D0*t377-184064.D0/3.D0*t380+19648*t382+259
     #84*t385+1692544.D0/3.D0*t387+159488*t389-408256*t391-420416*t393-2
     #38208.D0/3.D0*t395-5632*t397+150144*t399-69376*t401+605312.D0/3.D0
     #*t403+448*t646-64*t3
      t2151 = 1984.D0/3.D0-4672.D0/3.D0*t648+128.D0/3.D0*t650+t835-64*t8
     #36+768*t652-64448.D0/3.D0*t654+1546816.D0/3.D0*t592+93632.D0/3.D0*
     #t657-77120.D0/3.D0*t659-698368*t662-70976.D0/3.D0*t894+156800*t664
     #-610816.D0/3.D0*t667+157824*t757-396096*t669-596480.D0/3.D0*t671+1
     #591552.D0/3.D0*t673-1216.D0/3.D0*t675+832*t764
      t2166 = 2176.D0/3.D0+448.D0/3.D0*t5-4352.D0/3.D0*t18+7808*t359-t20
     #57-5504.D0/3.D0*l1-t2059+t2106+1408*t4+t1725-47232*t364
      t2178 = 10240*t623-3328.D0/3.D0*t366+13568*t369-5504.D0/3.D0*t371-
     #2944.D0/3.D0*t373-2176.D0/3.D0*t375+57856.D0/3.D0*t377-130496.D0/3
     #.D0*t380+76928.D0/3.D0*t382+54912*t387-86016*t391
      t2189 = -12480*t393+57216*t399-156544.D0/3.D0*t401+81152*t403-3008
     #*t648+896.D0/3.D0*t650-13632*t654+47616*t592+t2084-t2085-31488*t66
     #2
      t2199 = t2087-t2089-32000*t671+t2091-1280.D0/3.D0*t675+1024*t764+1
     #536*t677-3968.D0/3.D0*t466-19456.D0/3.D0*t683-6400.D0/3.D0*t685+24
     #064.D0/3.D0*t687+576*z
      t2211 = 2048*t596
      t2212 = 10240*t598
      t2213 = 20480*t600
      t2214 = 2048*t602
      t2218 = 20480*t606
      t2221 = 960+1280*t5-t2211-t2212+t2213+t2214-5632.D0/3.D0*t18-12044
     #8*t359+54528*t362-t2218-7936.D0/3.D0*l1+t2059-t914+1536*t4+t1596
      t2237 = 583936.D0/3.D0*t364-12736.D0/3.D0*t623+26752*t366-8832*t36
     #9-17728.D0/3.D0*t371-9920*t373+125440.D0/3.D0*t375-264896.D0/3.D0*
     #t377+58432*t380-31168.D0/3.D0*t382-4224*t385-813568.D0/3.D0*t387-2
     #7904*t389+238336*t391+429184.D0/3.D0*t393
      t2264 = 512.D0/3.D0*t369-512.D0/3.D0*t371-896.D0/3.D0*t373+1664*t3
     #75-23296.D0/3.D0*t377+29824.D0/3.D0*t380-14336.D0/3.D0*t382-50432.
     #D0/3.D0*t387+25088*t391+12416.D0/3.D0*t393-48896.D0/3.D0*t399
      t2285 = t1095-2048*t673+1664.D0/3.D0*t675+1152*t764+5888.D0/3.D0*t
     #677-4352.D0/3.D0*t466+384*t679-2560.D0/3.D0*t683-5888.D0/3.D0*t685
     #-3712*t687+t1842
      t2296 = 16000*t395+512*t397-287488.D0/3.D0*t399+149888.D0/3.D0*t40
     #1-156544*t403-t938+(1024-5888.D0/3.D0*t18-2304*t359+512*t362-8192.
     #D0/3.D0*l1+1792*t4+t1725+47104.D0/3.D0*t364-8320.D0/3.D0*t623+5632
     #.D0/3.D0*t366+t2264+13824*t401-23040*t403+1280.D0/3.D0*t648+11648.
     #D0/3.D0*t654-13824*t592-2304*t657+512*t659+9216*t662-t1090+3072*t6
     #69+t2285)*s2l1-448.D0/3.D0*t648+768*t650-t941+40576.D0/3.D0*t654-2
     #11968*t592-16768*t657+9472*t659+228352*t662
      t2312 = 10240*t969
      t2313 = -49152*t664+24960*t667-11264*t757+101888*t669+95872*t671-1
     #05472*t673+832.D0/3.D0*t675+960*t764+3904.D0/3.D0*t677-3712.D0/3.D
     #0*t466+1600.D0/3.D0*t679-3328.D0/3.D0*t683-1472*t685-299456.D0/3.D
     #0*t687+1984.D0/3.D0*z+t2312
      t2325 = 1088.D0/3.D0*t5+t2211+t2212-t2213-t2214-2816*t18+270848.D0
     #/3.D0*t359-47872*t362+t2218-3584*l1-t2059+1088.D0/3.D0*t368+2560*t
     #4+2368.D0/3.D0*t1
      t2340 = -391168.D0/3.D0*t364+40384.D0/3.D0*t623-1280.D0/3.D0*t366+
     #14848*t369-11456.D0/3.D0*t371-5824.D0/3.D0*t373-5824*t375+119360.D
     #0/3.D0*t377-196288.D0/3.D0*t380+95168.D0/3.D0*t382+645632.D0/3.D0*
     #t387+35840.D0/3.D0*t389-737792.D0/3.D0*t391-257536.D0/3.D0*t393
      t2356 = -10240.D0/3.D0*t395+389632.D0/3.D0*t399-254464.D0/3.D0*t40
     #1+493696.D0/3.D0*t403-3072*t648+2176.D0/3.D0*t650-75776.D0/3.D0*t6
     #54+648448.D0/3.D0*t592+70912.D0/3.D0*t657-11008*t659-600832.D0/3.D
     #0*t662+53248*t664-47104.D0/3.D0*t667+9216*t757
      t2369 = 4160.D0/3.D0-102912*t669-342272.D0/3.D0*t671+99328*t673-10
     #88*t675+1920*t764+2240*t677-2304*t466+64*t679-6400*t683-3136.D0/3.
     #D0*t685+121792.D0/3.D0*t687+1088*z-t2312
      t2375 = 1088*t677-3328.D0/3.D0*t466+192*t679-2368*t906+10112*t908+
     #11264*t910-128*t958+12288*t960+(t2166+t2178+t2189+t2199)*t693-1024
     #0.D0/3.D0*t683-3904.D0/3.D0*t685-77888.D0/3.D0*t964+97216.D0/3.D0*
     #t966+216064*t687+512*z-67584*t969+(t2221+t2237+t2296+t2313)*s1l1+(
     #t2325+t2340+t2356+t2369)*s2l1-30720*t971-2048*t973
      s1 = -t1095+20480*t673-10240.D0/3.D0*t906+35840.D0/3.D0*t908-t9*(2
     #90240.D0/3.D0*t654-78912*t373+1054336.D0/3.D0*t375+351296.D0/3.D0*
     #t659+t1141+t1165-3072*t677-2496*t679+32384.D0/3.D0*t906+44800.D0/3
     #.D0*t683-7110656.D0/3.D0*t592+t1573-36608*t623+256*t3+t1118+512*t9
     #58+4128704.D0/3.D0*t364+258112.D0/3.D0*t366+256*t836+886144.D0/3.D
     #0*t380-267584.D0/3.D0*t382-174656.D0/3.D0*t369+2624*t466-8512*t371
     #-424768.D0/3.D0*t657-145024.D0/3.D0*t908-1348096.D0/3.D0*t377-5708
     #8*t960+(t1192+t1212+t1232+t1252)*t693+2944.D0/3.D0*t650+7680*t648-
     #52800*t910-2816*t652-t1147)/6+t8*(-306112.D0/3.D0*t654+112832.D0/3
     #.D0*t373-634688.D0/3.D0*t375+t1691+t2031+12992.D0/3.D0*t677+1152*t
     #675+512*t679-23168.D0/3.D0*t906+2241536*t592+41920*t623-320*t3+192
     #0*t646-640*t958-3195584.D0/3.D0*t364-128192.D0/3.D0*t366-320*t836-
     #231360*t380+270784.D0/3.D0*t382+46144*t369-1984*t466+t1869+14848.D
     #0/3.D0*t371+t1597+426880.D0/3.D0*t657+113600.D0/3.D0*t908+288576*t
     #377+55808*t960+4160.D0/3.D0*t764-4288.D0/3.D0*t650-26240.D0/3.D0*t
     #648+51392*t910+3008*t652+t1147)/2
      t2380 = s1-5120.D0/3.D0*t683-5120.D0/3.D0*t685-2048*t964+9216*t966
     #+465920.D0/3.D0*t687+t494*(t2110+t2131+t2151+t2375)/24
      t2384 = (l1+z-t18)**2
      taux9 = (t379+t405+t1094+t2380)/t2384*t21
      t1 = l1**2
      t2 = t1**2
      t3 = t2*l1
      t4 = z*t3
      t6 = z**2
      t7 = t6**2
      t8 = t7*z
      t9 = t8*t3
      t11 = t1*l1
      t12 = t2*t11
      t13 = t8*t12
      t15 = t2**2
      t16 = t8*t15
      t18 = t6*z
      t19 = t7*t18
      t20 = t19*t3
      t22 = t8*t1
      t24 = t7*t3
      t26 = t18*t3
      t33 = 20992*t22
      t34 = t15*t11
      t35 = t34*t18
      t36 = 645120*t35
      t37 = t34*t7
      t38 = 645120*t37
      t39 = t34*t8
      t40 = 387072*t39
      t44 = 18432*t34
      t45 = t6*t11
      t47 = t34*t19
      t48 = 18432*t47
      t49 = t6*t3
      t51 = 1536*t6
      t52 = t7*t6
      t53 = t52*t12
      t55 = t6*t15
      t57 = t8*l1
      t58 = 4096*t57
      t59 = t18*l1
      t60 = 88576.D0/3.D0*t59
      t61 = z*t1
      t62 = 24576*t61
      t63 = 9814016.D0/3.D0*t26+t44-1456384.D0/3.D0*t45-t48-2563328*t49+
     #t51+2423296.D0/3.D0*t53-7305984*t55+t58+t60-t62
      t65 = t6*t1
      t67 = t6*l1
      t68 = 56320.D0/3.D0*t67
      t69 = z*t11
      t71 = t7*t11
      t73 = t15*l1
      t74 = t73*t19
      t76 = t2*t1
      t77 = t19*t76
      t79 = t15*t1
      t81 = t7*t12
      t83 = t52*t15
      t85 = t79*t8
      t88 = t79*t19
      t90 = t79*t52
      t92 = t7*t15
      t94 = t7*t1
      t96 = z*t15
      t98 = 2048*t7
      t99 = z*t12
      t101 = t18*t11
      t103 = t34*t6
      t104 = 387072*t103
      t105 = t19*t15
      t107 = z*t2
      t109 = 101376*t88-728064*t90-11147264*t92+107008*t94+2424320*t96+t
     #98-2891776.D0/3.D0*t99+2170880.D0/3.D0*t101+t104+207360*t105-15749
     #12.D0/3.D0*t107
      t112 = z*t79
      t114 = t8*t2
      t116 = t73*t18
      t119 = t7*t2
      t121 = t6*t12
      t123 = t8*t11
      t125 = t52*t2
      t127 = t7*l1
      t128 = 17920*t127
      t129 = t6*t2
      t134 = t19*t12
      t136 = t18*t1
      t138 = 2048.D0/3.D0*t1
      t139 = 3584*t18
      t141 = t18*t2
      t143 = t18*t15
      t145 = t79*t6
      t147 = t79*t7
      t149 = -332288*t15+220160.D0/3.D0*t12-67584*t134-484864.D0/3.D0*t1
     #36+t138-t139-180224*t3-6347776.D0/3.D0*t141+11838976*t143-2405376*
     #t145-3824640*t147
      t151 = t79*t18
      t153 = 16640.D0/3.D0*t11
      t155 = t52*t1
      t156 = 2048*t155
      t157 = t73*t6
      t159 = t73*t52
      t162 = t76*t7
      t164 = t52*t76
      t166 = t76*z
      t168 = t18*t76
      t170 = 3916800*t151-t153+174592.D0/3.D0*t2+t156+5995008*t157+16317
     #44*t159+561664.D0/3.D0*t76+438784.D0/3.D0*t162+205312.D0/3.D0*t164
     #-2280448.D0/3.D0*t166-2084864.D0/3.D0*t168
      t171 = t76*t8
      t173 = t52*t3
      t175 = t52*t11
      t176 = 6656*t175
      t177 = t76*t6
      t179 = t34*t52
      t180 = 129024*t179
      t181 = z*l1
      t182 = 3072*t181
      t183 = t73*t7
      t185 = t73*z
      t187 = t73*t8
      t189 = t18*t12
      t191 = t34*z
      t192 = 129024*t191
      t193 = -44544*t171-192768*t173-t176+1120768*t177+t180+t182+9210880
     #*t183-2056192*t185-5231616*t187-6523136*t189-t192
      t199 = 12*t3
      t200 = 10*t1
      t201 = 35*t11
      t202 = 36*t2
      t205 = 2*l6-1
      t206 = t205**2
      t209 = l5-1
      t211 = -1+l1
      t213 = (1-l1+t181)**2
      t214 = 1/t213
      t219 = dsqrt(z)
      t221 = dsqrt(l5)
      t223 = dsqrt(-t209)
      t226 = dsqrt(l1)
      t230 = l5*z*l1-l1+1-l5+l5*l1
      t233 = dsqrt(-t211)
      t238 = t230**2
      t241 = dsqrt(4*l1*t206*z*l5*t209*t211*t214-4*t205*t219*t221*t223*t
     #226*t230*t214*t233+t238*t214)
      t244 = dlog(z)
      t245 = t244**2
      t262 = 1774080*t35
      t263 = 1774080*t37
      t264 = 1064448*t39
      t266 = 2048+3686272.D0/3.D0*t4+9783168*t9+92898176.D0/3.D0*t13-296
     #43392*t16+309120*t20-46464*t22+t262-t263+t264-50508928.D0/3.D0*t24
      t268 = 50688*t34
      t270 = 50688*t47
      t272 = 2048*t6
      t278 = 45638528.D0/3.D0*t26-t268-364928*t45+t270-20461952.D0/3.D0*
     #t49+t272-10365952*t53+23489408*t55+4352*t57+104192.D0/3.D0*t59-682
     #24.D0/3.D0*t61
      t288 = t19*t2
      t292 = 116864.D0/3.D0*t65-15232.D0/3.D0*t67+101888.D0/3.D0*t69-331
     #7888.D0/3.D0*t71+925440*t74-868992*t77+278784*t79-147740672.D0/3.D
     #0*t81-49920*t288+10173568*t83-6614784*t85
      t300 = 2048*z
      t302 = 1064448*t103
      t304 = -329472*t88+2255616*t90+46949760*t92+461312.D0/3.D0*t94-676
     #4928*t96+2176*t7+17800832.D0/3.D0*t99+t300+1145728*t101-t302-14657
     #28*t105
      t313 = 20480.D0/3.D0*l1
      t318 = -241024*t107-2002176*t112-2771712*t114+27831552*t116-639488
     #*t73+16072832.D0/3.D0*t119+t313-68690432.D0/3.D0*t121+449664*t123+
     #2086016.D0/3.D0*t125-19328*t127
      t324 = 44032.D0/3.D0*t1
      t330 = 6020096.D0/3.D0*t129+787072*t15-1671680.D0/3.D0*t12+1428864
     #*t134-630016.D0/3.D0*t136-t324-3200*t18-148352.D0/3.D0*t3-5009536*
     #t141-43525760*t143+6158592*t145
      t343 = 10771200*t147-10517760*t151+2432.D0/3.D0*t11+12032*t2+2176*
     #t155-15725056*t157-6309376*t159+673408.D0/3.D0*t76+35004032*t162+2
     #0481664.D0/3.D0*t164-3407232*t166
      t349 = 354816*t179
      t350 = 20480.D0/3.D0*t181
      t355 = 354816*t191
      t356 = -94402688.D0/3.D0*t168-64479104.D0/3.D0*t171-8510080.D0/3.D
     #0*t173-261248.D0/3.D0*t175+45543296.D0/3.D0*t177-t349+t350-2926336
     #0*t183+4878592*t185+18301696*t187+134215040.D0/3.D0*t189+t355
      t362 = t52*l1
      t367 = 2048+798848.D0/3.D0*t4-1698176*t9-4843648.D0/3.D0*t13+73344
     #0*t16-710528.D0/3.D0*t20+117376*t22-3272192.D0/3.D0*t24+9542528.D0
     #/3.D0*t26-283136.D0/3.D0*t45-5390848.D0/3.D0*t49+(t266+t278+t292+t
     #304+t318+t330+t343+t356)*s2l1+t272+5218816.D0/3.D0*t53+6656*t362+1
     #291776*t55-22528*t57+30464.D0/3.D0*t59
      t375 = t19*t11
      t383 = t19*t1
      t388 = -49792.D0/3.D0*t61-13696.D0/3.D0*t65+8960.D0/3.D0*t67+19712
     #.D0/3.D0*t69+41984.D0/3.D0*t71-45568*t74+1368320.D0/3.D0*t77-15104
     #*t375-4859392.D0/3.D0*t81+198656.D0/3.D0*t288-827904*t83+2360320.D
     #0/3.D0*t92-143872.D0/3.D0*t94-912128.D0/3.D0*t96+3328*t383-384*t7+
     #1799296.D0/3.D0*t99+t300+333952*t101
      t407 = 685312.D0/3.D0*t105-130432.D0/3.D0*t107+2947456.D0/3.D0*t11
     #4+383488*t116+966400.D0/3.D0*t119+3328*t8+t313-7980032.D0/3.D0*t12
     #1-1228160.D0/3.D0*t123-587392*t125+15360*t127+1846912.D0/3.D0*t129
     #-456576*t134-56320*t136-t324+128*t18-2048*t3-3993088.D0/3.D0*t141-
     #1908480*t143
      t419 = 2048+2474240*t4+20415232.D0/3.D0*t9+23233024*t13-23742080*t
     #16+49920*t20-67712*t22+t262-t263+t264-14572544*t24
      t429 = 46508288.D0/3.D0*t26-t268-435712*t45+t270-8710400*t49+3456*
     #t6-18558848.D0/3.D0*t53+29224832*t55+13312*t57+45824*t59-111232.D0
     #/3.D0*t61
      t442 = 78080*t65-39808.D0/3.D0*t67+368128.D0/3.D0*t69-1154304*t71+
     #646656*t74-259200*t77+329472*t79-133551232.D0/3.D0*t81+6921344*t83
     #-6158592*t85-278784*t88
      t452 = 2002176*t90+43573120*t92+246656*t94-9963904*t96+6656*t7+296
     #30464.D0/3.D0*t99+t300+1038080*t101-t302-819072*t105-1847680.D0/3.
     #D0*t107
      t465 = -2255616*t112-5969152.D0/3.D0*t114+29225472*t116-918272*t73
     #+14239744.D0/3.D0*t119+t313-29895552*t121+391424*t123+243712*t125-
     #44928*t127+7153024.D0/3.D0*t129
      t477 = 1426560*t15-4032896.D0/3.D0*t12+609792*t134-212992*t136-398
     #08.D0/3.D0*t1-2944*t18-273920*t3-14317952.D0/3.D0*t141-46620800*t1
     #43+6614784*t145+10517760*t147
      t490 = -10771200*t151-15232.D0/3.D0*t11+178304.D0/3.D0*t2+6656*t15
     #5-18234112*t157-4915456*t159+781696*t76+92506496.D0/3.D0*t162+1077
     #8240.D0/3.D0*t164-18719360.D0/3.D0*t166-100259456.D0/3.D0*t168
      t500 = -15339648*t171-1284352*t173-25728*t175+20048512*t177-t349+2
     #8928.D0/3.D0*t181-27869440*t183+6272512*t185+15792640*t187+1446707
     #20.D0/3.D0*t189+t355
      t518 = 6656.D0/3.D0*t11+6144*t2-37248*t155-259072*t157+(t419+t429+
     #t442+t452+t465+t477+t490+t500)*s1l1+164864*t159+5337088.D0/3.D0*t1
     #62-5972992.D0/3.D0*t164-577792*t166-13859072.D0/3.D0*t168+2050304*
     #t171+4102144.D0/3.D0*t173+495104.D0/3.D0*t175+2903040*t177+t350-15
     #9744*t183+60928*t185-144896*t187+12034688.D0/3.D0*t189
      t525 = t245**2
      t536 = 256*t6
      t539 = 1024*t57
      t541 = 256+35584.D0/3.D0*t4+1199360.D0/3.D0*t9+479232*t13-217088*t
     #16-5120*t22-2337280.D0/3.D0*t24+616448*t26-1536*t45-555008.D0/3.D0
     #*t49+t536-99840*t53+192512*t55+t539+2304*t59
      t542 = 2304*t61
      t543 = 1280.D0/3.D0*t65
      t544 = 1024.D0/3.D0*t67
      t552 = 512*t7
      t554 = 256*z
      t557 = -t542-t543+t544+4352.D0/3.D0*t69-345088.D0/3.D0*t71-901120*
     #t81+45056*t83+417792*t92+25344*t94-36864*t96+t552+62464*t99+t554+1
     #78688.D0/3.D0*t101+1792.D0/3.D0*t107
      t560 = 81920*t116
      t562 = 2560.D0/3.D0*l1
      t569 = 2048*t1
      t572 = 1024*t11
      t573 = 512*t155
      t574 = -503296.D0/3.D0*t114+t560+1092352.D0/3.D0*t119+t562-367104*
     #t121+39680*t123+16896*t125-3072*t127+43264*t129-27392.D0/3.D0*t136
     #-t569-746752.D0/3.D0*t141-401408*t143+t572+t573
      t575 = 40960*t157
      t576 = 8192*t159
      t585 = 2560.D0/3.D0*t181
      t586 = 81920*t183
      t587 = 8192*t185
      t588 = 40960*t187
      t590 = -t575-t576+3217408.D0/3.D0*t162+111616*t164-139264.D0/3.D0*
     #t166-925696*t168-1712128.D0/3.D0*t171-64000*t173-2048*t175+1076224
     #.D0/3.D0*t177+t585-t586+t587+t588+826368*t189
      t607 = 256+111616.D0/3.D0*t4+(t541+t557+t574+t590)*s1l1-418304*t9-
     #1338112.D0/3.D0*t13+206080*t16-40448*t20+14336*t22+52736.D0/3.D0*t
     #24+1166336.D0/3.D0*t26-47104.D0/3.D0*t45-247808*t49+t536+1034240.D
     #0/3.D0*t53+512*t362+163328*t55-1792*t57
      t623 = 4864.D0/3.D0*t59-t542-256*t65+t544+1792*t69+8192*t71-8192*t
     #74+242944.D0/3.D0*t77-1280*t375+7168*t81+9472*t288-164352*t83-5120
     #.D0/3.D0*t92-11008.D0/3.D0*t94-122624.D0/3.D0*t96+256*t383+242432.
     #D0/3.D0*t99+t554
      t640 = 256+179456.D0/3.D0*t4+188672*t9+1103360.D0/3.D0*t13-192512*
     #t16-256*t22-1831424.D0/3.D0*t24+2277376.D0/3.D0*t26-92672.D0/3.D0*
     #t45-383488*t49+t536-188416.D0/3.D0*t53+217088*t55+7936.D0/3.D0*t59
      t654 = -2560*t61+1280*t65-512.D0/3.D0*t67+2304*t69-156160.D0/3.D0*
     #t71-2477056.D0/3.D0*t81+36864*t83+401408*t92+12032.D0/3.D0*t94-450
     #56*t96+297472.D0/3.D0*t99+t554+306176.D0/3.D0*t101-41984.D0/3.D0*t
     #107
      t666 = -44800*t114+t560+240896*t119+t562-1430528.D0/3.D0*t121+4864
     #.D0/3.D0*t123-512*t127+457984.D0/3.D0*t129-64768.D0/3.D0*t136-t569
     #-256*t18-1028864.D0/3.D0*t141-417792*t143+t572
      t675 = -t575-t576+2768896.D0/3.D0*t162+143360.D0/3.D0*t164-326656.
     #D0/3.D0*t166-3184640.D0/3.D0*t168-361472*t171-40960.D0/3.D0*t173+5
     #61152*t177+t585-t586+t587+t588+2695168.D0/3.D0*t189
      t690 = 57344*t101+123136.D0/3.D0*t105-7424*t107+(t640+t654+t666+t6
     #75)*s2l1+210176*t114+40960*t116-42752.D0/3.D0*t119+256*t8+t562-101
     #8880.D0/3.D0*t121-209408.D0/3.D0*t123-103936*t125+1024*t127+275200
     #.D0/3.D0*t129-245504.D0/3.D0*t134-34304.D0/3.D0*t136-t569
      t693 = 3072*t155
      t705 = -557312.D0/3.D0*t141-203520*t143+t572-t693-32768*t157+32768
     #*t159-14336*t162-392192*t164-235264.D0/3.D0*t166-521984*t168+16360
     #96.D0/3.D0*t171+788992.D0/3.D0*t173+22528*t175+1141760.D0/3.D0*t17
     #7+t585+t587-t588+1304320.D0/3.D0*t189
      t710 = -197120.D0/3.D0*t4-5235712.D0/3.D0*t9-6376448.D0/3.D0*t13+1
     #003520*t16-183040.D0/3.D0*t20+64000.D0/3.D0*t22+6899200.D0/3.D0*t2
     #4-1531136*t26+(3223040.D0/3.D0*t4+2870272.D0/3.D0*t9-9670912.D0/3.
     #D0*t13+6081024*t16+15360*t20-7130624.D0/3.D0*t24+t170-t33-t36+t38-
     #t40+1510656*t129+t63+323072.D0/3.D0*t65-t68+300032*t73-1460480.D0/
     #3.D0*t71+136960*t69-23040*t77-215040*t74+18722560.D0/3.D0*t81-1198
     #08*t79-1766144*t83+t193+2239488*t85+t109+820224*t112-9634816*t116-
     #1470464.D0/3.D0*t114+1502208*t119+163840.D0/3.D0*t125+116224*t123+
     #10968832.D0/3.D0*t121-t128+t149)*s1l1+51200.D0/3.D0*t45+256.D0/3.D
     #0*(t199-t200+t201-t202)*t241+t245*(t367+t388+t407+t518)/2+1518080.
     #D0/3.D0*t49+830976*t53+t525*(t607+t623+t690+t705)/24
      t713 = 13312*t74
      t725 = t245*t244
      t740 = 483840*t35
      t741 = 483840*t37
      t742 = 290304*t39
      t744 = -1280-465664*t4-3108864*t9-28082432.D0/3.D0*t13+8571264*t16
     #-74880*t20+11008*t22-t740+t741-t742+6200832*t24
      t746 = 13824*t34
      t748 = 13824*t47
      t750 = 1280*t6
      t755 = -18492800.D0/3.D0*t26+t746+543232.D0/3.D0*t45-t748+8533120.
     #D0/3.D0*t49-t750+2939776*t53-6928896*t55-t539-49664.D0/3.D0*t59+38
     #144.D0/3.D0*t61
      t768 = -40192.D0/3.D0*t65+4864.D0/3.D0*t67-38656.D0/3.D0*t69+13315
     #84.D0/3.D0*t71-251136*t74+220032*t77-76032*t79+15592704*t81+11520*
     #t288-8548864.D0/3.D0*t83+1804032*t85
      t775 = 1280*z
      t777 = 290304*t103
      t779 = 89856*t88-615168*t90-13855744*t92-161792.D0/3.D0*t94+192678
     #4*t96-t552-1814784*t99-t775-590336*t101+t777+392832*t105
      t788 = 12800.D0/3.D0*l1
      t793 = 94208*t107+546048*t112+2536064.D0/3.D0*t114-7800576*t116+17
     #2288*t73-6307328.D0/3.D0*t119-t788+22252928.D0/3.D0*t121-349312.D0
     #/3.D0*t123-172800*t125+6400*t127
      t799 = 28160.D0/3.D0*t1
      t805 = -2840320.D0/3.D0*t129-615424.D0/3.D0*t15+406784.D0/3.D0*t12
     #-374400*t134+115712*t136+t799+1536*t18+9472*t3+6868864.D0/3.D0*t14
     #1+38845568.D0/3.D0*t143-1679616*t145
      t817 = -2937600*t147+2868480*t151-5632.D0/3.D0*t11-4096*t2-t573+43
     #85024*t157+1738496*t159-144640.D0/3.D0*t76-11967488*t162-5764864.D
     #0/3.D0*t164+1159296*t166
      t823 = 96768*t179
      t824 = 12800.D0/3.D0*t181
      t829 = 96768*t191
      t830 = 34099328.D0/3.D0*t168+20225920.D0/3.D0*t171+2290816.D0/3.D0
     #*t173+63232.D0/3.D0*t175-16651264.D0/3.D0*t177+t823-t824+8195840*t
     #183-1343744*t185-5096192*t187-43607168.D0/3.D0*t189-t829
      t842 = -1280-555008.D0/3.D0*t4+4934912.D0/3.D0*t9+1705728*t13-7846
     #40*t16+182400*t20-69888*t22+307968*t24-6150784.D0/3.D0*t26+(t744+t
     #755+t768+t779+t793+t805+t817+t830)*s2l1+77824*t45+1236096*t49-t750
     #-4416640.D0/3.D0*t53-3072*t362-839680*t55+10496*t57-24320.D0/3.D0*
     #t59
      t859 = 32000.D0/3.D0*t61+t543-5120.D0/3.D0*t67-17920.D0/3.D0*t69-2
     #7776.D0/3.D0*t71+36352*t74-1082624.D0/3.D0*t77+7424*t375+1441280.D
     #0/3.D0*t81-135424.D0/3.D0*t288+701440*t83-711680.D0/3.D0*t92+52480
     #.D0/3.D0*t94+613120.D0/3.D0*t96-1536*t383-1212160.D0/3.D0*t99-t775
     #-269568*t101
      t877 = -1280-1826176.D0/3.D0*t4-2935168*t9-7591808*t13+7036800*t16
     #-11520*t20+31232*t22-t740+t741-t742+18196480.D0/3.D0*t24
      t886 = -17331968.D0/3.D0*t26+t746+92672*t45-t748+8102272.D0/3.D0*t
     #49-t51+1924480*t53-8411904*t55-6144*t57-16640*t59+42752.D0/3.D0*t6
     #1
      t899 = -35584.D0/3.D0*t65+3328.D0/3.D0*t67-76288.D0/3.D0*t69+62502
     #4*t71-175104*t74+63360*t77-89856*t79+43560704.D0/3.D0*t81-5970688.
     #D0/3.D0*t83+1679616*t85+76032*t88
      t909 = -546048*t90-12982016*t92-132608*t94+2768000*t96-3072*t7-825
     #0752.D0/3.D0*t99-t775-1270912.D0/3.D0*t101+t777+217728*t105+371200
     #.D0/3.D0*t107
      t922 = 615168*t112+3000448.D0/3.D0*t114-8180736*t116+248320*t73-22
     #80960*t119-t788+8963712*t121-213504*t123-110592*t125+19456*t127-68
     #0704*t129
      t928 = 512*t18
      t934 = -1132288.D0/3.D0*t15+340736*t12-156672*t134+75776*t136+2739
     #2.D0/3.D0*t1+t928+173056.D0/3.D0*t3+5756032.D0/3.D0*t141+41217152.
     #D0/3.D0*t143-1804032*t145-2868480*t147
      t946 = 2937600*t151-2560.D0/3.D0*t11-12544*t2-t693+5069312*t157+13
     #58336*t159-183808*t76-33782528.D0/3.D0*t162-3766016.D0/3.D0*t164+5
     #068928.D0/3.D0*t166+34363520.D0/3.D0*t168
      t956 = 17101952.D0/3.D0*t171+513408*t173+12288*t175-18624512.D0/3.
     #D0*t177+t823-14336.D0/3.D0*t181+7815680*t183-1723904*t185-4411904*
     #t187-15250432*t189-t829
      t965 = 1024*t3
      t967 = -546560.D0/3.D0*t105+34048*t107-2570752.D0/3.D0*t114-228864
     #*t116-253952.D0/3.D0*t119-1536*t8-t788+1739648*t121+303744*t123+13
     #73440.D0/3.D0*t125+(t877+t886+t899+t909+t922+t934+t946+t956)*s1l1-
     #5888*t127-1348352.D0/3.D0*t129+1090688.D0/3.D0*t134+156416.D0/3.D0
     #*t136+t799+t965+2796032.D0/3.D0*t141
      t970 = 3072*t2
      t986 = 1137920*t143-6400.D0/3.D0*t11-t970+17920*t155+168448*t157-1
     #39776*t159-517632*t162+1677568*t164+1176320.D0/3.D0*t166+2844416*t
     #168-6314752.D0/3.D0*t171-3396224.D0/3.D0*t173-322816.D0/3.D0*t175-
     #5792000.D0/3.D0*t177-t824+49152*t183-40960*t185+155648*t187-723929
     #6.D0/3.D0*t189
      t991 = -419840*t55+704000.D0/3.D0*t71-t713+386048.D0/3.D0*t77+2854
     #912*t81+33280.D0/3.D0*t288-399360*t83-1351680*t92-25600*t94+71680*
     #t96-426496.D0/3.D0*t99-99840*t101+66560*t105+35840.D0/3.D0*t107-t7
     #25*(t842+t859+t967+t986)/6
      t1013 = 1290240*t35
      t1014 = 1290240*t37
      t1015 = 774144*t39
      t1017 = -1024-3775360.D0/3.D0*t4-24872320.D0/3.D0*t9-71995904.D0/3
     #.D0*t13+22216448*t16-337920*t20+150784.D0/3.D0*t22-t1013+t1014-t10
     #15+39864448.D0/3.D0*t24
      t1019 = 36864*t34
      t1021 = 36864*t47
      t1023 = 1024*t6
      t1029 = -35483392.D0/3.D0*t26+t1019+839168.D0/3.D0*t45-t1021+56455
     #68*t49-t1023+25040512.D0/3.D0*t53-18047744*t55-5120*t57-25472*t59+
     #53248.D0/3.D0*t61
      t1042 = -118912.D0/3.D0*t65+18944.D0/3.D0*t67-46464*t69+2512768.D0
     #/3.D0*t71-688128*t74+835584*t77-202752*t79+112702592.D0/3.D0*t81+6
     #1440*t288-23256832.D0/3.D0*t83+4810752*t85
      t1050 = 1024*z
      t1052 = 774144*t103
      t1054 = 239616*t88-1640448*t90-35003136*t92-115968*t94+5438208*t96
     #-2560*t7-15686912.D0/3.D0*t99-t1050-2339072.D0/3.D0*t101+t1052+114
     #8928*t105
      t1063 = 10240.D0/3.D0*l1
      t1068 = 817280.D0/3.D0*t107+1456128*t112+7445120.D0/3.D0*t114-2046
     #5664*t116+490496*t73-12591488.D0/3.D0*t119-t1063+55110272.D0/3.D0*
     #t121-439168*t123-716288*t125+52352.D0/3.D0*t127
      t1074 = 27136.D0/3.D0*t1
      t1080 = -4779008.D0/3.D0*t129-2060032.D0/3.D0*t15+601984*t12-12226
     #56*t134+135552*t136+t1074+2560*t18+87168*t3+11176192.D0/3.D0*t141+
     #98058752.D0/3.D0*t143-4478976*t145
      t1093 = -7833600*t147+7649280*t151-13568.D0/3.D0*t11-24448.D0/3.D0
     #*t2-2560*t155+11653120*t157+4652032*t159-958336.D0/3.D0*t76-818919
     #68.D0/3.D0*t162-5928192*t164+9827456.D0/3.D0*t166
      t1099 = 258048*t179
      t1100 = 10240.D0/3.D0*t181
      t1105 = 258048*t191
      t1106 = 73864192.D0/3.D0*t168+17378944*t171+8068352.D0/3.D0*t173+2
     #80192.D0/3.D0*t175-12566784*t177+t1099-t1100+21452800*t183-3664384
     #*t185-13430272*t187-103308544.D0/3.D0*t189-t1105
      t1118 = -1024-34176*t4+2528896.D0/3.D0*t9+2150656.D0/3.D0*t13-3238
     #40*t16+233984.D0/3.D0*t20-97024*t22-314368*t24-1053184.D0/3.D0*t26
     #+(t1017+t1029+t1042+t1054+t1068+t1080+t1093+t1106)*s2l1-131200.D0/
     #3.D0*t45+258432*t49-t1023-1686016.D0/3.D0*t53-6144*t362-291328*t55
     #+21760*t57+4480*t59
      t1136 = 28672.D0/3.D0*t61+18176*t65-7168.D0/3.D0*t67-9856.D0/3.D0*
     #t69-764672.D0/3.D0*t71+t713-408320.D0/3.D0*t77+12544*t375+26368*t8
     #1-94720.D0/3.D0*t288+266752*t83-102400.D0/3.D0*t92+312448.D0/3.D0*
     #t94+214784.D0/3.D0*t96-3072*t383+1536*t7-409088.D0/3.D0*t99-t1050+
     #88832*t101
      t1158 = -1024-9183104.D0/3.D0*t4-5211008*t9-53535488.D0/3.D0*t13+1
     #7807616*t16-61440*t20+60928*t22-t1013+t1014-t1015+34517504.D0/3.D0
     #*t24
      t1168 = -40774144.D0/3.D0*t26+t1019+2141056.D0/3.D0*t45-t1021+2679
     #7952.D0/3.D0*t49-3584*t6+14935424.D0/3.D0*t53-22454016*t55-12288*t
     #57-173440.D0/3.D0*t59+48384*t61
      t1181 = -463616.D0/3.D0*t65+28800*t67-664960.D0/3.D0*t69+2677760.D
     #0/3.D0*t71-485376*t74+276480*t77-239616*t79+102158720.D0/3.D0*t81-
     #15942400.D0/3.D0*t83+4478976*t85+202752*t88
      t1191 = -1456128*t90-32532224*t92-209536*t94+7890176*t96-6144*t7-2
     #5980416.D0/3.D0*t99-t1050-3375104.D0/3.D0*t101+t1052+663552*t105+2
     #903168.D0/3.D0*t107
      t1204 = 1640448*t112+4343168.D0/3.D0*t114-21479424*t116+693248*t73
     #-3702016*t119-t1063+73541504.D0/3.D0*t121-852608.D0/3.D0*t123-1964
     #80*t125+45056*t127-2965376*t129
      t1216 = -3531520.D0/3.D0*t15+3865984.D0/3.D0*t12-559104*t134+26432
     #0*t136+19456.D0/3.D0*t1+5632*t18+1219840.D0/3.D0*t3+4546048*t141+1
     #05348608.D0/3.D0*t143-4810752*t145-7649280*t147
      t1229 = 7833600*t151+20096.D0/3.D0*t11-289024.D0/3.D0*t2-6144*t155
     #+13477888*t157+3638272*t159-2764288.D0/3.D0*t76-72364544.D0/3.D0*t
     #162-8995712.D0/3.D0*t164+19024256.D0/3.D0*t166+27542528*t168
      t1239 = 35802752.D0/3.D0*t171+3225472.D0/3.D0*t173+21504*t175-5415
     #9488.D0/3.D0*t177+t1099-25600.D0/3.D0*t181+20439040*t183-4678144*t
     #185-11605504*t187-37769472*t189-t1105
      t1245 = -199936.D0/3.D0*t105-15488.D0/3.D0*t107-1680512.D0/3.D0*t1
     #14-75776*t116+1123840.D0/3.D0*t119-3072*t8-t1063+1736960.D0/3.D0*t
     #121+836480.D0/3.D0*t123+252288*t125-70528.D0/3.D0*t127-8064*t129+4
     #00640.D0/3.D0*t134-110720.D0/3.D0*t136+t1074-t928+(t1158+t1168+t11
     #81+t1191+t1204+t1216+t1229+t1239)*s1l1-t965-28928.D0/3.D0*t141
      t1260 = 14336*t185
      t1263 = 377600*t143-21248.D0/3.D0*t11+t970+32768*t155+58368*t157-5
     #3248*t159+339712.D0/3.D0*t162+656000*t164+348544.D0/3.D0*t166+7587
     #84*t168-2820736.D0/3.D0*t171-1451648.D0/3.D0*t173-102528*t175-1703
     #552.D0/3.D0*t177-t1100+7168*t183-t1260+64512*t187-2272256.D0/3.D0*
     #t189
      t1282 = -637696.D0/3.D0*t4+46592.D0/3.D0*t9+240640.D0/3.D0*t13-453
     #12*t16+t33+369920.D0/3.D0*t24-403712*t26+220928*t45+1331968.D0/3.D
     #0*t49-t51-25088*t53-78336*t55-t58-t60+t62-85504*t65
      t1294 = 2048*t116
      t1295 = 2048*t73
      t1297 = t68-247040.D0/3.D0*t69+432640.D0/3.D0*t71-37632*t81+13824*
     #t83+36864*t92-68608*t94+49920*t96-t98-128768*t99-266240*t101+48179
     #2.D0/3.D0*t107+37888*t114-t1294+t1295-178688*t119
      t1310 = 664832.D0/3.D0*t121-44544*t123-17920.D0/3.D0*t125+t128-112
     #1792.D0/3.D0*t129-11264*t15+27136*t12+346624.D0/3.D0*t136-t138+t13
     #9+34304*t3+1187072.D0/3.D0*t141+34304*t143+t153-18944*t2-t156
      t1311 = 12288*t157
      t1312 = 3072*t159
      t1321 = 11264*t183
      t1322 = 8704*t185
      t1323 = 10752*t187
      t1325 = t1311-t1312-114688.D0/3.D0*t76-73216.D0/3.D0*t162+66560.D0
     #/3.D0*t164+600832.D0/3.D0*t166+869888.D0/3.D0*t168-71424*t171-1945
     #6.D0/3.D0*t173+t176-378368*t177-t182-t1321-t1322+t1323-137472*t189
      t1343 = 9728.D0/3.D0*t59
      t1344 = 13312.D0/3.D0*t61
      t1345 = 39424.D0/3.D0*t65
      t1346 = 2533376.D0/3.D0*t26+t44-226048.D0/3.D0*t45-t48-994048*t49+
     #6231808.D0/3.D0*t53-5504256*t55+t539+t1343-t1344+t1345
      t1348 = 2560*t67
      t1349 = 21248*t69
      t1368 = 119808*t88-820224*t90-12532736*t92-31744.D0/3.D0*t94+14105
     #60*t96+t552-610816.D0/3.D0*t99+114432*t101+t104+422400*t105-354304
     #.D0/3.D0*t107
      t1379 = 5120.D0/3.D0*t127
      t1386 = 2560.D0/3.D0*t1
      t1392 = -133632*t15-180736.D0/3.D0*t12-274944*t134-57344.D0/3.D0*t
     #136-t1386-t928-53248*t3-1530880.D0/3.D0*t141+11011584*t143-2239488
     #*t145-3916800*t147
      t1395 = 4352.D0/3.D0*t11
      t1396 = 11776.D0/3.D0*t2
      t1413 = 1111040*t168+2167296*t171-21248*t173-10752*t175+540672*t17
     #7+t180+9717760*t183-1549312*t185-6144000*t187-6601472*t189-t192
      t1419 = 2252800.D0/3.D0*t114-205824*t116-970240*t119+2559488.D0/3.
     #D0*t121-561920.D0/3.D0*t123-212480*t125-t244*(t1118+t1136+t1245+t1
     #263)-440320.D0/3.D0*t129-396032.D0/3.D0*t134+5120*t136+558592*t141
     #+(t1282+t1297+t1310+t1325)*s1l1+1029120*t143+(-15360*t288+1264640.
     #D0/3.D0*t4+3824640*t151+t1346-582656.D0/3.D0*t9-18006784.D0/3.D0*t
     #13-t1348+8219136*t16+5082624*t157+t1349+2138624*t159-8256512.D0/3.
     #D0*t162-661504*t164+10240.D0/3.D0*t22+38400*t20-1730560.D0/3.D0*t1
     #66-111104.D0/3.D0*t24+t1368-t36+t38-t40-t1379+t1392+1195520.D0/3.D
     #0*t129+t1395+198656*t73-12288*t71+44544*t77+126976*t76-316416*t74+
     #26021632.D0/3.D0*t81-101376*t79-2893056*t83+t1396+2405376*t85+t141
     #3+t573+728064*t112-9127936*t116-64512*t114+222208*t119+184064.D0/3
     #.D0*t125-12544.D0/3.D0*t123+7174144.D0/3.D0*t121)*s2l1+83968*t157
      t1470 = -220160.D0/3.D0*t4+557312.D0/3.D0*t9+468736.D0/3.D0*t13-62
     #208*t16+2560.D0/3.D0*t22-304384*t24+144640*t26+98048.D0/3.D0*t45+1
     #81504.D0/3.D0*t49-38912*t53-44544*t55-t539-t1343+t1344-t1345
      t1482 = t1348-t1349-160256.D0/3.D0*t71-212736*t81+16896*t83+75776*
     #t92+16384.D0/3.D0*t94+37120*t96-t552-71168*t99+20224.D0/3.D0*t101+
     #48128*t107-269824.D0/3.D0*t114-t1294+t1295+512000.D0/3.D0*t119
      t1494 = 208640.D0/3.D0*t121+19968*t123+10752*t125+t1379-51200*t129
     #-9216*t15+17920*t12+7168*t136+t1386+t928+41984.D0/3.D0*t3-67328*t1
     #41-13824*t143-t1395-t1396
      t1505 = -t573+t1311-t1312-60928.D0/3.D0*t76+328704*t162+143872.D0/
     #3.D0*t164+84736*t166-152064*t168-220416*t171-32768*t173-1024.D0/3.
     #D0*t175-68608*t177-t1321-t1322+t1323+79104*t189
      t1510 = 79872*t159-3309568*t162-927744*t164+415744.D0/3.D0*t166+23
     #82848*t168+7442432.D0/3.D0*t171+1788928.D0/3.D0*t173+97280.D0/3.D0
     #*t175-2680832.D0/3.D0*t177+270336*t183-t1260-200704*t187-256.D0/3.
     #D0*(t199-t200+t201-t202-t725*(-t202-15-50*l1+t199-25*t11+110*t1)/6
     #+t245*(72*t2-172*t1-24*t3+24+80*l1+26*t11)/2-t244*(106*t1-83*t11-t
     #199-40*l1+t202-12)+t525*(3+12*t11-24*t1+10*l1)/24)*t241+(t1470+t14
     #82+t1494+t1505)*s2l1-2139392*t189
      t1514 = (l1+z-t181)**2
      taux10 = (t710+t991+t1419+t1510)/t1514*t214

CCC Coefficients

      tcoff15 = taux9
      tcoff16 = taux9
      tcoff17 = taux10
      tcoff18 = -taux10

CCC Jet functions

      tFjet15 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff15)
      tFjet16 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff16)
      tFjet17 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff17)
      tFjet18 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff18)

CCC Sum

      RVnum5e4 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet18
     #*tcoff18
      RETURN
      END

      double precision function RVnum6e0(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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
      t2 = t1*z
      t4 = l1**2
      t5 = 128.D0/3.D0*t4
      t6 = 128.D0/3.D0*t1
      t7 = z*l1
      t8 = 64*t7
      t9 = t4*l1
      t10 = t4**2
      t11 = t10*t9
      t12 = t1*t11
      t14 = t1**2
      t15 = t14*z
      t16 = t15*t11
      t18 = t15*t4
      t20 = 128.D0/3.D0*l1
      t21 = t1*l1
      t23 = t2*l1
      t25 = z*t4
      t27 = t1*t4
      t29 = t2*t4
      t31 = z*t9
      t33 = t1*t9
      t35 = t14*t9
      t37 = t2*t9
      t39 = t14*t4
      t41 = -128.D0/3.D0-320.D0/3.D0*t2+t5+t6-t8-256*t12+64*t16-128*t18-
     #t20-1408.D0/3.D0*t21+2368.D0/3.D0*t23-448.D0/3.D0*t25+5504.D0/3.D0
     #*t27-8512.D0/3.D0*t29+768*t31-3968*t33-9472.D0/3.D0*t35+17600.D0/3
     #.D0*t37+1280*t39
      t42 = t15*t9
      t44 = t15*t10
      t46 = t14*t10
      t48 = z*t10
      t50 = t2*t10
      t52 = t1*t10
      t54 = t14*l1
      t56 = t10*t4
      t57 = t56*t14
      t59 = t56*z
      t61 = t2*t56
      t63 = t56*t15
      t65 = t56*t1
      t67 = t10*l1
      t68 = t1*t67
      t70 = t2*t67
      t72 = t14*t67
      t74 = t15*t67
      t76 = z*t67
      t78 = t2*t11
      t80 = t14*t11
      t82 = z*t11
      t84 = 512*t42-2624.D0/3.D0*t44+12928.D0/3.D0*t46-1216*t48-7296*t50
     #+15232.D0/3.D0*t52-256*t54+1408*t57-384*t59-2176*t61-1024.D0/3.D0*
     #t63+4480.D0/3.D0*t65-11264.D0/3.D0*t68+16192.D0/3.D0*t70-10112.D0/
     #3.D0*t72+768*t74+960*t76+384*t78-256*t80+64*t82
      t86 = s1l1**2
      t88 = t14*t1
      t89 = t88*t56
      t94 = 2*l6-1
      t95 = t94**2
      t98 = l5-1
      t100 = -1+l1
      t102 = (1-l1+t7)**2
      t103 = 1/t102
      t108 = dsqrt(z)
      t110 = dsqrt(l5)
      t112 = dsqrt(-t98)
      t115 = dsqrt(l1)
      t116 = l5*z
      t117 = t116*l1
      t118 = l5*l1
      t119 = t117-l1+1-l5+t118
      t121 = -t100
      t122 = dsqrt(t121)
      t127 = t119**2
      t130 = dsqrt(4*l1*t95*z*l5*t98*t100*t103-4*t94*t108*t110*t112*t115
     #*t119*t103*t122+t127*t103)
      t133 = 64.D0/3.D0*t1
      t137 = t15*l1
      t139 = t88*t10
      t142 = -448.D0/3.D0-896*t89+128.D0/3.D0*t2-64*z+256.D0/3.D0*t130+6
     #4*t4+t133-t8+3712.D0/3.D0*t12-4288.D0/3.D0*t16+1536*t18-256*t137-3
     #968.D0/3.D0*t139-64.D0/3.D0*l1
      t146 = t88*t4
      t158 = -704.D0/3.D0*t21-1600.D0/3.D0*t23+128*t27-128*t146+7040.D0/
     #3.D0*t29+128*t31+3136.D0/3.D0*t33+8768*t35-18688.D0/3.D0*t37-11072
     #.D0/3.D0*t39-13120.D0/3.D0*t42+21440.D0/3.D0*t44-13312*t46+320.D0/
     #3.D0*t48
      t174 = 31232.D0/3.D0*t50-9088.D0/3.D0*t52+896*t54-7872*t57+512*t59
     #+7104*t61+12736.D0/3.D0*t63-9280.D0/3.D0*t65+12352.D0/3.D0*t68-109
     #44*t70+38848.D0/3.D0*t72-21248.D0/3.D0*t74-1472.D0/3.D0*t76-7808.D
     #0/3.D0*t78
      t176 = t130*t4
      t178 = t130*l1
      t180 = t10**2
      t181 = t2*t180
      t183 = t14*t180
      t185 = t15*t180
      t187 = 128.D0/3.D0*z
      t202 = -128-t187+t5+t6-128.D0/3.D0*t7-512.D0/3.D0*t12+128.D0/3.D0*
     #t16-1024.D0/3.D0*t21+128*t23-256.D0/3.D0*t25+3200.D0/3.D0*t27-1024
     #*t29+512*t31-6784.D0/3.D0*t33-1152*t35+8192.D0/3.D0*t37+256*t39+12
     #8*t42
      t221 = -384*t44+6400.D0/3.D0*t46-2432.D0/3.D0*t48-12032.D0/3.D0*t5
     #0+3072*t52+896*t57-256*t59-1408*t61-640.D0/3.D0*t63+2944.D0/3.D0*t
     #65-7168.D0/3.D0*t68+9856.D0/3.D0*t70-5888.D0/3.D0*t72+1280.D0/3.D0
     #*t74+640*t76+256*t78-512.D0/3.D0*t80+128.D0/3.D0*t82
      t224 = t88*t67
      t226 = t88*t9
      t228 = z*t180
      t230 = t1*t180
      t232 = t88*t11
      t234 = t88*t180
      t238 = -128*t14-128.D0/3.D0*t176+64.D0/3.D0*t178+1280.D0/3.D0*t181
     #-1280.D0/3.D0*t183+640.D0/3.D0*t185+(t202+t221)*s2l1+4352.D0/3.D0*
     #t224+640*t226+128.D0/3.D0*t228-640.D0/3.D0*t230+896.D0/3.D0*t232-1
     #28.D0/3.D0*t234+8192.D0/3.D0*t80-704.D0/3.D0*t82
      t242 = 64*t2
      t243 = 64.D0/3.D0*t7
      t244 = 256.D0/3.D0*t12
      t245 = 64.D0/3.D0*t16
      t246 = 64*t18
      t257 = 256*t42
      t258 = -256.D0/3.D0-t242-t187+t243-t244+t245-t246+t20-256*t21+384*
     #t23-128*t25+896*t27-3904.D0/3.D0*t29+384*t31-5248.D0/3.D0*t33-4352
     #.D0/3.D0*t35+7616.D0/3.D0*t37+640*t39+t257
      t259 = 1216.D0/3.D0*t44
      t264 = 128*t54
      t268 = 128*t63
      t273 = 320*t74
      t275 = 128*t78
      t276 = 256.D0/3.D0*t80
      t277 = 64.D0/3.D0*t82
      t278 = -t259+5504.D0/3.D0*t46-1472.D0/3.D0*t48-2944*t50+6016.D0/3.
     #D0*t52-t264+512*t57-128*t59-768*t61-t268+512*t65-4096.D0/3.D0*t68+
     #6080.D0/3.D0*t70-3968.D0/3.D0*t72+t273+1024.D0/3.D0*t76+t275-t276+
     #t277
      t280 = s2l1**2
      t297 = 64*l1-832.D0/3.D0*t21+256*t23-320.D0/3.D0*t25+2944.D0/3.D0*
     #t27-3520.D0/3.D0*t29+448*t31-5824.D0/3.D0*t33-4160.D0/3.D0*t35+256
     #0*t37+576*t39
      t313 = -t268+576*t65-1600*t68+2240*t70-4160.D0/3.D0*t72+t273+1280.
     #D0/3.D0*t76+t275-64*t178-t276+t277
      t318 = t110*l5
      t319 = t108*z*t318
      t320 = t115*l1
      t322 = t112*t320*t122
      t325 = t108*t318
      t326 = t325*t112
      t328 = t320*t122*l6
      t332 = t108*t15
      t333 = t332*t318
      t334 = t333*t112
      t336 = t115*t122*l6
      t339 = l5**2
      t340 = t339*t14
      t343 = t339*t2
      t345 = t88*t339
      t348 = t14*t2
      t349 = t348*t339
      t353 = t108*t110
      t354 = t353*t112
      t358 = 7+8*t319*t322-24*t326*t328-4870*t89+8*t334*t336+4*t340*t9-t
     #343*t9-6*t345*l1+t349*t4-8*t333*t322+24*t354*t328+2*t339
      t361 = t180*t4
      t364 = t339*t15
      t367 = l5*t4
      t375 = t88*l5
      t376 = l1*l6
      t381 = t108*t14*t318
      t382 = t381*t112
      t386 = t112*t115*t122
      t389 = 2*t2+3*z+224*t361*t348-6*t364*t9-2*t367*t14+t364*l1-2*t345*
     #t4+6*t116*t9-16*t375*t376-7*t130+8*t382*t336-4*t333*t386
      t394 = t130*t1
      t396 = t339*z
      t407 = t180*l1
      t412 = -4*t381*t386+10*t4+t1+t394*t10-12*t396*l1+4*t118*t14-4*l5-8
     #*t325*t386+2*t340-2*t14*l5-29*t407*z-2055*t407*t15
      t419 = t180*t9
      t420 = t14**2
      t423 = t420*z
      t436 = 1648*t407*t14+1634*t407*t88+234*t407*t1+8*t419*t420-t419*t4
     #23+12*t420*t9-84*t420*t10+3*t423*t10-8*t367+4*t7+753*t12-6508*t16+
     #201*t18
      t440 = t115*t4
      t442 = t440*t122*l6
      t463 = -31*t137+8*t326*t442-1589*t139-t339*t9+46*t423*t56+55*t423*
     #t180-65*t423*t11-468*t420*t180+604*t420*t11-504*t420*t56-18*t423*t
     #67+268*t420*t67
      t470 = t319*t112
      t475 = t332*t110
      t480 = t112*t440*t122
      t483 = t130*z
      t486 = t15*l5
      t493 = -144*t348*t9+18*t348*t4+12*t325*t322-16*t470*t328-16*l1+12*
     #t117-8*t475*t112*t336-4*t325*t480+18*t483*t4-2*t486*l1+4*t486*t4+4
     #*t345*t9
      t495 = t339*l1
      t502 = t339*t1
      t516 = -3*t495*t2-6*t340*t4+3*t343*t4-t502*t9-6*t396*t9+4*t364*t4+
     #8*t470*t442+16*t326*t336+11*t21-12*t23-7*t25-42*t27
      t535 = -7*t394*t4-2*t483*t9-108*t146+4*t339*t4-5*t495-28*t419*t2+5
     #6*t419*t14-70*t419*t15+t343+94*t29-16*t31+137*t33+708*t35
      t551 = -386*t37-164*t39-832*t42+2343*t44-1995*t46+50*t48+1047*t50-
     #335*t52+22*t54-5389*t57+101*t59+2724*t61
      t567 = 6524*t63-785*t65+598*t68-1996*t70+3891*t72-4636*t74-83*t76-
     #2657*t78-t419*z+8*t419*t1-448*t361*t88-12*t353*t322
      t586 = t108*t88*t318
      t595 = -8*t382*t328-4*t319*t480+16*t334*t328-16*t354*t336+561*t348
     #*t10-1347*t348*t67-64*t361*t420+228*t407*t420+4*t586*t322+1713*t34
     #8*t180-2341*t348*t11+2*t394*t9
      t599 = l6**2
      t600 = l1*t599
      t605 = t4*l6
      t610 = t4*t599
      t625 = t130*t10-8*t354*t442+32*t116*t600-32*t396*t600-48*t396*t605
     #+32*t396*t376+48*t396*t610-32*t116*t376+48*t116*t605-2*t483*t10-18
     #*t116*t4+2*t502*t4+18*t396*t4
      t634 = t9*t599
      t651 = -t349*t9+6*t375*l1-2*t375*t4-10*t176-16*t396*t634+16*t116*t
     #634+4*t381*t322+16*t345*t376+16*t375*t600-16*t345*t600-48*t116*t61
     #0-64*t361*t1
      t668 = 2153*t348*t56+2*l5*t9+4*t475*t386+56*t419*t88-28*t419*t348+
     #10*t118+16*t178+1813*t181-3630*t183+4500*t185+4*t15-t10
      t692 = -448*t361*t14+560*t361*t15+12*t88*l1-8*t586*t112*t328+224*t
     #361*t2+t364-2*t486+8*t423*t361-28*t423*t407-823*t407*t2+8*z*t361+8
     #*t353*t386
      t705 = t9*l6
      t712 = 3323*t224+520*t226+63*t228-513*t230+5005*t232-3533*t234+530
     #2*t80-14*t483*l1-93*t82-809*t407*t348-16*t116*t705+16*t396*t705+4*
     #t353*t480
      t720 = zqe0(t121,z,l3)
      t723 = (l1+z-t7)**2
      taux11 = -((t41+t84)*t86+(t142+t158+t174+t238)*s1l1+(t258+t278)*t2
     #80+(-192-t242-320.D0/3.D0*z+128*t130+64.D0/3.D0*t4-t133+t243-t244+
     #t245-t246+t297+t257-t259+5632.D0/3.D0*t46-1856.D0/3.D0*t48-3200*t5
     #0+7040.D0/3.D0*t52-t264+1600.D0/3.D0*t57-448.D0/3.D0*t59-832*t61+t
     #313)*s2l1-64.D0/3.D0*(t358+t389+t412+t436+t463+t493+t516+t535+t551
     #+t567+t595+t625+t651+t668+t692+t712)*t103)*t720/t723*t103
      t2 = l1**2
      t3 = 3-3*l1+t2
      t4 = z**2
      t5 = t4**2
      t6 = t5*t4
      t7 = t6*t2
      t9 = t2**2
      t10 = t2*t9
      t11 = t10*t5
      t13 = t6*t10
      t15 = z*t10
      t17 = z*t4
      t18 = t17*t10
      t20 = t5*z
      t21 = t10*t20
      t23 = t9*l1
      t24 = t6*t23
      t26 = t2*l1
      t27 = t6*t26
      t29 = t10*t4
      t31 = t4*t23
      t33 = t17*t23
      t35 = t5*t23
      t37 = t20*t23
      t39 = z*t23
      t41 = t20*t2
      t43 = 2*z
      t44 = t6*t9
      t46 = t17*l1
      t48 = 2+4*t7+40*t11+4*t13-4*t15-40*t18-20*t21-16*t24-16*t27+20*t29
     #-72*t31+148*t33-152*t35+78*t37+14*t39-40*t41+t43+24*t44+18*t46
      t49 = z*t2
      t50 = 2*t49
      t51 = t2*t4
      t52 = 26*t51
      t53 = t4*l1
      t54 = 4*t53
      t55 = t17*t2
      t57 = z*t26
      t58 = 10*t57
      t59 = t4*t26
      t61 = t5*t26
      t63 = t17*t26
      t65 = t5*t2
      t67 = t20*l1
      t69 = t20*t26
      t71 = t20*t9
      t73 = t5*t9
      t75 = z*t9
      t77 = t17*t9
      t80 = t4*t9
      t82 = t5*l1
      t84 = 2*t4
      t85 = -t50+t52-t54-78*t55+t58-72*t59-192*t61+176*t63+90*t65+8*t67+
     #94*t69-120*t71+234*t73-18*t75-222*t77+4*t5+102*t80-24*t82+t84
      t99 = 2-8*t31+12*t33-8*t35+2*t37+2*t39-2*t41+t43+14*t46-4*t49+t52-
     #8*t53-38*t55
      t112 = t58-40*t59-32*t61+56*t63+18*t65+6*t69-6*t71+26*t73-8*t75-42
     #*t77+30*t80-2*t17-4*t82+t84
      t116 = t9*t26
      t123 = t5*t17
      t142 = 2-30*t17*t116-2*z*t116+40*t5*t116+10*t123*t10+12*t6*t116+12
     #*t116*t4-30*t20*t116+20*t123*t9-20*t123*t23-2*t123*t116-10*t123*t2
     #6+2*t123*t2
      t156 = -24*t7-160*t11-56*t13+6*t15+110*t18+130*t21+108*t24+68*t27-
     #40*t29+44*t31-156*t33+264*t35-236*t37
      t168 = -4*t39+58*t41+t43-112*t44+6*t46-t50+14*t51-t54-6*t55+6*t57-
     #16*t59+128*t61-28*t63
      t179 = 2*l6-1
      t180 = t179**2
      t183 = l5-1
      t185 = -1+l1
      t186 = z*l1
      t188 = (1-l1+t186)**2
      t189 = 1/t188
      t194 = dsqrt(z)
      t196 = dsqrt(l5)
      t198 = dsqrt(-t183)
      t201 = dsqrt(l1)
      t205 = l5*z*l1-l1+1-l5+l5*l1
      t207 = -t185
      t208 = dsqrt(t207)
      t213 = t205**2
      t216 = dsqrt(4*l1*t180*z*l5*t183*t185*t189-4*t179*t194*t196*t198*t
     #201*t205*t189*t208+t213*t189)
      t221 = -42*t65-14*t67-148*t69+238*t71-238*t73-4*t75+106*t77+2*t20-
     #10*t80-2*t216+4*t6*l1+8*t82+t84
      t227 = zqe0(t207,z,l3)
      t230 = (l1+z-t186)**2
      taux12 = -(128.D0/3.D0*t3*(t48+t85)*s1l1+128.D0/3.D0*t3*(t99+t112)
     #*s2l1+128.D0/3.D0*t3*(t142+t156+t168+t221))*t227/t230*t189

CCC Coefficients

      tcoff19 = taux11
      tcoff20 = taux11
      tcoff21 = taux12
      tcoff22 = -taux12

CCC Jet functions

      tFjet19 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff19)
      tFjet20 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff20)
      tFjet21 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff21)
      tFjet22 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff22)

CCC Sum

      RVnum6e0 = tFjet19*tcoff19+tFjet20*tcoff20+tFjet21*tcoff21+tFjet22
     #*tcoff22
      RETURN
      END

      double precision function RVnum6e1(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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
      t2 = t1**2
      t3 = t2*z
      t4 = l1**2
      t5 = t4**2
      t6 = t5*l1
      t7 = t3*t6
      t9 = t1*z
      t11 = t2*t6
      t13 = 128.D0/3.D0*t4
      t14 = 128.D0/3.D0*l1
      t15 = z*t4
      t17 = t1*t4
      t19 = t1*l1
      t21 = t9*t4
      t23 = t4*l1
      t24 = z*t23
      t26 = t1*t23
      t28 = t2*t23
      t30 = t9*t23
      t32 = t2*t4
      t34 = t3*t23
      t36 = t3*t5
      t38 = t5*t23
      t39 = t9*t38
      t41 = 128.D0/3.D0*t1
      t42 = -128.D0/3.D0+768*t7-320.D0/3.D0*t9-10112.D0/3.D0*t11+t13-t14
     #-448.D0/3.D0*t15+5504.D0/3.D0*t17-1408.D0/3.D0*t19-8512.D0/3.D0*t2
     #1+768*t24-3968*t26-9472.D0/3.D0*t28+17600.D0/3.D0*t30+1280*t32+512
     #*t34-2624.D0/3.D0*t36+384*t39+t41
      t43 = z*t6
      t45 = z*t38
      t47 = t2*t38
      t49 = z*l1
      t50 = 64*t49
      t51 = t1*t38
      t53 = t3*t38
      t55 = t3*t4
      t57 = t2*t5
      t59 = z*t5
      t61 = t9*t5
      t63 = t1*t5
      t65 = t2*l1
      t67 = t5*t4
      t68 = t67*t2
      t70 = t67*z
      t72 = t9*t67
      t74 = t67*t3
      t76 = t67*t1
      t78 = t1*t6
      t80 = t9*t6
      t82 = t9*l1
      t84 = 960*t43+64*t45-256*t47-t50-256*t51+64*t53-128*t55+12928.D0/3
     #.D0*t57-1216*t59-7296*t61+15232.D0/3.D0*t63-256*t65+1408*t68-384*t
     #70-2176*t72-1024.D0/3.D0*t74+4480.D0/3.D0*t76-11264.D0/3.D0*t78+16
     #192.D0/3.D0*t80+2368.D0/3.D0*t82
      t86 = s1l1**2
      t88 = t2*t9
      t89 = t88*t67
      t91 = t3*l1
      t93 = dsqrt(z)
      t94 = t93*t3
      t95 = dsqrt(l5)
      t96 = t94*t95
      t97 = 1-l5
      t98 = dsqrt(t97)
      t100 = dsqrt(l1)
      t101 = 1-l1
      t102 = dsqrt(t101)
      t104 = t100*t102*l6
      t105 = t96*t98*t104
      t107 = t2*t1
      t108 = t107*t38
      t111 = t95*l5
      t112 = t93*t2*t111
      t113 = t100*l1
      t115 = t98*t113*t102
      t116 = t112*t115
      t118 = l5**2
      t119 = t118*t1
      t120 = t119*t4
      t122 = t118*z
      t123 = t122*t4
      t126 = 2*l6-1
      t127 = t126**2
      t134 = (1-l1+t49)**2
      t135 = 1/t134
      t143 = l5*z
      t144 = t143*l1
      t145 = l5*l1
      t146 = t144-l1+1-l5+t145
      t152 = t146**2
      t155 = dsqrt(4*l1*t127*z*l5*t97*t101*t135-4*t126*t93*t95*t98*t100*
     #t146*t135*t102+t152*t135)
      t156 = t155*z
      t157 = t156*t23
      t161 = t93*t95
      t162 = t161*t98
      t164 = t113*t102*l6
      t165 = t162*t164
      t168 = t98*t100*t102
      t169 = t161*t168
      t171 = 7+2153*t89-31*t91-8*t105+5005*t108+4*t116+2*t120+18*t123-2*
     #t157-2*t156*t5+24*t165+8*t169
      t173 = t93*z*t111
      t174 = t173*t98
      t175 = t100*t4
      t177 = t175*t102*l6
      t178 = t174*t177
      t180 = t5**2
      t181 = t180*t23
      t182 = t181*t1
      t185 = t98*t175*t102
      t186 = t173*t185
      t189 = t180*t4
      t190 = t2**2
      t191 = t189*t190
      t193 = t88*t23
      t195 = t88*t4
      t197 = t107*t5
      t199 = t107*t4
      t201 = t107*t67
      t203 = t107*t6
      t205 = t107*t23
      t207 = 8*t178+8*t182-4*t186-4636*t7-64*t191-144*t193+18*t195-1589*
     #t197-108*t199-4870*t201+3323*t203+520*t205
      t210 = t88*t38
      t212 = t93*t111
      t213 = t212*t115
      t216 = t162*t177
      t218 = t88*t5
      t220 = t88*t6
      t222 = t118*t23
      t223 = t190*z
      t224 = t223*t6
      t226 = t190*t6
      t228 = t189*t2
      t230 = t189*t1
      t232 = 2*t9-2341*t210+12*t213+3*z-8*t216+561*t218-1347*t220-t222-1
     #8*t224+268*t226-448*t228-64*t230
      t234 = t180*l1
      t235 = t234*z
      t237 = t212*t98
      t238 = t237*t104
      t250 = 3891*t11-29*t235+16*t238+10*t4-16*l1-7*t15-42*t17+11*t19+94
     #*t21-16*t24+137*t26+708*t28-386*t30
      t256 = t174*t164
      t258 = t3*l5
      t259 = t258*l1
      t261 = t258*t4
      t263 = t143*t23
      t266 = t107*l5
      t267 = l1*l6
      t268 = t266*t267
      t270 = t107*t118
      t271 = t270*t267
      t273 = l6**2
      t274 = l1*t273
      t275 = t266*t274
      t277 = -164*t32-832*t34+2343*t36-16*t256-2*t259+4*t261+6*t263-2657
     #*t39+t1-16*t268+16*t271+16*t275
      t278 = l5*t23
      t280 = t155*t1
      t281 = t280*t4
      t283 = t4*l6
      t284 = t143*t283
      t286 = t234*t190
      t288 = z*t180
      t291 = t189*t107
      t294 = t189*t88
      t296 = t156*l1
      t298 = t94*t111
      t299 = t298*t98
      t300 = t299*t164
      t303 = 2*t278-7*t281+48*t284+228*t286+63*t288+10*t145-448*t291-83*
     #t43+224*t294-14*t296+16*t300-93*t45
      t310 = t118*t2
      t311 = t310*t4
      t313 = t118*t9
      t314 = t313*t4
      t316 = t313*t23
      t317 = t119*t23
      t318 = t122*t23
      t320 = t310*t23
      t322 = t161*t185
      t324 = 5302*t47+4*t49+753*t51-6508*t53+201*t55-6*t311+3*t314-t316-
     #t317-6*t318+4*t320+4*t322
      t325 = t96*t168
      t327 = t234*t9
      t329 = t234*t3
      t331 = t234*t2
      t333 = t118*l1
      t335 = t118*t4
      t338 = t2*t180
      t340 = t270*t23
      t342 = t88*t118
      t343 = t342*t23
      t344 = t266*l1
      t346 = t266*t4
      t349 = 4*t325-823*t327-2055*t329+1648*t331-5*t333+4*t335+2*t118-36
     #30*t338+4*t340-t343+6*t344-2*t346+t280*t5
      t353 = t3*t118
      t355 = z*t189
      t357 = t189*t9
      t359 = t112*t98
      t360 = t359*t104
      t362 = t359*t164
      t364 = t234*t1
      t366 = t237*t177
      t368 = t237*t164
      t370 = t270*t274
      t372 = t4*t273
      t373 = t143*t372
      t375 = t122*t274
      t377 = t353-2*t258+8*t355+224*t357+8*t360-8*t362+234*t364+8*t366-2
     #4*t368-16*t370-48*t373-32*t375
      t378 = t143*t274
      t380 = t143*t267
      t382 = t122*t372
      t384 = t122*t267
      t386 = t122*t283
      t389 = t212*t168
      t391 = t234*t107
      t394 = t181*t223
      t395 = t190*t23
      t397 = t190*t5
      t399 = 32*t378-32*t380+48*t382+32*t384-48*t386+4*t3-8*t389+1634*t3
     #91-7*t155-t394+12*t395-84*t397
      t401 = t223*t5
      t403 = t223*t67
      t405 = t223*t180
      t407 = t223*t38
      t409 = t190*t180
      t411 = t190*t38
      t413 = t190*t67
      t415 = t173*t115
      t417 = t280*t23
      t419 = t3*t180
      t421 = t1*t180
      t423 = t181*t107
      t425 = 3*t401+46*t403+55*t405-65*t407-468*t409+604*t411-504*t413+8
     #*t415+2*t417+4500*t419-513*t421+56*t423
      t426 = t181*t88
      t428 = t181*t190
      t430 = t107*t180
      t433 = t2*l5
      t435 = t88*t180
      t437 = t181*z
      t438 = t23*l6
      t439 = t122*t438
      t441 = t143*t438
      t443 = t23*t273
      t444 = t122*t443
      t446 = t143*t443
      t448 = -28*t426+8*t428+t313-3533*t430+2*t310-2*t433+1713*t435-t437
     #-t5+16*t439-16*t441-16*t444+16*t446
      t451 = t298*t115
      t453 = t155*l1
      t455 = t181*t9
      t457 = t181*t2
      t459 = t189*t3
      t461 = t161*t115
      t464 = t93*t107*t111
      t465 = t464*t115
      t467 = t156*t4
      t473 = -8*t451+16*t453-28*t455+56*t457+560*t459-12*t461+4*t465+18*
     #t467-1995*t57+50*t59+1047*t61-335*t63
      t483 = t464*t98*t164
      t486 = t270*l1
      t488 = t342*t4
      t489 = 22*t65-5389*t68+101*t70+2724*t72+6524*t74-785*t76+598*t78-1
     #996*t80-8*t483+12*t144-6*t486+t488
      t491 = t122*l1
      t493 = t155*t4
      t495 = t223*t189
      t497 = t223*t234
      t499 = t181*t3
      t501 = t298*t168
      t503 = t107*l1
      t505 = t9*t180
      t507 = t112*t168
      t509 = t234*t88
      t511 = t162*t104
      t513 = l5*t4
      t515 = -12*t491-10*t493+8*t495-28*t497-70*t499-4*t501+12*t503+1813
     #*t505-4*t507-809*t509-16*t511-8*t513
      t517 = t270*t4
      t519 = t143*t4
      t521 = t145*t2
      t523 = t353*t4
      t525 = t353*t23
      t527 = t513*t2
      t529 = t353*l1
      t530 = t333*t9
      t534 = t212*t185
      t536 = t299*t104
      t538 = -4*l5-2*t517-18*t519+4*t521+4*t523-6*t525-2*t527+t529-3*t53
     #0+t155*t5-12*t82-4*t534+8*t536
      t545 = 256*t91
      t549 = 128*t199
      t553 = 128.D0/3.D0*t9
      t556 = 128.D0/3.D0*z
      t570 = -128+1280.D0/3.D0*t7-t556-5888.D0/3.D0*t11+t13-256.D0/3.D0*
     #t15+3200.D0/3.D0*t17-1024.D0/3.D0*t19-1024*t21+512*t24-6784.D0/3.D
     #0*t26-1152*t28+8192.D0/3.D0*t30+256*t32+128*t34-384*t36+256*t39+t4
     #1
      t589 = 640*t43+128.D0/3.D0*t45-512.D0/3.D0*t47-128.D0/3.D0*t49-512
     #.D0/3.D0*t51+128.D0/3.D0*t53+6400.D0/3.D0*t57-2432.D0/3.D0*t59-120
     #32.D0/3.D0*t61+3072*t63+896*t68-256*t70-1408*t72-640.D0/3.D0*t74+2
     #944.D0/3.D0*t76-7168.D0/3.D0*t78+9856.D0/3.D0*t80+128*t82
      t595 = -t545+896.D0/3.D0*t108-21248.D0/3.D0*t7-3968.D0/3.D0*t197-t
     #549-896*t201+4352.D0/3.D0*t203+640*t205+t553-64*z+(t570+t589)*s2l1
     #+38848.D0/3.D0*t11+64*t4-64.D0/3.D0*l1
      t607 = 64.D0/3.D0*t1
      t610 = 128*t17-704.D0/3.D0*t19+7040.D0/3.D0*t21+128*t24+3136.D0/3.
     #D0*t26+8768*t28-18688.D0/3.D0*t30-11072.D0/3.D0*t32-13120.D0/3.D0*
     #t34+21440.D0/3.D0*t36-7808.D0/3.D0*t39+t607+128.D0/3.D0*t288-1472.
     #D0/3.D0*t43
      t617 = 1280.D0/3.D0*t338
      t618 = 128*t2
      t624 = -448.D0/3.D0-704.D0/3.D0*t45+8192.D0/3.D0*t47-t50+3712.D0/3
     #.D0*t51-4288.D0/3.D0*t53+1536*t55-t617-t618+256.D0/3.D0*t155+640.D
     #0/3.D0*t419-640.D0/3.D0*t421-128.D0/3.D0*t430+64.D0/3.D0*t453
      t637 = 128.D0/3.D0*t493
      t638 = 1280.D0/3.D0*t505
      t640 = -13312*t57+320.D0/3.D0*t59+31232.D0/3.D0*t61-9088.D0/3.D0*t
     #63+896*t65-7872*t68+512*t70+7104*t72+12736.D0/3.D0*t74-9280.D0/3.D
     #0*t76+12352.D0/3.D0*t78-10944*t80-t637+t638-1600.D0/3.D0*t82
      t644 = 320*t7
      t645 = 64*t9
      t656 = 256*t34
      t657 = 1216.D0/3.D0*t36
      t658 = 128*t39
      t660 = -256.D0/3.D0+t644-t645-t556-3968.D0/3.D0*t11+t14-128*t15+89
     #6*t17-256*t19-3904.D0/3.D0*t21+384*t24-5248.D0/3.D0*t26-4352.D0/3.
     #D0*t28+7616.D0/3.D0*t30+640*t32+t656-t657+t658+1024.D0/3.D0*t43
      t661 = 64.D0/3.D0*t45
      t662 = 256.D0/3.D0*t47
      t663 = 64.D0/3.D0*t49
      t664 = 256.D0/3.D0*t51
      t665 = 64.D0/3.D0*t53
      t666 = 64*t55
      t671 = 128*t65
      t675 = 128*t74
      t680 = t661-t662+t663-t664+t665-t666+5504.D0/3.D0*t57-1472.D0/3.D0
     #*t59-2944*t61+6016.D0/3.D0*t63-t671+512*t68-128*t70-768*t72-t675+5
     #12*t76-4096.D0/3.D0*t78+6080.D0/3.D0*t80+384*t82
      t682 = s2l1**2
      t699 = -3520.D0/3.D0*t21+448*t24-5824.D0/3.D0*t26-4160.D0/3.D0*t28
     #+2560*t30+576*t32+t656-t657+t658-t607+1280.D0/3.D0*t43
      t715 = -3200*t61+7040.D0/3.D0*t63-t671+1600.D0/3.D0*t68-448.D0/3.D
     #0*t70-832*t72-t675+576*t76-1600*t78+2240*t80+256*t82
      t720 = zqe0(t101,z,l3)
      t721 = dlog(l1)
      t723 = sqe1(t101)
      t724 = zqe1(t101,z,l3)
      t725 = sbqe0(t101)
      t726 = dlog(z)
      t741 = -2560*t7+4352.D0/3.D0*t6+1024.D0/3.D0*t9+128*z+37504.D0/3.D
     #0*t11-512.D0/3.D0*t4-t14+5632.D0/3.D0*t15-21760.D0/3.D0*t17+5120.D
     #0/3.D0*t19+26752.D0/3.D0*t21-17536.D0/3.D0*t24
      t750 = 640.D0/3.D0*t288
      t755 = 16896*t26+9728*t28-59648.D0/3.D0*t30-3840*t32-1536*t34+2688
     #*t36-10496.D0/3.D0*t39-256.D0/3.D0*t1+t750-27136.D0/3.D0*t43+640*t
     #23+896.D0/3.D0*t38-4736.D0/3.D0*t45
      t760 = 384*t53
      t761 = 384*t55
      t762 = 640.D0/3.D0*t338
      t763 = 128.D0/3.D0*t180
      t764 = 128.D0/3.D0*t419
      t765 = 1280.D0/3.D0*t421
      t769 = 640.D0/3.D0+5504.D0/3.D0*t47-256*t49+3328*t51-t760+t761-t76
     #2-t763+t764-t765-3968.D0/3.D0*t5-14208*t57+28544.D0/3.D0*t59
      t772 = 768*t65
      t782 = 82048.D0/3.D0*t61-72064.D0/3.D0*t63+t772-19456.D0/3.D0*t68+
     #15232.D0/3.D0*t70+12160*t72+4096.D0/3.D0*t74-33664.D0/3.D0*t76+211
     #20*t78-70400.D0/3.D0*t80+t638-7040.D0/3.D0*t82-896*t67
      t796 = 640.D0/3.D0*l1
      t798 = 768*t91-2240.D0/3.D0*t108+18368*t7+3520*t197+384*t199+6784.
     #D0/3.D0*t201-11200.D0/3.D0*t203-1792*t205+t553+384*z-103808.D0/3.D
     #0*t11-t13-t796-128.D0/3.D0*t15
      t814 = 3328.D0/3.D0*t17-128.D0/3.D0*t19-20864.D0/3.D0*t21+320*t24-
     #5376*t26-23232*t28+18560*t30+9536*t32+11584*t34-18688*t36+20480.D0
     #/3.D0*t39+640.D0/3.D0*t1-320.D0/3.D0*t288+6080.D0/3.D0*t43+640*t45
      t829 = 1664.D0/3.D0-7040*t47-128*t49-9920.D0/3.D0*t51+10880.D0/3.D
     #0*t53-4224*t55+3200.D0/3.D0*t338+384*t2-896.D0/3.D0*t155-1600.D0/3
     #.D0*t419+1600.D0/3.D0*t421+320.D0/3.D0*t430+256.D0/3.D0*t453+35648
     #*t57
      t853 = 512-3328*t7+4864.D0/3.D0*t6+256*t9+1024.D0/3.D0*z+47872.D0/
     #3.D0*t11-256*l1+768*t15-4608*t17+2048.D0/3.D0*t19+19840.D0/3.D0*t2
     #1-3968*t24
      t867 = 13824*t26+27392.D0/3.D0*t28-17792*t30-3072*t32-1408*t34+908
     #8.D0/3.D0*t36-5888*t39+512.D0/3.D0*t1+1280.D0/3.D0*t288-32000.D0/3
     #.D0*t43+256*t23+512*t38-2688*t45
      t870 = 512.D0/3.D0*t49
      t880 = 3072*t47-t870+5632*t51-640*t53+256*t55-t617-256.D0/3.D0*t18
     #0+256.D0/3.D0*t419-2560.D0/3.D0*t421-1024*t5-47104.D0/3.D0*t57+883
     #2*t59
      t894 = 29184*t61-24320*t63+512*t65-9472*t68+7296*t70+17664*t72+601
     #6.D0/3.D0*t74-48640.D0/3.D0*t76+77312.D0/3.D0*t78-88064.D0/3.D0*t8
     #0+2560.D0/3.D0*t505-1408*t82-1280*t67
      t899 = -1280*t59-30400*t61+11200*t63-2304*t65+20736*t68-1600*t70-5
     #8240.D0/3.D0*t72-32704.D0/3.D0*t74+26752.D0/3.D0*t76-13120*t78+931
     #84.D0/3.D0*t80+t637-3200.D0/3.D0*t505+(t853+t867+t880+t894)*s2l1+4
     #160.D0/3.D0*t82
      t908 = 1280.D0/3.D0*t15
      t909 = 11392.D0/3.D0*t17
      t913 = -8192.D0/3.D0*t7+384*t6+384*t9+640.D0/3.D0*z+11776*t11+t13-
     #t796+t908-t909+640*t19+7424*t21-1792*t24
      t925 = 896.D0/3.D0+9728*t26+9600*t28-47872.D0/3.D0*t30-3968*t32-16
     #64*t34+2944*t36-9472.D0/3.D0*t39+128*t1+t750-13312.D0/3.D0*t43+640
     #.D0/3.D0*t38-3712.D0/3.D0*t45
      t932 = 5248.D0/3.D0*t47-t870+2816*t51-t760+t761-t762-t763+t764-t76
     #5-128*t5-40832.D0/3.D0*t57+3712*t59
      t944 = 65408.D0/3.D0*t61-14720*t63+t772-18304.D0/3.D0*t68+9472.D0/
     #3.D0*t70+30848.D0/3.D0*t72+1408*t74-8320*t76+42112.D0/3.D0*t78-570
     #88.D0/3.D0*t80+t638-2176*t82-1280.D0/3.D0*t67
      t958 = 640+t545+256.D0/3.D0*t108-5824.D0/3.D0*t7+1664.D0/3.D0*t197
     #+t549-64*t201-704.D0/3.D0*t203-448*t205+512*t9+1408.D0/3.D0*z+3020
     #8.D0/3.D0*t11+t13-384*l1
      t971 = t908-t909+1664.D0/3.D0*t19+20672.D0/3.D0*t21-1664*t24+9280*
     #t26+21760.D0/3.D0*t28-14656*t30-7936.D0/3.D0*t32+640.D0/3.D0*t34+1
     #088*t36-6016.D0/3.D0*t39+896.D0/3.D0*t1+64.D0/3.D0*t288
      t985 = -7616.D0/3.D0*t43-832.D0/3.D0*t45+4864.D0/3.D0*t47-896.D0/3
     #.D0*t49+3584.D0/3.D0*t51-1856.D0/3.D0*t53-640*t55-t762+t618-384*t1
     #55+320.D0/3.D0*t419-320.D0/3.D0*t421-64.D0/3.D0*t430+256*t453
      t1000 = -33344.D0/3.D0*t57+2880*t59+19520*t61-12928*t63+320*t65-16
     #256.D0/3.D0*t68+3584.D0/3.D0*t70+7680*t72+1536*t74-4928*t76+10560*
     #t78-47744.D0/3.D0*t80-t637+640.D0/3.D0*t505-2048*t82
      t1015 = -20-7778*t89+116*t91+16*t105-19085*t108-8*t116-4*t120-36*t
     #123+20*t157-48*t165-16*t169-16*t178
      t1028 = -32*t182+8*t186+18996*t7+244*t191+540*t193-72*t195+6006*t1
     #97+408*t199+18500*t201-12598*t203-1954*t205-16*t9
      t1042 = 8522*t210-24*t213-12*z+16*t216-2046*t218+4866*t220+2*t222+
     #66*t224-955*t226+1820*t228+268*t230-17558*t11
      t1055 = 132*t235-32*t238-40*t4+50*l1+12*t15+184*t17-28*t19-564*t21
     #+68*t24-695*t26-3472*t28+2124*t30
      t1070 = 853*t32+3432*t34-9682*t36+32*t256+4*t259-8*t261-12*t263+12
     #110*t39-8*t1+32*t268-32*t271-32*t275
      t1084 = -4*t278+20*t281-96*t284-837*t286-308*t288-20*t145-10*t155*
     #t23+1764*t291+424*t43-868*t294+10*t23+40*t296
      t1098 = -32*t300+482*t45-22660*t47-2*t49-3665*t51+26192*t53-806*t5
     #5+12*t311-6*t314+2*t316+2*t317+12*t318
      t1112 = -8*t320-8*t322-8*t325+3502*t327+8220*t329-6795*t331+10*t33
     #3-8*t335-4*t118+15215*t338-8*t340+2*t343-12*t344
      t1128 = 4*t346-2*t353+4*t258-34*t355-924*t357-16*t360+16*t362-1029
     #*t364-16*t366+48*t368+32*t370+96*t373
      t1141 = 64*t375-64*t378+64*t380-96*t382-64*t384+96*t386-14*t3+16*t
     #389-6347*t391+20*t155+4*t394-48*t395
      t1155 = 312*t397-12*t401-160*t403-190*t405+222*t407+1669*t409-2118
     #*t411+1765*t413-16*t415-10*t417-18030*t419+2375*t421
      t1168 = -224*t423+112*t426-32*t428-2*t313+13573*t430-4*t310+4*t433
     #-6330*t435+4*t437-32*t439+32*t441+32*t444
      t1183 = -32*t446+16*t451-50*t453+112*t455-224*t457-2240*t459+24*t4
     #61-8*t465-60*t467+9339*t57-232*t59-5444*t61
      t1196 = 1759*t63-119*t65+23597*t68-534*t70-12910*t72-26466*t74+398
     #6*t76-3123*t78+9882*t80+16*t483-24*t144+12*t486
      t1210 = -2*t488+24*t491+40*t493-30*t495+100*t497+280*t499+8*t501-4
     #8*t503-7974*t505+8*t507+3054*t509+32*t511
      t1224 = 16*t513+8*l5+4*t517+36*t519-8*t521-8*t523+12*t525+4*t527-2
     #*t529+6*t530+96*t82+8*t534-16*t536
      t1235 = (l1+z-t49)**2
      taux11 = (((t42+t84)*t86-64.D0/3.D0*(t171+t207+t232+t250+t277+t303
     #+t324+t349+t377+t399+t425+t448+t473+t489+t515+t538)*t135+(t595+t61
     #0+t624+t640)*s1l1+(t660+t680)*t682+(-192+t644-t645-320.D0/3.D0*z-4
     #160.D0/3.D0*t11+64.D0/3.D0*t4+64*l1-320.D0/3.D0*t15+2944.D0/3.D0*t
     #17-832.D0/3.D0*t19+t699+t661-t662+t663-t664+t665-t666+128*t155-64*
     #t453+5632.D0/3.D0*t57-1856.D0/3.D0*t59+t715)*s2l1)*(t720*t721+t723
     #-t724+t725*t726)-((t741+t755+t769+t782)*t86+(t798+t814+t829+t899)*
     #s1l1+(t913+t925+t932+t944)*t682+(t958+t971+t985+t1000)*s2l1-64.D0/
     #3.D0*(t1015+t1028+t1042+t1055+t1070+t1084+t1098+t1112+t1128+t1141+
     #t1155+t1168+t1183+t1196+t1210+t1224)*t135)*t720)/t1235*t135
      t2 = l1**2
      t3 = 3-3*l1+t2
      t4 = z**2
      t5 = t4*z
      t6 = t5*l1
      t8 = z*t2
      t9 = 2*t8
      t10 = t4*t2
      t11 = 26*t10
      t12 = t4*l1
      t13 = 4*t12
      t14 = t5*t2
      t16 = t2*l1
      t17 = z*t16
      t18 = 10*t17
      t19 = t4*t16
      t21 = t4**2
      t22 = t21*t16
      t24 = t21*t2
      t26 = t5*t16
      t28 = 2*t4
      t29 = 4*t21
      t30 = t21*z
      t31 = t30*l1
      t32 = 8*t31
      t33 = t30*t16
      t35 = t2**2
      t36 = t30*t35
      t38 = t21*t35
      t40 = z*t35
      t42 = t5*t35
      t44 = 2+18*t6-t9+t11-t13-78*t14+t18-72*t19-192*t22+90*t24+176*t26+
     #t28+t29+t32+94*t33-120*t36+234*t38-18*t40-222*t42
      t45 = t4*t35
      t47 = t21*l1
      t49 = t21*t4
      t50 = t49*t35
      t52 = t49*t2
      t53 = 4*t52
      t54 = t35*t2
      t55 = t54*t21
      t57 = t49*t54
      t59 = t54*z
      t61 = t5*t54
      t63 = t54*t30
      t65 = t35*l1
      t66 = t49*t65
      t68 = t49*t16
      t70 = t54*t4
      t72 = t4*t65
      t74 = t5*t65
      t76 = 2*z
      t77 = t21*t65
      t79 = t30*t65
      t81 = z*t65
      t83 = t30*t2
      t85 = 102*t45-24*t47+24*t50+t53+40*t55+4*t57-4*t59-40*t61-20*t63-1
     #6*t66-16*t68+20*t70-72*t72+148*t74+t76-152*t77+78*t79+14*t81-40*t8
     #3
      t86 = t44+t85
      t98 = 2+14*t6-4*t8+t11-8*t12-38*t14+t18-40*t19-32*t22+18*t24-2*t5+
     #56*t26+t28
      t112 = 6*t33-6*t36+26*t38-8*t40-42*t42+30*t45-4*t47-8*t72+12*t74+t
     #76-8*t77+2*t79+2*t81-2*t83
      t113 = t98+t112
      t116 = 6*t6
      t117 = 14*t10
      t119 = 6*t17
      t122 = t35*t16
      t123 = t5*t122
      t124 = 30*t123
      t125 = z*t122
      t126 = 2*t125
      t127 = t21*t122
      t128 = 40*t127
      t129 = t21*t5
      t130 = t129*t54
      t132 = 2+t116-t9+t117-t13-6*t14+t119-16*t19+128*t22-t124-t126+t128
     #+10*t130
      t140 = 4*t40
      t145 = -42*t24-28*t26+t28+2*t30-14*t31-148*t33+238*t36-238*t38-t14
     #0+106*t42-10*t45+8*t47-112*t50
      t148 = t49*l1
      t150 = t49*t122
      t151 = 12*t150
      t152 = t129*t2
      t154 = 160*t55
      t156 = 6*t59
      t163 = -24*t52+4*t148+t151+2*t152-t154-56*t57+t156+110*t61+130*t63
     #+108*t66+68*t68-40*t70+44*t72
      t165 = t129*t16
      t167 = t4*t122
      t168 = 12*t167
      t170 = 2*t129*t122
      t171 = t30*t122
      t172 = 30*t171
      t173 = t129*t35
      t175 = t129*t65
      t179 = 4*t81
      t182 = 2*l6-1
      t183 = t182**2
      t186 = l5-1
      t188 = -1+l1
      t189 = z*l1
      t191 = (1-l1+t189)**2
      t192 = 1/t191
      t197 = dsqrt(z)
      t199 = dsqrt(l5)
      t201 = dsqrt(-t186)
      t204 = dsqrt(l1)
      t208 = l5*z*l1-l1+1-l5+l5*l1
      t210 = -t188
      t211 = dsqrt(t210)
      t216 = t208**2
      t219 = dsqrt(4*l1*t183*z*l5*t186*t188*t192-4*t182*t197*t199*t201*t
     #204*t208*t192*t211+t216*t192)
      t220 = 2*t219
      t221 = -156*t74-10*t165+t168+t76-t170-t172+20*t173-20*t175+264*t77
     #-236*t79-t179+58*t83-t220
      t223 = t132+t145+t163+t221
      t227 = zqe0(t210,z,l3)
      t228 = dlog(l1)
      t230 = sqe1(t210)
      t231 = zqe1(t210,z,l3)
      t232 = sbqe0(t210)
      t233 = dlog(z)
      t237 = 2*t122
      t246 = 40*t123
      t247 = 12*t125
      t248 = 30*t127
      t249 = -2+t237-58*t6+30*t8-128*t10+26*t12+242*t14-94*t17+326*t19+4
     #84*t22-t246-t247+t248
      t252 = 4*t5
      t263 = -218*t24-2*t2+t252-552*t26-4*t4-8*t21+10*t16-16*t31-206*t33
     #+280*t36-634*t38+154*t40+742*t42
      t270 = 2*t150
      t278 = -472*t45+56*t47-50*t50-8*t52-20*t35+t270-10*t54-190*t55-14*
     #t57+64*t59+240*t61+80*t63+38*t66
      t283 = 30*t167
      t285 = 12*t171
      t291 = 32*t68-170*t70+390*t72-580*t74+t283-t76+20*t65-t285-4*t189+
     #480*t77-210*t79-138*t81+84*t83
      t297 = -12+12*l1-4*t2
      t309 = -2+t237-40*t6+10*t8-74*t10+14*t12+158*t14-44*t17+196*t19+31
     #6*t22-t246-t247
      t318 = t248-138*t24+t252-360*t26-t28-t29+2*t16-t32-128*t33+188*t36
     #-434*t38+86*t40+502*t42
      t329 = -302*t45+34*t47-32*t50-t53-8*t35+t270-8*t54-t154-12*t57+52*
     #t59+200*t61+68*t63
      t340 = 28*t66+18*t68-140*t70+276*t72-424*t74+t283-t76+12*t65-t285+
     #356*t77-156*t79-92*t81+48*t83
      t350 = -2-t116+t9-t117+t13-t119+14*t19-200*t22+t124+t126-t128-12*t
     #130
      t362 = 72*t24+52*t26-t28-4*t30+26*t31+236*t33-346*t36+334*t38+t140
     #-146*t42+16*t45-14*t47+170*t50
      t375 = 44*t52-8*t148-t151-4*t152+180*t55+66*t57-t156-120*t61-150*t
     #63-146*t66-114*t68+42*t70-50*t72
      t383 = 188*t74+18*t165-t168-t76+t170+t172-32*t173+28*t175-332*t77+
     #308*t79+t179-100*t83+t220
      t393 = (l1+z-t189)**2
      taux12 = ((128.D0/3.D0*t3*t86*s1l1+128.D0/3.D0*t3*t113*s2l1+128.D0
     #/3.D0*t3*t223)*(t227*t228+t230-t231+t232*t233)-(128.D0/3.D0*(t3*(t
     #249+t263+t278+t291)+t297*t86)*s1l1+128.D0/3.D0*(t3*(t309+t318+t329
     #+t340)+t297*t113)*s2l1+128.D0/3.D0*t3*(t350+t362+t375+t383)+128.D0
     #/3.D0*t297*t223)*t227)/t393*t192

CCC Coefficients

      tcoff19 = taux11
      tcoff20 = taux11
      tcoff21 = taux12
      tcoff22 = -taux12

CCC Jet functions

      tFjet19 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff19)
      tFjet20 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff20)
      tFjet21 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff21)
      tFjet22 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff22)

CCC Sum

      RVnum6e1 = tFjet19*tcoff19+tFjet20*tcoff20+tFjet21*tcoff21+tFjet22
     #*tcoff22
      RETURN
      END

      double precision function RVnum6e2(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t1 = l1**2
      t2 = t1*l1
      t3 = t1**2
      t4 = t3*t2
      t6 = t3**2
      t7 = 128.D0/3.D0*t6
      t8 = z**2
      t9 = t8*t4
      t10 = 3328*t9
      t11 = t8**2
      t12 = t11*t6
      t13 = 640.D0/3.D0*t12
      t14 = t8*z
      t15 = t14*t6
      t16 = 1280.D0/3.D0*t15
      t17 = t14*t4
      t19 = z*t4
      t21 = z*t6
      t22 = 640.D0/3.D0*t21
      t24 = t11*t4
      t26 = t8*t6
      t27 = 1280.D0/3.D0*t26
      t28 = t8*l1
      t30 = 896.D0/3.D0*t4-t7+t10-t13+t16-10496.D0/3.D0*t17-4736.D0/3.D0
     #*t19+t22-3968.D0/3.D0*t3+5504.D0/3.D0*t24-t27+5120.D0/3.D0*t28
      t31 = t14*t1
      t33 = t3*t1
      t35 = t14*t33
      t37 = t11*z
      t38 = t33*t37
      t40 = t3*l1
      t41 = t8*t40
      t43 = t14*t40
      t45 = t11*t40
      t47 = t37*t40
      t49 = t33*z
      t51 = t8*t3
      t54 = t14*t3
      t56 = 640.D0/3.D0+26752.D0/3.D0*t31-896*t33+12160*t35+4096.D0/3.D0
     #*t38+21120*t41-70400.D0/3.D0*t43+37504.D0/3.D0*t45-2560*t47+15232.
     #D0/3.D0*t49-72064.D0/3.D0*t51+4352.D0/3.D0*t40+82048.D0/3.D0*t54
      t58 = t11*t1
      t60 = z*t3
      t62 = 256.D0/3.D0*t8
      t63 = t14*t2
      t65 = 128*z
      t66 = t11*t2
      t68 = 128.D0/3.D0*l1
      t69 = 1024.D0/3.D0*t14
      t70 = t8*t2
      t73 = z*t2
      t75 = t8*t1
      t77 = t14*l1
      t79 = -3840*t58+28544.D0/3.D0*t60-t62-59648.D0/3.D0*t63+t65+9728*t
     #66-t68+t69+16896*t70-512.D0/3.D0*t1-17536.D0/3.D0*t73-21760.D0/3.D
     #0*t75-7040.D0/3.D0*t77
      t80 = z*t1
      t82 = t11*t3
      t84 = t33*t11
      t86 = t33*t8
      t88 = z*t40
      t90 = z*l1
      t92 = t37*t4
      t93 = 384*t92
      t94 = t37*t6
      t95 = 128.D0/3.D0*t94
      t96 = t37*t3
      t98 = t37*t2
      t100 = t37*t1
      t101 = 384*t100
      t102 = t11*l1
      t103 = 768*t102
      t105 = 5632.D0/3.D0*t80-14208*t82-19456.D0/3.D0*t84-33664.D0/3.D0*
     #t86-27136.D0/3.D0*t88-256*t90-t93+t95+2688*t96-1536*t98+t101+t103+
     #640*t2
      t108 = s1l1**2
      t113 = 1280.D0/3.D0*t12
      t114 = 2560.D0/3.D0*t15
      t121 = 512+512*t4-256.D0/3.D0*t6+5632*t9-t113+t114-5888*t17-2688*t
     #19+1280.D0/3.D0*t21-1024*t3+3072*t24-2560.D0/3.D0*t26
      t135 = 2048.D0/3.D0*t28+19840.D0/3.D0*t31-1280*t33+17664*t35+6016.
     #D0/3.D0*t38+77312.D0/3.D0*t41-88064.D0/3.D0*t43+47872.D0/3.D0*t45-
     #3328*t47+7296*t49-24320*t51+4864.D0/3.D0*t40+29184*t54
      t143 = 256*l1
      t148 = 1408*t77
      t149 = -3072*t58+8832*t60+512.D0/3.D0*t8-17792*t63+1024.D0/3.D0*z+
     #27392.D0/3.D0*t66-t143+256*t14+13824*t70-3968*t73-4608*t75-t148
      t155 = 512.D0/3.D0*t90
      t160 = 256*t100
      t161 = 512*t102
      t163 = 768*t80-47104.D0/3.D0*t82-9472*t84-48640.D0/3.D0*t86-32000.
     #D0/3.D0*t88-t155-640*t92+256.D0/3.D0*t94+9088.D0/3.D0*t96-1408*t98
     #+t160+t161+256*t2
      t168 = 2*l6-1
      t169 = t168**2
      t172 = l5-1
      t174 = -1+l1
      t176 = (1-l1+t90)**2
      t177 = 1/t176
      t182 = dsqrt(z)
      t184 = dsqrt(l5)
      t186 = dsqrt(-t172)
      t189 = dsqrt(l1)
      t190 = l5*z
      t191 = t190*l1
      t192 = l5*l1
      t193 = t191-l1+1-l5+t192
      t195 = -t174
      t196 = dsqrt(t195)
      t201 = t193**2
      t204 = dsqrt(4*l1*t169*z*l5*t172*t174*t177-4*t168*t182*t184*t186*t
     #189*t193*t177*t196+t201*t177)
      t205 = t204*t1
      t206 = 128.D0/3.D0*t205
      t207 = t11*t8
      t208 = t207*t4
      t210 = 384*t11
      t217 = t37*l1
      t218 = 768*t217
      t221 = t207*t6
      t223 = (t121+t135+t149+t163)*s2l1+t206-2240.D0/3.D0*t208+t210-9920
     #.D0/3.D0*t9+3200.D0/3.D0*t12-3200.D0/3.D0*t15+20480.D0/3.D0*t17+64
     #0*t19-320.D0/3.D0*t21+t218-7040*t24+1600.D0/3.D0*t26+320.D0/3.D0*t
     #221
      t226 = t207*t3
      t229 = t207*t40
      t231 = t207*t1
      t232 = 384*t231
      t242 = -128.D0/3.D0*t28-20864.D0/3.D0*t31+3520*t226-896.D0/3.D0*t2
     #04-11200.D0/3.D0*t229+t232-58240.D0/3.D0*t35-32704.D0/3.D0*t38-131
     #20*t41+93184.D0/3.D0*t43-103808.D0/3.D0*t45+18368*t47-1600*t49+112
     #00*t51-30400*t54
      t246 = 640.D0/3.D0*t8
      t250 = 640.D0/3.D0*l1
      t251 = 128.D0/3.D0*t14
      t253 = 128.D0/3.D0*t1
      t257 = 1664.D0/3.D0+9536*t58-1280*t60+t246+18560*t63+384*z-23232*t
     #66-t250+t251-5376*t70-t253+320*t73+3328.D0/3.D0*t75+4160.D0/3.D0*t
     #77
      t260 = t207*t33
      t262 = t207*t2
      t270 = t204*l1
      t276 = -128.D0/3.D0*t80+35648*t82+6784.D0/3.D0*t260-1792*t262+2073
     #6*t84+26752.D0/3.D0*t86+6080.D0/3.D0*t88-128*t90+10880.D0/3.D0*t92
     #-1600.D0/3.D0*t94+256.D0/3.D0*t270-18688*t96+11584*t98-4224*t100-2
     #304*t102
      t286 = 896.D0/3.D0+640.D0/3.D0*t4-t7+2816*t9-t13+t16-9472.D0/3.D0*
     #t17-3712.D0/3.D0*t19+t22-128*t3+5248.D0/3.D0*t24-t27
      t300 = 640*t28+7424*t31-1280.D0/3.D0*t33+30848.D0/3.D0*t35+1408*t3
     #8+42112.D0/3.D0*t41-57088.D0/3.D0*t43+11776*t45-8192.D0/3.D0*t47+9
     #472.D0/3.D0*t49-14720*t51+384*t40+65408.D0/3.D0*t54
      t304 = 128*t8
      t308 = 384*t14
      t311 = 11392.D0/3.D0*t75
      t312 = -3968*t58+3712*t60+t304-47872.D0/3.D0*t63+640.D0/3.D0*z+960
     #0*t66-t250+t308+9728*t70+t253-1792*t73-t311
      t314 = 1280.D0/3.D0*t80
      t321 = -2176*t77+t314-40832.D0/3.D0*t82-18304.D0/3.D0*t84-8320*t86
     #-13312.D0/3.D0*t88-t155-t93+t95+2944*t96-1664*t98+t101+t103
      t324 = s2l1**2
      t327 = 128*t11
      t333 = 256*t217
      t336 = 64.D0/3.D0*t221
      t337 = 640-t206+256.D0/3.D0*t208+t327+3584.D0/3.D0*t9-t13+640.D0/3
     #.D0*t15-6016.D0/3.D0*t17-832.D0/3.D0*t19+64.D0/3.D0*t21+t333+4864.
     #D0/3.D0*t24-320.D0/3.D0*t26-t336
      t343 = 128*t231
      t352 = 1664.D0/3.D0*t28+20672.D0/3.D0*t31+1664.D0/3.D0*t226-384*t2
     #04-704.D0/3.D0*t229+t343+7680*t35+1536*t38+10560*t41-47744.D0/3.D0
     #*t43+30208.D0/3.D0*t45-5824.D0/3.D0*t47+3584.D0/3.D0*t49-12928*t51
      t366 = 19520*t54-7936.D0/3.D0*t58+2880*t60+896.D0/3.D0*t8-14656*t6
     #3+1408.D0/3.D0*z+21760.D0/3.D0*t66-384*l1+512*t14+9280*t70+t253-16
     #64*t73-t311-2048*t77
      t381 = t314-33344.D0/3.D0*t82-64*t260-448*t262-16256.D0/3.D0*t84-4
     #928*t86-7616.D0/3.D0*t88-896.D0/3.D0*t90-1856.D0/3.D0*t92+320.D0/3
     #.D0*t94+256*t270+1088*t96+640.D0/3.D0*t98-640*t100+320*t102
      t385 = t207*l5
      t386 = t385*l1
      t388 = l5**2
      t389 = t388*t14
      t390 = t389*t1
      t392 = t388*t11
      t393 = t392*t1
      t395 = t388*l1
      t396 = t395*t14
      t398 = t388*z
      t399 = t398*t1
      t401 = t11*t14
      t402 = t401*t1
      t404 = t6*t1
      t405 = t404*t8
      t407 = z*t404
      t409 = t6*l1
      t410 = t409*t14
      t413 = t11*l5
      t415 = -20-12*t386-6*t390+12*t393+6*t396-36*t399-72*t402+268*t405-
     #34*t407+3502*t410-4*t392+4*t413
      t416 = t37*l5
      t417 = t416*l1
      t419 = t207*t388
      t420 = t419*l1
      t422 = t401*t388
      t423 = t422*t1
      t425 = t416*t1
      t427 = t388*t2
      t429 = t182*t37
      t430 = t184*l5
      t431 = t429*t430
      t432 = t431*t186
      t434 = t189*t196*l6
      t435 = t432*t434
      t443 = 4*t417+12*t420-2*t423-8*t425+2*t427-16*t435+40*t205-19085*t
     #208-3665*t9+15215*t12-7974*t15+12110*t17
      t447 = t204*t8
      t448 = t447*t2
      t451 = t182*t207*t430
      t452 = t189*l1
      t454 = t186*t452*t196
      t455 = t451*t454
      t457 = t182*t184
      t459 = t186*t189*t196
      t460 = t457*t459
      t462 = t11**2
      t463 = t462*t33
      t465 = t462*z
      t466 = t465*t40
      t468 = t462*t40
      t470 = t401*t33
      t472 = t404*t401
      t474 = t404*t207
      t476 = t404*t37
      t478 = 482*t19-308*t21-10*t448-8*t455-16*t460+1765*t463+66*t466-95
     #5*t468-7778*t470-868*t472+1764*t474-2240*t476
      t479 = t431*t459
      t481 = t457*t186
      t482 = t481*t434
      t484 = t204*t2
      t486 = t37*t388
      t490 = l5*t1
      t491 = t490*t11
      t493 = t192*t11
      t495 = t486*t1
      t497 = t486*t2
      t499 = t486*l1
      t501 = t419*t1
      t503 = 8*t479+32*t482-10*t484-2*t486+4*t416-2*t389+4*t491-8*t493-8
     #*t495+12*t497-2*t499+4*t501
      t506 = t189*t1
      t508 = t506*t196*l6
      t509 = t481*t508
      t511 = t409*t207
      t513 = t409*t8
      t515 = t409*z
      t517 = t404*t462
      t519 = t409*t462
      t522 = t388*t1
      t525 = t452*t196*l6
      t526 = t481*t525
      t528 = t462*t4
      t530 = t388*t8
      t531 = t530*t1
      t533 = t389*t2
      t535 = 16*t509-6347*t511-1029*t513+132*t515+244*t517-837*t519+10*t
     #395-8*t522-48*t526-2118*t528-4*t531+2*t533
      t536 = t392*t2
      t541 = t190*t2
      t543 = t422*t2
      t545 = t457*t454
      t549 = t182*z*t430
      t550 = t549*t186
      t551 = t550*t508
      t554 = t204*z
      t555 = t554*l1
      t557 = t429*t184
      t558 = t557*t459
      t560 = -8*t536+116*t217-22660*t24+2375*t26-12*t541+2*t543+24*t545+
     #13573*t221-16*t551-14*t37+40*t555-8*t558
      t564 = t462*t6
      t566 = t431*t454
      t568 = t447*t1
      t570 = t465*t4
      t572 = t465*t6
      t574 = t550*t525
      t576 = t419*t2
      t578 = t404*t14
      t580 = l6**2
      t581 = t1*t580
      t582 = t190*t581
      t584 = l1*l6
      t585 = t398*t584
      t587 = -28*t28-564*t31+1669*t564+16*t566+20*t568+222*t570-190*t572
     #+32*t574-8*t576-924*t578+96*t582-64*t585
      t588 = t6*t2
      t589 = t588*t401
      t591 = t588*t462
      t595 = t182*t430
      t596 = t595*t459
      t606 = 112*t589-32*t591+6006*t226+20*t204+16*t596-12598*t229+408*t
     #231-12910*t35-26466*t38-3123*t41+9882*t43-17558*t45+18996*t47
      t611 = t588*t465
      t613 = t409*t37
      t615 = t409*t11
      t619 = t182*t11*t430
      t620 = t619*t454
      t622 = t595*t186
      t623 = t622*t508
      t625 = t465*t33
      t628 = t2*l6
      t629 = t190*t628
      t631 = t398*t628
      t633 = t2*t580
      t634 = t190*t633
      t636 = -534*t49+4*t611+8220*t613-6795*t615+1759*t51-8*t620-16*t623
     #-160*t625+16*t490+32*t629-32*t631-32*t634
      t637 = t398*t633
      t639 = t385*t1
      t642 = t465*t409
      t644 = t554*t1
      t646 = t190*t1
      t648 = t398*t2
      t650 = t530*t2
      t653 = t622*t525
      t657 = 32*t637+4*t639-5444*t54+100*t642-60*t644+36*t646+12*t648+2*
     #t650+853*t58+48*t653-232*t60-8*t8
      t661 = t619*t459
      t663 = t465*t3
      t668 = t186*t506*t196
      t669 = t549*t668
      t673 = t1*l6
      t674 = t398*t673
      t676 = l1*t580
      t677 = t398*t676
      t679 = t190*t676
      t681 = 8*l5+2124*t63+8*t661-12*t663-12*z-3472*t66+8*t669+50*l1-24*
     #t191+96*t674+64*t677-64*t679
      t682 = t190*t673
      t684 = t190*t584
      t686 = t398*t581
      t688 = t207*l1
      t690 = t622*t434
      t692 = t465*t404
      t694 = t588*t207
      t698 = t462*t2
      t700 = t462*t3
      t703 = -96*t682+64*t684-96*t686-48*t688-32*t690-30*t692-224*t694-1
     #6*t14-695*t70-48*t698+312*t700-40*t1
      t707 = t595*t454
      t709 = t588*z
      t711 = t404*t11
      t716 = t588*t14
      t718 = t588*t11
      t720 = t588*t37
      t722 = t595*t668
      t724 = t432*t525
      t726 = 68*t73-24*t707+4*t709+1820*t711+184*t75+96*t77+12*t80+112*t
     #716-224*t718+280*t720+8*t722-32*t724
      t727 = t401*t40
      t731 = t385*t676
      t733 = t419*t676
      t735 = t385*t584
      t737 = t419*t584
      t739 = t401*t6
      t741 = l5*t2
      t743 = t554*t2
      t745 = t398*l1
      t748 = 4866*t727+9339*t82-20*t192-32*t731+32*t733+32*t735-32*t737-
     #6330*t739-4*t741+20*t743+24*t745+18500*t260
      t752 = t451*t186*t525
      t754 = t549*t454
      t757 = t557*t186*t434
      t766 = t409*t401
      t768 = -1954*t262+16*t752-16*t754+16*t757+23597*t84+3986*t86+424*t
     #88-2*t90+26192*t92-18030*t94-4*t388+3054*t766
      t770 = t401*t3
      t773 = t588*t8
      t775 = t619*t186
      t776 = t775*t525
      t779 = t457*t668
      t781 = t775*t434
      t785 = t401*t4
      t788 = t401*t2
      t790 = -50*t270-2046*t770-9682*t96-32*t773+16*t776+3432*t98-8*t779
     #-16*t781-806*t100-119*t102+8522*t785+10*t2+540*t788
      t798 = zqe0(t195,z,l3)
      t799 = dlog(l1)
      t801 = sqe1(t195)
      t802 = zqe1(t195,z,l3)
      t803 = sbqe0(t195)
      t804 = dlog(z)
      t815 = 1024.D0/3.D0*t38
      t825 = 128.D0/3.D0*t8
      t826 = -128.D0/3.D0-256*t9+384*t17+64*t19-256*t24-1408.D0/3.D0*t28
     #-8512.D0/3.D0*t31-2176*t35-t815-11264.D0/3.D0*t41+16192.D0/3.D0*t4
     #3-10112.D0/3.D0*t45+768*t47-384*t49+15232.D0/3.D0*t51-7296*t54+128
     #0*t58-1216*t60+t825
      t831 = 768*t73
      t839 = 64*t90
      t845 = 17600.D0/3.D0*t63-9472.D0/3.D0*t66-t68-320.D0/3.D0*t14-3968
     #*t70+t253+t831+5504.D0/3.D0*t75+2368.D0/3.D0*t77-448.D0/3.D0*t80+1
     #2928.D0/3.D0*t82+1408*t84+4480.D0/3.D0*t86+960*t88-t839+64*t92-262
     #4.D0/3.D0*t96+512*t98-128*t100-256*t102
      t848 = 6*t386
      t849 = 3*t390
      t850 = 6*t393
      t851 = 3*t396
      t852 = 18*t399
      t857 = 2*t392
      t858 = 2*t413
      t859 = 7+t848+t849-t850-t851+t852+18*t402-64*t405+8*t407-823*t410+
     #t857-t858
      t860 = 2*t417
      t861 = 6*t420
      t862 = 4*t425
      t863 = 8*t435
      t870 = -t860-t861+t423+t862-t427+t863-10*t205+5005*t208+753*t9-363
     #0*t12+1813*t15-2657*t17
      t875 = 4*t455
      t877 = 8*t460
      t884 = -93*t19+63*t21+2*t448+t875+t204*t3+t877-504*t463-18*t466+26
     #8*t468+2153*t470+224*t472-448*t474
      t886 = 4*t479
      t887 = 16*t482
      t888 = 2*t416
      t889 = 2*t491
      t890 = 4*t493
      t891 = 4*t495
      t892 = 6*t497
      t893 = 2*t501
      t894 = 560*t476-t886-t887+t486-t888+t389-t889+t890+t891-t892+t499-
     #t893-t3
      t897 = 8*t509
      t903 = 5*t395
      t904 = 4*t522
      t905 = 24*t526
      t907 = 2*t531
      t908 = -t897+1634*t511+234*t513-29*t515-64*t517+228*t519-t903+t904
     #+t905+604*t528+t907-t533
      t909 = 4*t536
      t913 = 6*t541
      t914 = 12*t545
      t916 = 8*t551
      t919 = 4*t558
      t920 = t909-31*t217+5302*t24-513*t26+t913-t543-t914-3533*t221+t916
     #+4*t37-14*t555+t919
      t925 = 8*t566
      t929 = 16*t574
      t930 = 4*t576
      t932 = 48*t582
      t933 = 32*t585
      t934 = 11*t28+94*t31-468*t564-t925-7*t568-65*t570+55*t572-t929+t93
     #0+224*t578-t932+t933
      t939 = 8*t596
      t948 = -28*t589+8*t591-1589*t226-7*t204-t939+3323*t229-108*t231+27
     #24*t35+6524*t38+598*t41-1996*t43+3891*t45-4636*t47
      t956 = 4*t620
      t957 = 8*t623
      t959 = 8*t490
      t960 = 16*t629
      t961 = 16*t631
      t962 = 16*t634
      t963 = 101*t49-t611-2055*t613+1648*t615-335*t51+t956+t957+46*t625-
     #t959-t960+t961+t962
      t964 = 16*t637
      t965 = 2*t639
      t969 = 18*t646
      t970 = 6*t648
      t972 = 24*t653
      t974 = -t964-t965+1047*t54-28*t642+18*t644-t969-t970-t650-164*t58-
     #t972+50*t60+t8
      t976 = 4*l5
      t978 = 4*t661
      t982 = 4*t669
      t985 = 12*t191
      t986 = 48*t674
      t987 = 32*t677
      t988 = -t976-386*t63-t978+3*t663+3*z+708*t66-t982+t447*t3-16*l1+t9
     #85-t986-t987
      t989 = 32*t679
      t990 = 48*t682
      t991 = 32*t684
      t992 = 48*t686
      t994 = 16*t690
      t1002 = t989+t990-t991+t992+12*t688+t994+8*t692+56*t694+2*t14+137*
     #t70+12*t698-84*t700+10*t1
      t1006 = 12*t707
      t1014 = 4*t722
      t1015 = 16*t724
      t1016 = -16*t73+t1006-t709-448*t711-42*t75-12*t77-7*t80-28*t716+56
     #*t718-70*t720-t1014+t1015
      t1019 = 10*t192
      t1020 = 16*t731
      t1021 = 16*t733
      t1022 = 16*t735
      t1023 = 16*t737
      t1025 = 2*t741
      t1027 = 12*t745
      t1029 = -1347*t727-1995*t82+t1019+t1020-t1021-t1022+t1023+1713*t73
     #9+t1025-2*t743-t1027-4870*t260
      t1032 = 8*t752
      t1033 = 8*t754
      t1034 = 8*t757
      t1041 = 2*t388
      t1043 = 520*t262-t1032+t1033-t1034-5389*t84-785*t86-83*t88+4*t90-6
     #508*t92+4500*t94+t1041-809*t766
      t1048 = 8*t776
      t1050 = 4*t779
      t1051 = 8*t781
      t1058 = 16*t270+561*t770+2343*t96+8*t773-t1048-832*t98+t1050+t1051
     #-2*t554*t3+201*t100+22*t102-2341*t785-144*t788
      t1069 = 128.D0/3.D0*t21
      t1071 = 640.D0/3.D0*t26
      t1072 = 128.D0/3.D0*t221
      t1074 = -t206+896.D0/3.D0*t208-t327+3712.D0/3.D0*t9-t113+t16-7808.
     #D0/3.D0*t17-704.D0/3.D0*t19+t1069-t333+8192.D0/3.D0*t24-t1071-t107
     #2-704.D0/3.D0*t28
      t1077 = 256.D0/3.D0*t204
      t1087 = -448.D0/3.D0+7040.D0/3.D0*t31-3968.D0/3.D0*t226+t1077+4352
     #.D0/3.D0*t229-t343+7104*t35+12736.D0/3.D0*t38+12352.D0/3.D0*t41-10
     #944*t43+38848.D0/3.D0*t45-21248.D0/3.D0*t47+512*t49-9088.D0/3.D0*t
     #51
      t1092 = 64.D0/3.D0*t8
      t1101 = 1024.D0/3.D0*t28
      t1114 = -128-512.D0/3.D0*t9+256*t17+128.D0/3.D0*t19-512.D0/3.D0*t2
     #4-t1101-1024*t31-1408*t35-640.D0/3.D0*t38-7168.D0/3.D0*t41+9856.D0
     #/3.D0*t43-5888.D0/3.D0*t45+1280.D0/3.D0*t47-256*t49+3072*t51-12032
     #.D0/3.D0*t54+256*t58-2432.D0/3.D0*t60
      t1116 = 128.D0/3.D0*z
      t1131 = t825+8192.D0/3.D0*t63-t1116-1152*t66-6784.D0/3.D0*t70+t253
     #+512*t73+3200.D0/3.D0*t75+128*t77-256.D0/3.D0*t80+6400.D0/3.D0*t82
     #+896*t84+2944.D0/3.D0*t86+640*t88-128.D0/3.D0*t90+128.D0/3.D0*t92-
     #384*t96+128*t98
      t1138 = 31232.D0/3.D0*t54-11072.D0/3.D0*t58+320.D0/3.D0*t60+t1092-
     #18688.D0/3.D0*t63-64*z+8768*t66-64.D0/3.D0*l1+(t1114+t1131)*s2l1+t
     #251+3136.D0/3.D0*t70+64*t1+128*t73+128*t75
      t1147 = 640.D0/3.D0*t94
      t1148 = 64.D0/3.D0*t270
      t1153 = -1600.D0/3.D0*t77-13312*t82-896*t260+640*t262-7872*t84-928
     #0.D0/3.D0*t86-1472.D0/3.D0*t88-t839-4288.D0/3.D0*t92+t1147+t1148+2
     #1440.D0/3.D0*t96-13120.D0/3.D0*t98+1536*t100+896*t102
      t1157 = 256.D0/3.D0*t9
      t1158 = 128*t17
      t1159 = 64.D0/3.D0*t19
      t1160 = 256.D0/3.D0*t24
      t1164 = 128*t38
      t1168 = 320*t47
      t1169 = 128*t49
      t1175 = -256.D0/3.D0-t1157+t1158+t1159-t1160-256*t28-3904.D0/3.D0*
     #t31-768*t35-t1164-4096.D0/3.D0*t41+6080.D0/3.D0*t43-3968.D0/3.D0*t
     #45+t1168-t1169+6016.D0/3.D0*t51-2944*t54+640*t58-1472.D0/3.D0*t60+
     #7616.D0/3.D0*t63
      t1177 = 64*t14
      t1186 = 1024.D0/3.D0*t88
      t1187 = 64.D0/3.D0*t90
      t1188 = 64.D0/3.D0*t92
      t1189 = 1216.D0/3.D0*t96
      t1190 = 256*t98
      t1191 = 64*t100
      t1192 = 128*t102
      t1193 = -t1116-4352.D0/3.D0*t66+t68-t1177-5248.D0/3.D0*t70+384*t73
     #+896*t75+384*t77-128*t80+5504.D0/3.D0*t82+512*t84+512*t86+t1186+t1
     #187+t1188-t1189+t1190-t1191-t1192
      t1196 = 832.D0/3.D0*t28
      t1206 = 3200*t54
      t1210 = -1600*t41+2240*t43-4160.D0/3.D0*t45+t1168-448.D0/3.D0*t49+
     #7040.D0/3.D0*t51-t1206+576*t58-1856.D0/3.D0*t60-t1092+2560*t63
      t1215 = 5824.D0/3.D0*t70
      t1216 = 64.D0/3.D0*t1
      t1227 = 5632.D0/3.D0*t82+1600.D0/3.D0*t84+576*t86+1280.D0/3.D0*t88
     #+t1187+t1188-64*t270-t1189+t1190-t1191-t1192
      t1232 = t804**2
      t1235 = sbqe1(t195)
      t1237 = t799**2
      t1242 = sqe2(t195)
      t1243 = zqe2(t195,z,l3)
      t1247 = 64.D0/3.D0*t6
      t1250 = 320*t12
      t1251 = 1920*t17
      t1253 = 128*t21
      t1256 = 320*t26
      t1257 = -128.D0/3.D0+448.D0/3.D0*t4-t1247+320.D0/3.D0*t208+1600*t9
     #-t1250+t16-t1251-2240.D0/3.D0*t19+t1253-1984.D0/3.D0*t3+4160.D0/3.
     #D0*t24-t1256
      t1262 = 9088.D0/3.D0*t35
      t1265 = 1216.D0/3.D0*t43
      t1270 = -t336-t1196-3776*t31-448*t33-64*t226+192*t229+t1262+2816.D
     #0/3.D0*t38+6976.D0/3.D0*t41-t1265-128*t45-1024.D0/3.D0*t47+5504.D0
     #/3.D0*t49+512*t51
      t1273 = 13184.D0/3.D0*t54
      t1276 = 17984.D0/3.D0*t63
      t1282 = 2176.D0/3.D0*t40-t1273+2176*t58+5312.D0/3.D0*t60-t1092+t12
     #76-12416.D0/3.D0*t66-t68-704.D0/3.D0*t14-t1215-t1216-2816.D0/3.D0*
     #t73+3712.D0/3.D0*t75+3904.D0/3.D0*t77
      t1291 = 128*t94
      t1293 = 768*t98
      t1295 = 1472.D0/3.D0*t80+10496.D0/3.D0*t82-640.D0/3.D0*t260-6080.D
     #0/3.D0*t84-9344.D0/3.D0*t86-2368*t88-192*t90-576*t92+t1291-1984.D0
     #/3.D0*t96+t1293-t160-t161+320*t2
      t1299 = 256*t208
      t1307 = -t206+t1299-256*t11+768*t9-t113+t16-1792*t17-128*t19+t1069
     #-512*t217+2048*t24-t1071-t1072-2432.D0/3.D0*t28
      t1320 = -10624.D0/3.D0*t31-4096.D0/3.D0*t226+t1077+3712.D0/3.D0*t2
     #29-256*t231+7552.D0/3.D0*t35+8704.D0/3.D0*t38-2816*t41+1024.D0/3.D
     #0*t43+13696.D0/3.D0*t45-13312.D0/3.D0*t47-t1169+18944.D0/3.D0*t51-
     #13952.D0/3.D0*t54
      t1332 = -256+256*t4-t7+t1299+t10-640*t12+t114-4352*t17-1408*t19+25
     #6*t21-512*t3+3328*t24-640*t26-t1072
      t1347 = -3328.D0/3.D0*t28-25600.D0/3.D0*t31-640*t33-512*t226+2432.
     #D0/3.D0*t229+8192*t35+3072*t38+12800.D0/3.D0*t41-12800.D0/3.D0*t43
     #+12800.D0/3.D0*t45-2944*t47+3072*t49+9856.D0/3.D0*t51+2432.D0/3.D0
     #*t40
      t1353 = 512.D0/3.D0*z
      t1361 = -20480.D0/3.D0*t54+4480*t58+640*t60-t62+36992.D0/3.D0*t63-
     #t1353-7040*t66+128*l1-640*t14-7040*t70+t831+4480*t75+3072*t77-384*
     #t80
      t1374 = 9856.D0/3.D0*t82-640*t260+128*t262-6528*t84-6528*t86-2944*
     #t88+256.D0/3.D0*t90-1408*t92+256*t94+640*t96+t1293-t101-1152*t102+
     #128*t2
      t1383 = -640.D0/3.D0-3712.D0/3.D0*t58-1920*t60-t825+(t1332+t1347+t
     #1361+t1374)*s2l1+17152.D0/3.D0*t63-t65+1920*t66+t68-t69-6144*t70+t
     #253+1408*t73+9728.D0/3.D0*t75
      t1397 = 3968.D0/3.D0*t77-1024.D0/3.D0*t80-3200*t82-2176.D0/3.D0*t2
     #60+896*t262-4224*t84-1024.D0/3.D0*t86+3328.D0/3.D0*t88-1152*t92+t1
     #147+t1148+14464.D0/3.D0*t96-11264.D0/3.D0*t98+1920*t100+640*t102
      t1409 = 320.D0/3.D0*t4-t1247+448.D0/3.D0*t208-64*t11+4160.D0/3.D0*
     #t9-t1250+t16-t1251-576*t19+t1253-64*t3-128*t217+1600*t24-t1256
      t1422 = -t336-384*t28-11456.D0/3.D0*t31-640.D0/3.D0*t33-1984.D0/3.
     #D0*t226+2176.D0/3.D0*t229-64*t231+t1262+5504.D0/3.D0*t38-128*t41-t
     #1265+6976.D0/3.D0*t45-2368*t47+2816.D0/3.D0*t49+10496.D0/3.D0*t51
      t1432 = -640.D0/3.D0+192*t40-t1273+3584.D0/3.D0*t58-1984.D0/3.D0*t
     #60-t304+t1276-t1353-5824.D0/3.D0*t66+512.D0/3.D0*l1-320*t14-12416.
     #D0/3.D0*t70-t253+t831+6400.D0/3.D0*t75
      t1433 = 896.D0/3.D0*t80
      t1445 = t148-t1433+512*t82-448*t260+320*t262-9344.D0/3.D0*t84-6080
     #.D0/3.D0*t86-t1186+448.D0/3.D0*t90-2240.D0/3.D0*t92+t1291+5312.D0/
     #3.D0*t96-2816.D0/3.D0*t98+448*t100-192*t102
      t1454 = -384+t206-128*t208-t210-1152*t9+t113-t16+2048*t17+256*t19-
     #t1069-t218-1792*t24+t1071+t1072
      t1466 = -t1101-1280*t31-1920*t226+256*t204+3328.D0/3.D0*t229-t232-
     #4224*t35-t815-13312.D0/3.D0*t41+13696.D0/3.D0*t43+1024.D0/3.D0*t45
     #-2816*t47-2176.D0/3.D0*t49+14464.D0/3.D0*t51
      t1477 = -t1206-3584*t58-4096.D0/3.D0*t60-t246+1920*t63-896.D0/3.D0
     #*z+17152.D0/3.D0*t66+t143-t308-11264.D0/3.D0*t70-t253+896*t73+5632
     #.D0/3.D0*t75+768*t77
      t1491 = -t1433-13952.D0/3.D0*t82-128*t260+1408*t262+7552.D0/3.D0*t
     #84+8704.D0/3.D0*t86+3712.D0/3.D0*t88+640.D0/3.D0*t90+768*t92-t1147
     #-192*t270+18944.D0/3.D0*t96-6144*t98+3200*t100+1408*t102
      t1499 = 10+t848+t849-t850-t851+t852+90*t402-334*t405+43*t407-4229*
     #t410+t857-t858
      t1506 = -t860-t861+t423+t862-t427+t863-20*t205+19155*t208-3*t11+45
     #77*t9-17375*t12+9513*t15
      t1518 = -14376*t17-630*t19+399*t21+5*t448+t875+t877-1608*t463-69*t
     #466+929*t468+7313*t470+994*t472-2058*t474
      t1521 = 2660*t476-t886-t887+5*t484+t486-t888+t389-t889+t890+t891-t
     #892+t499
      t1530 = -t893-t897+6964*t511+1280*t513-169*t515-274*t517+854*t519-
     #t903+t904+t905+1911*t528+t907
      t1537 = -t533+t909-142*t217+25477*t24-2955*t26+t913-t543-t914-1412
     #9*t221+t916+16*t37-20*t555
      t1546 = t919+32*t28+744*t31-1573*t564-t925-10*t568-192*t570+170*t5
     #72-t929+t930+1134*t578-t932
      t1558 = t933-140*t589+40*t591-6180*t226-10*t204-t939+12572*t229-47
     #4*t231+15405*t35+27828*t38+3980*t41-11979*t43+19915*t45
      t1569 = -20096*t47+700*t49-5*t611-9325*t613+7958*t615-2269*t51+t95
     #6+t957+146*t625-t959-t960+t961
      t1575 = t962-t964-t965+6759*t54-98*t642+30*t644-t969-t970-t650-105
     #6*t58-t972+299*t60
      t1583 = 4*t8-t976-2713*t63-t978+15*t663+6*z+4155*t66-t982-25*l1+t9
     #85-t986-t987
      t1592 = t989+t990-t991+t992+60*t688+t994+33*t692+280*t694+17*t14+9
     #01*t70+60*t698-339*t700+20*t1
      t1604 = -94*t73+t1006-5*t709-2198*t711-233*t75-132*t77+3*t80-140*t
     #716+280*t718-350*t720-t1014+t1015
      t1610 = -4699*t727-10836*t82+t1019+t1020-t1021-t1022+t1023+6275*t7
     #39+t1025-10*t743-t1027-18299*t260
      t1619 = 2114*t262-t1032+t1033-t1034-26466*t84-5019*t86-553*t88+t90
     #-27830*t92+19675*t94+t1041-3235*t766
      t1630 = 25*t270+2100*t770+10451*t96+40*t773-t1048-3817*t98+t1050+t
     #1051+931*t100+151*t102-8092*t785-5*t2-606*t788
      t1641 = (l1+z-t90)**2
      s2 = ((t30+t56+t79+t105)*t108+(t223+t242+t257+t276)*s1l1+(t286+t30
     #0+t312+t321)*t324+(t337+t352+t366+t381)*s2l1-64.D0/3.D0*(t415+t443
     #+t478+t503+t535+t560+t587+t606+t636+t657+t681+t703+t726+t748+t768+
     #t790)*t177)*(t798*t799+t801-t802+t803*t804)
      s3 = ((t826+t845)*t108-64.D0/3.D0*(t859+t870+t884+t894+t908+t920+t
     #934+t948+t963+t974+t988+t1002+t1016+t1029+t1043+t1058)*t177+(t1074
     #+t1087+t1138+t1153)*s1l1+(t1175+t1193)*t324+(-192-t1157+t1158+t115
     #9-t1160-t1196-3520.D0/3.D0*t31+128*t204-832*t35-t1164+t1210-320.D0
     #/3.D0*z-4160.D0/3.D0*t66+64*l1-t1177-t1215+t1216+448*t73+2944.D0/3
     #.D0*t75+256*t77-320.D0/3.D0*t80+t1227)*s2l1)*(-t803*t1232/2+t1235*
     #t804-t798*t1237/2-(t801-t802)*t799+t1242-t1243)-((t1257+t1270+t128
     #2+t1295)*t108+(t1307+t1320+t1383+t1397)*s1l1+(t1409+t1422+t1432+t1
     #445)*t324+(t1454+t1466+t1477+t1491)*s2l1-64.D0/3.D0*(t1499+t1506+t
     #1518+t1521+t1530+t1537+t1546+t1558+t1569+t1575+t1583+t1592+t1604+t
     #1610+t1619+t1630)*t177)*t798
      s1 = s2+s3
      s2 = 1/t1641*t177
      taux11 = s1*s2
      t2 = l1**2
      t3 = 3-3*l1+t2
      t4 = t2**2
      t5 = t4*t2
      t7 = z*t2
      t9 = z**2
      t10 = t9*t2
      t12 = t9**2
      t13 = t12*t9
      t14 = t2*l1
      t15 = t4*t14
      t16 = t13*t15
      t17 = 2*t16
      t18 = z*t14
      t20 = t9*t14
      t22 = t12*t14
      t24 = t9*z
      t25 = t24*t14
      t27 = t12*t2
      t29 = z*t4
      t31 = t24*t4
      t33 = t9*t4
      t35 = -2-10*t5+30*t7-128*t10+t17-94*t18+326*t20+484*t22-552*t25-21
     #8*t27+154*t29+742*t31-472*t33
      t36 = t12*l1
      t38 = t12*z
      t39 = t38*t2
      t41 = t38*t14
      t43 = t38*t4
      t45 = t12*t4
      t47 = t38*l1
      t48 = 16*t47
      t49 = t24*t5
      t51 = t5*t38
      t53 = t4*l1
      t54 = t9*t53
      t56 = t24*t53
      t58 = t13*t53
      t60 = t12*t53
      t62 = t38*t53
      t64 = 56*t36+84*t39-206*t41+280*t43-634*t45-t48+240*t49+80*t51+390
     #*t54-580*t56+38*t58+480*t60-210*t62
      t66 = t5*z
      t68 = t13*t4
      t70 = t13*t2
      t72 = t5*t12
      t74 = t13*t5
      t76 = t5*t9
      t78 = z*t53
      t80 = t13*t14
      t82 = z*l1
      t84 = t38*t15
      t85 = 12*t84
      t86 = t9*t15
      t87 = 30*t86
      t88 = t24*t15
      t89 = 40*t88
      t90 = z*t15
      t91 = 12*t90
      t92 = 64*t66-50*t68-8*t70-190*t72-14*t74-170*t76-138*t78+32*t80-4*
     #t82-t85+t87-t89-t91
      t93 = t12*t15
      t94 = 30*t93
      t96 = 4*t24
      t98 = t24*t2
      t102 = t24*l1
      t104 = t9*l1
      t106 = 2*z
      t109 = 2*t15
      t110 = t94-4*t9+t96-2*t2+242*t98+10*t14-8*t12-58*t102+26*t104-t106
     #-20*t4+20*t53+t109
      t112 = t35+t64+t92+t110
      t116 = -12+12*l1-4*t2
      t117 = 2*t7
      t118 = 26*t10
      t119 = 10*t18
      t132 = 8*t47
      t135 = 2-t117+t118+t119-72*t20-192*t22+176*t25+90*t27-18*t29-222*t
     #31+102*t33-24*t36-40*t39+94*t41-120*t43+234*t45+t132-40*t49-20*t51
      t143 = 4*t70
      t145 = 4*t74
      t149 = 2*t9
      t151 = 4*t12
      t153 = 4*t104
      t154 = -72*t54+148*t56-16*t58-152*t60+78*t62-4*t66+24*t68+t143+40*
     #t72+t145+20*t76+14*t78-16*t80+t149-78*t98+t151+18*t102-t153+t106
      t155 = t135+t154
      t159 = 8*t5
      t169 = -2-t159+10*t7-74*t10+t17-44*t18+196*t20+316*t22-360*t25-138
     #*t27+86*t29+502*t31
      t182 = -302*t33+34*t36+48*t39-128*t41+188*t43-434*t45-t132+200*t49
     #+68*t51+276*t54-424*t56+28*t58+356*t60
      t187 = 160*t72
      t192 = -156*t62+52*t66-32*t68-t143-t187-12*t74-140*t76-92*t78+18*t
     #80-t85+t87-t89
      t199 = -t91+t94-t149+t96+158*t98+2*t14-t151-40*t102+14*t104-t106-8
     #*t4+12*t53+t109
      t201 = t169+t182+t192+t199
      t206 = 56*t25
      t213 = 2-4*t7+t118+t119-40*t20-32*t22+t206+18*t27-8*t29-42*t31+30*
     #t33-4*t36-2*t39
      t220 = 2*t62
      t226 = 6*t41-6*t43+26*t45-8*t54+12*t56-8*t60+t220+2*t78+t149-2*t24
     #-38*t98+14*t102-8*t104+t106
      t227 = t213+t226
      t231 = 14*t10
      t232 = 12*t16
      t233 = 6*t18
      t238 = 4*t29
      t241 = -2+t117-t231-t232-t233+14*t20-200*t22+52*t25+72*t27+t238-14
     #6*t31+16*t33
      t255 = -14*t36-100*t39+236*t41-346*t43+334*t45+26*t47-120*t49-150*
     #t51-50*t54+188*t56-146*t58-332*t60+308*t62
      t257 = 6*t66
      t263 = 4*t78
      t265 = 30*t84
      t266 = 12*t86
      t267 = 30*t88
      t268 = 2*t90
      t269 = 40*t93
      t270 = -t257+170*t68+44*t70+180*t72+66*t74+42*t76+t263-114*t80+t26
     #5-t266+t267+t268-t269
      t271 = t12*t24
      t272 = t271*t5
      t274 = t271*t14
      t276 = 6*t102
      t278 = 2*l6-1
      t279 = t278**2
      t282 = l5-1
      t284 = -1+l1
      t286 = (1-l1+t82)**2
      t287 = 1/t286
      t292 = dsqrt(z)
      t294 = dsqrt(l5)
      t296 = dsqrt(-t282)
      t299 = dsqrt(l1)
      t303 = l5*z*l1-l1+1-l5+l5*l1
      t305 = -t284
      t306 = dsqrt(t305)
      t311 = t303**2
      t314 = dsqrt(4*l1*t279*z*l5*t282*t284*t287-4*t278*t292*t294*t296*t
     #299*t303*t287*t306+t311*t287)
      t315 = 2*t314
      t317 = t271*t2
      t319 = t13*l1
      t322 = 2*t271*t15
      t323 = t271*t4
      t325 = t271*t53
      t327 = -t149-12*t272+18*t274-t276+t153+t315-4*t38-t106-4*t317-8*t3
     #19+t322-32*t323+28*t325
      t329 = t241+t255+t270+t327
      t338 = 2-t117+t231+t232+t233-16*t20+128*t22-28*t25-42*t27-t238+106
     #*t31-10*t33+8*t36
      t344 = 110*t49
      t351 = 58*t39-148*t41+238*t43-238*t45-14*t47+t344+130*t51+44*t54-1
     #56*t56+108*t58+264*t60-236*t62+t257
      t358 = -112*t68-24*t70-t187-56*t74-40*t76-t263+68*t80-t265+t266-t2
     #67-t268+t269+t149
      t362 = 2*t38
      t363 = 2*t317
      t364 = 4*t319
      t367 = -6*t98+10*t272-10*t274+t276-t153-t315+t362+t106+t363+t364-t
     #322+20*t323-20*t325
      t369 = t338+t351+t358+t367
      t373 = zqe0(t305,z,l3)
      t374 = dlog(l1)
      t376 = sqe1(t305)
      t377 = zqe1(t305,z,l3)
      t378 = sbqe0(t305)
      t379 = dlog(z)
      t393 = -t159-27*t7+80*t10+t17+38*t18-93*t20-111*t22+132*t25+91*t27
     #+3*t29+10*t31-14*t33
      t406 = -31*t36-33*t39+50*t41-13*t43+9*t45+t132+140*t49+38*t51+131*
     #t54-164*t56+5*t58+111*t60-38*t62
      t416 = 46*t66+4*t68+t143-100*t72-6*t74-110*t76-54*t78-9*t80+6*t82-
     #t85+t87-t89
      t425 = -t91+t94+3*t9-7*t24+3*t2-118*t98-7*t14+t151+44*t102-27*t104
     #+t4+9*t53+t109
      t430 = -t116
      t444 = -6*t5-7*t7+29*t10+t17+22*t18-65*t20-83*t22+100*t25+40*t27-1
     #9*t29-74*t31
      t456 = 47*t33-9*t36-13*t39+34*t41-35*t43+70*t45+2*t47+100*t49+26*t
     #51+41*t54-44*t56-t58
      t465 = 21*t60-t220+34*t66+7*t68+t70-70*t72-t145-80*t76-18*t78-5*t8
     #0-t85+t87
      t468 = 12*t102
      t472 = -t89-t91+t94-t24-50*t98-3*t14+t12+t468-5*t104+4*t4+3*t53+t1
     #09
      t488 = t7-13*t10-t232-4*t18+22*t20+84*t22-t206-66*t27+t29-14*t31+3
     #*t33+22*t36
      t500 = 53*t39-74*t41+9*t43+10*t45-t48-t344-100*t51-48*t54+116*t56-
     #32*t58-144*t60+96*t62-8*t66
      t509 = -13*t68-19*t70+140*t72+38*t74+46*t76+8*t78+34*t80+t265-t266
     #+t267+t268-t269
      t517 = t24+42*t98-6*t272-3*t12-6*t274-t468+2*t104+t362+t363+t364+t
     #322+4*t323+4*t325
      t533 = t379**2
      t536 = sbqe1(t305)
      t538 = t374**2
      t543 = sqe2(t305)
      t544 = zqe2(t305,z,l3)
      t549 = (l1+z-t82)**2
      taux12 = ((128.D0/3.D0*(t3*t112+t116*t155)*s1l1+128.D0/3.D0*(t3*t2
     #01+t116*t227)*s2l1+128.D0/3.D0*t3*t329+128.D0/3.D0*t116*t369)*(t37
     #3*t374+t376-t377+t378*t379)-(128.D0/3.D0*(t3*(t393+t406+t416+t425)
     #+t116*t112+t430*t155)*s1l1+128.D0/3.D0*(t3*(t444+t456+t465+t472)+t
     #116*t201+t430*t227)*s2l1+128.D0/3.D0*t3*(t488+t500+t509+t517)+128.
     #D0/3.D0*t116*t329+128.D0/3.D0*t430*t369)*t373+(128.D0/3.D0*t3*t155
     #*s1l1+128.D0/3.D0*t3*t227*s2l1+128.D0/3.D0*t3*t369)*(-t378*t533/2+
     #t536*t379-t373*t538/2-(t376-t377)*t374+t543-t544))/t549*t287

CCC Coefficients

      tcoff19 = taux11
      tcoff20 = taux11
      tcoff21 = taux12
      tcoff22 = -taux12

CCC Jet functions

      tFjet19 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff19)
      tFjet20 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff20)
      tFjet21 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff21)
      tFjet22 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff22)

CCC Sum

      RVnum6e2 = tFjet19*tcoff19+tFjet20*tcoff20+tFjet21*tcoff21+tFjet22
     #*tcoff22
      RETURN
      END

      double precision function RVnum6e3(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t1 = l1**2
      t2 = t1**2
      t3 = t2*t1
      t5 = t2*l1
      t7 = z*t5
      t9 = t1*l1
      t11 = z**2
      t12 = t11**2
      t13 = t12*l1
      t14 = 768*t13
      t15 = t11*l1
      t17 = t11*z
      t18 = t17*t1
      t20 = t17*l1
      t22 = z*t1
      t24 = t11*t1
      t26 = z*t9
      t28 = t11*t9
      t30 = -896*t3+4352.D0/3.D0*t5-27136.D0/3.D0*t7+640*t9+t14+5120.D0/
     #3.D0*t15+26752.D0/3.D0*t18-7040.D0/3.D0*t20+5632.D0/3.D0*t22-21760
     #.D0/3.D0*t24-17536.D0/3.D0*t26+16896*t28
      t31 = t12*t9
      t33 = t17*t9
      t35 = t12*t1
      t37 = z*t2
      t39 = t17*t2
      t42 = t12*z
      t43 = t42*t9
      t44 = 1536*t43
      t45 = t42*t2
      t46 = 2688*t45
      t47 = t12*t2
      t49 = t17*t3
      t51 = t3*t42
      t53 = t11*t5
      t55 = 640.D0/3.D0+9728*t31-59648.D0/3.D0*t33-3840*t35+28544.D0/3.D
     #0*t37+82048.D0/3.D0*t39-3968.D0/3.D0*t2-t44+t46-14208*t47+12160*t4
     #9+4096.D0/3.D0*t51+21120*t53
      t57 = t17*t5
      t59 = t12*t5
      t61 = t42*t5
      t63 = t3*z
      t65 = t3*t12
      t67 = t42*t1
      t68 = 384*t67
      t69 = t3*t11
      t71 = 256.D0/3.D0*t11
      t72 = 1024.D0/3.D0*t17
      t74 = z*l1
      t76 = t2*t9
      t77 = t42*t76
      t78 = 384*t77
      t79 = t11*t76
      t80 = 3328*t79
      t81 = -70400.D0/3.D0*t57+37504.D0/3.D0*t59-2560*t61+15232.D0/3.D0*
     #t63-19456.D0/3.D0*t65+t68-33664.D0/3.D0*t69-t71+t72-512.D0/3.D0*t1
     #-256*t74-t78+t80
      t82 = t17*t76
      t84 = z*t76
      t86 = t12*t76
      t88 = t2**2
      t89 = t11*t88
      t90 = 1280.D0/3.D0*t89
      t91 = z*t88
      t92 = 640.D0/3.D0*t91
      t93 = t42*t88
      t94 = 128.D0/3.D0*t93
      t95 = t17*t88
      t96 = 1280.D0/3.D0*t95
      t97 = t12*t88
      t98 = 640.D0/3.D0*t97
      t100 = 128.D0/3.D0*t88
      t101 = 128*z
      t102 = 128.D0/3.D0*l1
      t103 = t11*t2
      t105 = -10496.D0/3.D0*t82-4736.D0/3.D0*t84+5504.D0/3.D0*t86-t90+t9
     #2+t94+t96-t98+896.D0/3.D0*t76-t100+t101-t102-72064.D0/3.D0*t103
      t108 = s1l1**2
      t110 = t12*t11
      t111 = t110*t88
      t126 = 320.D0/3.D0*t111+6080.D0/3.D0*t7-2304*t13-128.D0/3.D0*t15-2
     #0864.D0/3.D0*t18+4160.D0/3.D0*t20-128.D0/3.D0*t22+3328.D0/3.D0*t24
     #+320*t26-5376*t28-23232*t31+18560*t33+9536*t35-1280*t37
      t129 = 2*l6-1
      t130 = t129**2
      t133 = l5-1
      t135 = -1+l1
      t137 = (1-l1+t74)**2
      t138 = 1/t137
      t143 = dsqrt(z)
      t145 = dsqrt(l5)
      t147 = dsqrt(-t133)
      t150 = dsqrt(l1)
      t151 = l5*z
      t152 = t151*l1
      t153 = l5*l1
      t154 = t152-l1+1-l5+t153
      t156 = -t135
      t157 = dsqrt(t156)
      t162 = t154**2
      t165 = dsqrt(4*l1*t130*z*l5*t133*t135*t138-4*t129*t143*t145*t147*t
     #150*t154*t138*t157+t162*t138)
      t166 = t165*t1
      t167 = 128.D0/3.D0*t166
      t168 = 384*t12
      t180 = 1664.D0/3.D0-30400*t39+t167+t168+11584*t43-18688*t45+35648*
     #t47-58240.D0/3.D0*t49-32704.D0/3.D0*t51-13120*t53+93184.D0/3.D0*t5
     #7-103808.D0/3.D0*t59+18368*t61-1600*t63+20736*t65
      t186 = 512*t13
      t189 = 1408*t20
      t190 = 768*t22
      t193 = 512-1280*t3+4864.D0/3.D0*t5-32000.D0/3.D0*t7+256*t9+t186+20
     #48.D0/3.D0*t15+19840.D0/3.D0*t18-t189+t190-4608*t24-3968*t26
      t207 = 13824*t28+27392.D0/3.D0*t31-17792*t33-3072*t35+8832*t37+291
     #84*t39-1024*t2-1408*t43+9088.D0/3.D0*t45-47104.D0/3.D0*t47+17664*t
     #49+6016.D0/3.D0*t51+77312.D0/3.D0*t53
      t214 = 256*t67
      t218 = 512.D0/3.D0*t74
      t221 = -88064.D0/3.D0*t57+47872.D0/3.D0*t59-3328*t61+7296*t63-9472
     #*t65+t214-48640.D0/3.D0*t69+512.D0/3.D0*t11+256*t17-t218-640*t77+5
     #632*t79
      t228 = 2560.D0/3.D0*t95
      t229 = 1280.D0/3.D0*t97
      t233 = 256*l1
      t235 = -5888*t82-2688*t84+3072*t86-2560.D0/3.D0*t89+1280.D0/3.D0*t
     #91+256.D0/3.D0*t93+t228-t229+512*t76-256.D0/3.D0*t88+1024.D0/3.D0*
     #z-t233-24320*t103
      t241 = 640.D0/3.D0*t11
      t242 = 128.D0/3.D0*t17
      t243 = 128.D0/3.D0*t1
      t244 = 128*t74
      t252 = (t193+t207+t221+t235)*s2l1-4224*t67+26752.D0/3.D0*t69+t241+
     #t242-t243-t244+10880.D0/3.D0*t77-9920.D0/3.D0*t79+20480.D0/3.D0*t8
     #2+640*t84-7040*t86+1600.D0/3.D0*t89-320.D0/3.D0*t91
      t256 = t110*t76
      t258 = t42*l1
      t259 = 768*t258
      t260 = t110*t5
      t262 = t110*t2
      t264 = t110*t1
      t265 = 384*t264
      t266 = t110*t3
      t268 = t110*t9
      t270 = t165*l1
      t274 = 640.D0/3.D0*l1
      t276 = -1600.D0/3.D0*t93-3200.D0/3.D0*t95+3200.D0/3.D0*t97-2240.D0
     #/3.D0*t256+t259-11200.D0/3.D0*t260+3520*t262+t265+6784.D0/3.D0*t26
     #6-1792*t268+256.D0/3.D0*t270-896.D0/3.D0*t165+384*z-t274+11200*t10
     #3
      t286 = 1280.D0/3.D0*t22
      t287 = 11392.D0/3.D0*t24
      t291 = -1280.D0/3.D0*t3+384*t5-13312.D0/3.D0*t7+t14+640*t15+7424*t
     #18-2176*t20+t286-t287-1792*t26+9728*t28+9600*t31
      t304 = 896.D0/3.D0-47872.D0/3.D0*t33-3968*t35+3712*t37+65408.D0/3.
     #D0*t39-128*t2-1664*t43+2944*t45-40832.D0/3.D0*t47+30848.D0/3.D0*t4
     #9+1408*t51+42112.D0/3.D0*t53-57088.D0/3.D0*t57
      t311 = 128*t11
      t312 = 384*t17
      t314 = 11776*t59-8192.D0/3.D0*t61+9472.D0/3.D0*t63-18304.D0/3.D0*t
     #65+t68-8320*t69+t311+t312+t243-t218-t78+2816*t79
      t321 = -9472.D0/3.D0*t82-3712.D0/3.D0*t84+5248.D0/3.D0*t86-t90+t92
     #+t94+t96-t98+640.D0/3.D0*t76-t100+640.D0/3.D0*z-t274-14720*t103
      t324 = s2l1**2
      t326 = 64.D0/3.D0*t111
      t337 = 640-t326-7616.D0/3.D0*t7+320*t13+1664.D0/3.D0*t15+20672.D0/
     #3.D0*t18-2048*t20+t286-t287-1664*t26+9280*t28+21760.D0/3.D0*t31-14
     #656*t33-7936.D0/3.D0*t35
      t340 = 128*t12
      t345 = 1536*t51
      t351 = 2880*t37+19520*t39-t167+t340+640.D0/3.D0*t43+1088*t45-33344
     #.D0/3.D0*t47+7680*t49+t345+10560*t53-47744.D0/3.D0*t57+30208.D0/3.
     #D0*t59-5824.D0/3.D0*t61+3584.D0/3.D0*t63
      t366 = -16256.D0/3.D0*t65-640*t67-4928*t69+896.D0/3.D0*t11+512*t17
     #+t243-896.D0/3.D0*t74-1856.D0/3.D0*t77+3584.D0/3.D0*t79-6016.D0/3.
     #D0*t82-832.D0/3.D0*t84+4864.D0/3.D0*t86-320.D0/3.D0*t89+64.D0/3.D0
     #*t91
      t370 = 256*t258
      t373 = 128*t264
      t381 = 320.D0/3.D0*t93+640.D0/3.D0*t95-t98+256.D0/3.D0*t256+t370-7
     #04.D0/3.D0*t260+1664.D0/3.D0*t262+t373-64*t266-448*t268+256*t270-3
     #84*t165+1408.D0/3.D0*z-384*l1-12928*t103
      t386 = t9*l6
      t387 = t151*t386
      t389 = l5**2
      t390 = t389*z
      t391 = t390*t386
      t393 = t110*l5
      t394 = l1*l6
      t395 = t393*t394
      t397 = t12*t17
      t398 = t397*t76
      t400 = t151*t394
      t402 = l6**2
      t403 = t1*t402
      t404 = t390*t403
      t406 = t390*t394
      t408 = t1*l6
      t409 = t390*t408
      t411 = t9*t402
      t412 = t151*t411
      t414 = t390*t411
      t416 = -20+13573*t111+32*t387-32*t391+32*t395+8522*t398+64*t400-96
     #*t404-64*t406+96*t409-32*t412+32*t414
      t417 = l1*t402
      t418 = t390*t417
      t420 = t151*t417
      t422 = t151*t403
      t424 = t110*t389
      t425 = t424*t417
      t427 = t393*t417
      t429 = t424*t394
      t431 = t390*l1
      t433 = t151*t408
      t435 = t397*t389
      t436 = t435*t9
      t438 = t393*t1
      t440 = t393*l1
      t443 = t145*l5
      t444 = t143*z*t443
      t445 = t444*t147
      t446 = t150*l1
      t448 = t446*t157*l6
      t449 = t445*t448
      t451 = 64*t418-64*t420+96*t422+32*t425-32*t427-32*t429+24*t431-96*
     #t433+2*t436+4*t438-12*t440+32*t449
      t453 = t165*z
      t454 = t453*t1
      t456 = t143*t42
      t457 = t456*t145
      t459 = t147*t150*t157
      t460 = t457*t459
      t462 = t143*t443
      t464 = t147*t446*t157
      t465 = t462*t464
      t468 = t143*t12*t443
      t469 = t468*t147
      t470 = t469*t448
      t472 = t151*t9
      t474 = t389*l1
      t475 = t474*t17
      t477 = t390*t1
      t479 = t389*t11
      t480 = t479*t1
      t482 = t151*t1
      t484 = t165*t11
      t485 = t484*t1
      t487 = t389*t17
      t488 = t487*t9
      t490 = t389*t12
      t491 = t490*t9
      t493 = -60*t454-8*t460-24*t465+16*t470-12*t472+6*t475-36*t477-4*t4
     #80+36*t482+20*t485+2*t488-8*t491
      t494 = t390*t9
      t496 = t479*t9
      t498 = t487*t1
      t500 = t490*t1
      t502 = t143*t145
      t503 = t502*t147
      t504 = t150*t1
      t506 = t504*t157*l6
      t507 = t503*t506
      t509 = t453*t9
      t511 = t484*t9
      t513 = t462*t147
      t514 = t513*t448
      t516 = t165*t9
      t518 = t456*t443
      t519 = t518*t147
      t521 = t150*t157*l6
      t522 = t519*t521
      t524 = t453*l1
      t526 = t435*t1
      t528 = 12*t494+2*t496-6*t498+12*t500+16*t507+20*t509-10*t511+48*t5
     #14-10*t516-16*t522+40*t524-2*t526
      t531 = t424*l1
      t533 = t424*t9
      t535 = t42*l5
      t536 = t535*t1
      t538 = t535*l1
      t540 = t424*t1
      t542 = t42*t389
      t543 = t542*l1
      t545 = l5*t1
      t546 = t545*t12
      t548 = t542*t9
      t550 = t542*t1
      t552 = t153*t12
      t555 = t469*t521
      t557 = 12*t531-8*t533-8*t536+4*t538+4*t540-2*t543+4*t546+12*t548-8
     #*t550-8*t552+4*t535-16*t555
      t559 = t445*t506
      t561 = t518*t464
      t564 = t513*t521
      t567 = t88*t9
      t568 = t567*t12
      t570 = t567*t42
      t572 = t389*t1
      t574 = t88*l1
      t575 = t574*t397
      t577 = t567*t17
      t580 = -2*t542-16*t559+16*t561+424*t7-32*t564+10*t9-224*t568+280*t
     #570-8*t572+3054*t575+112*t577+10*t474
      t582 = t389*t9
      t595 = 2*t582-119*t13-28*t15-564*t18+96*t20+12*t22+184*t24+68*t26-
     #695*t28-3472*t31+2124*t33+853*t35
      t598 = t397*t88
      t600 = t468*t459
      t603 = t567*t110
      t605 = t397*t3
      t607 = t503*t448
      t610 = t457*t147*t521
      t614 = t147*t504*t157
      t615 = t462*t614
      t617 = t468*t464
      t619 = t574*t42
      t621 = -232*t37-5444*t39-6330*t598+8*t600-4*t389-224*t603-7778*t60
     #5-48*t607+16*t610-14*t42+8*t615-8*t617+8220*t619
      t625 = t574*t12
      t627 = t513*t506
      t629 = t12*l5
      t631 = t574*z
      t635 = t12**2
      t636 = t635*t5
      t638 = t503*t521
      t640 = t574*t11
      t642 = t574*t17
      t644 = t574*t110
      t646 = t635*z
      t647 = t646*t5
      t649 = -6795*t625-16*t627+4*t629+132*t631-2*t487-4*t490-955*t636+3
     #2*t638-1029*t640+3502*t642-6347*t644+66*t647
      t650 = t635*t3
      t652 = t88*t1
      t653 = t652*t635
      t655 = t635*t76
      t657 = t635*t88
      t659 = t110*l1
      t661 = t444*t614
      t663 = t646*t76
      t666 = t143*t110*t443
      t667 = t666*t464
      t669 = t397*t1
      t671 = t502*t614
      t673 = t397*t9
      t675 = t646*t3
      t677 = 1765*t650+244*t653-2118*t655+1669*t657-48*t659+8*t661+222*t
     #663-8*t667-72*t669-8*t671+540*t673-160*t675
      t679 = t646*t88
      t681 = t646*t2
      t683 = t635*t2
      t686 = t635*t9
      t689 = t652*t11
      t691 = t567*t646
      t693 = t462*t459
      t695 = t652*t12
      t699 = -190*t679-12*t681+312*t683+40*t166-48*t686+8*l5+268*t689+4*
     #t691+16*t693+1820*t695+3432*t43-9682*t45
      t709 = t652*t17
      t711 = t652*t110
      t713 = z*t652
      t715 = 9339*t47-12910*t49-26466*t51-3123*t53+9882*t57-17558*t59+18
     #996*t61-534*t63+23597*t65-924*t709+1764*t711-34*t713
      t718 = t502*t464
      t722 = t518*t459
      t725 = t666*t147*t448
      t728 = t652*t397
      t731 = t652*t42
      t733 = t502*t459
      t737 = 24*t718-20*t153-806*t67+8*t722+16*t725+3986*t69-868*t728+16
     #*t545-2240*t731-16*t733-8*t11-16*t17
      t739 = l5*t9
      t751 = -40*t1-4*t739-2*t74+26192*t77-3665*t79+12110*t82+482*t84-22
     #660*t86+2375*t89-308*t91-18030*t93-7974*t95
      t761 = t646*t574
      t764 = t444*t464
      t766 = t567*t11
      t768 = 15215*t97-19085*t256+116*t258-12598*t260+6006*t262+408*t264
     #+18500*t266-1954*t268+100*t761-50*t270-16*t764-32*t766
      t769 = t519*t448
      t771 = t567*t397
      t774 = t567*t635
      t776 = t567*z
      t779 = t646*t652
      t781 = t397*t2
      t783 = t574*t635
      t785 = t397*t5
      t790 = -32*t769+112*t771+20*t165-32*t774+4*t776-24*t152-30*t779-20
     #46*t781-837*t783+4866*t785-12*z+50*l1+1759*t103
      t798 = sbqe0(t156)
      t799 = dlog(z)
      t800 = t799**2
      t803 = sbqe1(t156)
      t805 = zqe0(t156,z,l3)
      t806 = dlog(l1)
      t807 = t806**2
      t810 = sqe1(t156)
      t811 = zqe1(t156,z,l3)
      t812 = t810-t811
      t814 = sqe2(t156)
      t815 = zqe2(t156,z,l3)
      t818 = 448*t3
      t819 = 2176.D0/3.D0*t5
      t821 = 320*t9
      t822 = 832.D0/3.D0*t15
      t828 = -128.D0/3.D0-t326-t818+t819-2368*t7+t821-t186-t822-3776*t18
     #+3904.D0/3.D0*t20+1472.D0/3.D0*t22+3712.D0/3.D0*t24-2816.D0/3.D0*t
     #26
      t829 = 5824.D0/3.D0*t28
      t831 = 17984.D0/3.D0*t33
      t834 = 13184.D0/3.D0*t39
      t835 = 1984.D0/3.D0*t2
      t836 = 768*t43
      t839 = 9088.D0/3.D0*t49
      t842 = 1216.D0/3.D0*t57
      t843 = -t829-12416.D0/3.D0*t31+t831+2176*t35+5312.D0/3.D0*t37-t834
     #-t835+t836-1984.D0/3.D0*t45+10496.D0/3.D0*t47+t839+2816.D0/3.D0*t5
     #1+6976.D0/3.D0*t53-t842
      t850 = 64.D0/3.D0*t11
      t852 = 64.D0/3.D0*t1
      t856 = 1920*t82
      t858 = -128*t59-1024.D0/3.D0*t61+5504.D0/3.D0*t63-6080.D0/3.D0*t65
     #-t214-9344.D0/3.D0*t69-t850-704.D0/3.D0*t17-t852-192*t74-576*t77+1
     #600*t79-t856-2240.D0/3.D0*t84
      t860 = 320*t89
      t861 = 128*t91
      t862 = 128*t93
      t863 = 320*t97
      t864 = 320.D0/3.D0*t256
      t865 = 192*t260
      t866 = 64*t262
      t867 = 640.D0/3.D0*t266
      t868 = 448.D0/3.D0*t76
      t869 = 64.D0/3.D0*t88
      t871 = 4160.D0/3.D0*t86-t860+t861+t862+t96-t863+t864+t865-t866-t86
     #7+t868-t869-t102+512*t103
      t875 = 128.D0/3.D0*t111
      t885 = 1920*t31
      t889 = -t875+3328.D0/3.D0*t7+640*t13-2432.D0/3.D0*t15-10624.D0/3.D
     #0*t18+3968.D0/3.D0*t20-1024.D0/3.D0*t22+9728.D0/3.D0*t24+1408*t26-
     #6144*t28+t885+17152.D0/3.D0*t33-3712.D0/3.D0*t35-1920*t37
      t891 = 640*t3
      t892 = 2432.D0/3.D0*t5
      t894 = 128*t9
      t901 = 768*t26
      t903 = -256-t875-t891+t892-2944*t7+t894-1152*t13-3328.D0/3.D0*t15-
     #25600.D0/3.D0*t18+3072*t20-384*t22+4480*t24+t901-7040*t28
      t909 = 512*t2
      t917 = -7040*t31+36992.D0/3.D0*t33+4480*t35+640*t37-20480.D0/3.D0*
     #t39-t909+t836+640*t45+9856.D0/3.D0*t47+8192*t49+3072*t51+12800.D0/
     #3.D0*t53-12800.D0/3.D0*t57+12800.D0/3.D0*t59
      t929 = 640*t89
      t930 = -2944*t61+3072*t63-6528*t65-t68-6528*t69-t71-640*t17+256.D0
     #/3.D0*t74-1408*t77+t80-4352*t82-1408*t84+3328*t86-t929
      t931 = 256*t91
      t932 = 256*t93
      t933 = 640*t97
      t934 = 256*t256
      t935 = 2432.D0/3.D0*t260
      t936 = 512*t262
      t937 = 640*t266
      t938 = 128*t268
      t939 = 256*t76
      t940 = 512.D0/3.D0*z
      t943 = t931+t932+t228-t933+t934+t935-t936-t937+t938+t939-t100-t940
     #+128*l1+9856.D0/3.D0*t103
      t957 = -640.D0/3.D0-13952.D0/3.D0*t39+(t903+t917+t930+t943)*s2l1-t
     #167-256*t12-11264.D0/3.D0*t43+14464.D0/3.D0*t45-3200*t47+7552.D0/3
     #.D0*t49+8704.D0/3.D0*t51-2816*t53+1024.D0/3.D0*t57+13696.D0/3.D0*t
     #59-13312.D0/3.D0*t61
      t959 = 128*t63
      t963 = 128.D0/3.D0*t11
      t969 = 640.D0/3.D0*t89
      t970 = 128.D0/3.D0*t91
      t971 = -t959-4224*t65+1920*t67-1024.D0/3.D0*t69-t963-t72+t243-1152
     #*t77+768*t79-1792*t82-128*t84+2048*t86-t969+t970
      t972 = 640.D0/3.D0*t93
      t979 = 64.D0/3.D0*t270
      t980 = 256.D0/3.D0*t165
      t982 = t972+t96-t229+t934-512*t258+3712.D0/3.D0*t260-4096.D0/3.D0*
     #t262-256*t264-2176.D0/3.D0*t266+896*t268+t979+t980-t101+t102+18944
     #.D0/3.D0*t103
      t986 = 640.D0/3.D0*t3
      t987 = 192*t5
      t988 = 1024.D0/3.D0*t7
      t992 = 896.D0/3.D0*t22
      t996 = -t326-t986+t987-t988-192*t13-384*t15-11456.D0/3.D0*t18+t189
     #-t992+6400.D0/3.D0*t24+t901-12416.D0/3.D0*t28-5824.D0/3.D0*t31+t83
     #1
      t999 = 64*t12
      t1000 = 64*t2
      t1008 = -640.D0/3.D0+3584.D0/3.D0*t35-1984.D0/3.D0*t37-t834-t999-t
     #1000-2816.D0/3.D0*t43+5312.D0/3.D0*t45+512*t47+t839+5504.D0/3.D0*t
     #51-128*t53-t842+6976.D0/3.D0*t59-2368*t61
      t1020 = 2816.D0/3.D0*t63-9344.D0/3.D0*t65+448*t67-6080.D0/3.D0*t69
     #-t311-320*t17-t243+448.D0/3.D0*t74-2240.D0/3.D0*t77+4160.D0/3.D0*t
     #79-t856-576*t84+1600*t86-t860+t861
      t1021 = 448.D0/3.D0*t256
      t1022 = 128*t258
      t1023 = 2176.D0/3.D0*t260
      t1024 = 1984.D0/3.D0*t262
      t1025 = 64*t264
      t1026 = 448*t266
      t1027 = 320*t268
      t1028 = 320.D0/3.D0*t76
      t1031 = t862+t96-t863+t1021-t1022+t1023-t1024-t1025-t1026+t1027+t1
     #028-t869-t940+512.D0/3.D0*l1+10496.D0/3.D0*t103
      t1037 = 1024.D0/3.D0*t15
      t1039 = 768*t20
      t1040 = 5632.D0/3.D0*t24
      t1046 = -384+t875+3712.D0/3.D0*t7+1408*t13-t1037-1280*t18+t1039-t9
     #92+t1040+896*t26-11264.D0/3.D0*t28+17152.D0/3.D0*t31+1920*t33-3584
     #*t35
      t1048 = 3200*t39
      t1053 = 1024.D0/3.D0*t51
      t1059 = -4096.D0/3.D0*t37-t1048+t167-t168-6144*t43+18944.D0/3.D0*t
     #45-13952.D0/3.D0*t47-4224*t49-t1053-13312.D0/3.D0*t53+13696.D0/3.D
     #0*t57+1024.D0/3.D0*t59-2816*t61-2176.D0/3.D0*t63
      t1070 = 7552.D0/3.D0*t65+3200*t67+8704.D0/3.D0*t69-t241-t312-t243+
     #640.D0/3.D0*t74+768*t77-1152*t79+2048*t82+256*t84-1792*t86+t969-t9
     #70
      t1080 = -t972-t96+t229-128*t256-t259+3328.D0/3.D0*t260-1920*t262-t
     #265-128*t266+1408*t268-192*t270+256*t165-896.D0/3.D0*z+t233+14464.
     #D0/3.D0*t103
      t1085 = 16*t387
      t1086 = 16*t391
      t1087 = 16*t395
      t1089 = 32*t400
      t1090 = 48*t404
      t1091 = 32*t406
      t1092 = 48*t409
      t1093 = 16*t412
      t1094 = 16*t414
      t1095 = 10-14129*t111-t1085+t1086-t1087-8092*t398-t1089+t1090+t109
     #1-t1092+t1093-t1094
      t1107 = -32*t418+32*t420-48*t422-16*t425+16*t427+16*t429-12*t431+4
     #8*t433-t436-2*t438+6*t440-16*t449
      t1110 = 4*t460
      t1111 = 12*t465
      t1112 = 8*t470
      t1113 = 6*t472
      t1114 = 3*t475
      t1115 = 18*t477
      t1116 = 2*t480
      t1117 = 18*t482
      t1119 = 4*t491
      t1120 = 30*t454+t1110+t1111-t1112+t1113-t1114+t1115+t1116-t1117-10
     #*t485-t488+t1119
      t1121 = 6*t494
      t1122 = 3*t498
      t1123 = 6*t500
      t1124 = 8*t507
      t1127 = 24*t514
      t1129 = 8*t522
      t1131 = -t1121-t496+t1122-t1123-t1124-10*t509+5*t511-t1127+5*t516+
     #t1129-20*t524+t526
      t1134 = 6*t531
      t1135 = 4*t533
      t1136 = 4*t536
      t1137 = 2*t538
      t1138 = 2*t540
      t1139 = 2*t546
      t1140 = 6*t548
      t1141 = 4*t550
      t1142 = 4*t552
      t1143 = 2*t535
      t1144 = 8*t555
      t1145 = -t1134+t1135+t1136-t1137-t1138+t543-t1139-t1140+t1141+t114
     #2-t1143+t1144
      t1146 = 8*t559
      t1147 = 8*t561
      t1149 = 16*t564
      t1153 = 4*t572
      t1156 = 5*t474
      t1157 = t542+t1146-t1147-553*t7+t1149-5*t9+280*t568-350*t570+t1153
     #-3235*t575-140*t577-t1156
      t1170 = -t582+151*t13+32*t15+744*t18-132*t20+3*t22-233*t24-94*t26+
     #901*t28+4155*t31-2713*t33-1056*t35
      t1174 = 4*t600
      t1175 = 2*t389
      t1178 = 24*t607
      t1179 = 8*t610
      t1181 = 4*t615
      t1182 = 4*t617
      t1184 = 299*t37+6759*t39+6275*t598-t1174+t1175+280*t603+7313*t605+
     #t1178-t1179+16*t42-t1181+t1182-9325*t619
      t1189 = 8*t627
      t1190 = 2*t629
      t1192 = 2*t490
      t1194 = 16*t638
      t1199 = 7958*t625+t1189-t1190-169*t631+t487+t1192+929*t636-t1194+1
     #280*t640-4229*t642+6964*t644-69*t647
      t1205 = 4*t661
      t1207 = 4*t667
      t1209 = 4*t671
      t1212 = -1608*t650-274*t653+1911*t655-1573*t657+60*t659-t1205-192*
     #t663+t1207+90*t669+t1209-606*t673+146*t675
      t1220 = 4*l5
      t1223 = 8*t693
      t1226 = 170*t679+15*t681-339*t683-20*t166-3*t12+60*t686-t1220-334*
     #t689-5*t691-t1223-2198*t695-3817*t43
      t1240 = 10451*t45-10836*t47+15405*t49+27828*t51+3980*t53-11979*t57
     #+19915*t59-20096*t61+700*t63-26466*t65+1134*t709-2058*t711+43*t713
      t1243 = 12*t718
      t1244 = 10*t153
      t1246 = 4*t722
      t1247 = 8*t725
      t1250 = 8*t545
      t1252 = 8*t733
      t1255 = -t1243+t1244+931*t67-t1246-t1247-5019*t69+994*t728-t1250+2
     #660*t731+t1252+4*t11+17*t17
      t1257 = 2*t739
      t1267 = 20*t1+t1257+t74-27830*t77+4577*t79-14376*t82-630*t84+25477
     #*t86-2955*t89+399*t91+19675*t93+9513*t95
      t1279 = 8*t764
      t1281 = -17375*t97+19155*t256-142*t258+12572*t260-6180*t262-474*t2
     #64-18299*t266+2114*t268-98*t761+25*t270+t1279+40*t766
      t1282 = 16*t769
      t1287 = 12*t152
      t1295 = t1282-140*t771-10*t165+40*t774-5*t776+t1287+33*t779+2100*t
     #781+854*t783-4699*t785+6*z-25*l1-2269*t103
      t1324 = -128.D0/3.D0+960*t7-256*t13-1408.D0/3.D0*t15-8512.D0/3.D0*
     #t18+2368.D0/3.D0*t20-448.D0/3.D0*t22+5504.D0/3.D0*t24+t901-3968*t2
     #8-9472.D0/3.D0*t31+17600.D0/3.D0*t33+1280*t35-1216*t37-7296*t39+51
     #2*t43-2624.D0/3.D0*t45+12928.D0/3.D0*t47-2176*t49
      t1331 = 128*t67
      t1334 = 64*t74
      t1341 = -t1053-11264.D0/3.D0*t53+16192.D0/3.D0*t57-10112.D0/3.D0*t
     #59+768*t61-384*t63+1408*t65-t1331+4480.D0/3.D0*t69+t963-320.D0/3.D
     #0*t17+t243-t1334+64*t77-256*t79+384*t82+64*t84-256*t86-t102+15232.
     #D0/3.D0*t103
      t1346 = 7-3533*t111-t1085+t1086-t1087-2341*t398-t1089+t1090+t1091-
     #t1092+t1093-t1094
      t1350 = 18*t454+t1110+t1111-t1112+t1113-t1114+t1115+t1116-t1117-7*
     #t485-t488+t1119
      t1354 = -t1121-t496+t1122-t1123-t1124-2*t509+2*t511-t1127+t1129-14
     #*t524+t526-t1134+t1135
      t1357 = t1136-t1137-t1138+t543-t1139-t1140+t1141+t1142-t1143+t1144
     #+t542+t1146
      t1365 = -t1147-83*t7+t1149+56*t568-70*t570+t1153-809*t575-28*t577-
     #t1156-t582+22*t13+11*t15
      t1379 = 94*t18-12*t20-7*t22-42*t24-16*t26+137*t28+708*t31-386*t33-
     #164*t35+50*t37+1047*t39+1713*t598
      t1387 = -t1174+t1175+t165*t2+t484*t2+56*t603+2153*t605+t1178-t1179
     #+4*t42-t1181+t1182-2055*t619+1648*t625
      t1398 = t1189-t1190-29*t631+t487+t1192+268*t636-t1194+234*t640-823
     #*t642+1634*t644-18*t647-504*t650
      t1409 = -64*t653+604*t655-2*t453*t2-468*t657+12*t659-t1205-65*t663
     #+t1207+18*t669+t1209-144*t673+46*t675
      t1419 = 55*t679+3*t681-84*t683-10*t166-t2+12*t686-t1220-64*t689-t6
     #91-t1223-448*t695-832*t43
      t1433 = 2343*t45-1995*t47+2724*t49+6524*t51+598*t53-1996*t57+3891*
     #t59-4636*t61+101*t63-5389*t65+224*t709-448*t711+8*t713
      t1441 = -t1243+t1244+201*t67-t1246-t1247-785*t69+224*t728-t1250+56
     #0*t731+t1252+t11+2*t17
      t1453 = 10*t1+t1257+4*t74-6508*t77+753*t79-2657*t82-93*t84+5302*t8
     #6-513*t89+63*t91+4500*t93+1813*t95
      t1466 = -3630*t97+5005*t256-31*t258+3323*t260-1589*t262-108*t264-4
     #870*t266+520*t268-28*t761+16*t270+t1279+8*t766
      t1477 = t1282-28*t771-7*t165+8*t774-t776+t1287+8*t779+561*t781+228
     #*t783-1347*t785+3*z-16*l1-335*t103
      t1490 = 128*t26
      t1497 = -t875-1472.D0/3.D0*t7+896*t13-704.D0/3.D0*t15+7040.D0/3.D0
     #*t18-1600.D0/3.D0*t20+128*t24+t1490+3136.D0/3.D0*t28+8768*t31-1868
     #8.D0/3.D0*t33-11072.D0/3.D0*t35+320.D0/3.D0*t37+31232.D0/3.D0*t39
      t1509 = -448.D0/3.D0-t167-t340-13120.D0/3.D0*t43+21440.D0/3.D0*t45
     #-13312*t47+7104*t49+12736.D0/3.D0*t51+12352.D0/3.D0*t53-10944*t57+
     #38848.D0/3.D0*t59-21248.D0/3.D0*t61+512*t63-7872*t65
      t1513 = 640*t7
      t1518 = 512*t26
      t1522 = 256*t35
      t1523 = 2432.D0/3.D0*t37
      t1529 = -128+t1513-t1037-1024*t18+128*t20-256.D0/3.D0*t22+3200.D0/
     #3.D0*t24+t1518-6784.D0/3.D0*t28-1152*t31+8192.D0/3.D0*t33+t1522-t1
     #523-12032.D0/3.D0*t39+128*t43-384*t45+6400.D0/3.D0*t47-1408*t49
      t1535 = 256*t63
      t1542 = 128.D0/3.D0*t84
      t1544 = 128.D0/3.D0*z
      t1546 = -640.D0/3.D0*t51-7168.D0/3.D0*t53+9856.D0/3.D0*t57-5888.D0
     #/3.D0*t59+1280.D0/3.D0*t61-t1535+896*t65+2944.D0/3.D0*t69+t963+t24
     #3-128.D0/3.D0*t74+128.D0/3.D0*t77-512.D0/3.D0*t79+256*t82+t1542-51
     #2.D0/3.D0*t86-t1544+3072*t103
      t1549 = 64*t1
      t1555 = 1536*t67-9280.D0/3.D0*t69+(t1529+t1546)*s2l1+t850+t242+t15
     #49-t1334-4288.D0/3.D0*t77+3712.D0/3.D0*t79-7808.D0/3.D0*t82-704.D0
     #/3.D0*t84+8192.D0/3.D0*t86-t969+t970
      t1564 = t972+t96-t229+896.D0/3.D0*t256-t370+4352.D0/3.D0*t260-3968
     #.D0/3.D0*t262-t373-896*t266+640*t268+t979+t980-64*z-64.D0/3.D0*l1-
     #9088.D0/3.D0*t103
      t1568 = 128*t13
      t1569 = 256*t15
      t1571 = 384*t20
      t1572 = 128*t22
      t1581 = 256*t43
      t1582 = 1216.D0/3.D0*t45
      t1585 = -256.D0/3.D0+t988-t1568-t1569-3904.D0/3.D0*t18+t1571-t1572
     #+896*t24+384*t26-5248.D0/3.D0*t28-4352.D0/3.D0*t31+7616.D0/3.D0*t3
     #3+640*t35-1472.D0/3.D0*t37-2944*t39+t1581-t1582+5504.D0/3.D0*t47-7
     #68*t49
      t1586 = 128*t51
      t1590 = 320*t61
      t1592 = 64*t67
      t1594 = 64*t17
      t1595 = 64.D0/3.D0*t74
      t1596 = 64.D0/3.D0*t77
      t1597 = 256.D0/3.D0*t79
      t1598 = 128*t82
      t1599 = 64.D0/3.D0*t84
      t1600 = 256.D0/3.D0*t86
      t1602 = -t1586-4096.D0/3.D0*t53+6080.D0/3.D0*t57-3968.D0/3.D0*t59+
     #t1590-t959+512*t65-t1592+512*t69-t1594+t1595+t1596-t1597+t1598+t15
     #99-t1600-t1544+t102+6016.D0/3.D0*t103
      t1605 = 1280.D0/3.D0*t7
      t1607 = 256*t20
      t1619 = -4160.D0/3.D0*t31+2560*t33+576*t35-1856.D0/3.D0*t37-t1048+
     #t1581-t1582+5632.D0/3.D0*t47-832*t49-t1586-1600*t53
      t1632 = t1595+t1596-t1597+t1598+t1599-t1600-64*t270+128*t165-320.D
     #0/3.D0*z+64*l1+7040.D0/3.D0*t103
      t1642 = sbqe2(t156)
      t1655 = 1536*t18
      t1658 = -t326-t818+t819-3968*t7+t821+448*t15+t1655-t1607+t190-t104
     #0-6656.D0/3.D0*t26+4928*t28
      t1659 = 14080.D0/3.D0*t33
      t1662 = 25856.D0/3.D0*t39
      t1665 = 20224.D0/3.D0*t49
      t1667 = 29440.D0/3.D0*t57
      t1668 = t885-t1659-384*t35+11392.D0/3.D0*t37+t1662-t835-t1581+1024
     #*t45-4480*t47+t1665+t345+26176.D0/3.D0*t53-t1667
      t1678 = 2560*t82
      t1680 = 5888*t59-5248.D0/3.D0*t61+7424.D0/3.D0*t63-13376.D0/3.D0*t
     #65-5632*t69-64*t11-t1549-t244-2048.D0/3.D0*t77+6080.D0/3.D0*t79-t1
     #678-2560.D0/3.D0*t84
      t1683 = 5440.D0/3.D0*t86-t860+t861+t862+t96-t863+t864+t865-t866-t8
     #67+t868-t869-24704.D0/3.D0*t103
      t1700 = -t875-t891+t892-15232.D0/3.D0*t7+t894-t1568-128*t15+640.D0
     #/3.D0*t18+t1571+t1572-256*t24-1152*t26
      t1712 = 9856.D0/3.D0*t28+9856.D0/3.D0*t31-13568.D0/3.D0*t33-t1522+
     #10496.D0/3.D0*t37+38144.D0/3.D0*t39-t909-1152*t43+10496.D0/3.D0*t4
     #5-9344*t47+13056*t49+11648.D0/3.D0*t51+39040.D0/3.D0*t53
      t1720 = 128*t17
      t1726 = -52480.D0/3.D0*t57+39040.D0/3.D0*t59-15232.D0/3.D0*t61+116
     #48.D0/3.D0*t63-29312.D0/3.D0*t65+t1331-29312.D0/3.D0*t69-t1720-153
     #6*t77+3840*t79-5120*t82-1536*t84+3840*t86
      t1728 = -t929+t931+t932+t228-t933+t934+t935-t936-t937+t938+t939-t1
     #00-9344*t103
      t1735 = t1513-384*t13-t1569-8192.D0/3.D0*t18+t1039-t1572+1280*t24+
     #t1518-9088.D0/3.D0*t28-4736*t31+17536.D0/3.D0*t33+1920*t35-t1523-2
     #2784.D0/3.D0*t39+(t1700+t1712+t1726+t1728)*s2l1+t44-t46+6656*t47-7
     #168.D0/3.D0*t49-3328.D0/3.D0*t51
      t1742 = 640.D0/3.D0*t77
      t1743 = 640.D0/3.D0*t79
      t1744 = 1280.D0/3.D0*t82
      t1745 = 1280.D0/3.D0*t86
      t1746 = 128.D0/3.D0*t256
      t1750 = -3072*t53+17408.D0/3.D0*t57-5376*t59+2432*t61-t1535+2304*t
     #65-t68+3712.D0/3.D0*t69-t1720+t1742-t1743+t1744+t1542-t1745-t1746-
     #1280.D0/3.D0*t260+384*t262+t867-t938+12160.D0/3.D0*t103
      t1755 = 384*t24
      t1759 = -t326-t986+t987-5248.D0/3.D0*t7+448*t13+t1655-t1607-t1755-
     #256*t26+1920*t28+4928*t31-t1659
      t1768 = -5632.D0/3.D0*t35+1024*t37+t1662-t999-t1000-6656.D0/3.D0*t
     #43+11392.D0/3.D0*t45-24704.D0/3.D0*t47+t1665+7424.D0/3.D0*t51+5888
     #*t53-t1667+26176.D0/3.D0*t59
      t1779 = -3968*t61+1536*t63-5632*t65+768*t67-13376.D0/3.D0*t69-2560
     #.D0/3.D0*t77+5440.D0/3.D0*t79-t1678-2048.D0/3.D0*t84+6080.D0/3.D0*
     #t86-t860+t861
      t1781 = t862+t96-t863+t1021-t1022+t1023-t1024-t1025-t1026+t1027+t1
     #028-t869-4480*t103
      t1798 = -t1605+t14+1920*t18-t1571-t1755-t1490+1536*t28+17536.D0/3.
     #D0*t31-4736*t33-8192.D0/3.D0*t35+384*t37+6656*t39-t340-9088.D0/3.D
     #0*t43+12160.D0/3.D0*t45-22784.D0/3.D0*t47+2304*t49+3712.D0/3.D0*t5
     #1+2432*t53
      t1811 = -5376*t57+17408.D0/3.D0*t59-3072*t61+640.D0/3.D0*t63-7168.
     #D0/3.D0*t65+1280*t67-3328.D0/3.D0*t69-t1742+t1743-t1744-t1542+t174
     #5+t1746-t370+640*t260-2432.D0/3.D0*t262-t373-256*t266+512*t268-268
     #8*t103
      t1825 = 3443*t111+1770*t398+120*t7-112*t568+140*t570+908*t575+56*t
     #577-147*t13-12*t15-398*t18+75*t20
      t1838 = -6*t22+93*t24+33*t26-331*t28-2159*t31+1214*t33+722*t35-80*
     #t37-2407*t39-1477*t598-112*t603-1784*t605
      t1852 = -6*t42+2700*t619-2326*t625+50*t631-275*t636-378*t640+1244*
     #t642-1990*t644+21*t647+372*t650+92*t653-389*t655
      t1865 = 351*t657-24*t659+35*t663-36*t669+258*t673-32*t675-35*t679-
     #6*t681+123*t683+15*t12-24*t686+116*t689
      t1880 = 2*t691+756*t695+2043*t43-4313*t45+4284*t47-3646*t49-7057*t
     #51-1000*t53+3369*t57-6045*t59+6335*t61-137*t63
      t1893 = 6574*t65-392*t709+700*t711-15*t713-620*t67+1099*t69-336*t7
     #28-910*t731-6*t17+6478*t77-1007*t79+3282*t82
      t1907 = 131*t84-5927*t86+737*t89-98*t91-4895*t93-2391*t95+4365*t97
     #-4373*t256+105*t258-3964*t260+2498*t262+246*t264
      t1920 = 4611*t266-1034*t268+26*t761-16*t766+56*t771-16*t774+2*t776
     #-11*t779-798*t781-234*t783+1439*t785+699*t103
      t1930 = (l1+z-t74)**2
      s2 = ((t30+t55+t81+t105)*t108+(t126+t180+t252+t276)*s1l1+(t291+t30
     #4+t314+t321)*t324+(t337+t351+t366+t381)*s2l1-64.D0/3.D0*(t416+t451
     #+t493+t528+t557+t580+t595+t621+t649+t677+t699+t715+t737+t751+t768+
     #t790)*t138)*(-t798*t800/2+t803*t799-t805*t807/2-t812*t806+t814-t81
     #5)+((t828+t843+t858+t871)*t108+(t889+t957+t971+t982)*s1l1+(t996+t1
     #008+t1020+t1031)*t324+(t1046+t1059+t1070+t1080)*s2l1-64.D0/3.D0*(t
     #1095+t1107+t1120+t1131+t1145+t1157+t1170+t1184+t1199+t1212+t1226+t
     #1240+t1255+t1267+t1281+t1295)*t138)*(t805*t806+t810-t811+t798*t799
     #)
      s1 = s2+((t1324+t1341)*t108-64.D0/3.D0*(t1346+t1107+t1350+t1354+t1
     #357+t1365+t1379+t1387+t1398+t1409+t1419+t1433+t1441+t1453+t1466+t1
     #477)*t138+(t1497+t1509+t1555+t1564)*s1l1+(t1585+t1602)*t324+(-192+
     #t1605-t1568-t822-3520.D0/3.D0*t18+t1607-320.D0/3.D0*t22+2944.D0/3.
     #D0*t24+448*t26-t829+t1619+2240*t57-4160.D0/3.D0*t59+t1590-448.D0/3
     #.D0*t63+1600.D0/3.D0*t65-t1592+576*t69-t850-t1594+t852+t1632)*s2l1
     #)*(t798*t800*t799/6-t803*t800/2+t1642*t799+t805*t807*t806/6+t812*t
     #807/2-(t814-t815)*t806)-((t1658+t1668+t1680+t1683)*t108+(t1735+t17
     #50)*s1l1+(t1759+t1768+t1779+t1781)*t324+(t1798+t1811)*s2l1-64.D0/3
     #.D0*(t1825+t1838+t1852+t1865+t1880+t1893+t1907+t1920)*t138)*t805
      s2 = 1/t1930*t138
      taux11 = s1*s2
      t2 = l1**2
      t3 = 3-3*l1+t2
      t4 = z**2
      t5 = t2*l1
      t6 = t4*t5
      t8 = z*t5
      t10 = t4**2
      t11 = t10*t4
      t12 = t2**2
      t13 = t12*t5
      t14 = t11*t13
      t15 = 2*t14
      t16 = t4*t2
      t18 = z*t2
      t20 = t4*z
      t21 = t20*l1
      t23 = t20*t2
      t25 = t4*l1
      t27 = z*t12
      t29 = t10*t2
      t31 = t20*t5
      t33 = t10*t5
      t35 = -93*t6+38*t8+t15+80*t16-27*t18+44*t21-118*t23-27*t25+3*t27+9
     #1*t29+132*t31-111*t33
      t36 = t12*l1
      t37 = t10*t36
      t39 = t11*t36
      t41 = t20*t36
      t43 = t4*t36
      t45 = t12*t2
      t46 = t10*z
      t47 = t45*t46
      t49 = t20*t45
      t50 = 140*t49
      t51 = t46*l1
      t52 = 8*t51
      t53 = t10*t12
      t55 = t46*t12
      t57 = t46*t5
      t59 = t46*t2
      t61 = t10*l1
      t63 = t4*t12
      t65 = 111*t37+5*t39-164*t41+131*t43+38*t47+t50+t52+9*t53-13*t55+50
     #*t57-33*t59-31*t61-14*t63
      t67 = t20*t12
      t69 = t11*t5
      t71 = z*t36
      t73 = t45*t4
      t75 = t11*t45
      t77 = t45*t10
      t79 = t11*t2
      t80 = 4*t79
      t81 = t11*t12
      t83 = t45*z
      t85 = t46*t36
      t87 = t10*t13
      t88 = 30*t87
      t89 = z*t13
      t90 = 12*t89
      t91 = 10*t67-9*t69-54*t71-110*t73-6*t75-100*t77+t80+4*t81+46*t83-3
     #8*t85+t88-t90
      t92 = t20*t13
      t93 = 40*t92
      t94 = t4*t13
      t95 = 30*t94
      t96 = t46*t13
      t97 = 12*t96
      t98 = z*l1
      t103 = 2*t13
      t105 = 4*t10
      t106 = 8*t45
      t108 = -t93+t95-t97+6*t98+3*t2+3*t4-7*t20+t103-7*t5+t105+t12-t106+
     #9*t36
      t110 = t35+t65+t91+t108
      t114 = -12+12*l1-4*t2
      t126 = -2+326*t6-94*t8+t15-128*t16+30*t18-58*t21+242*t23+26*t25+15
     #4*t27-218*t29-552*t31+484*t33
      t133 = 16*t51
      t140 = 480*t37+38*t39-580*t41+390*t43+80*t47+240*t49-t133-634*t53+
     #280*t55-206*t57+84*t59+56*t61-472*t63
      t152 = 742*t67+32*t69-138*t71-170*t73-14*t75-190*t77-8*t79-50*t81+
     #64*t83-210*t85+t88-t90-t93
      t154 = 2*z
      t157 = 4*t20
      t163 = t95-t97-4*t98-t154-2*t2-4*t4+t157+t103+10*t5-8*t10-20*t12-1
     #0*t45+20*t36
      t165 = t126+t140+t152+t163
      t167 = -t114
      t169 = 10*t8
      t170 = 26*t16
      t171 = 2*t18
      t174 = 4*t25
      t177 = 176*t31
      t180 = 16*t39
      t185 = 2-72*t6+t169+t170-t171+18*t21-78*t23-t174-18*t27+90*t29+t17
     #7-192*t33-152*t37-t180+148*t41-72*t43-20*t47-40*t49+t52
      t196 = 4*t75
      t201 = 2*t4
      t202 = 234*t53-120*t55+94*t57-40*t59-24*t61+102*t63-222*t67-16*t69
     #+14*t71+20*t73+t196+40*t77+t80+24*t81-4*t83+78*t85+t154+t201+t105
      t203 = t185+t202
      t211 = 12*t21
      t217 = -65*t6+22*t8+t15+29*t16-7*t18+t211-50*t23-5*t25-19*t27+40*t
     #29+100*t31
      t224 = 2*t51
      t229 = -83*t33+21*t37-t39-44*t41+41*t43+26*t47+100*t49+t224+70*t53
     #-35*t55+34*t57-13*t59
      t240 = 2*t85
      t241 = -9*t61+47*t63-74*t67-5*t69-18*t71-80*t73-t196-70*t77+t79+7*
     #t81+34*t83-t240
      t242 = 3*t5
      t246 = t88-t90-t93+t95-t97-t20+t103-t242+t10+4*t12-6*t45+3*t36
      t248 = t217+t229+t241+t246
      t260 = -2+196*t6-44*t8+t15-74*t16+10*t18-40*t21+158*t23+14*t25+86*
     #t27-138*t29-360*t31
      t273 = 316*t33+356*t37+28*t39-424*t41+276*t43+68*t47+200*t49-t52-4
     #34*t53+188*t55-128*t57+48*t59+34*t61
      t281 = 160*t77
      t285 = -302*t63+502*t67+18*t69-92*t71-140*t73-12*t75-t281-t80-32*t
     #81+52*t83-156*t85+t88
      t289 = -t90-t93+t95-t97-t154-t201+t157+t103+2*t5-t105-8*t12-t106+1
     #2*t36
      t291 = t260+t273+t285+t289
      t300 = 56*t31
      t303 = 2-40*t6+t169+t170-4*t18+14*t21-38*t23-8*t25-8*t27+18*t29+t3
     #00-32*t33-8*t37
      t309 = 2*t59
      t310 = 4*t61
      t315 = 12*t41-8*t43+26*t53-6*t55+6*t57-t309-t310+30*t63-42*t67+2*t
     #71+t240+t154+t201-2*t20
      t316 = t303+t315
      t322 = 12*t14
      t328 = 22*t6-4*t8-t322-13*t16+t18-t211+42*t23+2*t25+t27-66*t29-t30
     #0+84*t33
      t334 = 110*t49
      t341 = -144*t37-32*t39+116*t41-48*t43-100*t47-t334-t133+10*t53+9*t
     #55-74*t57+53*t59+22*t61+3*t63
      t350 = 13*t81
      t351 = 8*t83
      t353 = 40*t87
      t354 = 2*t89
      t355 = -14*t67+34*t69+8*t71+46*t73+38*t75+140*t77-19*t79-t350-t351
     #+96*t85-t353+t354
      t356 = 30*t92
      t357 = 12*t94
      t358 = 30*t96
      t359 = t10*t20
      t360 = t359*t45
      t362 = t359*t5
      t364 = t359*t2
      t365 = 2*t364
      t367 = 2*t46
      t368 = t359*t36
      t370 = t359*t12
      t373 = 2*t359*t13
      t374 = t11*l1
      t375 = 4*t374
      t376 = t356-t357+t358-6*t360-6*t362+t20+t365-3*t10+t367+4*t368+4*t
     #370+t373+t375
      t378 = t328+t341+t355+t376
      t381 = 6*t8
      t382 = 14*t16
      t383 = 6*t21
      t384 = 4*t27
      t388 = -2+14*t6-t381-t322-t382+t171-t383+t174+t384+72*t29+52*t31-2
     #00*t33
      t394 = 120*t49
      t402 = -332*t37-146*t39+188*t41-50*t43-150*t47-t394+26*t51+334*t53
     #-346*t55+236*t57-100*t59-14*t61+16*t63
      t406 = 4*t71
      t412 = 6*t83
      t414 = -146*t67-114*t69+t406+42*t73+66*t75+180*t77+44*t79+170*t81-
     #t412+308*t85-t353+t354+t356
      t420 = 2*l6-1
      t421 = t420**2
      t424 = l5-1
      t426 = -1+l1
      t428 = (1-l1+t98)**2
      t429 = 1/t428
      t434 = dsqrt(z)
      t436 = dsqrt(l5)
      t438 = dsqrt(-t424)
      t441 = dsqrt(l1)
      t445 = l5*z*l1-l1+1-l5+l5*l1
      t447 = -t426
      t448 = dsqrt(t447)
      t453 = t445**2
      t456 = dsqrt(4*l1*t421*z*l5*t424*t426*t429-4*t420*t434*t436*t438*t
     #441*t445*t429*t448+t453*t429)
      t457 = 2*t456
      t461 = -t357+t358-t154-12*t360-t201+18*t362-4*t364-4*t46+t457+28*t
     #368-32*t370+t373-8*t374
      t463 = t388+t402+t414+t461
      t470 = 2-16*t6+t381+t322+t382-t171+t383-6*t23-t174-t384-42*t29-28*
     #t31+128*t33
      t483 = 264*t37+108*t39-156*t41+44*t43+130*t47+t334-14*t51-238*t53+
     #238*t55-148*t57+58*t59+8*t61-10*t63
      t492 = 106*t67+68*t69-t406-40*t73-56*t75-t281-24*t79-112*t81+t412-
     #236*t85+t353-t354-t356
      t497 = t357-t358+t154+10*t360+t201-10*t362+t365+t367-t457-20*t368+
     #20*t370-t373+t375
      t499 = t470+t483+t492+t497
      t503 = zqe0(t447,z,l3)
      t504 = dlog(l1)
      t506 = sqe1(t447)
      t507 = zqe1(t447,z,l3)
      t508 = sbqe0(t447)
      t509 = dlog(z)
      t520 = t509**2
      t524 = sbqe1(t447)
      t527 = sbqe2(t447)
      t529 = t504**2
      t533 = t506-t507
      t536 = sqe2(t447)
      t537 = zqe2(t447,z,l3)
      t561 = -t177+127*t33+251*t37+t180-334*t41+246*t43+54*t47+180*t49-2
     #41*t53+90*t55-41*t57
      t572 = 8*t59+t310-231*t63+324*t67+4*t69-95*t71-135*t73-9*t75-135*t
     #77-t350+54*t83
      t577 = -99*t85+t88-t90-t93+t95-t97+t103+t242-11*t12-9*t45+15*t36
      t600 = 33*t6-7*t8+t15-5*t16-t21+12*t23+30*t27-9*t29-t300+42*t33+14
     #6*t37+9*t39-194*t41+141*t43+42*t47+t50-105*t53+38*t55-13*t57+t309
      t613 = t61-95*t63+140*t67+t69-53*t71-105*t73-7*t75-105*t77-5*t81+4
     #2*t83-57*t85+t88-t90-t93+t95-t97+t103-3*t12-7*t45+8*t36
      t637 = -254*t37-75*t39+186*t41-71*t43-120*t47-t394+t224+212*t53-15
     #8*t55+73*t57-18*t59
      t647 = -3*t61+50*t63-148*t67-23*t69+11*t71+48*t73+48*t75+t281+t80+
     #58*t81-t351
      t653 = 191*t85-t353+t354+t356-t357+t358-8*t360+2*t362+12*t368-8*t3
     #70+t373
      t684 = (l1+z-t98)**2
      s2 = (128.D0/3.D0*(t3*t110+t114*t165+t167*t203)*s1l1+128.D0/3.D0*(
     #t3*t248+t114*t291+t167*t316)*s2l1+128.D0/3.D0*t3*t378+128.D0/3.D0*
     #t114*t463+128.D0/3.D0*t167*t499)*(t503*t504+t506-t507+t508*t509)+(
     #128.D0/3.D0*t3*t203*s1l1+128.D0/3.D0*t3*t316*s2l1+128.D0/3.D0*t3*t
     #499)*(t508*t520*t509/6-t524*t520/2+t527*t509+t503*t529*t504/6+t533
     #*t529/2-(t536-t537)*t504)
      s1 = s2-(128.D0/3.D0*(t3*(118*t6-35*t8+t15-31*t16+6*t18-7*t21+52*t
     #23+3*t25+82*t27-35*t29+t561+t572+t577)+t114*t110+t167*t165)*s1l1+1
     #28.D0/3.D0*(t3*(t600+t613)+t114*t248+t167*t291)*s2l1+128.D0/3.D0*t
     #3*(-17*t6+t8-t322+2*t16+t21-14*t23-6*t27+26*t29+64*t31-100*t33+t63
     #7+t647+t653)+128.D0/3.D0*t114*t378+128.D0/3.D0*t167*t463)*t503+(12
     #8.D0/3.D0*(t3*t165+t114*t203)*s1l1+128.D0/3.D0*(t3*t291+t114*t316)
     #*s2l1+128.D0/3.D0*t3*t463+128.D0/3.D0*t114*t499)*(-t508*t520/2+t52
     #4*t509-t503*t529/2-t533*t504+t536-t537)
      s2 = 1/t684*t429
      taux12 = s1*s2

CCC Coefficients

      tcoff19 = taux11
      tcoff20 = taux11
      tcoff21 = taux12
      tcoff22 = -taux12

CCC Jet functions

      tFjet19 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff19)
      tFjet20 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff20)
      tFjet21 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff21)
      tFjet22 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff22)

CCC Sum

      RVnum6e3 = tFjet19*tcoff19+tFjet20*tcoff20+tFjet21*tcoff21+tFjet22
     #*tcoff22
      RETURN
      END

      double precision function RVnum6e4(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t1 = l1**2
      t2 = t1**2
      t3 = t2**2
      t4 = t3*l1
      t5 = t4*z
      t7 = z**2
      t8 = t7**2
      t9 = t8*t7
      t10 = t4*t9
      t12 = t4*t7
      t14 = t7*z
      t15 = t4*t14
      t17 = t8*t14
      t18 = t2*l1
      t19 = t17*t18
      t21 = t1*l1
      t22 = t2*t21
      t23 = t17*t22
      t25 = t17*t2
      t27 = t3*t1
      t28 = t27*t14
      t30 = z*t27
      t32 = t9*t3
      t34 = t14*t1
      t36 = t8**2
      t37 = t36*z
      t38 = t37*t27
      t40 = t37*t4
      t42 = t14*l1
      t44 = t7*t1
      t46 = t9*t22
      t48 = t2*t1
      t49 = t17*t48
      t51 = t27*t8
      t53 = t9*t21
      t55 = -16*t5+756*t10+124*t12-420*t15-846*t19-1624*t23+264*t25+56*t
     #28+2*t30-2256*t32+72*t34+2*t38-12*t40-6*t42-12*t44+3800*t46+1516*t
     #49-112*t51+240*t53
      t56 = z*t18
      t58 = z*t21
      t60 = t7*t21
      t62 = t8*t21
      t64 = t14*t21
      t66 = t8*t1
      t68 = t17*t21
      t70 = z*t2
      t72 = t14*t2
      t74 = t36*t18
      t76 = t36*t2
      t78 = t4*t8
      t80 = t4*t17
      t82 = t7*t3
      t84 = t8*z
      t85 = t84*t22
      t87 = t7*t22
      t89 = t14*t22
      t91 = z*t22
      t93 = t27*t36
      t95 = -92*t56-6*t58+96*t60+688*t62-392*t64-132*t66-36*t68+36*t70+1
     #234*t72+132*t74-24*t76+812*t78-364*t80-424*t82-5420*t85+824*t87-26
     #64*t89-110*t91-16*t93
      t97 = t4*t36
      t99 = t8*t22
      t101 = z*t3
      t103 = t27*t7
      t105 = t37*t48
      t107 = t37*t3
      t109 = t37*t22
      t111 = t36*t3
      t113 = t36*t22
      t115 = t36*t48
      t117 = t37*t18
      t119 = t7*t2
      t121 = t8*l1
      t123 = t84*t1
      t125 = t4*t84
      t127 = t27*t17
      t129 = t27*t9
      t131 = t84*t21
      t133 = t84*t2
      t135 = 100*t97+4848*t99+56*t101-16*t103+24*t105+30*t107-38*t109-26
     #8*t111+384*t113-308*t115-6*t117-352*t119+12*t121+96*t123-980*t125+
     #56*t127-112*t129-590*t131+2022*t133
      t136 = t8*t2
      t138 = t84*l1
      t140 = t14*t48
      t142 = t48*t84
      t144 = t7*t18
      t146 = t14*t18
      t148 = t9*t18
      t150 = t14*t3
      t152 = t8*t3
      t154 = t84*t3
      t156 = t27*t84
      t158 = t8*t18
      t160 = t84*t18
      t162 = t17*t3
      t164 = t48*z
      t166 = t9*t2
      t168 = t9*t1
      t170 = t48*t8
      t172 = t9*t48
      t174 = t48*t7
      t176 = -2140*t136-6*t138+3176*t140+6018*t142+748*t144-2454*t146+25
     #68*t148+1398*t150-2620*t152+3050*t154+140*t156+4280*t158-4330*t160
     #+1034*t162+130*t164-1040*t166-24*t168-5636*t170-3932*t172-988*t174
      t179 = z*l1
      t181 = (1-l1+t179)**2
      t182 = 1/t181
      t184 = 1-l1
      t185 = zqe0(t184,z,l3)
      t191 = t7*l1
      t194 = 128.D0/3.D0*t3
      t195 = 1024.D0/3.D0*t14
      t199 = 256.D0/3.D0*t7
      t201 = 26752.D0/3.D0*t34-7040.D0/3.D0*t42-21760.D0/3.D0*t44+5120.D
     #0/3.D0*t191-27136.D0/3.D0*t56-t194+t195-512.D0/3.D0*t1+640*t21+896
     #.D0/3.D0*t22-t199-17536.D0/3.D0*t58
      t212 = 1280.D0/3.D0*t82
      t213 = 384*t85
      t214 = 3328*t87
      t215 = 16896*t60+9728*t62-59648.D0/3.D0*t64-3840*t66+28544.D0/3.D0
     #*t70+82048.D0/3.D0*t72-3968.D0/3.D0*t2+4352.D0/3.D0*t18-256*t179-8
     #96*t48-t212-t213+t214
      t220 = 640.D0/3.D0*t101
      t222 = 768*t121
      t223 = 384*t123
      t224 = 128.D0/3.D0*l1
      t225 = 128*z
      t226 = z*t1
      t228 = 1536*t131
      t229 = 2688*t133
      t230 = 640.D0/3.D0-10496.D0/3.D0*t89-4736.D0/3.D0*t91+5504.D0/3.D0
     #*t99+t220-72064.D0/3.D0*t119+t222+t223-t224+t225+5632.D0/3.D0*t226
     #-t228+t229
      t236 = 1280.D0/3.D0*t150
      t237 = 640.D0/3.D0*t152
      t238 = 128.D0/3.D0*t154
      t244 = -14208*t136+12160*t140+4096.D0/3.D0*t142+21120*t144-70400.D
     #0/3.D0*t146+t236-t237+t238+37504.D0/3.D0*t158-2560*t160+15232.D0/3
     #.D0*t164-19456.D0/3.D0*t170-33664.D0/3.D0*t174
      t247 = s1l1**2
      t255 = 2*l6-1
      t256 = t255**2
      t259 = l5-1
      t266 = dsqrt(z)
      t268 = dsqrt(l5)
      t270 = dsqrt(-t259)
      t273 = dsqrt(l1)
      t274 = l5*z
      t275 = t274*l1
      t276 = l5*l1
      t277 = t275-l1+1-l5+t276
      t279 = dsqrt(t184)
      t284 = t277**2
      t287 = dsqrt(-4*l1*t256*z*l5*t259*t184*t182-4*t255*t266*t268*t270*
     #t273*t277*t182*t279+t284*t182)
      t288 = t287*l1
      t293 = 128.D0/3.D0*t14
      t294 = 128.D0/3.D0*t1
      t295 = 640.D0/3.D0*t7
      t298 = 320.D0/3.D0*t32-20864.D0/3.D0*t34+4160.D0/3.D0*t42+3328.D0/
     #3.D0*t44-2240.D0/3.D0*t46+256.D0/3.D0*t288-128.D0/3.D0*t191-1792*t
     #53+6080.D0/3.D0*t56+t293-t294+t295+320*t58-5376*t60
      t304 = 384*t8
      t305 = 128*t179
      t307 = 1408*t42
      t317 = 512+19840.D0/3.D0*t34-t307-4608*t44+2048.D0/3.D0*t191-32000
     #.D0/3.D0*t56-256.D0/3.D0*t3+256*t14+256*t21+512*t22+512.D0/3.D0*t7
     #-3968*t58
      t326 = 512.D0/3.D0*t179
      t331 = 13824*t60+27392.D0/3.D0*t62-17792*t64-3072*t66+8832*t70+291
     #84*t72-1024*t2+4864.D0/3.D0*t18-t326-1280*t48-2560.D0/3.D0*t82-640
     #*t85+5632*t87
      t338 = 512*t121
      t339 = 256*t123
      t340 = 256*l1
      t342 = 768*t226
      t345 = -5888*t89-2688*t91+3072*t99+1280.D0/3.D0*t101-24320*t119+t3
     #38+t339-t340+1024.D0/3.D0*z+t342-1408*t131+9088.D0/3.D0*t133
      t351 = 2560.D0/3.D0*t150
      t352 = 1280.D0/3.D0*t152
      t359 = -47104.D0/3.D0*t136+17664*t140+6016.D0/3.D0*t142+77312.D0/3
     #.D0*t144-88064.D0/3.D0*t146+t351-t352+256.D0/3.D0*t154+47872.D0/3.
     #D0*t158-3328*t160+7296*t164-9472*t170-48640.D0/3.D0*t174
      t370 = -23232*t62+18560*t64+9536*t66-1280*t70-30400*t72+t304-t305+
     #(t317+t331+t345+t359)*s2l1+1600.D0/3.D0*t82+10880.D0/3.D0*t85-9920
     #.D0/3.D0*t87+20480.D0/3.D0*t89+640*t91-7040*t99-320.D0/3.D0*t101
      t375 = 640.D0/3.D0*l1
      t378 = t287*t1
      t379 = 128.D0/3.D0*t378
      t384 = 768*t138
      t386 = 1664.D0/3.D0+11200*t119-2304*t121-4224*t123-t375-896.D0/3.D
     #0*t287+384*z+t379-128.D0/3.D0*t226+11584*t131-18688*t133+35648*t13
     #6+t384-58240.D0/3.D0*t140
      t398 = 384*t168
      t402 = -32704.D0/3.D0*t142-13120*t144+93184.D0/3.D0*t146-11200.D0/
     #3.D0*t148-3200.D0/3.D0*t150+3200.D0/3.D0*t152-1600.D0/3.D0*t154-10
     #3808.D0/3.D0*t158+18368*t160-1600*t164+3520*t166+t398+20736*t170+6
     #784.D0/3.D0*t172+26752.D0/3.D0*t174
      t408 = 11392.D0/3.D0*t44
      t411 = 384*t14
      t413 = 128*t7
      t416 = 7424*t34-2176*t42-t408+640*t191-13312.D0/3.D0*t56-t194+t411
     #+t294+640.D0/3.D0*t22+t413-1792*t58+9728*t60
      t426 = 896.D0/3.D0+9600*t62-47872.D0/3.D0*t64-3968*t66+3712*t70+65
     #408.D0/3.D0*t72-128*t2+384*t18-t326-1280.D0/3.D0*t48-t212-t213+281
     #6*t87
      t433 = 1280.D0/3.D0*t226
      t436 = -9472.D0/3.D0*t89-3712.D0/3.D0*t91+5248.D0/3.D0*t99+t220-14
     #720*t119+t222+t223-t375+640.D0/3.D0*z+t433-1664*t131+2944*t133
      t447 = -40832.D0/3.D0*t136+30848.D0/3.D0*t140+1408*t142+42112.D0/3
     #.D0*t144-57088.D0/3.D0*t146+t236-t237+t238+11776*t158-8192.D0/3.D0
     #*t160+9472.D0/3.D0*t164-18304.D0/3.D0*t170-8320*t174
      t450 = s2l1**2
      t452 = 64.D0/3.D0*t32
      t463 = 640-t452+20672.D0/3.D0*t34-2048*t42-t408+256.D0/3.D0*t46+25
     #6*t288+1664.D0/3.D0*t191-448*t53-7616.D0/3.D0*t56+512*t14+t294+896
     #.D0/3.D0*t7-1664*t58
      t470 = 128*t8
      t478 = 9280*t60+21760.D0/3.D0*t62-14656*t64-7936.D0/3.D0*t66+2880*
     #t70+19520*t72+t470-896.D0/3.D0*t179-320.D0/3.D0*t82-1856.D0/3.D0*t
     #85+3584.D0/3.D0*t87-6016.D0/3.D0*t89-832.D0/3.D0*t91+4864.D0/3.D0*
     #t99
      t490 = 256*t138
      t492 = 64.D0/3.D0*t101-12928*t119+320*t121-640*t123-384*l1-384*t28
     #7+1408.D0/3.D0*z-t379+t433+640.D0/3.D0*t131+1088*t133-33344.D0/3.D
     #0*t136+t490+7680*t140
      t493 = 1536*t142
      t503 = 128*t168
      t507 = t493+10560*t144-47744.D0/3.D0*t146-704.D0/3.D0*t148+640.D0/
     #3.D0*t150-t237+320.D0/3.D0*t154+30208.D0/3.D0*t158-5824.D0/3.D0*t1
     #60+3584.D0/3.D0*t164+1664.D0/3.D0*t166+t503-16256.D0/3.D0*t170-64*
     #t172-4928*t174
      t511 = t268*l5
      t512 = t266*t511
      t513 = t273*l1
      t515 = t270*t513*t279
      t516 = t512*t515
      t522 = t273*t1
      t524 = t270*t522*t279
      t525 = t512*t524
      t532 = t266*z*t511
      t533 = t532*t270
      t535 = t513*t279*l6
      t536 = t533*t535
      t538 = -20-24*t516+132*t5-6347*t10-1029*t12+3502*t15+8*t525+4866*t
     #19+8522*t23-2046*t25-924*t28+32*t536
      t539 = t287*t7
      t540 = t539*t1
      t542 = t84*l5
      t543 = t542*l1
      t546 = t270*t273*t279
      t547 = t512*t546
      t550 = l5**2
      t551 = t550*t1
      t553 = t550*l1
      t555 = t266*t268
      t556 = t555*t270
      t557 = t556*t535
      t559 = l5*t1
      t565 = 20*t540+4*t543+16*t547-34*t30-8*t551+10*t553-48*t557+16*t55
     #9+13573*t32-564*t34-30*t38+100*t40
      t570 = t532*t524
      t573 = t550*t21
      t577 = t266*t84
      t578 = t577*t511
      t579 = t578*t515
      t581 = t577*t268
      t584 = t273*t279*l6
      t585 = t581*t270*t584
      t589 = 96*t42+184*t44-19085*t46+8*t570-7778*t49+2*t573+1820*t51-50
     #*t288+16*t579+16*t585-28*t191-1954*t53
      t590 = t550*z
      t591 = t590*l1
      t593 = t274*t21
      t595 = t555*t515
      t597 = t3*t21
      t598 = t597*t37
      t601 = t578*t270
      t602 = t601*t535
      t604 = t276*t8
      t606 = t550*t84
      t607 = t606*t1
      t609 = t606*t21
      t611 = t559*t8
      t614 = l5*t21
      t616 = 24*t591-12*t593+24*t595+4*t598+424*t56-32*t602-8*t604-8*t60
     #7+12*t609+4*t611-16*t14-4*t614
      t621 = t512*t270
      t622 = t621*t535
      t633 = -40*t1+10*t21+48*t622-8*t7+68*t58-695*t60-3472*t62+2124*t64
     #+853*t66+540*t68-232*t70-5444*t72
      t636 = t532*t515
      t639 = t266*t8*t511
      t640 = t639*t546
      t642 = t550*t14
      t644 = t597*z
      t646 = t597*t7
      t650 = t9*t550
      t651 = t650*l1
      t653 = t17*t550
      t654 = t653*t1
      t656 = t274*t1
      t658 = -955*t74+312*t76-16*t636+8*t640-2*t642+4*t644-32*t646-2*t17
     #9-6795*t78+12*t651-2*t654+36*t656
      t660 = t553*t14
      t662 = t550*t8
      t663 = t662*t1
      t665 = t642*t1
      t667 = t550*t7
      t668 = t667*t1
      t670 = t590*t1
      t679 = 6*t660+12*t663-6*t665-4*t668-36*t670+3054*t80+2375*t82+2619
     #2*t85-3665*t87+12110*t89+482*t91+244*t93
      t683 = t621*t584
      t686 = t639*t515
      t695 = -837*t97-22660*t99-308*t101-32*t683+268*t103-8*t686-160*t10
     #5-190*t107+222*t109+1669*t111-2118*t113+1765*t115+66*t117
      t703 = t597*t14
      t705 = t597*t8
      t707 = t597*t84
      t710 = t578*t546
      t712 = t287*t21
      t714 = t597*t36
      t716 = t9*l1
      t718 = 1759*t119-119*t121-806*t123+8220*t125+112*t703-224*t705+280
     #*t707-4*t550+8*t710-10*t712-32*t714-48*t716
      t720 = t522*t279*l6
      t721 = t556*t720
      t724 = t9*l5
      t725 = t724*l1
      t727 = t724*t1
      t729 = t287*z
      t730 = t729*t1
      t735 = l1*l6
      t736 = t724*t735
      t738 = t650*t735
      t740 = l6**2
      t741 = t1*t740
      t742 = t274*t741
      t744 = t667*t21
      t746 = 16*t721+50*l1-12*t725+4*t727-60*t730+20*t287-868*t127-12*z+
     #32*t736-32*t738+96*t742+2*t744
      t748 = t556*t584
      t753 = t266*t9*t511
      t755 = t753*t270*t535
      t757 = t533*t720
      t761 = t639*t270
      t762 = t761*t584
      t766 = l1*t740
      t767 = t274*t766
      t769 = t590*t735
      t771 = 32*t748+40*t378+8*l5+16*t755-16*t757-14*t84-24*t275-16*t762
     #+1764*t129+12*t226-64*t767-64*t769
      t772 = t1*l6
      t773 = t590*t772
      t775 = t590*t766
      t777 = t590*t741
      t779 = t274*t735
      t781 = t274*t772
      t790 = 96*t773+64*t775-96*t777+64*t779-96*t781+3432*t131-9682*t133
     #+9339*t136+116*t138-12910*t140-26466*t142-3123*t144
      t795 = t581*t546
      t797 = t36*t21
      t799 = t542*t1
      t801 = t650*t21
      t803 = t729*l1
      t805 = t653*t21
      t807 = t539*t21
      t809 = t729*t21
      t811 = t761*t535
      t814 = 9882*t146-12598*t148-8*t795-48*t797-8*t799-8*t801+40*t803+2
     #*t805-10*t807+20*t809+16*t811-4*t662
      t815 = t8*l5
      t819 = t606*l1
      t821 = t650*t1
      t823 = t753*t515
      t828 = t621*t720
      t831 = t17*t1
      t833 = 4*t815-2*t606+4*t542-2*t819+4*t821-8*t823-7974*t150+15215*t
     #152-18030*t154-16*t828-20*t276-72*t831
      t835 = t601*t584
      t837 = t590*t21
      t839 = t662*t21
      t841 = t642*t21
      t844 = t555*t546
      t846 = t37*t2
      t848 = t555*t524
      t850 = t724*t766
      t852 = t650*t766
      t854 = t21*t740
      t855 = t274*t854
      t857 = t590*t854
      t859 = -16*t835+12*t837-8*t839+2*t841-2240*t156-16*t844-12*t846-8*
     #t848-32*t850+32*t852-32*t855+32*t857
      t860 = t21*l6
      t861 = t590*t860
      t863 = t274*t860
      t865 = t597*t9
      t867 = t597*t17
      t878 = -32*t861+32*t863-224*t865+112*t867-17558*t158+18996*t160-63
     #30*t162-534*t164+6006*t166+408*t168+23597*t170+18500*t172+3986*t17
     #4
      t886 = sbqe0(t184)
      t887 = dlog(z)
      t888 = t887**2
      t889 = t888*t887
      t892 = sbqe1(t184)
      t895 = sbqe2(t184)
      t897 = dlog(l1)
      t898 = t897**2
      t899 = t898*t897
      t902 = sqe1(t184)
      t903 = zqe1(t184,z,l3)
      t904 = t902-t903
      t907 = sqe2(t184)
      t908 = zqe2(t184,z,l3)
      t909 = t907-t908
      t919 = 128.D0/3.D0*t7
      t920 = 768*t58
      t927 = 64*t179
      t930 = -128.D0/3.D0-8512.D0/3.D0*t34+2368.D0/3.D0*t42+5504.D0/3.D0
     #*t44-1408.D0/3.D0*t191+960*t56-320.D0/3.D0*t14+t294+t919+t920-3968
     #*t60-9472.D0/3.D0*t62+17600.D0/3.D0*t64+1280*t66-1216*t70-7296*t72
     #-t927+64*t85-256*t87
      t936 = 128*t123
      t942 = 1024.D0/3.D0*t142
      t950 = 384*t89+64*t91-256*t99+15232.D0/3.D0*t119-256*t121-t936-t22
     #4-448.D0/3.D0*t226+512*t131-2624.D0/3.D0*t133+12928.D0/3.D0*t136-2
     #176*t140-t942-11264.D0/3.D0*t144+16192.D0/3.D0*t146-10112.D0/3.D0*
     #t158+768*t160-384*t164+1408*t170+4480.D0/3.D0*t174
      t953 = 12*t516
      t958 = 4*t525
      t963 = 16*t536
      t964 = 7+t953-29*t5+1634*t10+234*t12-823*t15-t958-1347*t19-2341*t2
     #3+561*t25+224*t28-t963
      t966 = 2*t543
      t967 = 8*t547
      t971 = 4*t551
      t972 = 5*t553
      t973 = 24*t557
      t974 = 8*t559
      t978 = -7*t540-t966-t967-2*t729*t2+8*t30+t971-t972+t973-t974-3533*
     #t32+94*t34+8*t38
      t984 = 4*t570
      t988 = 8*t579
      t989 = 8*t585
      t991 = -28*t40-12*t42-42*t44+5005*t46-t984+2153*t49-t573-448*t51+1
     #6*t288-t988-t989+11*t191
      t993 = 12*t591
      t994 = 6*t593
      t995 = 12*t595
      t997 = 16*t602
      t998 = 4*t604
      t999 = 4*t607
      t1000 = 6*t609
      t1001 = 2*t611
      t1003 = 2*t614
      t1004 = 520*t53-t993+t994-t995-t598-83*t56+t997+t998+t999-t1000-t1
     #001+2*t14+t1003
      t1008 = 24*t622
      t1018 = 10*t1-t1008+t7-16*t58+137*t60+708*t62-386*t64-164*t66-144*
     #t68+50*t70+1047*t72+268*t74
      t1020 = 8*t636
      t1021 = 4*t640
      t1026 = 6*t651
      t1027 = -84*t76+t1020-t1021-t2+t642-t644+8*t646+4*t179+t539*t2+164
     #8*t78-t1026+t654
      t1029 = 18*t656
      t1030 = 3*t660
      t1031 = 6*t663
      t1032 = 3*t665
      t1033 = 2*t668
      t1034 = 18*t670
      t1041 = -t1029-t1030-t1031+t1032+t1033+t1034-809*t80-513*t82-6508*
     #t85+753*t87-2657*t89-93*t91
      t1046 = 16*t683
      t1048 = 4*t686
      t1055 = -64*t93+228*t97+5302*t99+63*t101+t1046-64*t103+t1048+46*t1
     #05+55*t107-65*t109-468*t111+604*t113-504*t115
      t1067 = 2*t550
      t1068 = 4*t710
      t1071 = -18*t117-335*t119+22*t121+201*t123-2055*t125-28*t703+56*t7
     #05-70*t707+t1067-t1068+8*t714+12*t716
      t1072 = 8*t721
      t1074 = 6*t725
      t1075 = 2*t727
      t1080 = 16*t736
      t1081 = 16*t738
      t1082 = 48*t742
      t1083 = -t1072-16*l1+t1074-t1075+18*t730-7*t287+224*t127+3*z-t1080
     #+t1081-t1082-t744
      t1085 = 16*t748
      t1087 = 4*l5
      t1088 = 8*t755
      t1089 = 8*t757
      t1091 = 12*t275
      t1092 = 8*t762
      t1095 = 32*t767
      t1096 = 32*t769
      t1097 = -t1085-10*t378-t1087-t1088+t1089+4*t84+t1091+t1092-448*t12
     #9-7*t226+t1095+t1096
      t1098 = 48*t773
      t1099 = 32*t775
      t1100 = 48*t777
      t1101 = 32*t779
      t1102 = 48*t781
      t1111 = -t1098-t1099+t1100-t1101+t1102-832*t131+2343*t133-1995*t13
     #6-31*t138+2724*t140+6524*t142+598*t144-1996*t146
      t1115 = 4*t795
      t1117 = 4*t799
      t1118 = 4*t801
      t1122 = 8*t811
      t1123 = 2*t662
      t1124 = 2*t815
      t1125 = 3323*t148+t1115+12*t797+t1117+t1118-14*t803-t805+2*t807-2*
     #t809-t1122+t1123-t1124
      t1126 = 2*t542
      t1127 = 2*t821
      t1128 = 4*t823
      t1132 = 8*t828
      t1133 = 10*t276
      t1135 = 8*t835
      t1136 = t606-t1126+t819-t1127+t1128+1813*t150-3630*t152+4500*t154+
     #t1132+t1133+18*t831+t1135
      t1138 = 6*t837
      t1139 = 4*t839
      t1141 = 8*t844
      t1143 = 4*t848
      t1144 = 16*t850
      t1145 = 16*t852
      t1146 = 16*t855
      t1147 = 16*t857
      t1148 = 16*t861
      t1149 = -t1138+t1139-t841+560*t156+t1141+3*t846+t1143+t1144-t1145+
     #t1146-t1147+t1148
      t1150 = 16*t863
      t1163 = -t1150+56*t865-28*t867+3891*t158-4636*t160+1713*t162+101*t
     #164-1589*t166-108*t168-5389*t170-4870*t172-785*t174+t287*t2
      t1170 = 128.D0/3.D0*t32
      t1175 = 64.D0/3.D0*t288
      t1179 = 64*t1
      t1180 = 64.D0/3.D0*t7
      t1181 = 128*t58
      t1183 = -t1170+7040.D0/3.D0*t34-1600.D0/3.D0*t42+128*t44+896.D0/3.
     #D0*t46+t1175-704.D0/3.D0*t191+640*t53-1472.D0/3.D0*t56+t293+t1179+
     #t1180+t1181+3136.D0/3.D0*t60
      t1189 = 640.D0/3.D0*t82
      t1195 = 128.D0/3.D0*t101
      t1196 = 8768*t62-18688.D0/3.D0*t64-11072.D0/3.D0*t66+320.D0/3.D0*t
     #70+31232.D0/3.D0*t72-t470-t927-t1189-4288.D0/3.D0*t85+3712.D0/3.D0
     #*t87-7808.D0/3.D0*t89-704.D0/3.D0*t91+8192.D0/3.D0*t99+t1195
      t1202 = 256.D0/3.D0*t287
      t1207 = 1024.D0/3.D0*t191
      t1208 = 640*t56
      t1209 = 512*t58
      t1213 = 256*t66
      t1214 = 2432.D0/3.D0*t70
      t1219 = -128-1024*t34+128*t42+3200.D0/3.D0*t44-t1207+t1208+t294+t9
     #19+t1209-6784.D0/3.D0*t60-1152*t62+8192.D0/3.D0*t64+t1213-t1214-12
     #032.D0/3.D0*t72-128.D0/3.D0*t179+128.D0/3.D0*t85-512.D0/3.D0*t87
      t1221 = 128.D0/3.D0*t91
      t1224 = 128.D0/3.D0*z
      t1235 = 256*t164
      t1238 = 256*t89+t1221-512.D0/3.D0*t99+3072*t119-t1224-256.D0/3.D0*
     #t226+128*t131-384*t133+6400.D0/3.D0*t136-1408*t140-640.D0/3.D0*t14
     #2-7168.D0/3.D0*t144+9856.D0/3.D0*t146-5888.D0/3.D0*t158+1280.D0/3.
     #D0*t160-t1235+896*t170+2944.D0/3.D0*t174
      t1245 = -448.D0/3.D0-9088.D0/3.D0*t119+896*t121+1536*t123-64.D0/3.
     #D0*l1+t1202-64*z-t379+(t1219+t1238)*s2l1-13120.D0/3.D0*t131+21440.
     #D0/3.D0*t133-13312*t136-t490+7104*t140
      t1250 = 640.D0/3.D0*t154
      t1258 = 12736.D0/3.D0*t142+12352.D0/3.D0*t144-10944*t146+4352.D0/3
     #.D0*t148+t236-t352+t1250+38848.D0/3.D0*t158-21248.D0/3.D0*t160+512
     #*t164-3968.D0/3.D0*t166-t503-7872*t170-896*t172-9280.D0/3.D0*t174
      t1263 = 384*t42
      t1265 = 256*t191
      t1266 = 1024.D0/3.D0*t56
      t1267 = 64*t14
      t1275 = 64.D0/3.D0*t179
      t1276 = 64.D0/3.D0*t85
      t1277 = 256.D0/3.D0*t87
      t1278 = 128*t89
      t1279 = 64.D0/3.D0*t91
      t1280 = -256.D0/3.D0-3904.D0/3.D0*t34+t1263+896*t44-t1265+t1266-t1
     #267+384*t58-5248.D0/3.D0*t60-4352.D0/3.D0*t62+7616.D0/3.D0*t64+640
     #*t66-1472.D0/3.D0*t70-2944*t72+t1275+t1276-t1277+t1278+t1279
      t1281 = 256.D0/3.D0*t99
      t1283 = 128*t121
      t1284 = 64*t123
      t1285 = 128*t226
      t1286 = 256*t131
      t1287 = 1216.D0/3.D0*t133
      t1290 = 128*t142
      t1294 = 320*t160
      t1295 = 128*t164
      t1298 = -t1281+6016.D0/3.D0*t119-t1283-t1284+t224-t1224-t1285+t128
     #6-t1287+5504.D0/3.D0*t136-768*t140-t1290-4096.D0/3.D0*t144+6080.D0
     #/3.D0*t146-3968.D0/3.D0*t158+t1294-t1295+512*t170+512*t174
      t1302 = 256*t42
      t1305 = 832.D0/3.D0*t191
      t1306 = 1280.D0/3.D0*t56
      t1307 = 64.D0/3.D0*t1
      t1310 = 5824.D0/3.D0*t60
      t1315 = 3200*t72
      t1316 = 448*t58-t1310-4160.D0/3.D0*t62+2560*t64+576*t66-1856.D0/3.
     #D0*t70-t1315+t1275+t1276-t1277+t1278
      t1332 = -t1287+5632.D0/3.D0*t136-832*t140-t1290-1600*t144+2240*t14
     #6-4160.D0/3.D0*t158+t1294-448.D0/3.D0*t164+1600.D0/3.D0*t170+576*t
     #174
      t1350 = 320.D0/3.D0*t46
      t1352 = 64.D0/3.D0*t3
      t1354 = 320*t21
      t1355 = 448.D0/3.D0*t22
      t1356 = -128.D0/3.D0-t452-3776*t34+3904.D0/3.D0*t42+3712.D0/3.D0*t
     #44+t1350-t1305-2368*t56-t1352-704.D0/3.D0*t14-t1307+t1354+t1355
      t1359 = 17984.D0/3.D0*t64
      t1362 = 13184.D0/3.D0*t72
      t1363 = 1984.D0/3.D0*t2
      t1364 = 2176.D0/3.D0*t18
      t1366 = 448*t48
      t1367 = 320*t82
      t1369 = -t1180-2816.D0/3.D0*t58-t1310-12416.D0/3.D0*t62+t1359+2176
     #*t66+5312.D0/3.D0*t70-t1362-t1363+t1364-192*t179-t1366-t1367-576*t
     #85
      t1372 = 1920*t89
      t1375 = 128*t101
      t1378 = 768*t131
      t1381 = 9088.D0/3.D0*t140
      t1382 = 1600*t87-t1372-2240.D0/3.D0*t91+4160.D0/3.D0*t99+t1375+512
     #*t119-t338-t339-t224+1472.D0/3.D0*t226+t1378-1984.D0/3.D0*t133+104
     #96.D0/3.D0*t136+t1381
      t1385 = 1216.D0/3.D0*t146
      t1386 = 192*t148
      t1387 = 320*t152
      t1388 = 128*t154
      t1392 = 64*t166
      t1394 = 640.D0/3.D0*t172
      t1396 = 2816.D0/3.D0*t142+6976.D0/3.D0*t144-t1385+t1386+t236-t1387
     #+t1388-128*t158-1024.D0/3.D0*t160+5504.D0/3.D0*t164-t1392-6080.D0/
     #3.D0*t170-t1394-9344.D0/3.D0*t174
      t1403 = 256*t46
      t1409 = -t1170-10624.D0/3.D0*t34+3968.D0/3.D0*t42+9728.D0/3.D0*t44
     #+t1403+t1175-2432.D0/3.D0*t191+896*t53+3328.D0/3.D0*t56-t195+t294-
     #t919+1408*t58-6144*t60
      t1410 = 1920*t62
      t1422 = t1410+17152.D0/3.D0*t64-3712.D0/3.D0*t66-1920*t70-13952.D0
     #/3.D0*t72-256*t8-t1189-1152*t85+768*t87-1792*t89-128*t91+2048*t99+
     #t1195+18944.D0/3.D0*t119
      t1433 = -640.D0/3.D0+640*t121+1920*t123+t224+t1202-t225-t379-1024.
     #D0/3.D0*t226-11264.D0/3.D0*t131+14464.D0/3.D0*t133-3200*t136-512*t
     #138+7552.D0/3.D0*t140+8704.D0/3.D0*t142
      t1441 = 128*t53
      t1444 = 128*t21
      t1445 = 256*t22
      t1446 = -256-t1170-25600.D0/3.D0*t34+3072*t42+4480*t44+t1403-3328.
     #D0/3.D0*t191+t1441-2944*t56-t194-640*t14+t1444+t1445-t199
      t1453 = 512*t2
      t1454 = 2432.D0/3.D0*t18
      t1456 = 640*t48
      t1457 = 640*t82
      t1459 = t920-7040*t60-7040*t62+36992.D0/3.D0*t64+4480*t66+640*t70-
     #20480.D0/3.D0*t72-t1453+t1454+256.D0/3.D0*t179-t1456-t1457-1408*t8
     #5+t214
      t1464 = 256*t101
      t1468 = 512.D0/3.D0*z
      t1473 = -4352*t89-1408*t91+3328*t99+t1464+9856.D0/3.D0*t119-1152*t
     #121-t223+128*l1-t1468-384*t226+t1378+640*t133+9856.D0/3.D0*t136+81
     #92*t140
      t1477 = 2432.D0/3.D0*t148
      t1478 = 640*t152
      t1479 = 256*t154
      t1483 = 512*t166
      t1485 = 640*t172
      t1487 = 3072*t142+12800.D0/3.D0*t144-12800.D0/3.D0*t146+t1477+t351
     #-t1478+t1479+12800.D0/3.D0*t158-2944*t160+3072*t164-t1483-6528*t17
     #0-t1485-6528*t174
      t1498 = -2816*t144+1024.D0/3.D0*t146+3712.D0/3.D0*t148+t236-t352+t
     #1250+(t1446+t1459+t1473+t1487)*s2l1+13696.D0/3.D0*t158-13312.D0/3.
     #D0*t160-t1295-4096.D0/3.D0*t166-256*t168-4224*t170-2176.D0/3.D0*t1
     #72-1024.D0/3.D0*t174
      t1504 = 448.D0/3.D0*t46
      t1506 = 320*t53
      t1508 = 320.D0/3.D0*t22
      t1509 = -t452-11456.D0/3.D0*t34+t307+6400.D0/3.D0*t44+t1504-384*t1
     #91+t1506-t1266-t1352-320*t14-t294+t1508-t413+t920
      t1514 = 64*t8
      t1515 = 64*t2
      t1516 = 192*t18
      t1518 = 640.D0/3.D0*t48
      t1521 = -12416.D0/3.D0*t60-5824.D0/3.D0*t62+t1359+3584.D0/3.D0*t66
     #-1984.D0/3.D0*t70-t1362-t1514-t1515+t1516+448.D0/3.D0*t179-t1518-t
     #1367-2240.D0/3.D0*t85+4160.D0/3.D0*t87-t1372
      t1529 = 896.D0/3.D0*t226
      t1533 = 128*t138
      t1534 = -640.D0/3.D0-576*t91+1600*t99+t1375+10496.D0/3.D0*t119-192
     #*t121+448*t123+512.D0/3.D0*l1-t1468-t1529-2816.D0/3.D0*t131+5312.D
     #0/3.D0*t133+512*t136-t1533+t1381
      t1537 = 2176.D0/3.D0*t148
      t1541 = 1984.D0/3.D0*t166
      t1542 = 64*t168
      t1544 = 448*t172
      t1546 = 5504.D0/3.D0*t142-128*t144-t1385+t1537+t236-t1387+t1388+69
     #76.D0/3.D0*t158-2368*t160+2816.D0/3.D0*t164-t1541-t1542-9344.D0/3.
     #D0*t170-t1544-6080.D0/3.D0*t174
      t1551 = 768*t42
      t1552 = 5632.D0/3.D0*t44
      t1558 = -384+t1170-1280*t34+t1551+t1552-128*t46-192*t288-t1207+140
     #8*t53+3712.D0/3.D0*t56-t411-t294-t295+896*t58
      t1570 = -11264.D0/3.D0*t60+17152.D0/3.D0*t62+1920*t64-3584*t66-409
     #6.D0/3.D0*t70-t1315-t304+640.D0/3.D0*t179+t1189+768*t85-1152*t87+2
     #048*t89+256*t91-1792*t99
      t1581 = -t1195+14464.D0/3.D0*t119+1408*t121+3200*t123+t340+256*t28
     #7-896.D0/3.D0*z+t379-t1529-6144*t131+18944.D0/3.D0*t133-13952.D0/3
     #.D0*t136-t384-4224*t140
      t1592 = -t942-13312.D0/3.D0*t144+13696.D0/3.D0*t146+3328.D0/3.D0*t
     #148-t236+t352-t1250+1024.D0/3.D0*t158-2816*t160-2176.D0/3.D0*t164-
     #1920*t166-t398+7552.D0/3.D0*t170-128*t172+8704.D0/3.D0*t174
      t1604 = 10+t953-169*t5+6964*t10+1280*t12-4229*t15-t958-4699*t19-80
     #92*t23+2100*t25+1134*t28-t963
      t1611 = -10*t540-t966-t967+43*t30+t971-t972+t973-t974-14129*t32+74
     #4*t34+33*t38-98*t40
      t1621 = -132*t42-233*t44+19155*t46-t984+7313*t49-t573-2198*t51+25*
     #t288-t988-t989+32*t191+2114*t53
      t1625 = -t993+t994-t995-5*t598-553*t56+t997+t998+t999-t1000-t1001+
     #17*t14+t1003
      t1639 = 20*t1-5*t21-t1008+4*t7-94*t58+901*t60+4155*t62-2713*t64-10
     #56*t66-606*t68+299*t70+6759*t72
      t1646 = 929*t74-339*t76+t1020-t1021-3*t8+t642-5*t644+40*t646+t179+
     #7958*t78-t1026+t654
      t1654 = -t1029-t1030-t1031+t1032+t1033+t1034-3235*t80-2955*t82-278
     #30*t85+4577*t87-14376*t89-630*t91
      t1666 = -274*t93+854*t97+25477*t99+399*t101+t1046-334*t103+t1048+1
     #46*t105+170*t107-192*t109-1573*t111+1911*t113-1608*t115
      t1680 = -69*t117-2269*t119+151*t121+931*t123-9325*t125-140*t703+28
     #0*t705-350*t707+t1067-t1068+5*t712+40*t714
      t1687 = 60*t716-t1072-25*l1+t1074-t1075+30*t730-10*t287+994*t127+6
     #*z-t1080+t1081-t1082
      t1693 = -t744-t1085-20*t378-t1087-t1088+t1089+16*t84+t1091+t1092-2
     #058*t129+3*t226+t1095
      t1701 = t1096-t1098-t1099+t1100-t1101+t1102-3817*t131+10451*t133-1
     #0836*t136-142*t138+15405*t140+27828*t142+3980*t144
      t1710 = -11979*t146+12572*t148+t1115+60*t797+t1117+t1118-20*t803-t
     #805+5*t807-10*t809-t1122+t1123
      t1715 = -t1124+t606-t1126+t819-t1127+t1128+9513*t150-17375*t152+19
     #675*t154+t1132+t1133+90*t831
      t1719 = t1135-t1138+t1139-t841+2660*t156+t1141+15*t846+t1143+t1144
     #-t1145+t1146-t1147
      t1731 = t1148-t1150+280*t865-140*t867+19915*t158-20096*t160+6275*t
     #162+700*t164-6180*t166-474*t168-26466*t170-18299*t172-5019*t174
      t1749 = 128.D0/3.D0*t46
      t1754 = 128*t14
      t1755 = -t1170+640.D0/3.D0*t34+t1263-256*t44+t1403-128*t191+t1441-
     #15232.D0/3.D0*t56-t194-t1754+t1444+t1445
      t1764 = -1152*t58+9856.D0/3.D0*t60+9856.D0/3.D0*t62-13568.D0/3.D0*
     #t64-t1213+10496.D0/3.D0*t70+38144.D0/3.D0*t72-t1453+t1454-t1456-t1
     #457-1536*t85+3840*t87
      t1775 = -5120*t89-1536*t91+3840*t99+t1464-9344*t119-t1283+t936+t12
     #85-1152*t131+10496.D0/3.D0*t133-9344*t136+13056*t140+11648.D0/3.D0
     #*t142
      t1783 = 39040.D0/3.D0*t144-52480.D0/3.D0*t146+t1477+t351-t1478+t14
     #79+39040.D0/3.D0*t158-15232.D0/3.D0*t160+11648.D0/3.D0*t164-t1483-
     #29312.D0/3.D0*t170-t1485-29312.D0/3.D0*t174
      t1792 = 640.D0/3.D0*t85
      t1793 = 640.D0/3.D0*t87
      t1794 = 1280.D0/3.D0*t89
      t1795 = -8192.D0/3.D0*t34+t1551+1280*t44-t1749-t1265-t1441+(t1755+
     #t1764+t1775+t1783)*s2l1+t1208-t1754+t1209-9088.D0/3.D0*t60-4736*t6
     #2+17536.D0/3.D0*t64+1920*t66-t1214-22784.D0/3.D0*t72+t1792-t1793+t
     #1794+t1221
      t1796 = 1280.D0/3.D0*t99
      t1810 = -t1796+12160.D0/3.D0*t119-384*t121-t223-t1285+t228-t229+66
     #56*t136-7168.D0/3.D0*t140-3328.D0/3.D0*t142-3072*t144+17408.D0/3.D
     #0*t146-1280.D0/3.D0*t148-5376*t158+2432*t160-t1235+384*t166+2304*t
     #170+t1394+3712.D0/3.D0*t174
      t1813 = 1536*t34
      t1814 = 384*t44
      t1819 = -t452+t1813-t1302-t1814+t1504+t1506-5248.D0/3.D0*t56-t1352
     #+t1508-256*t58+1920*t60+4928*t62
      t1820 = 14080.D0/3.D0*t64
      t1823 = 25856.D0/3.D0*t72
      t1826 = 2560*t89
      t1828 = -t1820-5632.D0/3.D0*t66+1024*t70+t1823-t1514-t1515+t1516-t
     #1518-t1367-2560.D0/3.D0*t85+5440.D0/3.D0*t87-t1826-2048.D0/3.D0*t9
     #1
      t1837 = 20224.D0/3.D0*t140
      t1840 = 6080.D0/3.D0*t99+t1375-4480*t119+448*t121+768*t123-6656.D0
     #/3.D0*t131+11392.D0/3.D0*t133-24704.D0/3.D0*t136-t1533+t1837+7424.
     #D0/3.D0*t142+5888*t144
      t1841 = 29440.D0/3.D0*t146
      t1847 = -t1841+t1537+t236-t1387+t1388+26176.D0/3.D0*t158-3968*t160
     #+1536*t164-t1541-t1542-5632*t170-t1544-13376.D0/3.D0*t174
      t1859 = 1920*t34-t1263-t1814+t1749+512*t53-t1306-t1181+1536*t60+17
     #536.D0/3.D0*t62-4736*t64-8192.D0/3.D0*t66+384*t70+6656*t72-t470-t1
     #792+t1793-t1794-t1221+t1796
      t1877 = -2688*t119+t222+1280*t123-9088.D0/3.D0*t131+12160.D0/3.D0*
     #t133-22784.D0/3.D0*t136-t490+2304*t140+3712.D0/3.D0*t142+2432*t144
     #-5376*t146+640*t148+17408.D0/3.D0*t158-3072*t160+640.D0/3.D0*t164-
     #2432.D0/3.D0*t166-t503-7168.D0/3.D0*t170-256*t172-3328.D0/3.D0*t17
     #4
      t1883 = -t452+t1813-t1302-t1552+t1350+448*t191-3968*t56-t1352-t117
     #9+t1354+t1355-64*t7
      t1889 = -6656.D0/3.D0*t58+4928*t60+t1410-t1820-384*t66+11392.D0/3.
     #D0*t70+t1823-t1363+t1364-t305-t1366-t1367-2048.D0/3.D0*t85
      t1897 = 6080.D0/3.D0*t87-t1826-2560.D0/3.D0*t91+5440.D0/3.D0*t99+t
     #1375-24704.D0/3.D0*t119+t342-t1286+1024*t133-4480*t136+t1837+t493
      t1904 = 26176.D0/3.D0*t144-t1841+t1386+t236-t1387+t1388+5888*t158-
     #5248.D0/3.D0*t160+7424.D0/3.D0*t164-t1392-13376.D0/3.D0*t170-t1394
     #-5632*t174
      t1919 = 50*t5-1990*t10-378*t12+1244*t15+1439*t19+1770*t23-798*t25-
     #392*t28-15*t30+3443*t32-398*t34
      t1932 = -11*t38+26*t40+75*t42+93*t44-4373*t46-1784*t49+756*t51-12*
     #t191-1034*t53+2*t598+120*t56-6*t14
      t1946 = 33*t58-331*t60-2159*t62+1214*t64+722*t66+258*t68-80*t70-24
     #07*t72-275*t74+123*t76+15*t8+2*t644
      t1959 = -16*t646-2326*t78+908*t80+737*t82+6478*t85-1007*t87+3282*t
     #89+131*t91+92*t93-234*t97-5927*t99-98*t101
      t1974 = 116*t103-32*t105-35*t107+35*t109+351*t111-389*t113+372*t11
     #5+21*t117+699*t119-147*t121-620*t123+2700*t125
      t1987 = 56*t703-112*t705+140*t707-16*t714-24*t716-336*t127-6*t84+7
     #00*t129-6*t226+2043*t131-4313*t133+4284*t136
      t2001 = 105*t138-3646*t140-7057*t142-1000*t144+3369*t146-3964*t148
     #-24*t797-2391*t150+4365*t152-4895*t154-36*t831-910*t156
      t2014 = -6*t846-112*t865+56*t867-6045*t158+6335*t160-1477*t162-137
     #*t164+2498*t166+246*t168+6574*t170+4611*t172+1099*t174
      t2027 = (l1+z-t179)**2
      s2 = 64.D0/3.D0*(t55+t95+t135+t176)*t182*t185+((t201+t215+t230+t24
     #4)*t247+(t298+t370+t386+t402)*s1l1+(t416+t426+t436+t447)*t450+(t46
     #3+t478+t492+t507)*s2l1-64.D0/3.D0*(t538+t565+t589+t616+t633+t658+t
     #679+t695+t718+t746+t771+t790+t814+t833+t859+t878)*t182)*(t886*t889
     #/6-t892*t888/2+t895*t887+t185*t899/6+t904*t898/2-t909*t897)
      s3 = s2+((t930+t950)*t247-64.D0/3.D0*(t964+t978+t991+t1004+t1018+t
     #1027+t1041+t1055+t1071+t1083+t1097+t1111+t1125+t1136+t1149+t1163)*
     #t182+(t1183+t1196+t1245+t1258)*s1l1+(t1280+t1298)*t450+(-192-3520.
     #D0/3.D0*t34+t1302+2944.D0/3.D0*t44-64*t288-t1305+t1306-t1267+t1307
     #-t1180+t1316+t1279-t1281+7040.D0/3.D0*t119-t1283-t1284+64*l1+128*t
     #287-320.D0/3.D0*z-320.D0/3.D0*t226+t1286+t1332)*s2l1)*(t892*t889/6
     #-t895*t888/2-t904*t899/6+t909*t898/2)
      s1 = s3+((t1356+t1369+t1382+t1396)*t247+(t1409+t1422+t1433+t1498)*
     #s1l1+(t1509+t1521+t1534+t1546)*t450+(t1558+t1570+t1581+t1592)*s2l1
     #-64.D0/3.D0*(t1604+t1611+t1621+t1625+t1639+t1646+t1654+t1666+t1680
     #+t1687+t1693+t1701+t1710+t1715+t1719+t1731)*t182)*(-t886*t888/2+t8
     #92*t887-t185*t898/2-t904*t897+t907-t908)+((t1795+t1810)*s1l1+(t181
     #9+t1828+t1840+t1847)*t450+(t1859+t1877)*s2l1+(t1883+t1889+t1897+t1
     #904)*t247-64.D0/3.D0*(t1919+t1932+t1946+t1959+t1974+t1987+t2001+t2
     #014)*t182)*(t185*t897+t902-t903+t886*t887)
      s2 = 1/t2027*t182
      taux11 = s1*s2
      t2 = l1**2
      t3 = 3-3*l1+t2
      t4 = t2**2
      t5 = z*t4
      t7 = z**2
      t8 = t7*z
      t9 = t8*t4
      t11 = t7*t4
      t13 = t7**2
      t14 = t13*l1
      t15 = 4*t14
      t16 = t13*z
      t17 = t16*t2
      t19 = t2*l1
      t20 = t16*t19
      t22 = t16*t4
      t24 = t13*t4
      t26 = t4*t2
      t27 = t8*t26
      t29 = t26*t16
      t32 = t4*l1
      t33 = t7*t32
      t35 = t8*t32
      t37 = t13*t7
      t38 = t37*t32
      t39 = 16*t38
      t40 = t13*t32
      t42 = t16*t32
      t44 = t26*z
      t46 = t37*t4
      t47 = 13*t46
      t48 = t26*t13
      t50 = t37*t26
      t52 = t26*t7
      t54 = z*t32
      t56 = 246*t33-334*t35+t39+251*t40-99*t42+54*t44-t47-135*t48-9*t50-
     #135*t52-95*t54
      t58 = t37*t19
      t60 = t4*t19
      t61 = t8*t60
      t62 = 40*t61
      t63 = z*t60
      t64 = 12*t63
      t65 = t7*t60
      t66 = 30*t65
      t69 = t13*t60
      t70 = 30*t69
      t71 = 2*t60
      t72 = 3*t19
      t73 = t13*t19
      t75 = t7*l1
      t77 = 4*t58-t62-t64+t66-9*t26+15*t32+t70+t71+t72+127*t73+3*t75
      t78 = t8*t2
      t80 = t8*l1
      t82 = z*t2
      t84 = t7*t2
      t86 = t37*t60
      t87 = 2*t86
      t88 = z*t19
      t90 = t7*t19
      t92 = t8*t19
      t93 = 176*t92
      t94 = t16*t60
      t95 = 12*t94
      t96 = t13*t2
      t99 = 52*t78-7*t80+6*t82-31*t84+t87-35*t88+118*t90-t93-t95-35*t96-
     #11*t4
      t101 = 82*t5+324*t9-231*t11+t15+8*t17-41*t20+90*t22-241*t24+180*t2
     #7+54*t29+t56+t77+t99
      t105 = -12+12*l1-4*t2
      t114 = 140*t27
      t118 = 3*t5+10*t9-14*t11-31*t14-33*t17+50*t20-13*t22+9*t24+t114+38
     #*t29+131*t33-164*t35
      t129 = 5*t38+111*t40-38*t42+46*t44+4*t46-100*t48-6*t50-110*t52-54*
     #t54-9*t58-t62-t64+t66
      t131 = 8*t26
      t132 = t37*t2
      t133 = 4*t132
      t134 = z*l1
      t143 = -t131+t133+6*t134+9*t32+3*t2+3*t7-7*t8+t70+t71-7*t19-111*t7
     #3-27*t75
      t151 = 4*t13
      t152 = t16*l1
      t153 = 8*t152
      t155 = -118*t78+44*t80-27*t82+80*t84+t87+38*t88-93*t90+132*t92+t15
     #1+t153-t95+91*t96+t4
      t157 = t118+t129+t143+t155
      t159 = -t105
      t172 = -2+154*t5+742*t9-472*t11+56*t14+84*t17-206*t20+280*t22-634*
     #t24+240*t27+80*t29+390*t33-580*t35
      t183 = 38*t38+480*t40-210*t42+64*t44-50*t46-190*t48-14*t50-170*t52
     #-138*t54+32*t58-t62-t64+t66
      t191 = 4*t8
      t196 = -10*t26-8*t132-4*t134+20*t32-2*t2-4*t7+t191+t70+t71+10*t19+
     #484*t73+26*t75+242*t78
      t204 = 2*z
      t205 = 16*t152
      t208 = -58*t80+30*t82-128*t84+t87-94*t88+326*t90-552*t92-8*t13-t20
     #4-t205-t95-218*t96-20*t4
      t210 = t172+t183+t196+t208
      t217 = 2*t17
      t232 = 30*t5+140*t9-95*t11+t14+t217-13*t20+38*t22-105*t24+t114+42*
     #t29+141*t33-194*t35+9*t38+146*t40-57*t42+42*t44-5*t46-105*t48-7*t5
     #0-105*t52
      t241 = 56*t92
      t244 = -53*t54+t58-t62-t64+t66-7*t26+8*t32+t70+t71+42*t73+12*t78-t
     #80-5*t84+t87-7*t88+33*t90-t241-t95-9*t96-3*t4
      t245 = t232+t244
      t258 = -19*t5-74*t9+47*t11-9*t14-13*t17+34*t20-35*t22+70*t24+100*t
     #27+26*t29+41*t33
      t261 = 2*t42
      t265 = 4*t50
      t269 = -44*t35-t38+21*t40-t261+34*t44+7*t46-70*t48-t265-80*t52-18*
     #t54-5*t58-t62
      t276 = -t64+t66-6*t26+t132+3*t32-t8+t70+t71-t72-83*t73-5*t75-50*t7
     #8
      t277 = 12*t80
      t283 = 2*t152
      t286 = t277-7*t82+29*t84+t87+22*t88-65*t90+100*t92+t13+t283-t95+40
     #*t96+4*t4
      t288 = t258+t269+t276+t286
      t301 = -2+86*t5+502*t9-302*t11+34*t14+48*t17-128*t20+188*t22-434*t
     #24+200*t27+68*t29+276*t33
      t308 = 160*t48
      t313 = -424*t35+28*t38+356*t40-156*t42+52*t44-32*t46-t308-12*t50-1
     #40*t52-92*t54+18*t58-t62-t64
      t316 = 2*t7
      t321 = t66-t131-t133+12*t32-t316+t191+t70+t71+2*t19+316*t73+14*t75
     #+158*t78
      t330 = -40*t80+10*t82-74*t84+t87-44*t88+196*t90-360*t92-t151-t204-
     #t153-t95-138*t96-8*t4
      t332 = t301+t313+t321+t330
      t344 = 120*t27
      t352 = 8*t44
      t357 = -71*t33+186*t35-75*t38-254*t40+191*t42-t352+58*t46+t308+48*
     #t50+48*t52+11*t54
      t359 = t13*t8
      t360 = t359*t26
      t363 = 2*t359*t60
      t364 = t359*t4
      t366 = t359*t32
      t369 = t359*t19
      t371 = 30*t61
      t372 = 2*t63
      t373 = 12*t65
      t374 = 40*t69
      t375 = -8*t360+t363-8*t364+12*t366-23*t58+2*t369+t371+t372-t373+t1
     #33-t374
      t379 = 12*t86
      t382 = 30*t94
      t384 = -100*t73-14*t78+t80+2*t84-t379+t88-17*t90+64*t92+t283+t382+
     #26*t96
      t386 = -6*t5-148*t9+50*t11-3*t14-18*t17+73*t20-158*t22+212*t24-t34
     #4-120*t29+t357+t375+t384
      t388 = t37*l1
      t389 = 4*t388
      t397 = 110*t27
      t400 = t389+t5-14*t9+3*t11+22*t14+53*t17-74*t20+9*t22+10*t24-t397-
     #100*t29-48*t33
      t411 = 116*t35-32*t38-144*t40+96*t42-t352-t47+140*t48+38*t50+46*t5
     #2+8*t54-6*t360+t363+4*t364
      t415 = 2*t16
      t420 = 4*t366+34*t58+t415-6*t369+t371+t372-t373-19*t132+t8-t374+84
     #*t73+2*t75
      t427 = t359*t2
      t428 = 2*t427
      t429 = 42*t78-t277+t82-13*t84-t379-4*t88+22*t90-t241-3*t13-t205+t3
     #82-66*t96+t428
      t431 = t400+t411+t420+t429
      t434 = 4*t5
      t443 = -2-8*t388+t434-146*t9+16*t11-14*t14-100*t17+236*t20-346*t22
     #+334*t24-t344-150*t29
      t449 = 6*t44
      t454 = 4*t54
      t456 = -50*t33+188*t35-146*t38-332*t40+308*t42-t449+170*t46+180*t4
     #8+66*t50+42*t52+t454-12*t360+t363
      t465 = 2*l6-1
      t466 = t465**2
      t469 = l5-1
      t471 = -1+l1
      t473 = (1-l1+t134)**2
      t474 = 1/t473
      t479 = dsqrt(z)
      t481 = dsqrt(l5)
      t483 = dsqrt(-t469)
      t486 = dsqrt(l1)
      t490 = l5*z*l1-l1+1-l5+l5*l1
      t492 = -t471
      t493 = dsqrt(t492)
      t498 = t490**2
      t501 = dsqrt(4*l1*t466*z*l5*t469*t471*t474-4*t465*t479*t481*t483*t
     #486*t490*t474*t493+t498*t474)
      t502 = 2*t501
      t504 = -32*t364+28*t366-114*t58-4*t16+18*t369+t371+t372-t373+44*t1
     #32-t316+t502-t374-200*t73
      t505 = 4*t75
      t506 = 6*t80
      t507 = 2*t82
      t508 = 14*t84
      t509 = 6*t88
      t515 = t505-t506+t507-t508-t379-t509+14*t90+52*t92-t204+26*t152+t3
     #82+72*t96-4*t427
      t517 = t443+t456+t504+t515
      t521 = zqe0(t492,z,l3)
      t522 = dlog(l1)
      t524 = sqe1(t492)
      t525 = zqe1(t492,z,l3)
      t526 = sbqe0(t492)
      t527 = dlog(z)
      t550 = 2-18*t5-222*t9+102*t11-24*t14-40*t17+94*t20-120*t22+234*t24
     #-40*t27-20*t29-72*t33+148*t35-t39-152*t40+78*t42-4*t44+24*t46+40*t
     #48
      t557 = 26*t84
      t558 = 10*t88
      t561 = t265+20*t52+14*t54-16*t58+t133+t316-192*t73-t505-78*t78+18*
     #t80-t507+t557+t558-72*t90+t93+t151+t204+t153+90*t96
      t562 = t550+t561
      t577 = 2-8*t5-42*t9+30*t11-t15-t217+6*t20-6*t22+26*t24-8*t33+12*t3
     #5-8*t40+t261
      t587 = 2*t54+t316-2*t8-32*t73-8*t75-38*t78+14*t80-4*t82+t557+t558-
     #40*t90+t241+t204+18*t96
      t588 = t577+t587
      t603 = 2+t389-t434+106*t9-10*t11+8*t14+58*t17-148*t20+238*t22-238*
     #t24+t397+130*t29+44*t33
      t613 = -156*t35+108*t38+264*t40-236*t42+t449-112*t46-t308-56*t50-4
     #0*t52-t454+10*t360-t363+20*t364
      t620 = -20*t366+68*t58+t415-10*t369-t371-t372+t373-24*t132+t316-t5
     #02+t374+128*t73-t505
      t626 = -6*t78+t506-t507+t508+t379+t509-16*t90-28*t92+t204-14*t152-
     #t382-42*t96+t428
      t628 = t603+t613+t620+t626
      t632 = t527**2
      t635 = sbqe1(t492)
      t637 = t522**2
      t640 = t524-t525
      t642 = sqe2(t492)
      t643 = zqe2(t492,z,l3)
      t658 = t632*t527
      t663 = sbqe2(t492)
      t665 = t637*t522
      t670 = t642-t643
      t706 = (l1+z-t134)**2
      s2 = (128.D0/3.D0*(t3*t101+t105*t157+t159*t210)*s1l1+128.D0/3.D0*(
     #t3*t245+t105*t288+t159*t332)*s2l1+128.D0/3.D0*t3*t386+128.D0/3.D0*
     #t105*t431+128.D0/3.D0*t159*t517)*(t521*t522+t524-t525+t526*t527)+(
     #128.D0/3.D0*(t3*t157+t105*t210+t159*t562)*s1l1+128.D0/3.D0*(t3*t28
     #8+t105*t332+t159*t588)*s2l1+128.D0/3.D0*t3*t431+128.D0/3.D0*t105*t
     #517+128.D0/3.D0*t159*t628)*(-t526*t632/2+t635*t527-t521*t637/2-t64
     #0*t522+t642-t643)
      s1 = s2+(128.D0/3.D0*(t3*t210+t105*t562)*s1l1+128.D0/3.D0*(t3*t332
     #+t105*t588)*s2l1+128.D0/3.D0*t3*t517+128.D0/3.D0*t105*t628)*(t526*
     #t658/6-t635*t632/2+t663*t527+t521*t665/6+t640*t637/2-t670*t522)+(1
     #28.D0/3.D0*t3*t562*s1l1+128.D0/3.D0*t3*t588*s2l1+128.D0/3.D0*t3*t6
     #28)*(t635*t658/6-t663*t632/2-t640*t665/6+t670*t637/2)-(128.D0/3.D0
     #*(t105*t101+t159*t157)*s1l1+128.D0/3.D0*(t105*t245+t159*t288)*s2l1
     #+128.D0/3.D0*t105*t386+128.D0/3.D0*t159*t431)*t521
      s2 = 1/t706*t474
      taux12 = s1*s2

CCC Coefficients

      tcoff19 = taux11
      tcoff20 = taux11
      tcoff21 = taux12
      tcoff22 = -taux12

CCC Jet functions

      tFjet19 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff19)
      tFjet20 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff20)
      tFjet21 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff21)
      tFjet22 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff22)

CCC Sum

      RVnum6e4 = tFjet19*tcoff19+tFjet20*tcoff20+tFjet21*tcoff21+tFjet22
     #*tcoff22
      RETURN
      END

      double precision function RVnum7e0(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      taux13 = 0
      taux14 = 0

CCC Coefficients

      tcoff23 = taux13
      tcoff24 = taux13
      tcoff25 = taux14
      tcoff26 = -taux14

CCC Jet functions

      tFjet23 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff23)
      tFjet24 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff24)
      tFjet25 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff25)
      tFjet26 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff26)

CCC Sum

      RVnum7e0 = tFjet23*tcoff23+tFjet24*tcoff24+tFjet25*tcoff25+tFjet26
     #*tcoff26
      RETURN
      END

      double precision function RVnum7e1(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t1 = l1**2
      t2 = z*t1
      t3 = 192*t2
      t6 = z*l1
      t9 = s1l1**2
      t11 = z**2
      t12 = t11*l1
      t14 = t11*t1
      t16 = t1*l1
      t17 = z*t16
      t19 = t11*t16
      t25 = 2*l6-1
      t26 = t25**2
      t29 = l5-1
      t31 = -1+l1
      t33 = (1-l1+t6)**2
      t34 = 1/t33
      t39 = dsqrt(z)
      t41 = dsqrt(l5)
      t43 = dsqrt(-t29)
      t46 = dsqrt(l1)
      t47 = z*l5
      t48 = t47*l1
      t49 = l5*l1
      t50 = t48-l1+1-l5+t49
      t53 = dsqrt(-t31)
      t58 = t50**2
      t61 = dsqrt(4*l1*t26*z*l5*t29*t31*t34-4*t25*t39*t41*t43*t46*t50*t3
     #4*t53+t58*t34)
      t62 = 64*t61
      t66 = t61*l1
      t68 = -64+384*t12+t3-384*t14-128*t17+128*t19-192*l1-192*t6-128*t11
     #+t62+(-128*t2-128*l1+128*t6)*s2l1+128*t66
      t72 = s2l1**2
      t81 = t61*t11
      t84 = t61*z
      t88 = t11**2
      t91 = t88*z
      t95 = t1**2
      t100 = t41*l5
      t101 = t39*t100
      t102 = t101*t43
      t103 = t46*l1
      t105 = t103*t53*l6
      t112 = l6**2
      t113 = t1*t112
      t117 = t46*t53*l6
      t120 = -2+2*t81*t1-2*t84*t16+t81*t16-2*t88*l1-t91*t1+4*t91*t16-6*t
     #91*t95+26*t88*t95+8*t102*t105-5*t2+t6-4*t12+13*t14+3*l1+16*t47*t11
     #3-8*t102*t117
      t122 = l5**2
      t124 = t39*t41
      t125 = t124*t43
      t135 = t43*t46*t53
      t140 = t11*z
      t151 = t61*t16-t122*t1-t16-8*t125*t105+9*t17-27*t19+2*t61-z+2*l5+8
     #*t125*t117-t122+4*t101*t135+4*t84*l1-20*t140*t1+6*t140*l1-22*t88*t
     #16+37*t140*t16+10*t88*t1
      t156 = l1*t112
      t159 = t122*z
      t162 = t1*l6
      t165 = l1*l6
      t187 = t43*t103*t53
      t191 = t39*z*t100
      t194 = -3*t66-41*t140*t95-16*t47*t156+16*t159*t156+16*t159*t162-16
     #*t159*t165-16*t159*t113+16*t47*t165-16*t47*t162+6*t159*l1+6*t47*t1
     #-t122*t11*t1-6*t159*t1-4*t49-7*z*t95-6*t48+4*t124*t187-4*t191*t187
      t201 = t95*t1
      t214 = t95*l1
      t230 = -4*t124*t135+28*t11*t95+2*l5*t1-t140-t201*t91-6*t140*t201-4
     #*t101*t187-2*t84*t1+8*t191*t43*t105+2*t122*l1-16*t11*t214+24*t140*
     #t214-16*t88*t214+4*t91*t214-t201*z+4*t201*t88+4*t201*t11+4*z*t214
      t238 = dlog(1-l3+l3*l1)
      t239 = z*l3
      t242 = dlog(l1+z-t239-t6+t239*l1)
      t243 = dlog(z)
      taux13 = ((-t3-128*l1-128*z+256*t6)*t9+t68*s1l1+64*t31*(-t6+z)*t72
     #+(-64-64*t2-64*l1-64*z+128*t6+t62)*s2l1+64*(t120+t151+t194+t230)*t
     #34)*(t238-t242+t243)
      t1 = -1+l1
      t2 = z**2
      t4 = t2*l1
      t6 = l1**2
      t7 = z*t6
      t9 = t2*t6
      t11 = z*l1
      t20 = t2*z
      t29 = t6*l1
      t38 = 2*l6-1
      t39 = t38**2
      t42 = l5-1
      t45 = (1-l1+t11)**2
      t46 = 1/t45
      t51 = dsqrt(z)
      t53 = dsqrt(l5)
      t55 = dsqrt(-t42)
      t58 = dsqrt(l1)
      t62 = l5*z*l1-l1+1-l5+l5*l1
      t65 = dsqrt(-t1)
      t70 = t62**2
      t73 = dsqrt(4*l1*t39*z*l5*t42*t1*t46-4*t38*t51*t53*t55*t58*t62*t46
     #*t65+t70*t46)
      t75 = 2+2*t20+4*t4+6*t20*t6-6*t20*l1+2*t7-8*t9-2*z*t29+4*t2*t29-2*
     #t20*t29+2*t11-2*t73
      t81 = dlog(1-l3+l3*l1)
      t82 = z*l3
      t85 = dlog(l1+z-t82-t11+t82*l1)
      t86 = dlog(z)
      taux14 = (-128*t1*(2+4*t2-8*t4-4*t7+4*t9+6*t11)*s1l1-128*t1*(2-2*z
     #+2*t11)*s2l1-128*t1*t75)*(t81-t85+t86)

CCC Coefficients

      tcoff23 = taux13
      tcoff24 = taux13
      tcoff25 = taux14
      tcoff26 = -taux14

CCC Jet functions

      tFjet23 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff23)
      tFjet24 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff24)
      tFjet25 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff25)
      tFjet26 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff26)

CCC Sum

      RVnum7e1 = tFjet23*tcoff23+tFjet24*tcoff24+tFjet25*tcoff25+tFjet26
     #*tcoff26
      RETURN
      END

      double precision function RVnum7e2(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t1 = l1**2
      t2 = z*t1
      t3 = 192*t2
      t4 = 128*l1
      t5 = z*l1
      t7 = 128*z
      t9 = s1l1**2
      t11 = l5**2
      t12 = t11*z
      t13 = l6**2
      t14 = t1*t13
      t15 = t12*t14
      t17 = l1*l6
      t18 = t12*t17
      t20 = z*l5
      t21 = t1*l6
      t22 = t20*t21
      t24 = dsqrt(z)
      t25 = dsqrt(l5)
      t26 = t25*l5
      t27 = t24*t26
      t28 = 1-l5
      t29 = dsqrt(t28)
      t30 = t27*t29
      t31 = dsqrt(l1)
      t32 = 1-l1
      t33 = dsqrt(t32)
      t35 = t31*t33*l6
      t36 = t30*t35
      t39 = 2*l6-1
      t40 = t39**2
      t45 = -t32
      t47 = (1-l1+t5)**2
      t48 = 1/t47
      t56 = t20*l1
      t57 = l5*l1
      t58 = t56-l1+1-l5+t57
      t64 = t58**2
      t67 = dsqrt(-4*l1*t40*z*l5*t28*t45*t48-4*t39*t24*t25*t29*t31*t58*t
     #48*t33+t64*t48)
      t68 = t1*l1
      t70 = z**2
      t71 = t70*z
      t72 = t71*l1
      t75 = t29*t31*t33
      t76 = t27*t75
      t78 = t11*t1
      t79 = t1**2
      t80 = t70*t79
      t82 = t11*l1
      t84 = t31*l1
      t86 = t29*t84*t33
      t87 = t27*t86
      t89 = t79*l1
      t90 = t70*t89
      t92 = t71*t89
      t94 = t70**2
      t95 = t94*t89
      t97 = t94*z
      t98 = t97*t89
      t100 = t79*t1
      t101 = t100*z
      t102 = -2-16*t15-16*t18-16*t22-8*t36+t67*t68+6*t72+4*t76-t78+28*t8
     #0+2*t82-4*t87-16*t90+24*t92-16*t95+4*t98-t101
      t103 = t100*t70
      t105 = z*t89
      t107 = t100*t94
      t110 = t24*t25
      t111 = t110*t75
      t113 = t67*z
      t116 = t67*t70
      t118 = t100*t97
      t119 = t71*t100
      t121 = t110*t86
      t124 = t84*t33*l6
      t125 = t30*t124
      t128 = t67*l1
      t130 = t70*t1
      t132 = t70*l1
      t134 = t71*t1
      t136 = 4*t103+4*t105+4*t107-5*t2-4*t111-2*t113*t68+t116*t68-z-t118
     #-6*t119+4*t121+8*t125+t5+2*l5-3*t128+13*t130-4*t132-20*t134
      t138 = t94*l1
      t140 = t97*t1
      t141 = z*t79
      t143 = t12*l1
      t145 = t94*t79
      t147 = t97*t79
      t149 = z*t68
      t151 = t97*t68
      t153 = t71*t79
      t156 = t24*z*t26
      t157 = t156*t86
      t159 = t113*t1
      t161 = t94*t68
      t163 = t113*l1
      t165 = t70*t68
      t167 = t71*t68
      t169 = -2*t138-t140-7*t141+6*t143+26*t145-t11-6*t147-t71+9*t149-t6
     #8+4*t151-41*t153-4*t157-2*t159-22*t161+4*t163-27*t165+37*t167
      t172 = t11*t70*t1
      t173 = t12*t1
      t175 = t110*t29
      t176 = t175*t124
      t178 = t175*t35
      t180 = t116*t1
      t183 = t94*t1
      t185 = l5*t1
      t188 = t156*t29*t124
      t190 = t20*t14
      t192 = t12*t21
      t194 = l1*t13
      t195 = t12*t194
      t197 = t20*t194
      t199 = t20*t17
      t203 = t20*t1
      t205 = -4*t57-t172-6*t173-8*t176+8*t178+2*t180+3*l1+10*t183+2*t185
     #+8*t188+16*t190+16*t192+16*t195-16*t197+16*t199+2*t67-6*t56+6*t203
      t213 = 128*t128
      t214 = 128*t70
      t215 = 64*t67
      t218 = 128*t149
      t222 = -64-192*t5-192*l1-384*t130+t213-t214+t3+t215+128*t165+384*t
     #132-t218+(-128*l1+128*t5-128*t2)*s2l1
      t226 = s2l1**2
      t236 = dlog(z)
      t237 = z*l3
      t240 = dlog(l1+z-t237-t5+t237*l1)
      t242 = t240**2
      t244 = t236**2
      t248 = dlog(1-l3+l3*l1)
      t249 = t248**2
      t257 = 640*t2
      t263 = 256*t67
      t277 = t7-t257+448*t5+320*t149-t263+960*t130-320*t165+384*t70+t4-1
     #024*t132-t213+(256+256*z+256*l1-896*t5+256*t68+896*t2-512*t1-256*t
     #149)*s2l1+384
      t282 = 6*z
      t293 = 384+64*t165-64*t149-t4+t214+576*t2+t213-192*t132-1024*t5+51
     #2*z-t263
      t312 = 6+32*t15+32*t18+32*t22+16*t36-22*t72-8*t76+2*t78-125*t80-4*
     #t82+8*t87-6*t67*t1+68*t90-96*t92+60*t95-14*t98+4*t101
      t329 = -16*t103-18*t105-16*t107+30*t2+8*t111+t282+4*t118+24*t119-8
     #*t121-16*t125+6*t1-14*t5-4*l5+12*t128-71*t130+21*t132+76*t134
      t348 = 8*t138+4*t140+36*t141-12*t143-93*t145+2*t11+20*t147+4*t71-4
     #4*t149-14*t151+162*t153+8*t157+12*t159+77*t161-12*t163+127*t165-14
     #6*t167
      t367 = 8*t57+2*t172+12*t173+16*t176-16*t178-6*t180-12*l1-36*t183-4
     #*t185-16*t188-32*t190-32*t192-32*t195+32*t197-32*t199-6*t67+12*t56
     #-12*t203
      taux13 = ((-t3-t4+256*t5-t7)*t9+64*(t102+t136+t169+t205)*t48+t222*
     #s1l1+64*t45*(z-t5)*t226+(-64+t215-64*l1-64*t2+128*t5-64*z)*s2l1)*(
     #t236*t240-t242/2-t244/2+t249/2)+((128*t68-384*t1+512*l1+384*z+t257
     #-t218-768*t5)*t9+t277*s1l1+64*t45*(8*t5-2-2*t2+2*t1-t282)*t226+t29
     #3*s2l1+64*(t312+t329+t348+t367)*t48)*(t248-t240+t236)
      t1 = -1+l1
      t2 = z**2
      t3 = t2*l1
      t5 = z*l1
      t7 = l1**2
      t8 = z*t7
      t10 = t2*t7
      t12 = 4*t2
      t13 = 2-8*t3+6*t5-4*t8+4*t10+t12
      t17 = 2-2*z+2*t5
      t20 = 8*t10
      t22 = 2*l6-1
      t23 = t22**2
      t26 = l5-1
      t29 = (1-l1+t5)**2
      t30 = 1/t29
      t35 = dsqrt(z)
      t37 = dsqrt(l5)
      t39 = dsqrt(-t26)
      t42 = dsqrt(l1)
      t46 = l5*z*l1-l1+1-l5+l5*l1
      t48 = -t1
      t49 = dsqrt(t48)
      t54 = t46**2
      t57 = dsqrt(4*l1*t23*z*l5*t26*t1*t30-4*t22*t35*t37*t39*t42*t46*t30
     #*t49+t54*t30)
      t58 = 2*t57
      t59 = t2*z
      t60 = t59*l1
      t62 = t7*l1
      t63 = z*t62
      t64 = 2*t63
      t66 = 2*t8
      t67 = t2*t62
      t68 = 4*t67
      t69 = 2*t5
      t70 = t59*t7
      t74 = 2*t59*t62
      t75 = 2-t20-t58-6*t60-t64+4*t3+t66+t68+t69+6*t70+2*t59-t74
      t79 = dlog(z)
      t80 = z*l3
      t83 = dlog(l1+z-t80-t5+t80*l1)
      t85 = t83**2
      t87 = t79**2
      t91 = dlog(1-l3+l3*l1)
      t92 = t91**2
      t97 = 4*t63
      t100 = 4*z
      t101 = 2*t67
      t102 = 2*t62
      t105 = 10*t10
      t107 = -4+6*l1-t97-18*t5-6*t7+t100+t101+t102+16*t8-8*t2-t105+16*t3
      t109 = 4*t48
      t118 = -12*t5+12*t8+t100-t12+t102-4*t7+2*l1-2+10*t3+t101-t97-t20
      t127 = -2-t68+t74-t69-4*t59+10*t60+t58-t66-8*t70+t64-6*t3+t105
      taux14 = (-128*t1*t13*s1l1-128*t1*t17*s2l1-128*t1*t75)*(t79*t83-t8
     #5/2-t87/2+t92/2)+(-128*(t1*t107+t109*t13)*s1l1-128*(t1*t118+t109*t
     #17)*s2l1-128*t1*t127-128*t109*t75)*(t91-t83+t79)

CCC Coefficients

      tcoff23 = taux13
      tcoff24 = taux13
      tcoff25 = taux14
      tcoff26 = -taux14

CCC Jet functions

      tFjet23 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff23)
      tFjet24 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff24)
      tFjet25 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff25)
      tFjet26 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff26)

CCC Sum

      RVnum7e2 = tFjet23*tcoff23+tFjet24*tcoff24+tFjet25*tcoff25+tFjet26
     #*tcoff26
      RETURN
      END

      double precision function RVnum7e3(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t1 = 128*l1
      t2 = z*l1
      t4 = l1**2
      t5 = z*t4
      t6 = 192*t5
      t7 = 128*z
      t9 = s1l1**2
      t11 = dsqrt(z)
      t12 = dsqrt(l5)
      t13 = t12*l5
      t14 = t11*t13
      t15 = 1-l5
      t16 = dsqrt(t15)
      t17 = dsqrt(l1)
      t19 = 1-l1
      t20 = dsqrt(t19)
      t21 = t16*t17*t20
      t22 = t14*t21
      t23 = 4*t22
      t24 = t11*t12
      t25 = t24*t16
      t27 = t17*t20*l6
      t28 = t25*t27
      t29 = 8*t28
      t30 = z**2
      t31 = t30**2
      t32 = t31*z
      t33 = t4**2
      t34 = t32*t33
      t36 = t33*t4
      t37 = t36*t31
      t39 = t4*l1
      t40 = t31*t39
      t42 = t17*l1
      t44 = t42*t20*l6
      t45 = t25*t44
      t46 = 8*t45
      t47 = t33*l1
      t48 = z*t47
      t50 = t36*t30
      t53 = 2*l6-1
      t54 = t53**2
      t59 = -t19
      t61 = (1-l1+t2)**2
      t62 = 1/t61
      t70 = z*l5
      t71 = t70*l1
      t72 = l5*l1
      t73 = t71-l1+1-l5+t72
      t79 = t73**2
      t82 = dsqrt(-4*l1*t54*z*l5*t15*t59*t62-4*t53*t11*t12*t16*t17*t73*t
     #62*t20+t79*t62)
      t84 = t82*z
      t85 = t84*t4
      t87 = t70*t4
      t88 = 6*t87
      t89 = t30*t4
      t91 = t30*l1
      t93 = l6**2
      t94 = t4*t93
      t95 = t70*t94
      t96 = 16*t95
      t97 = l5**2
      t98 = t97*z
      t99 = l1*t93
      t100 = t98*t99
      t101 = 16*t100
      t102 = t70*t99
      t103 = 16*t102
      t104 = -2+t23+t29-6*t34+4*t37-22*t40-t46+4*t48+4*t50+t82*t39-2*t85
     #+t88+13*t89-4*t91+t96+t101-t103
      t105 = l1*l6
      t106 = t98*t105
      t107 = 16*t106
      t108 = t31*t4
      t110 = t82*l1
      t112 = t14*t16
      t113 = t112*t44
      t114 = 8*t113
      t115 = t30*z
      t116 = t115*t33
      t119 = t16*t42*t20
      t120 = t24*t119
      t121 = 4*t120
      t122 = t97*t4
      t123 = t97*l1
      t124 = 2*t123
      t126 = t11*z*t13
      t127 = t126*t119
      t128 = 4*t127
      t129 = t112*t27
      t130 = 8*t129
      t131 = z*t33
      t134 = t126*t16*t44
      t135 = 8*t134
      t136 = 6*t71
      t137 = t115*t4
      t139 = -t107+10*t108-3*t110+t114-41*t116+t121-t39-t122+t124-z-t128
     #-t130-7*t131-t115-t97+t135-t136-20*t137
      t141 = t31*t33
      t144 = t82*t30
      t145 = t144*t4
      t148 = t98*l1
      t149 = 6*t148
      t150 = t24*t21
      t151 = 4*t150
      t152 = t32*t39
      t154 = t32*t4
      t155 = t14*t119
      t156 = 4*t155
      t157 = t30*t33
      t162 = t115*l1
      t164 = z*t39
      t166 = t115*t39
      t168 = t30*t39
      t171 = 26*t141+2*t82+2*t145+3*l1+t149-t151+4*t152-t154-t156+28*t15
     #7-2*t84*t39+t144*t39+6*t162+9*t164+37*t166-27*t168+t2-5*t5
      t172 = t4*l6
      t173 = t98*t172
      t174 = 16*t173
      t175 = t70*t105
      t176 = 16*t175
      t177 = t98*t94
      t178 = 16*t177
      t179 = t70*t172
      t180 = 16*t179
      t181 = t30*t47
      t183 = t115*t47
      t185 = t31*t47
      t187 = t115*t36
      t189 = t36*t32
      t190 = t36*z
      t191 = t84*l1
      t193 = t32*t47
      t195 = 2*l5
      t196 = t31*l1
      t198 = 4*t72
      t199 = l5*t4
      t200 = 2*t199
      t202 = t97*t30*t4
      t203 = t98*t4
      t204 = 6*t203
      t205 = t174+t176-t178-t180-16*t181+24*t183-16*t185-6*t187-t189-t19
     #0+4*t191+4*t193+t195-2*t196-t198+t200-t202-t204
      t214 = 64*t82
      t215 = 128*t168
      t216 = 128*t110
      t219 = 384*t91
      t220 = 128*t164
      t221 = 128*t30
      t222 = -64-192*t2+(-128*t5-128*l1+128*t2)*s2l1+t214+t215+t216-192*
     #l1-384*t89+t219+t6-t220-t221
      t226 = s2l1**2
      t229 = 64*l1
      t230 = 64*t5
      t236 = dlog(z)
      t237 = z*l3
      t240 = dlog(l1+z-t237-t2+t237*l1)
      t241 = t240**2
      t244 = t236**2
      t253 = dlog(1-l3+l3*l1)
      t254 = t253**2
      t259 = 128*t39
      t260 = 640*t5
      t262 = 384*z
      t263 = 768*t2
      t268 = 384*t30
      t272 = 256*t164
      t276 = 256*z
      t281 = 256*t82
      t283 = 320*t164+384+t268+448*t2-320*t168+(-512*t4+256-t272+256*t39
     #+896*t5-896*t2+t276+256*l1)*s2l1-1024*t91-t216-t281+t1+960*t89+t7-
     #t260
      t286 = 2*t5
      t288 = 6*z
      t293 = 64*t168
      t299 = 384+t293-192*t91+512*z+t216-t281-t1+576*t5-1024*t2+t221-64*
     #t164
      t311 = t82*t4
      t318 = 6-8*t22-16*t28+20*t34-16*t37+77*t40+16*t45-18*t48-16*t50+12
     #*t85-12*t87-6*t311-71*t89+21*t91-32*t95-32*t100+32*t102
      t335 = 32*t106-36*t108+12*t110-16*t113+162*t116-8*t120+2*t122-4*t1
     #23+t288+8*t127+16*t129+36*t131+4*t115+2*t97-16*t134+12*t71+76*t137
      t354 = -93*t141-6*t82-6*t145-12*l1-12*t148+8*t150-14*t152+4*t154+8
     #*t155-125*t157-22*t162-44*t164-146*t166+127*t168-14*t2+30*t5-32*t1
     #73
      t373 = -32*t175+32*t177+32*t179+68*t181-96*t183+60*t185+24*t187+4*
     #t189+4*t190-12*t191-14*t193-4*l5+8*t196+8*t72+6*t4-4*t199+2*t202+1
     #2*t203
      t391 = 256*t89
      t395 = 128*t5
      t396 = -256*t4+t215+t263-t391-640*z+t1-128+t259+128*t91+t395-t272
      t398 = 512*t2
      t400 = t215-t262+t216-t220-t391+t214+t396*s2l1+t398+t219-t1-128-25
     #6*t30-t395
      t413 = -t215+t220-256-t262-t216+t1+192*t82+640*t91+t398-128*t89-t2
     #68-256*t5
      t424 = -3+t23+t29-16*t34+20*t37-69*t40-t46+23*t48+20*t50-6*t85+t88
     #+3*t311+85*t89-25*t91+t96+t101-t103
      t431 = -t107+37*t108-6*t110+t114-167*t116+t121-t122+t124-t288-t128
     #-t130+t30-47*t131-5*t115-t97+t135-t136-77*t137
      t446 = 84*t141+3*t82+3*t145+6*l1+t149-t151+13*t152-5*t154-t156+146
     #*t157+25*t162+58*t164+145*t166-147*t168+21*t2-44*t5+t174
      t457 = t176-t178-t180-82*t181+108*t183-62*t185-30*t187-5*t189-5*t1
     #90+6*t191+13*t193+t195-10*t196-t198-3*t4+t200-t202-t204
      taux13 = ((-t1+256*t2-t6-t7)*t9+64*(t104+t139+t171+t205)*t62+t222*
     #s1l1+64*t59*(z-t2)*t226+(-64-t229-t230+128*t2-64*z+t214)*s2l1)*(t2
     #36*t241/2-t244*t240/2-t241*t240/6+t244*t236/6+t254*t253/6)+((t259-
     #t220+t260-384*t4+t262-t263+512*l1)*t9+t283*s1l1+64*t59*(8*t2-2-t28
     #6+2*t4-t288)*t226+t299*s2l1+64*(t318+t335+t354+t373)*t62)*(t236*t2
     #40-t241/2-t244/2+t254/2)+((-t276+t229-64+384*t2-t220-64*t89+64*t39
     #-192*t4+t293-t230)*t9+t400*s1l1+64*t59*(-t286+2-2*t91+t89-3*t2+5*z
     #+t4+t30)*t226+t413*s2l1+64*(t424+t431+t446+t457)*t62)*(t253-t240+t
     #236)
      t1 = -1+l1
      t2 = z**2
      t3 = l1**2
      t4 = t2*t3
      t6 = t2*l1
      t8 = z*l1
      t9 = 6*t8
      t10 = z*t3
      t11 = 4*t10
      t12 = 4*t2
      t13 = 2+4*t4-8*t6+t9-t11+t12
      t17 = 2-2*z+2*t8
      t20 = t3*l1
      t21 = z*t20
      t22 = 2*t21
      t23 = t2*z
      t25 = 2*t23*t20
      t26 = t2*t20
      t27 = 4*t26
      t29 = 2*l6-1
      t30 = t29**2
      t33 = l5-1
      t36 = (1-l1+t8)**2
      t37 = 1/t36
      t42 = dsqrt(z)
      t44 = dsqrt(l5)
      t46 = dsqrt(-t33)
      t49 = dsqrt(l1)
      t53 = l5*z*l1-l1+1-l5+l5*l1
      t55 = -t1
      t56 = dsqrt(t55)
      t61 = t53**2
      t64 = dsqrt(4*l1*t30*z*l5*t33*t1*t37-4*t29*t42*t44*t46*t49*t53*t37
     #*t56+t61*t37)
      t65 = 2*t64
      t66 = 2*t8
      t67 = t23*l1
      t69 = 8*t4
      t70 = 4*t6
      t71 = t23*t3
      t73 = 2*t10
      t74 = 2*t23
      t75 = 2-t22-t25+t27-t65+t66-6*t67-t69+t70+6*t71+t73+t74
      t79 = dlog(z)
      t80 = z*l3
      t83 = dlog(l1+z-t80-t8+t80*l1)
      t84 = t83**2
      t87 = t79**2
      t96 = dlog(1-l3+l3*l1)
      t97 = t96**2
      t103 = 4*t21
      t106 = 10*t4
      t109 = 2*t26
      t110 = 4*z
      t111 = 2*t20
      t113 = 16*t10-t103+16*t6-8*t2-t106-6*t3-18*t8+t109+t110+t111+6*l1-
     #4
      t115 = 4*t55
      t119 = 4*t3
      t124 = -t119-12*t8-t69+t111-2-t12+12*t10-t103+t110+t109+2*l1+10*t6
      t133 = 10*t67-6*t6-t66+t106+t22-t27-2+t25+t65-8*t71-4*t23-t73
      t147 = -t119-t103+6*t10+t111-l1+t66-7*z+t109-2*t4+3-t6+t12
      t150 = -t115
      t167 = 6*t4-t66+z+t74-2*t67-t27+t25-t11-2*t71+t22-3*t2+t70
      taux14 = (-128*t1*t13*s1l1-128*t1*t17*s2l1-128*t1*t75)*(t79*t84/2-
     #t87*t83/2-t84*t83/6+t87*t79/6+t97*t96/6)+(-128*(t1*t113+t115*t13)*
     #s1l1-128*(t1*t124+t115*t17)*s2l1-128*t1*t133-128*t115*t75)*(t79*t8
     #3-t84/2-t87/2+t97/2)+(-128*(t1*t147+t115*t113+t150*t13)*s1l1-128*(
     #t1*(-2*t3+t2+t109-3*l1-3*t6+t73+t9+t111-t103-z)+t115*t124+t150*t17
     #)*s2l1-128*t1*t167-128*t115*t133-128*t150*t75)*(t96-t83+t79)

CCC Coefficients

      tcoff23 = taux13
      tcoff24 = taux13
      tcoff25 = taux14
      tcoff26 = -taux14

CCC Jet functions

      tFjet23 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff23)
      tFjet24 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff24)
      tFjet25 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff25)
      tFjet26 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff26)

CCC Sum

      RVnum7e3 = tFjet23*tcoff23+tFjet24*tcoff24+tFjet25*tcoff25+tFjet26
     #*tcoff26
      RETURN
      END

      double precision function RVnum7e4(var,ndim,z,ux,l1,
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
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
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

      t1 = 128*l1
      t2 = z*l1
      t3 = 256*t2
      t4 = l1**2
      t5 = z*t4
      t6 = 192*t5
      t7 = 128*z
      t9 = s1l1**2
      t12 = 2*l6-1
      t13 = t12**2
      t16 = l5-1
      t18 = -1+l1
      t20 = (1-l1+t2)**2
      t21 = 1/t20
      t26 = dsqrt(z)
      t28 = dsqrt(l5)
      t30 = dsqrt(-t16)
      t33 = dsqrt(l1)
      t34 = l5*z
      t35 = t34*l1
      t36 = l5*l1
      t37 = t35-l1+1-l5+t36
      t40 = dsqrt(-t18)
      t45 = t37**2
      t48 = dsqrt(4*l1*t13*z*l5*t16*t18*t21-4*t12*t26*t28*t30*t33*t37*t2
     #1*t40+t45*t21)
      t49 = t48*z
      t50 = t49*l1
      t52 = z**2
      t53 = t48*t52
      t54 = t53*t4
      t56 = l5**2
      t58 = t56*t52*t4
      t59 = t56*z
      t60 = t59*t4
      t61 = 6*t60
      t62 = t4*l1
      t66 = 2*l5
      t67 = t59*l1
      t68 = 6*t67
      t69 = t28*l5
      t70 = t26*t69
      t71 = t33*l1
      t73 = t30*t71*t40
      t74 = t70*t73
      t75 = 4*t74
      t76 = t26*t28
      t77 = t76*t30
      t79 = t71*t40*l6
      t80 = t77*t79
      t81 = 8*t80
      t82 = 4*t36
      t83 = t48*l1
      t85 = t52**2
      t86 = t85*t4
      t88 = t85*t62
      t90 = t85*z
      t91 = t4**2
      t92 = t90*t91
      t94 = t85*t91
      t96 = -2+4*t50+2*t54-t58-t61-2*t49*t62+t53*t62+t66+t68-t75-t81-t82
     #-3*t83+10*t86-22*t88-6*t92+26*t94
      t97 = z*t91
      t99 = t52*z
      t100 = t99*t91
      t103 = t30*t33*t40
      t104 = t70*t103
      t105 = 4*t104
      t106 = t91*l1
      t107 = t52*t106
      t109 = t99*t106
      t111 = t85*t106
      t112 = 16*t111
      t113 = t91*t4
      t114 = t99*t113
      t116 = t113*t90
      t117 = t90*t106
      t119 = t113*z
      t120 = t113*t85
      t122 = z*t106
      t124 = t113*t52
      t126 = z*t62
      t128 = t90*t4
      t129 = t85*l1
      t131 = t49*t4
      t133 = -7*t97-41*t100+t105-16*t107+24*t109-t112-6*t114-t116+4*t117
     #-t119+4*t120+4*t122+4*t124+t2+9*t126-t128-2*t129-2*t131
      t135 = t34*t4
      t136 = 6*t135
      t137 = t99*t62
      t139 = t52*t62
      t141 = t76*t73
      t142 = 4*t141
      t143 = t52*t91
      t145 = t90*t62
      t147 = l6**2
      t148 = l1*t147
      t149 = t34*t148
      t150 = 16*t149
      t151 = l1*l6
      t152 = t59*t151
      t153 = 16*t152
      t154 = t4*l6
      t155 = t59*t154
      t156 = 16*t155
      t157 = t34*t151
      t158 = 16*t157
      t159 = t4*t147
      t160 = t59*t159
      t161 = 16*t160
      t162 = t34*t154
      t163 = 16*t162
      t165 = t34*t159
      t166 = 16*t165
      t167 = t59*t148
      t168 = 16*t167
      t170 = t26*z*t69
      t171 = t170*t73
      t172 = 4*t171
      t173 = t56*t4
      t174 = t136+37*t137-27*t139+t142+28*t143-z+4*t145-t150-t153+t156+t
     #158-t161-t163+t48*t62+t166+t168-t172-t173
      t175 = t56*l1
      t176 = 2*t175
      t178 = l5*t4
      t179 = 2*t178
      t180 = t52*t4
      t182 = t52*l1
      t184 = t99*t4
      t186 = t70*t30
      t188 = t33*t40*l6
      t189 = t186*t188
      t190 = 8*t189
      t191 = t76*t103
      t192 = 4*t191
      t196 = t170*t30*t79
      t197 = 8*t196
      t198 = t99*l1
      t200 = t186*t79
      t201 = 8*t200
      t202 = t77*t188
      t203 = 8*t202
      t204 = 6*t35
      t205 = t176-5*t5+t179+13*t180-4*t182-20*t184-t190-t192-t99+2*t48+3
     #*l1+t197-t62+6*t198+t201+t203-t204-t56
      t214 = 64*t48
      t215 = 128*t139
      t216 = 128*t83
      t217 = 192*l1
      t219 = 384*t182
      t220 = 128*t126
      t221 = 128*t52
      t222 = -64-192*t2+(-128*t5-128*l1+128*t2)*s2l1+t214+t215+t216-t217
     #-384*t180+t219+t6-t220-t221
      t226 = s2l1**2
      t229 = 64*l1
      t230 = 64*t5
      t231 = 128*t2
      t238 = dlog(1-l3+l3*l1)
      t239 = t238**2
      t240 = t239**2
      t242 = dlog(z)
      t243 = t242**2
      t244 = t243**2
      t246 = z*l3
      t249 = dlog(l1+z-t246-t2+t246*l1)
      t250 = t249**2
      t251 = t250*t249
      t254 = t250**2
      t258 = t243*t242
      t263 = 128*t62
      t264 = 640*t5
      t266 = 384*z
      t267 = 768*t2
      t272 = 384*t52
      t276 = 256*t126
      t280 = 256*z
      t285 = 256*t48
      t287 = 320*t126+384+t272+448*t2-320*t139+(-512*t4+256-t276+256*t62
     #+896*t5-896*t2+t280+256*l1)*s2l1-1024*t182-t216-t285+t1+960*t180+t
     #7-t264
      t290 = 2*t5
      t292 = 6*z
      t297 = 64*t139
      t303 = 384+t297-192*t182+512*z+t216-t285-t1+576*t5-1024*t2+t221-64
     #*t126
      t321 = 6-12*t50-6*t54+2*t58+12*t60-4*l5-12*t67+8*t74+16*t80+8*t36+
     #12*t83-36*t86+77*t88+20*t92-93*t94+36*t97+162*t100
      t339 = -8*t104+68*t107-96*t109+60*t111+24*t114+4*t116-14*t117+4*t1
     #19-16*t120-18*t122-16*t124-14*t2-44*t126+4*t128+8*t129+12*t131-12*
     #t135
      t357 = -146*t137+127*t139-8*t141-125*t143+t292-14*t145+32*t149+32*
     #t152-32*t155-32*t157+32*t160+32*t162+6*t4-32*t165-32*t167+8*t171+2
     #*t173
      t362 = 21*t182
      t366 = t48*t4
      t377 = -4*t175+30*t5-4*t178-71*t180+t362+76*t184+16*t189+8*t191-6*
     #t366+4*t99-6*t48-12*l1-16*t196-22*t198-16*t200-16*t202+12*t35+2*t5
     #6
      t393 = 128*t182
      t394 = 128*t5
      t396 = 256*t180
      t397 = 256*t4
      t400 = 128*t180
      t404 = 2*t182
      t412 = 256*t5
      t413 = 64*t180
      t414 = 192*t4
      t415 = 64*t62
      t433 = -15*t86+17*t88+2*t92-14*t94+10*t97+32*t100+24*t107-30*t109+
     #t112+12*t114+2*t116-3*t117+2*t119-8*t120-7*t122-8*t124
      t449 = -7*t2-10*t126+2*t128+4*t129-35*t137+31*t139-30*t143+2*z-3*t
     #145+10*t5-33*t180+t362+36*t184-5*t52+2*t99-17*t198
      t460 = -t397+t215+t267-t396-640*z+t1-128+t263+t393+t394-t276
      t462 = 512*t2
      t464 = t215-t266+t216-t220-t396+t214+t460*s2l1+t462+t219-t1-128-25
     #6*t52-t394
      t474 = -t215+t220-256-t266-t216+t1+192*t48+640*t182+t462-t400-t272
     #-t412
      t485 = -3+6*t50+3*t54-t58-t61+t66+t68-t75-t81-t82-6*t83+37*t86-69*
     #t88-16*t92+84*t94-47*t97-167*t100
      t502 = t105-82*t107+108*t109-62*t111-30*t114-5*t116+13*t117-5*t119
     #+20*t120+23*t122+20*t124+21*t2+58*t126-5*t128-10*t129-6*t131+t136+
     #145*t137
      t508 = -147*t139+t142+146*t143-t292+13*t145-t150-t153+t156+t158-t1
     #61-t163-3*t4+t166+t168-t172-t173+t176
      t518 = -44*t5+t179+85*t180-25*t182-77*t184+t52-t190-t192+3*t366-5*
     #t99+3*t48+6*l1+t197+25*t198+t201+t203-t204-t56
      s1 = ((-t1+t3-t6-t7)*t9+64*(t96+t133+t174+t205)*t21+t222*s1l1+64*t
     #18*(z-t2)*t226+(-64-t229-t230+t231-64*z+t214)*s2l1)*(t240/24-t244/
     #24+t242*t251/6-t254/24-t243*t250/4+t258*t249/6)+((t263-t220+t264-3
     #84*t4+t266-t267+512*l1)*t9+t287*s1l1+64*t18*(8*t2-2-t290+2*t4-t292
     #)*t226+t303*s2l1+64*(t321+t339+t357+t377)*t21)*(t242*t250/2-t243*t
     #249/2-t251/6+t258/6+t239*t238/6)
      taux13 = s1+((t3-t393-t7-t394+(512*t5-t396-t231+t263-t397+t1+t393+
     #t215-t7-t276)*s2l1+t400)*s1l1+64*t18*(2*t2+t52-t404+t180-t290+t4)*
     #t226+(t394-t400-t221-t231+256*t182)*s2l1+(t412+t217-t231-t413-t414
     #+t297-64-t220+t415)*t9+64*(t433+t449)*t21)*(t238-t249+t242)+((-t28
     #0+t229-64+384*t2-t220-t413+t415-t414+t297-t230)*t9+t464*s1l1+64*t1
     #8*(-t290+2-t404+t180-3*t2+5*z+t4+t52)*t226+t474*s2l1+64*(t485+t502
     #+t508+t518)*t21)*(t242*t249-t250/2-t243/2+t239/2)
      t1 = -1+l1
      t2 = z*l1
      t3 = 6*t2
      t4 = z**2
      t5 = t4*l1
      t7 = l1**2
      t8 = z*t7
      t9 = 4*t8
      t10 = t4*t7
      t12 = 4*t4
      t13 = 2+t3-8*t5-t9+4*t10+t12
      t17 = 2-2*z+2*t2
      t21 = 2*l6-1
      t22 = t21**2
      t25 = l5-1
      t28 = (1-l1+t2)**2
      t29 = 1/t28
      t34 = dsqrt(z)
      t36 = dsqrt(l5)
      t38 = dsqrt(-t25)
      t41 = dsqrt(l1)
      t45 = l5*z*l1-l1+1-l5+l5*l1
      t47 = -t1
      t48 = dsqrt(t47)
      t53 = t45**2
      t56 = dsqrt(4*l1*t22*z*l5*t25*t1*t29-4*t21*t34*t36*t38*t41*t45*t29
     #*t48+t53*t29)
      t57 = 2*t56
      t58 = 8*t10
      t59 = z*t4
      t60 = t59*t7
      t62 = t59*l1
      t64 = t7*l1
      t65 = t4*t64
      t66 = 4*t65
      t67 = 4*t5
      t68 = z*t64
      t69 = 2*t68
      t70 = 2*t2
      t72 = 2*t59*t64
      t73 = 2*t59
      t74 = 2*t8
      t75 = 2-t57-t58+6*t60-6*t62+t66+t67-t69+t70-t72+t73+t74
      t81 = dlog(1-l3+l3*l1)
      t82 = t81**2
      t83 = t82**2
      t85 = dlog(z)
      t86 = t85**2
      t87 = t86**2
      t89 = z*l3
      t92 = dlog(l1+z-t89-t2+t89*l1)
      t93 = t92**2
      t94 = t93*t92
      t97 = t93**2
      t101 = t86*t85
      t108 = 2*t64
      t109 = 4*z
      t113 = 10*t10
      t114 = 4*t68
      t115 = 2*t65
      t117 = -4-18*t2-6*t7+t108+t109+16*t8+6*l1+16*t5-t113-t114+t115-8*t
     #4
      t119 = 4*t47
      t124 = 4*t7
      t128 = -2+2*l1-t124+t108-t58+t109-12*t2-t114+t115-t12+10*t5+12*t8
      t137 = -t66+t57-2+t69-6*t5-4*t59+10*t62-8*t60-t70+t72+t113-t74
      t154 = 3*l1
      t161 = 6*t8
      t162 = 3-l1-t124-7*z-2*t10+t161-t5+t108+t115-t114+t12+t70
      t164 = -t119
      t172 = 3*t5
      t174 = t3-t172+t115+t4-t154+t108-z+t74-2*t7-t114
      t180 = 2*t62
      t186 = z+t73-t70+t67+6*t10-3*t4-t66+t69+t72-2*t60-t180-t9
      s1 = (-128*t1*t13*s1l1-128*t1*t17*s2l1-128*t1*t75)*(t83/24-t87/24+
     #t85*t94/6-t97/24-t86*t93/4+t101*t92/6)+(-128*(t1*t117+t119*t13)*s1
     #l1-128*(t1*t128+t119*t17)*s2l1-128*t1*t137-128*t119*t75)*(t85*t93/
     #2-t86*t92/2-t94/6+t101/6+t82*t81/6)
      taux14 = s1+(-128*(t1*(t67+10*t8-5*t7+t154-5*t10+t108-7*t2-t114+t1
     #15)+t119*t162+t164*t117)*s1l1-128*(t1*(-t2+t161-3*t10-3*t7+t108+t5
     #-t114+t115)+t119*t174+t164*t128)*s2l1-128*t1*(-t9+t2+t58-4*t60+t18
     #0+t69-t172+t72-t66)-128*t119*t186-128*t164*t137)*(t81-t92+t85)+(-1
     #28*(t1*t162+t119*t117+t164*t13)*s1l1-128*(t1*t174+t119*t128+t164*t
     #17)*s2l1-128*t1*t186-128*t119*t137-128*t164*t75)*(t85*t92-t93/2-t8
     #6/2+t82/2)

CCC Coefficients

      tcoff23 = taux13
      tcoff24 = taux13
      tcoff25 = taux14
      tcoff26 = -taux14

CCC Jet functions

      tFjet23 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff23)
      tFjet24 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff24)
      tFjet25 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff25)
      tFjet26 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff26)

CCC Sum

      RVnum7e4 = tFjet23*tcoff23+tFjet24*tcoff24+tFjet25*tcoff25+tFjet26
     #*tcoff26
      RETURN
      END

