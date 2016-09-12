CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCC Numerator functions CCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRVnum1e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = t1**2
      t3 = t2*s1l1
      t6 = s2l1**2
      t7 = t6*t1
      t8 = l1**2
      t11 = s1l1**2
      t12 = t11*t2
      t15 = t1*z
      t16 = s1l1*t15
      t19 = s2l1*t15
      t20 = t8*s1l1
      t23 = s1l1*s2l1
      t24 = z*t8
      t27 = s2l1*t2
      t30 = t1*t8
      t33 = s2l1*t1
      t34 = t8*l1
      t35 = t34*s1l1
      t38 = s2l1*l1
      t41 = s2l1*z
      t44 = -64*t3*l1-576*t7*t8-64*t12*t8+64*t16*t8+896*t19*t20+384*t23*
     #t24-256*t27*t20-1024*t23*t30+384*t33*t35-512*t16*t38-128*t41*t35
      t52 = z*l1
      t59 = t1*l1
      t64 = t6*t15
      t68 = -320*t41*l1-128*t11*t15*l1-384*t35*t19-384*t23*t52+128*t27*t
     #35+128*t3*t38+768*t23*t59+64*z+192*t15+576*t64*t8-128*t2
      t75 = t6*t2
      t82 = s1l1*z
      t90 = -128*t1+128*t23*z+256*t41*t8-192*t75*t8-704*t33*t8+192*t35*t
     #1+64*t82*t8-576*t64*l1+128*t19-64*t75-192*t16
      t92 = t15*t8
      t95 = t6*z
      t105 = z*t34
      t110 = t11*t1
      t115 = -64*t27+320*t92*t11+192*t95*t8+192*t27*l1+192*t64+192*t75*l
     #1+640*t7*l1-64*t105*t6+320*t16*l1+448*t110*l1-576*t19*l1
      t118 = t15*t34
      t135 = t11*z
      t137 = -192*t118*t11-64*t105*s2l1-192*t27*t8-256*t95*l1-192*t52+64
     #*t82+704*t33*l1-128*t23*t1+128*t95-256*t7+192*t135
      t138 = t1*t34
      t143 = t2*t34
      t161 = 192*t138*t11-64*t82*l1+64*t143*s1l1-192*t33-192*t118*s1l1-6
     #4*t11*t34*z-448*t110*t8-320*t135*l1+64*t75*t34+192*t138*s2l1+128*t
     #41
      t170 = s1l1*t1
      t182 = 192*t135*t8+64*t12*t34-192*t19*t34-192*t110-192*t170*l1-320
     #*t2*t8+64*t27*t34+896*t92-384*t118-832*t30+192*t138*t6
      t201 = 640*t19*t8+320*t2*l1+128*t143+128*t170-192*t118*t6+256*t24-
     #704*t15*l1-64*t105*s1l1+384*t138-128*t170*t8-128*t105+576*t59
      taux1 = (t44+t68+t90+t115+t137+t161+t182+t201)/(1-l3*l1*(1-z))
      t2 = z**2
      t3 = s2l1*t2
      t6 = s1l1*z
      t7 = l1**2
      t10 = s2l1*z
      t15 = t2*z
      t16 = s2l1*t15
      t19 = s1l1*t15
      t25 = t2*s1l1
      t30 = 10*t3*l1+2*t6*t7+2*t10*t7-6*t10*l1-4*t16*l1-2*t19*l1+2*t19*t
     #7+2*t16+10*t25*l1-6*t3+2*t6
      t50 = -2*t2+2*t16*t7-4*z*l1-6*t25+2*z-2*t15-8*t6*l1+4*t10-4*t3*t7+
     #2*t15*l1-4*t25*t7+2*t2*l1
      t53 = -1+z
      taux2 = -128*(1-l1)*(t30+t50)*t53/(1+l3*l1*t53)

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*
     #tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoef
     #f*tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*
     #tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoef
     #f*tcoff4)

CCC Sum

      qgRVnum1e0 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcof
     #f4
      RETURN
      END

      double precision function qgRVnum1e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = t1*z
      t4 = z*l1
      t6 = s2l1**2
      t7 = t6*t2
      t8 = l1**2
      t11 = t6*z
      t15 = s2l1*t2
      t18 = t8*l1
      t19 = t2*t18
      t22 = z*t18
      t25 = s1l1*t1
      t30 = s1l1*z
      t33 = -512*t2+640*t4-768*t7*t8+384*t11*l1-384*t11-1472*t15*t8+256*
     #t19*t6+192*t22*s1l1+960*t25*t8+768*t7*l1-384*t30*t8
      t34 = s1l1*t18
      t37 = s2l1*t1
      t40 = t2*s1l1
      t43 = t6*t1
      t46 = s1l1**2
      t47 = t1**2
      t48 = t46*t47
      t52 = s2l1*t47
      t61 = t46*t1
      t66 = -576*t1*t34+1408*t37*t8-768*t40*t8+256*t43*t8+128*t48*t8-256
     #*t7+512*t52*t8+192*t22*s2l1+256*t19*t46+896*t15*l1-384*t61*l1-192*
     #t40*l1
      t72 = t6*t47
      t79 = t2*t8
      t82 = s1l1*s2l1
      t89 = t47*t18
      t94 = -256*t22*t6-768*t43*l1-384*t72*l1-448*t52*l1+256*t11*t8-384*
     #t79*t46+512*t82*t1-1088*t37*l1+576*t19*s1l1-192*t89*s1l1+576*t30*l
     #1
      t95 = t46*t18
      t100 = t46*z
      t103 = t1*t18
      t121 = -256*t95*z-128*t48*t18+640*t100*t8-576*t103*s2l1-128*t72*t1
     #8-128*t100*l1-384*t25*l1+576*t15*t18-192*t52*t18+128*t95-384*t46*t
     #8-128*t6*t8
      t124 = t1*l1
      t134 = z*t8
      t140 = t1*t8
      t142 = -1728*t124+384*t22+128*t6*t18+384*t46*l1-1152*t103+2048*t2*
     #l1-768*t134-960*t47*l1-384*t89+1152*t19+2496*t140
      t144 = t47*t8
      t147 = s2l1*z
      t151 = s2l1*l1
      t162 = t8*s1l1
      t165 = -2688*t79+960*t144+384*t61-384*t147+256*t37-256*t47*s1l1*t1
     #51-256*t100+384*t40+896*t82*t134+256*t82*t140+640*t40*t151-1152*t1
     #5*t162
      t186 = 512*t52*t162-768*t82*t124+512*t15*t34+384*t1+128*t72-256*z-
     #512*t162*s2l1+256*t34*s2l1+256*t82*l1+192*t144*s1l1+640*t147*l1
      t207 = 128*t46*t2*l1-512*t82*z-448*t147*t8+384*t72*t8-256*t52*t34+
     #128*t82*t4-512*t147*t34+128*t52-128*t46+384*t47-384*t30+512*t43
      taux1 = (t33+t66+t94+t121+t142+t165+t186+t207)/(1-l3*l1*(1-z))
      t1 = 1-l1
      t2 = z**2
      t3 = s1l1*t2
      t4 = t3*l1
      t6 = s1l1*z
      t7 = t6*l1
      t9 = s2l1*t2
      t10 = t9*l1
      t12 = t2*z
      t13 = s1l1*t12
      t14 = t13*l1
      t16 = z*l1
      t19 = s2l1*z
      t21 = s2l1*t12
      t22 = t21*l1
      t25 = t12*l1
      t27 = t2*l1
      t29 = l1**2
      t30 = t13*t29
      t32 = t9*t29
      t34 = -28*t4+24*t7-28*t10+6*t14+10*t16+16*t9-14*t19+12*t22+18*t3-6
     #*t25-4*t27-6*t30+12*t32
      t35 = t6*t29
      t37 = t3*t29
      t41 = t21*t29
      t47 = t19*t29
      t49 = t19*l1
      t55 = -6*t35+12*t37-2*s1l1*l1-6*t41-2*s2l1*l1-6*t21-8*t6-6*t47+18*
     #t49-6*z+2*s1l1+6*t12+4*t2
      t71 = 10*t10+2*t35+2*t47-6*t49-4*t22-2*t14+2*t30+2*t21+10*t4-6*t9+
     #2*t6
      t84 = -2*t2+2*t41-4*t16-6*t3+2*z-2*t12-8*t7+4*t19-4*t32+2*t25-4*t3
     #7+2*t27
      t88 = -1+z
      taux2 = -128*(t1*(t34+t55)-2*t1*(t71+t84))*t88/(1+l3*l1*t88)

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*
     #tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoef
     #f*tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*
     #tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoef
     #f*tcoff4)

CCC Sum

      qgRVnum1e1 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcof
     #f4
      RETURN
      END

      double precision function qgRVnum1e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = t1*z
      t3 = s1l1*t2
      t5 = s1l1*s2l1
      t6 = l1**2
      t7 = t1*t6
      t10 = s2l1*t2
      t11 = t6*s1l1
      t14 = t1**2
      t15 = s2l1*t14
      t18 = t1*l1
      t23 = s2l1*t1
      t26 = t6*l1
      t27 = t2*t26
      t30 = t14*t26
      t33 = s1l1*z
      t36 = s1l1**2
      t37 = t36*t1
      t40 = t36*t26
      t43 = -128*t3-768*t5*t7+640*t10*t11-256*t15*t11-256*t5*t18+128*t5*
     #t1+1280*t23*l1-384*t27*s1l1+128*t30*s1l1+384*t33*l1-128*t37*t6-256
     #*t40*z
      t44 = t36*t14
      t47 = t36*z
      t50 = t1*t26
      t53 = s2l1**2
      t54 = t53*t14
      t56 = t14*t6
      t59 = s2l1*z
      t73 = s1l1*t26
      t76 = 64*t44*t26+320*t47*t6+384*t50*s2l1-64*t54-256*t56*s1l1-640*t
     #59*l1+128*t36*t2*l1-128*t5*z+512*t59*t6-192*t54*t6+128*t5*l1+128*t
     #73*s2l1
      t78 = t14*s1l1
      t83 = s1l1*t1
      t91 = z*t26
      t96 = z*t6
      t98 = 128*t78*l1-128*t5*t2+128*t83-64*t37+448*t30+192*t36*l1+64*t4
     #0+1920*t18-448*t91+1344*t50-1856*t2*l1+896*t96
      t110 = t53*t2
      t121 = t53*t1
      t124 = 704*t14*l1-128*t33+1280*t10*t6+384*t50*t53-256*t27*t53-128*
     #t91*s1l1-384*t110*l1-128*t33*t6+384*t73*t1-1408*t23*t6+384*t3*t6-6
     #40*t121*t6
      t137 = t53*z
      t147 = -768*z*l1-128*t15-64*t44*t6+576*t2+320*z+64*t47+576*t110*t6
     #-384*t137*l1-384*t15*t6-128*t91*s2l1-192*t36*t6-64*t36
      t152 = t2*t6
      t167 = -832*t56-64*t53*t6-320*t14+2560*t152+128*t10-320*t121+64*t1
     #10+640*t5*t96-256*t27*t36-1024*t10*l1-192*t37*l1-256*t91*t53
      t191 = 576*t121*l1+192*t54*l1+384*t15*l1+320*t137*t6+64*t152*t36+3
     #84*t50*t36+256*t137-640*t1+64*t53*t26-512*t73*t10+64*t54*t26-128*t
     #47*l1
      t213 = -512*t83*l1-384*t10*t26+128*t15*t26-256*t11*s2l1-384*t23+12
     #8*t78*s2l1*l1+128*t15*t73+768*t23*t73-512*t59*t73-2624*t7+256*t59-
     #1344*t27
      taux1 = (t43+t76+t98+t124+t147+t167+t191+t213)/(1-l3*l1*(1-z))
      t1 = 1-l1
      t2 = s1l1*z
      t3 = 8*t2
      t4 = z**2
      t5 = s2l1*t4
      t7 = s2l1*z
      t8 = t7*l1
      t10 = l1**2
      t11 = t7*t10
      t14 = t10*l1
      t15 = z*t14
      t18 = t4*z
      t19 = s2l1*t18
      t20 = t19*l1
      t22 = s1l1*t18
      t23 = t22*l1
      t25 = t5*l1
      t30 = t2*l1
      t35 = s1l1*t4
      t36 = t35*l1
      t40 = t19*t10
      t44 = t35*t10
      t46 = t2*t10
      t48 = t3-8*t5-17*t8+23*t11+13*t7-9*t15*s2l1-8*t20-4*t23+18*t25-3*t
     #18*t14*s1l1-27*t30+9*t4*t14*s2l1+16*t36-3*t19*t14+7*t40-9*t15*s1l1
     #-19*t44+26*t46
      t49 = t14*s1l1
      t52 = t5*t10
      t54 = t22*t10
      t60 = 4*t4*l1
      t64 = t18*l1
      t74 = s2l1*l1
      t76 = s1l1*l1
      t82 = 9*t49*t4-22*t52+4*t54+3*z-7*s1l1+4*t19-t60+3*t49+3*s2l1*t14+
     #4*t64-3*z*t10-12*t35+6*t4*t10-3*t18*t10-4*t18+7*t74+15*t76-11*s1l1
     #*t10-8*s2l1*t10
      t99 = -t3+16*t5+18*t8-6*t11-14*t7+12*t20+6*t23-28*t25+24*t30-28*t3
     #6-6*t40+12*t44-6*t46
      t113 = 12*t52-6*t54-6*z+2*s1l1-6*t19-t60-6*t64+18*t35+6*t18+4*t4-2
     #*t74-2*t76+10*z*l1
      t117 = -1+z
      taux2 = -128*(t1*(t48+t82)-2*t1*(t99+t113))*t117/(1+l3*l1*t117)

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*
     #tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoef
     #f*tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*
     #tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoef
     #f*tcoff4)

CCC Sum

      qgRVnum1e2 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcof
     #f4
      RETURN
      END

      double precision function qgRVnum1e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = s2l1**2
      t2 = t1*z
      t5 = z**2
      t6 = t5*z
      t7 = t1*t6
      t8 = l1**2
      t13 = t1*t5
      t16 = t8*l1
      t17 = z*t16
      t20 = s2l1*t6
      t23 = t6*t16
      t24 = s1l1**2
      t28 = t5**2
      t31 = z*t8
      t35 = t5*t16
      t40 = -128*t2*l1+128*t7*t8+256*t2*t8+192*t13*l1-192*t17*t1-128*t20
     #*l1-64*t23*t24+960*t23-320*t28*t16-640*t31+256*t6*l1-960*t35+320*t
     #17-704*t5*l1
      t41 = t24*t16
      t45 = s1l1*t5
      t48 = t24*z
      t62 = s1l1*z
      t65 = s2l1*t5
      t71 = z*l1
      t75 = 64*t41+192*t24*l1-256*t45*l1-192*t48*l1+384*t48*t8-192*t41*z
     #-192*t24*t5*t8+64*t1*t16-64*t1*t8+256*t62*l1+256*t65*l1-192*t24*t8
     #-128*t62+448*t71+192*t35*t24
      t83 = t5*t8
      t87 = s1l1*t8
      t105 = s1l1*s2l1
      t108 = -64*t13-128*t65+320*t28*t8-1280*t6*t8+1600*t83+128*t20*t8-2
     #56*t87*s2l1-128*s1l1*t6*t8-256*t65*t8-128*t62*t8-64*t7*l1+256*t45*
     #t8-64*t23*t1+192*t35*t1+640*t105*t31
      t112 = s2l1*z
      t119 = s1l1*t16
      t138 = -128*z-320*t13*t8+128*t112*t8-128*t105*z-128*t112*l1+128*t1
     #19*s2l1+128*t105*l1+128*t20*t87-512*t105*t83-128*t105*t71-128*t119
     #*t20-384*t112*t119+384*t65*t119-64*t24+64*t5
      taux1 = (t40+t75+t108+t138)/(1-l3*l1*(1-z))
      t1 = 1-l1
      t2 = l1**2
      t3 = t2*l1
      t4 = s1l1*t3
      t5 = z**2
      t6 = t4*t5
      t8 = s2l1*t5
      t9 = t2*t8
      t11 = s2l1*t2
      t13 = s1l1*t2
      t15 = s1l1*l1
      t17 = s2l1*l1
      t22 = s2l1*z
      t23 = t22*l1
      t25 = t22*t2
      t27 = s2l1*t3
      t30 = t5*l1
      t32 = z*t2
      t34 = -15*t6+23*t9+13*t11+18*t13-19*t15-6*t17+2*z-7*z*l1+9*t23-31*
     #t25-5*t27-5*t4+7*t30+5*t32
      t35 = t5*t2
      t37 = t5*z
      t38 = t37*t2
      t41 = z*t3
      t42 = t41*s2l1
      t44 = t8*l1
      t47 = t37*t3*s1l1
      t50 = s1l1*z
      t51 = t50*l1
      t54 = t5*t3*s2l1
      t56 = s2l1*t37
      t57 = t56*t3
      t59 = t56*t2
      t61 = t41*s1l1
      t63 = s1l1*t5
      t64 = t63*t2
      t66 = t50*t2
      t68 = -10*t35+5*t38-2*t22+15*t42-3*t44+5*t47+6*s1l1+19*t51-15*t54+
     #5*t57-5*t59+15*t61+18*t64-36*t66
      t80 = s1l1*t37
      t94 = 9*t6-22*t9-8*t11-11*t13+15*t15+7*t17+3*z+4*t80*t2-17*t23+23*
     #t25-8*t56*l1+3*t27+3*t4-4*t30-3*t32+6*t35-3*t38+13*t22
      t117 = 4*t56-9*t42+18*t44-3*t47+4*t37*l1-12*t63+8*t50-4*t37-4*t80*
     #l1+16*t63*l1-8*t8-7*s1l1-27*t51+9*t54-3*t57+7*t59-9*t61-19*t64+26*
     #t66
      t121 = -1+z
      taux2 = -128*(t1*(t34+t68)-2*t1*(t94+t117))*t121/(1+l3*l1*t121)

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      tFjet1 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*
     #tcoff1)
      tFjet2 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoef
     #f*tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*
     #tcoff3)
      tFjet4 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoef
     #f*tcoff4)

CCC Sum

      qgRVnum1e3 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcof
     #f4
      RETURN
      END

      double precision function qgRVnum2e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = s2l1**2
      t3 = z**2
      t4 = t3**2
      t5 = t2*t4
      t7 = z*t3
      t9 = l1**2
      t10 = t9*l1
      t11 = s1l1*t10
      t12 = s2l1*t7
      t15 = t4*s1l1
      t16 = s2l1*l1
      t19 = s2l1*t4
      t22 = s1l1*s2l1
      t23 = z*l1
      t26 = s2l1*t3
      t29 = s2l1*z
      t33 = z*t9
      t36 = 64*z-64*t5+192*t7-384*t11*t12+128*t15*t16+128*t19*t11-384*t2
     #2*t23+384*t26*t11-128*t29*t11-128*t4+384*t22*t33
      t37 = t3*t9
      t40 = s1l1*t7
      t43 = t9*s1l1
      t48 = t3*l1
      t57 = t3*t10
      t60 = t7*t10
      t63 = z*t10
      t66 = -1024*t22*t37-512*t40*t16+896*t12*t43-256*t19*t43+768*t22*t4
     #8-192*t12*t10+64*t19*t10+640*t12*t9+192*t57*t2-192*t60*t2-64*t63*s
     #1l1
      t68 = s1l1*t3
      t71 = t2*t7
      t74 = s1l1*z
      t83 = t2*t3
      t86 = s1l1**2
      t87 = t86*t4
      t97 = -128*t68*t9-576*t71*l1+64*t74*t9+192*t11*t3-704*t26*t9+64*t4
     #0*t9-576*t83*t9-64*t87*t9-64*t15*l1-320*t29*l1-128*t86*t7*l1
      t105 = z*t2
      t118 = t86*t3
      t121 = 128*t22*z+256*t29*t9-192*t5*t9+128*t12-256*t105*l1+576*t71*
     #t9-192*t19*t9-64*t63*s2l1-192*t60*t86-576*t12*l1+448*t118*l1
      t136 = t7*t9
      t147 = 320*t40*l1-64*t63*t2+640*t83*l1+192*t5*l1+192*t19*l1+192*t1
     #05*t9+320*t136*t86+192*t57*t86-128*t22*t3+704*t26*l1-192*t60*s1l1
      t148 = t4*t10
      t160 = t86*z
      t172 = 64*t148*s1l1-64*t74*l1-448*t118*t9-64*t86*t10*z+64*t87*t10+
     #192*t160*t9+192*t57*s2l1+64*t5*t10-320*t160*l1-192*t68*l1-128*t3
      t187 = -192*t40+576*t48-128*t63+384*t57-704*t7*l1+256*t33+320*t4*l
     #1+128*t148+128*t68-384*t60-832*t37
      t201 = 896*t136-320*t4*t9-192*t118+128*t29-192*t26+128*t105-256*t8
     #3+192*t160+192*t71-192*t23+64*t74-64*t19
      taux3 = (t36+t66+t97+t121+t147+t172+t187+t201)/(1-l3*(1-l1)*(1-z))
      t1 = 1-l1
      t2 = z**2
      t3 = s2l1*t2
      t6 = s1l1*z
      t7 = l1**2
      t10 = s2l1*z
      t15 = t2*z
      t16 = s2l1*t15
      t19 = s1l1*t15
      t25 = s1l1*t2
      t30 = 10*t3*l1+2*t6*t7+2*t10*t7-6*t10*l1-4*t16*l1-2*t19*l1+2*t19*t
     #7+2*t16+10*t25*l1-6*t3+2*t6
      t50 = -2*t2+2*t16*t7-4*z*l1-6*t25+2*z-2*t15-8*t6*l1+4*t10-4*t3*t7+
     #2*t15*l1-4*t25*t7+2*t2*l1
      t53 = -1+z
      taux4 = -128*t1*(t30+t50)*t53/(1+l3*t1*t53)

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoef
     #f*tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*
     #tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoef
     #f*tcoff8)

CCC Sum

      qgRVnum2e0 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcof
     #f8
      RETURN
      END

      double precision function qgRVnum2e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = l1**2
      t3 = t1*t2
      t5 = s1l1*s2l1
      t6 = t2*z
      t9 = s2l1*z
      t11 = t1*z
      t12 = t2*l1
      t13 = t11*t12
      t15 = t1**2
      t16 = t15*t12
      t18 = s2l1**2
      t19 = t18*t1
      t21 = t12*s1l1
      t24 = t18*z
      t26 = s2l1*t11
      t29 = s2l1*t15
      t34 = 2496*t3+896*t5*t6-384*t9+1152*t13-384*t16+512*t19-512*t9*t21
     #-384*t24+576*t26*t12-192*t29*t12-1472*t26*t2
      t37 = z*t12
      t40 = s1l1*t1
      t48 = t1*t12
      t53 = t18*t11
      t58 = 256*t13*t18+192*t37*s1l1+960*t40*t2+256*t21*s2l1-960*t15*l1-
     #768*t6-1152*t48-512*t11+384*t37-256*z-768*t53*t2+384*t24*l1
      t62 = t1*l1
      t68 = s1l1*t11
      t69 = s2l1*l1
      t72 = s1l1*t2
      t79 = s1l1**2
      t84 = 512*t29*t2-1728*t62+2048*t11*l1+256*t5*t3+640*t68*t69-1152*t
     #26*t72+512*t29*t72-768*t5*t62-128*t79+384*t15-512*t72*s2l1
      t85 = t15*t2
      t103 = t79*t1
      t112 = 192*t85*s1l1+640*t9*l1+128*t79*t11*l1-512*t5*z-448*t9*t2+19
     #2*t37*s2l1+256*t13*t79+896*t26*l1-384*t103*l1-192*t68*l1-256*t37*t
     #18-768*t19*l1
      t124 = s1l1*z
      t127 = t79*t12
      t130 = t79*t15
      t133 = t79*z
      t138 = t18*t15
      t141 = -384*t79*t2-128*t18*t2+384*t79*l1+128*t18*t12+384*t68+576*t
     #124*l1-256*t127*z-128*t130*t12+640*t133*t2-576*t48*s2l1-128*t138*t
     #12
      t152 = z*l1
      t160 = -128*t133*l1-384*t40*l1+128*t127+128*t138+384*t1+128*t29-38
     #4*t124-256*t53+640*t152+768*t53*l1-384*t124*t2-576*t21*t1
      t162 = s2l1*t1
      t180 = t11*t2
      t185 = 1408*t162*t2-768*t68*t2+256*t19*t2+128*t130*t2+256*t5*l1+96
     #0*t85-384*t138*l1-448*t29*l1+256*t24*t2-384*t180*t79+512*t5*t1
      t207 = -1088*t162*l1+576*t13*s1l1-192*t16*s1l1+512*t21*t26-256*t15
     #*s1l1*t69-256*t29*t21+128*t5*t152+384*t138*t2+256*t162-256*t133+38
     #4*t103-2688*t180
      taux3 = (t34+t58+t84+t112+t141+t160+t185+t207)/(1-l3*(1-l1)*(1-z))
      t1 = 1-l1
      t2 = z**2
      t4 = z*l1
      t6 = s1l1*z
      t8 = t2*z
      t9 = s2l1*t8
      t15 = s1l1*t2
      t16 = l1**2
      t17 = t15*t16
      t19 = t6*t16
      t21 = s2l1*t2
      t22 = t21*t16
      t24 = s1l1*t8
      t25 = t24*t16
      t27 = s2l1*z
      t28 = t27*l1
      t30 = t27*t16
      t32 = t15*l1
      t34 = 4*t2+10*t4-8*t6-6*t9-2*s2l1*l1-2*s1l1*l1+12*t17-6*t19+12*t22
     #-6*t25+18*t28-6*t30-28*t32
      t35 = t9*t16
      t37 = t6*l1
      t39 = t21*l1
      t41 = t9*l1
      t43 = t24*l1
      t46 = t2*l1
      t48 = t8*l1
      t55 = -6*t35+24*t37-28*t39+12*t41+6*t43-14*t27-4*t46-6*t48+18*t15+
     #16*t21+6*t8-6*z+2*s1l1
      t71 = 10*t39+2*t19+2*t30-6*t28-4*t41-2*t43+2*t25+2*t9+10*t32-6*t21
     #+2*t6
      t84 = -2*t2+2*t35-4*t4-6*t15+2*z-2*t8-8*t37+4*t27-4*t22+2*t48-4*t1
     #7+2*t46
      t88 = -1+z
      taux4 = -128*(t1*(t34+t55)-2*t1*(t71+t84))*t88/(1+l3*t1*t88)

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoef
     #f*tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*
     #tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoef
     #f*tcoff8)

CCC Sum

      qgRVnum2e1 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcof
     #f8
      RETURN
      END

      double precision function qgRVnum2e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = s2l1**2
      t2 = t1*z
      t5 = z**2
      t6 = t5*z
      t7 = t1*t6
      t8 = l1**2
      t11 = t1*t5
      t14 = t8*l1
      t15 = z*t14
      t18 = s1l1**2
      t19 = t18*t5
      t22 = s2l1*t6
      t25 = t6*t14
      t29 = s1l1*t5
      t31 = t5**2
      t32 = t31*t14
      t36 = s2l1*t31
      t39 = -384*t2*l1+576*t7*t8+576*t11*l1-256*t15*t1-192*t19*l1-1024*t
     #22*l1-256*t25*t18-1344*t25+128*t29+448*t32-128*t15*s2l1-384*t36*t8
      t42 = z*t8
      t46 = t5*t14
      t49 = t5*l1
      t51 = s1l1*s2l1
      t56 = s2l1*t5
      t63 = t6*t8
      t66 = 704*t31*l1+896*t42-1856*t6*l1+1344*t46-448*t15+1920*t49-128*
     #t51*t6-384*t25*s1l1+1280*t56*l1+128*t51*t5+384*t46*t18+64*t63*t18
      t77 = s1l1*t6
      t79 = t1*t31
      t82 = s1l1*z
      t87 = t18*z
      t89 = -192*t18*t8+192*t18*l1-64*t1*t8+64*t1*t14+128*t22-128*t77-64
     #*t79-128*t36-128*t82-768*z*l1+64*t7+64*t87
      t99 = s2l1*z
      t102 = t31*t8
      t105 = t5*t8
      t107 = t18*t31
      t110 = 320*t2*t8+384*t36*l1-320*t11+256*t2-384*t56+192*t79*l1+256*
     #t99-64*t19-832*t102+2560*t63-2624*t105+64*t107*t14
      t113 = t18*t14
      t137 = -256*t113*z-128*t19*t8+384*t82*l1+128*t32*s1l1+64*t113+128*
     #t36*t14-384*t22*t14-512*t29*l1-128*t87*l1+64*t79*t14+384*t46*s2l1+
     #320*t87*t8
      t147 = t31*s1l1
      t158 = s1l1*t14
      t163 = -256*t25*t1+384*t46*t1+1280*t22*t8+320*z+128*t51*l1+128*t14
     #7*l1-64*t107*t8-640*t11*t8+384*t77*t8-1408*t56*t8+384*t158*t5-128*
     #t82*t8
      t180 = s1l1*t8
      t190 = -384*t7*l1-128*t15*s1l1+512*t99*t8-128*t51*z+128*t18*t6*l1-
     #640*t99*l1-256*t102*s1l1-256*t180*s2l1+128*t158*s2l1-192*t79*t8-76
     #8*t51*t105+576*t6
      t213 = 640*t51*t42-256*t51*t49-256*t36*t180+640*t22*t180-64*t18-32
     #0*t31-640*t5-512*t158*t22-512*t99*t158+768*t56*t158+128*t36*t158+1
     #28*t147*s2l1*l1
      taux3 = (t39+t66+t89+t110+t137+t163+t190+t213)/(1-l3*(1-l1)*(1-z))
      t1 = 1-l1
      t3 = z**2
      t4 = s1l1*t3
      t5 = l1**2
      t6 = t4*t5
      t8 = s1l1*z
      t9 = t5*t8
      t11 = t5*l1
      t12 = s1l1*t11
      t15 = s2l1*t3
      t16 = t15*t5
      t18 = t3*z
      t19 = s1l1*t18
      t20 = t19*t5
      t22 = t18*l1
      t24 = z*t11
      t27 = s2l1*t18
      t28 = t27*l1
      t30 = t19*l1
      t37 = s2l1*z
      t40 = 8*t8
      t42 = s2l1*l1
      t44 = -7*s1l1-19*t6+26*t9+9*t12*t3-22*t16+4*t20+4*t22-9*t24*s2l1-8
     #*t28-4*t30-12*t4+6*t3*t5-3*t18*t5+13*t37-8*t15+t40+4*t27+7*t42
      t45 = s1l1*l1
      t51 = t37*l1
      t53 = t37*t5
      t56 = t15*l1
      t61 = t8*l1
      t66 = t4*l1
      t70 = t27*t5
      t78 = 4*t3*l1
      t82 = 15*t45-11*s1l1*t5-8*s2l1*t5-17*t51+23*t53-4*t18+18*t56-3*t18
     #*t11*s1l1-27*t61+9*t3*t11*s2l1+16*t66-3*t27*t11+7*t70-9*t24*s1l1+3
     #*z-3*t5*z-t78+3*s2l1*t11+3*t12
      t99 = 2*s1l1+12*t6-6*t9+12*t16-6*t20-6*t22+12*t28+6*t30+18*t4-14*t
     #37+16*t15-t40-6*t27
      t113 = -2*t42-2*t45+18*t51-6*t53+10*z*l1+4*t3+6*t18-28*t56+24*t61-
     #28*t66-6*t70-6*z-t78
      t117 = -1+z
      taux4 = -128*(t1*(t44+t82)-2*t1*(t99+t113))*t117/(1+l3*t1*t117)

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoef
     #f*tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*
     #tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoef
     #f*tcoff8)

CCC Sum

      qgRVnum2e2 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcof
     #f8
      RETURN
      END

      double precision function qgRVnum2e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = l1**2
      t2 = t1*l1
      t3 = s1l1*t2
      t4 = z**2
      t5 = t4*z
      t6 = s2l1*t5
      t9 = s1l1*s2l1
      t10 = z*l1
      t13 = s2l1*t4
      t16 = s2l1*z
      t19 = t5*t2
      t21 = t4**2
      t24 = s1l1**2
      t29 = z*t2
      t30 = s2l1**2
      t33 = t30*t4
      t36 = t30*z
      t39 = t4*t2
      t44 = s1l1*z
      t47 = -128*t3*t6-128*t9*t10+384*t13*t3-384*t16*t3+960*t19-320*t21*
     #t2-64*t19*t24-128*t6*l1-192*t29*t30+192*t33*l1+256*t36*t1+192*t39*
     #t24+256*t13*l1+256*t44*l1
      t51 = t24*t2
      t54 = t24*z
      t59 = s1l1*t4
      t70 = t30*t5
      t84 = -192*t24*t4*t1-192*t51*z+384*t54*t1-192*t54*l1-256*t59*l1+12
     #8*t1*t6+192*t39*t30-64*t19*t30+256*t59*t1-64*t70*l1-128*t44*t1-256
     #*t13*t1-128*s1l1*t5*t1-320*t33*t1+128*t9*l1
      t89 = s1l1*t1
      t103 = t4*t1
      t105 = t1*z
      t115 = 320*t29+128*s2l1*t3-256*t89*s2l1-128*t16*l1-128*t9*z+128*t1
     #6*t1+320*t21*t1+64*t4+128*t70*t1+1600*t103+640*t9*t105-512*t9*t103
     #+128*t6*t89-128*t36*l1+448*t10
      t138 = -128*t44-640*t105-128*z+256*t5*l1-704*t4*l1+64*t30*t2-64*t3
     #0*t1+192*t24*l1-192*t24*t1+64*t51-64*t33-128*t13-64*t24-1280*t1*t5
     #-960*t39
      taux3 = (t47+t84+t115+t138)/(1-l3*(1-l1)*(1-z))
      t1 = 1-l1
      t2 = z**2
      t3 = l1**2
      t4 = t3*l1
      t6 = t2*t4*s2l1
      t8 = s1l1*z
      t9 = t8*l1
      t11 = t2*z
      t13 = t11*t4*s1l1
      t15 = s2l1*t2
      t16 = t15*l1
      t18 = z*t4
      t19 = t18*s2l1
      t21 = t8*t3
      t23 = s1l1*t2
      t24 = t23*t3
      t26 = t18*s1l1
      t28 = s2l1*t11
      t29 = t28*t3
      t31 = t11*t3
      t33 = t2*t3
      t35 = z*t3
      t37 = s1l1*t4
      t39 = s2l1*t4
      t41 = -15*t6+19*t9+5*t13-3*t16+15*t19-36*t21+18*t24+15*t26-5*t29+5
     #*t31-10*t33+5*t35-5*t37-5*t39
      t42 = t2*l1
      t44 = t28*t4
      t46 = t15*t3
      t48 = s2l1*t3
      t50 = t37*t2
      t52 = s1l1*t3
      t54 = s1l1*l1
      t56 = s2l1*l1
      t60 = s2l1*z
      t62 = t60*t3
      t64 = t60*l1
      t68 = 7*t42+5*t44+23*t46+13*t48-15*t50+18*t52-19*t54-6*t56-7*z*l1-
     #2*t60-31*t62+9*t64+2*z+6*s1l1
      t92 = 9*t6-27*t9-3*t13+18*t16-9*t19+26*t21-19*t24-9*t26+7*t29-3*t3
     #1+6*t33-3*t35+3*t37+3*t39-4*t42-3*t44-12*t23+4*t11*l1
      t107 = s1l1*t11
      t117 = -22*t46-8*t48+9*t50-11*t52+15*t54+7*t56+13*t60+23*t62+4*t28
     #-17*t64+8*t8-8*t15+16*t23*l1-4*t107*l1-8*t28*l1+3*z+4*t107*t3-7*s1
     #l1-4*t11
      t121 = -1+z
      taux4 = -128*(t1*(t41+t68)-2*t1*(t92+t117))*t121/(1+l3*t1*t121)

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoef
     #f*tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff*
     #tcoff7)
      tFjet8 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,1
     #.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoef
     #f*tcoff8)

CCC Sum

      qgRVnum2e3 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcof
     #f8
      RETURN
      END

      double precision function qgRVnum3e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = s2l1*z
      t4 = s1l1*s2l1
      t7 = l1**2
      t10 = z**2
      t12 = z*t7
      t15 = t10*t7
      t18 = z*l1
      t20 = t10*z
      t21 = s2l1*t20
      t24 = t7*l1
      t25 = t10*t24
      t26 = s2l1**2
      t29 = t20*t24
      t32 = z*t24
      t35 = s1l1*t10
      t38 = t26*t20
      t41 = s1l1*z
      t44 = s1l1*t24
      t47 = s2l1*t10
      t50 = -640*t1*l1+256*z*t4+512*t1*t7-128*t10+768*t4*t12-1280*t4*t15
     #-384*t18+384*t21*t7+256*t25*t26-128*t29*t26-128*t32*s1l1-128*t35*t
     #7-384*t38*l1+128*t41*t7+256*t44*t10-896*t47*t7
      t51 = t26*t10
      t58 = t20*t7
      t62 = t26*z
      t78 = s1l1*t20
      t82 = -768*t51*t7+256*t20-128*t47+128*t38+256*t1+640*t58+384*t38*t
     #7-512*t62*l1-128*t32*s2l1+768*t47*l1-256*t44*t21-768*t4*t18+512*t4
     #7*t44-256*t1*t44-1152*t15-256*t78*s2l1*l1
      t87 = t10*l1
      t98 = s1l1**2
      t99 = t98*z
      t106 = t98*t10
      t114 = 512*t21*t7*s1l1+768*t4*t87+768*t87-256*t32+512*t25-640*t20*
     #l1+512*t12+384*t35+128*t41+384*t99-256*t51-128*t29*s1l1-128*t41*l1
     #-512*t106*t7-128*t98*t24*z+384*t99*t7
      t145 = 256*t25*s2l1-640*t99*l1-512*t35*l1-128*t21*t24+256*t62+128*
     #t21-128*t29*t98-384*t21*l1+256*t106*l1+128*t78*l1-128*t32*t26+768*
     #t51*l1+384*t62*t7+128*t58*t98+256*t25*t98-256*t29+128*z
      taux5 = t50+t82+t114+t145
      t1 = -1+l1
      t5 = s2l1*l1
      t6 = s1l1*l1
      t8 = 1-z
      t9 = t8**2
      t13 = 2*t5
      t16 = l1**2
      taux6 = 512*t1**2*(t5+t6-s2l1+1)*t9*t8+256*t1*(8+t13-4*t6-8*l1-2*s
     #1l1*t16-2*s2l1*t16+6*s1l1)*t9-256*t1*(-6*t6+8+10*s1l1-6*l1-t13+2*s
     #2l1)*t8+256*t1*(2+4*s1l1)

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux5
      tcoff11 = taux6
      tcoff12 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*
     #tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff12)

CCC Sum

      qgRVnum3e0 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12
     #*tcoff12
      RETURN
      END

      double precision function qgRVnum3e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = s2l1**2
      t2 = t1*z
      t5 = z**2
      t6 = t5*z
      t7 = t1*t6
      t8 = l1**2
      t11 = s1l1*s2l1
      t12 = z*t8
      t15 = t5*t8
      t18 = s1l1*t6
      t22 = s2l1*t6
      t23 = t8*s1l1
      t26 = t5*l1
      t30 = t8*l1
      t31 = z*t30
      t33 = s1l1**2
      t34 = t33*t30
      t38 = t5*t30
      t43 = s1l1*t5
      t45 = t6*t30
      t48 = t6*t8
      t50 = 768*t2*l1-768*t7*t8+768*t11*t12+1280*t11*t15+512*t18*s2l1*l1
     #-1024*t22*t23-768*t11*t26-2176*t26+768*t31+256*t34+768*t33*l1-1536
     #*t38+1920*t6*l1-1536*t12-768*t43+768*t45+3456*t15-1920*t48
      t51 = s2l1*z
      t53 = s2l1*t5
      t56 = t1*t5
      t58 = t33*z
      t61 = z*l1
      t63 = s1l1*z
      t81 = t33*t5
      t84 = -768*t51-256*t53-768*t2+256*t56-768*t58-256*t7+1280*t61-768*
     #t63-256*t22+256*t1*t30-768*t33*t8-256*t1*t8-768*t6-512*z-256*t33+3
     #84*t31*s2l1+256*t45*t33+896*t22*l1-256*t81*l1
      t122 = -256*t31*t1-768*t56*l1+256*t2*t8-256*t48*t33-256*t38*t33-89
     #6*t53*l1+384*t45*s1l1+1152*t63*l1+512*t81*t8-256*t34*z+512*t58*t8-
     #768*t38*s2l1+512*t58*l1+384*t43*l1+384*t22*t30-1024*t22*t8-256*t38
     #*t1+256*t45*t1
      t131 = s1l1*t30
      t161 = 384*t31*s1l1+1152*t43*t8+768*t7*l1-768*t63*t8-768*t131*t5+1
     #920*t53*t8-384*t18*t8+768*t56*t8+512*t11*l1+512*t131*s2l1-1024*t23
     #*s2l1+1280*t51*l1-1024*t11*z-896*t51*t8+256*t5+512*t131*t22+768*t1
     #1*t61-512*t53*t131-512*t51*t131
      taux5 = t50+t84+t122+t161
      t1 = -1+l1
      t2 = 2*t1
      t5 = -4*t1
      t9 = s2l1*l1
      t10 = s1l1*l1
      t13 = 1-z
      t14 = t13**2
      t17 = l1**2
      t19 = 2*s2l1*t17
      t21 = 2*s1l1*t17
      t22 = 6*l1
      t23 = 6*s1l1
      t25 = 2*s2l1
      t26 = 2*t10
      t29 = 2*t9
      taux6 = -256*(t1*t2-t5*t2)*(t9+t10-s2l1+1)*t14*t13+256*(t1*(t19+t2
     #1+t22-t23-4*t9-6+t25+t26)+t5*(8+t23+t29-8*l1-t19-t21-4*t10))*t14-2
     #56*(t1*(4*s2l1-4+4*l1-t29-8*s1l1+t26)+t5*(-6*t10+8+10*s1l1-t22-t29
     #+t25))*t13+256*t1*(4*s2l1-4*s1l1)+256*t5*(2+4*s1l1)

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux5
      tcoff11 = taux6
      tcoff12 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*
     #tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff12)

CCC Sum

      qgRVnum3e1 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12
     #*tcoff12
      RETURN
      END

      double precision function qgRVnum3e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = l1**2
      t3 = t1*t2
      t5 = t2*l1
      t6 = t1*t5
      t9 = s1l1**2
      t10 = t9*z
      t13 = s1l1*t1
      t16 = t1*z
      t17 = s2l1*t16
      t22 = s2l1**2
      t25 = t16*t5
      t28 = t5*z
      t33 = t22*t16
      t36 = s1l1*z
      t39 = s1l1*t5
      t42 = s2l1*t1
      t45 = t1*l1
      t50 = 128-3456*t3+512*s2l1*t6+128*t10*l1-256*t13*l1-256*t17*t5+768
     #*t17*t2+384*t6*t22-128*t25*t22-256*t28*s1l1-256*t13*t2-384*t33*l1-
     #256*t36*t2+512*t39*t1-1792*t42*t2+2304*t45+128*t22-896*t25+640*t16
      t65 = t22*t1
      t68 = t22*z
      t71 = t16*t2
      t82 = t9*t5
      t89 = 256*s2l1+256*t13+768*z-256*t28*s2l1-128*t25*t9-768*t17*l1-25
     #6*t9*t1*l1-384*t28*t22+384*t65*l1+512*t68*t2+128*t71*t9+384*t6*t9+
     #1280*t42*l1-256*t25*s1l1+768*t36*l1-384*t82*z+256*t10*t2+256*s1l1-
     #512*t1
      t96 = s1l1*s2l1
      t97 = z*l1
      t102 = s2l1*z
      t105 = s1l1*t16
      t109 = t2*s1l1
      t117 = dsqrt(z)
      t118 = dsqrt(l5)
      t119 = t118*l5
      t120 = t117*t119
      t122 = dsqrt(1-l5)
      t123 = t120*t122
      t124 = dsqrt(l1)
      t125 = t124*l1
      t127 = dsqrt(1-l1)
      t129 = t125*t127*l6
      t132 = t117*t118
      t133 = t132*t122
      t135 = t124*t127*l6
      t139 = t122*t125*t127
      t144 = t117*z*t119
      t148 = l5**2
      t149 = t148*z
      t150 = l6**2
      t151 = l1*t150
      t154 = l5*z
      t155 = t2*t150
      t172 = -1+8*t123*t129+8*t133*t135-t2+4*t132*t139+2*l1+8*t144*t122*
     #t129+16*t149*t151+16*t154*t155-16*t149*t155-16*t154*t151-8*t133*t1
     #29+6*t149*l1-6*t149*t2+6*t154*t2-t148*t1*t2
      t180 = t122*t124*t127
      t192 = l1*l6
      t195 = t2*l6
      t205 = -8*t123*t135-6*t154*l1-4*t144*t139-t148-4*t132*t180-4*l5*l1
     #-t148*t2+2*t148*l1+2*l5*t2-4*t120*t139-16*t149*t192+16*t149*t195+1
     #6*t154*t192-16*t154*t195+2*l5+4*t120*t180
      t209 = (1-l1+t97)**2
      t224 = 1664*t71+256*t17+640*t68-256*t39*t17+256*t96*t97+768*t42*t3
     #9-768*t102*t39-256*t105*s2l1*l1+512*t17*t109-256*t96*t45+128*t33-1
     #536*t97-128*(-1+z)*(t172+t205)/t209+256*t96*t1-256*t105*l1-1280*t1
     #02*l1+1024*t102*t2+128*t10-768*t68*l1
      t238 = z*t2
      t258 = 384*t33*t2+512*t105*t2-768*t65*t2+256*t96*l1+256*t39*s2l1-5
     #12*t109*s2l1+128*t82+1792*t238-1408*t16*l1+256*t96+384*t9*l1-128*t
     #22*t2-896*t28+768*t102+1792*t6-384*t9*t2+768*t96*t238-768*t96*t3+1
     #28*t22*t5
      taux5 = t50+t89+t224+t258
      t1 = -1+l1
      t2 = l1**2
      t6 = -4*t1
      t7 = 2*t1
      t9 = -t6
      t13 = s2l1*l1
      t14 = s1l1*l1
      t17 = 1-z
      t18 = t17**2
      t21 = t2*s1l1
      t24 = s2l1*t2
      t26 = 2*t13
      t27 = 4*l1
      t30 = 2*t24
      t31 = 2*t21
      t32 = 6*l1
      t33 = 6*s1l1
      t35 = 2*s2l1
      t36 = 2*t14
      taux6 = -256*(3*t1*t2+t6*t7-t9*t7)*(t13+t14-s2l1+1)*t18*t17+256*(t
     #1*(11*t21+3*t2+5*t24-t26+t27)+t6*(t30+t31+t32-t33-4*t13-6+t35+t36)
     #+t9*(8+t33+t26-8*l1-t30-t31-4*t14))*t18-256*(t1*(-s2l1+t13+1+t27+1
     #1*t14)+t6*(4*s2l1-4+t27-t26-8*s1l1+t36)+t9*(-6*t14+8+10*s1l1-t32-t
     #26+t35))*t17+256*t1*(-s2l1+3*s1l1+1)+256*t6*(4*s2l1-4*s1l1)+256*t9
     #*(2+4*s1l1)

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux5
      tcoff11 = taux6
      tcoff12 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*
     #tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff12)

CCC Sum

      qgRVnum3e2 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12
     #*tcoff12
      RETURN
      END

      double precision function qgRVnum3e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z*l1
      t3 = s1l1**2
      t4 = l1**2
      t7 = t4*l1
      t8 = t3*t7
      t10 = s2l1**2
      t11 = t10*z
      t13 = s1l1*s2l1
      t15 = z**2
      t16 = t10*t15
      t19 = t15*z
      t23 = s2l1*t15
      t26 = s2l1*z
      t29 = 128+896*t1-384*t3*t4+128*t8+128*t11+256*t13-256*t16*t4-640*t
     #19*t4-128*z-256*t23*t4+256*t26*t4
      t30 = s1l1*z
      t33 = t15*t7
      t39 = t15*t4
      t42 = dsqrt(z)
      t43 = dsqrt(l5)
      t44 = t43*l5
      t45 = t42*t44
      t47 = dsqrt(1-l5)
      t48 = t45*t47
      t49 = dsqrt(l1)
      t51 = dsqrt(1-l1)
      t53 = t49*t51*l6
      t57 = t42*z*t44
      t59 = t49*l1
      t61 = t59*t51*l6
      t64 = l5*z
      t65 = l6**2
      t66 = l1*t65
      t71 = l5**2
      t76 = t42*t43
      t78 = t47*t59*t51
      t83 = t76*t47
      t95 = -1-8*t48*t53+8*t57*t47*t61-16*t64*t66+8*t48*t61-t71*t15*t4-6
     #*t64*l1+4*t76*t78+2*l5*t4+8*t83*t53-t4-4*t45*t78+2*t71*l1-t71*t4-4
     #*l5*l1-8*t83*t61
      t97 = t47*t49*t51
      t103 = t71*z
      t110 = t4*t65
      t117 = l1*l6
      t120 = t4*l6
      t130 = -4*t76*t97-4*t57*t78+2*l5-t71+6*t103*l1-6*t103*t4+6*t64*t4+
     #16*t64*t110-16*t103*t110+4*t45*t97-16*t103*t117+16*t103*t120+16*t6
     #4*t117-16*t64*t120+16*t103*t66+2*l1
      t134 = (1-l1+t1)**2
      t143 = z*t7
      t147 = -256*t30*t4+128*t33*t10+256*s1l1*t15*t4+1920*t39-128*(-1+z)
     #*(t95+t130)/t134+256*t13*l1+128*t10+256*t26-1280*t33+640*t143-512*
     #t15*l1
      t159 = t7*s1l1
      t172 = 384*t3*l1+128*t10*t7-256*t11*l1-256*t26*l1+128*t16*l1-512*t
     #26*t159-256*t13*t39+384*t11*t4-256*t143*t10+128*t33*t3+512*t30*l1
      t182 = z*t4
      t197 = 256*t23*l1+384*t3*z*t4-256*t8*z+640*t19*t7-1280*t182+768*t1
     #3*t182+256*t23*t159-128*t10*t4+256*t159*s2l1-512*s2l1*t4*s1l1+256*
     #s2l1+256*s1l1
      taux5 = t29+t147+t172+t197
      t1 = -1+l1
      t2 = l1**2
      t5 = -4*t1
      t8 = -t5
      t12 = s2l1*l1
      t13 = s1l1*l1
      t16 = 1-z
      t17 = t16**2
      t20 = t2*s1l1
      t22 = s2l1*t2
      t23 = 2*t22
      t29 = 2*t12
      t30 = 4*l1
      t38 = 2*t13
      taux6 = -256*(t1*t2+3*t5*t2+2*t8*t1)*(t12+t13-s2l1+1)*t17*t16+256*
     #(t1*(4*t20-t12+l1+t2+t23)+t5*(11*t20+3*t2+5*t22-t29+t30)+t8*(t23+2
     #*t20+6*l1-6*s1l1-4*t12-6+2*s2l1+t38))*t17-256*(t1*(-t12+3*t13+l1)+
     #t5*(-s2l1+t12+1+t30+11*t13)+t8*(4*s2l1-4+t30-t29-8*s1l1+t38))*t16+
     #256*t5*(-s2l1+3*s1l1+1)+256*t8*(4*s2l1-4*s1l1)

CCC Coefficients

      tcoff9 = taux5
      tcoff10 = taux5
      tcoff11 = taux6
      tcoff12 = taux6

CCC Jet functions

      tFjet9 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.D
     #0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff*
     #tcoff9)
      tFjet10 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff11)
      tFjet12 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff12)

CCC Sum

      qgRVnum3e3 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12
     #*tcoff12
      RETURN
      END

      double precision function qgRVnum4e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = s2l1*z
      t2 = l1**2
      t3 = t2*l1
      t4 = t3*s1l1
      t7 = s2l1**2
      t8 = t7*z
      t10 = z**2
      t11 = t10*l1
      t14 = t10*z
      t15 = t14*t2
      t16 = s1l1**2
      t21 = z*t3
      t24 = t10*t3
      t27 = s1l1*s2l1
      t30 = s2l1*t14
      t34 = s2l1*t10
      t37 = z*l1
      t45 = -128*t1*t4+128*t8+384*t11+128*t1+64*t15*t16+192*t8*t2-64*t21
     #*t7+128*t24*t16+384*t27*t11+256*t30*t2*s1l1+256*t34*t4-384*t27*t37
     #-128*t21+256*t24-320*t14*l1-64*t34
      t46 = t7*t10
      t52 = t10*t2
      t54 = t14*t3
      t56 = s1l1*t10
      t58 = s1l1*z
      t60 = z*t2
      t62 = t16*z
      t64 = t7*t14
      t74 = -128*t46+64*z-192*t37+64*t30+320*t15-576*t52-128*t54+192*t56
     #+64*t58+256*t60+192*t62+64*t64-384*t46*t2-448*t34*t2+192*t30*t2-64
     #*t30*t3
      t108 = 192*t64*t2+256*t1*t2+128*t4*t10+64*t58*t2-64*t54*t7+128*t24
     #*t7-64*t21*s1l1-64*t56*t2-192*t64*l1-64*t58*l1-64*t54*s1l1+384*t34
     #*l1-256*t56*l1-320*t62*l1+128*t24*s2l1+192*t62*t2
      t112 = t16*t10
      t123 = s1l1*t14
      t145 = -64*t16*t3*z-256*t112*t2+128*t14-64*t10-320*t1*l1+128*t27*z
     #-256*t8*l1+64*t123*l1+128*t112*l1-64*t21*s2l1+384*t46*l1-64*t54*t1
     #6-192*t30*l1-128*t4*t30-640*t27*t52+384*t27*t60-128*t123*s2l1*l1
      taux7 = t45+t74+t108+t145
      t1 = -1+l1
      t5 = s2l1*l1
      t6 = s1l1*l1
      t8 = 1-z
      t9 = t8**2
      t15 = 2*t5
      t16 = l1**2
      taux8 = 256*t1**2*(t5+t6-s2l1+1)*t9*t8+128*t1*(8+6*s1l1-4*t6+t15-2
     #*s1l1*t16-2*s2l1*t16-8*l1)*t9-128*t1*(-6*t6+8+10*s1l1-6*l1-t15+2*s
     #2l1)*t8+128*t1*(2+4*s1l1)

CCC Coefficients

      tcoff13 = taux7
      tcoff14 = taux7
      tcoff15 = taux8
      tcoff16 = taux8

CCC Jet functions

      tFjet13 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff16)

CCC Sum

      qgRVnum4e0 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFjet
     #16*tcoff16
      RETURN
      END

      double precision function qgRVnum4e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = t1*l1
      t4 = t1*z
      t5 = s2l1*t4
      t7 = s1l1**2
      t8 = l1**2
      t11 = s2l1**2
      t14 = t11*z
      t16 = z*l1
      t18 = t4*t8
      t20 = t1*t8
      t22 = t8*l1
      t23 = t4*t22
      t27 = t7*t22
      t29 = s1l1*t1
      t31 = s1l1*z
      t33 = z*t8
      t35 = t7*z
      t37 = t11*t4
      t39 = s1l1*s2l1
      t42 = t11*t1
      t45 = -1088*t2-128*t5-384*t7*t8-128*t11*t8-384*t14+640*t16-960*t18
     #+1728*t20+384*t23+128*t11*t22+128*t27-384*t29-384*t31-768*t33-384*
     #t35-128*t37+256*t39*l1+384*t42*t8
      t46 = s1l1*t4
      t49 = s2l1*t1
      t63 = s2l1*z
      t66 = t7*t1
      t73 = z*t22
      t83 = s1l1*t22
      t86 = -192*t46*t8+960*t49*t8-512*t5*t8+192*t5*t22+960*t4*l1-128*t4
     #9+384*t14*l1-512*t39*z+640*t63*l1-128*t66*l1-384*t42*l1+128*t23*t7
     #+192*t73*s2l1+448*t5*l1+128*t42-256*z-128*t7+192*t23*s1l1+256*t83*
     #s2l1
      t94 = t1*t22
      t123 = -448*t49*l1+192*t29*l1+256*t35*l1-384*t94*s2l1+256*t35*t8-1
     #28*t27*z+128*t1-384*t4+256*t66*t8-256*t63*t83+256*t83*t5+640*t39*t
     #20+384*t39*t33+256*t46*s2l1*l1+384*t73-384*t37*t8-448*t63*t8-384*t
     #83*t1
      t161 = -384*t31*t8+128*t23*t11-128*t94*t11+192*t73*s1l1+576*t29*t8
     #+384*t37*l1+576*t31*l1-384*t39*t2-512*t5*t8*s1l1-256*t49*t83+384*t
     #39*t16-768*t94-128*t18*t7+128*t14*t8-128*t73*t11-128*t94*t7-512*t3
     #9*t8-384*t63+384*t7*l1
      taux7 = t45+t86+t123+t161
      t1 = -1+l1
      t2 = 2*t1
      t5 = -4*t1
      t9 = s2l1*l1
      t10 = s1l1*l1
      t13 = 1-z
      t14 = t13**2
      t17 = 2*s2l1
      t18 = l1**2
      t20 = 2*s2l1*t18
      t21 = 6*l1
      t23 = 6*s1l1
      t25 = 2*s1l1*t18
      t26 = 2*t10
      t30 = 2*t9
      taux8 = -128*(t1*t2-t5*t2)*(t9+t10-s2l1+1)*t14*t13+128*(t1*(t17-6+
     #t20+t21-4*t9-t23+t25+t26)+t5*(8+t23-4*t10+t30-t25-t20-8*l1))*t14-1
     #28*(t1*(4*s2l1-4+4*l1-t30-8*s1l1+t26)+t5*(-6*t10+8+10*s1l1-t21-t30
     #+t17))*t13+128*t1*(4*s2l1-4*s1l1)+128*t5*(2+4*s1l1)

CCC Coefficients

      tcoff13 = taux7
      tcoff14 = taux7
      tcoff15 = taux8
      tcoff16 = taux8

CCC Jet functions

      tFjet13 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff16)

CCC Sum

      qgRVnum4e1 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFjet
     #16*tcoff16
      RETURN
      END

      double precision function qgRVnum4e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = s2l1*z
      t4 = s2l1**2
      t5 = t4*z
      t8 = s1l1**2
      t9 = t8*z
      t13 = l5**2
      t14 = t13*z
      t17 = l1**2
      t20 = l5*z
      t23 = dsqrt(z)
      t24 = dsqrt(l5)
      t25 = t23*t24
      t27 = dsqrt(1-l5)
      t28 = dsqrt(l1)
      t31 = dsqrt(1-l1)
      t32 = t27*t28*t31
      t35 = t24*l5
      t36 = t23*t35
      t39 = z**2
      t43 = t23*z*t35
      t44 = t28*l1
      t46 = t27*t44*t31
      t55 = t44*t31*l6
      t62 = t36*t27
      t64 = t28*t31*l6
      t70 = -1+6*t14*l1-6*t14*t17+6*t20*t17-4*t25*t32+4*t36*t32-t13*t39*
     #t17-4*t43*t46-6*t20*l1-4*t36*t46+8*t43*t27*t55+2*t13*l1+2*l5*t17-8
     #*t62*t64+4*t25*t46-t13*t17
      t74 = t17*l6
      t77 = l1*l6
      t82 = l6**2
      t83 = l1*t82
      t86 = t17*t82
      t95 = t25*t27
      t103 = -4*l5*l1+2*l1+16*t14*t74+16*t20*t77-16*t14*t77-16*t20*t83-1
     #6*t14*t86+16*t20*t86+16*t14*t83-16*t20*t74-t17+8*t95*t64-8*t95*t55
     #-t13+2*l5+8*t62*t55
      t106 = z*l1
      t108 = (1-l1+t106)**2
      t112 = s2l1*t39
      t115 = s1l1*t39
      t118 = s1l1*s2l1
      t123 = t17*l1
      t124 = t39*t123
      t127 = z*t123
      t132 = t39*z
      t133 = t132*t17
      t136 = s2l1*t132
      t141 = t132*t123
      t144 = t4*t39
      t150 = s1l1*t132
      t153 = 64-640*t1*l1-384*t5*l1+64*t9*l1-64*(-1+z)*(t70+t103)/t108+6
     #40*t112*l1-128*t115*l1+128*t118*l1-256*t118*t17+192*t124*t8-192*t1
     #27*t4+256*t5*t17+64*t133*t8-384*t136*l1-128*t127*s2l1-64*t141*t8+1
     #92*t144*l1-128*t8*t39*l1-128*t150*l1
      t177 = s1l1*z
      t180 = t39*l1
      t187 = s1l1*t123
      t190 = 384*t136*t17-896*t112*t17+256*t150*t17-384*t144*t17-128*t11
     #5*t17-128*t127*s1l1+192*t124*t4-64*t141*t4-64*t4*t17-192*t8*t17+12
     #8*t136+128*t118*t39-128*t177*t17+1152*t180-448*t127+896*t124-704*t
     #132*l1+384*t1+256*t187*t39
      t194 = t4*t132
      t208 = t8*t123
      t215 = t39*t17
      t225 = 512*t1*t17+192*t194*t17-128*t136*t123+128*t187*s2l1+64*t194
     #+192*t8*l1-128*t141*s1l1-384*t1*t187+64*t208+64*t4*t123+384*t177*l
     #1-448*t141-1728*t215+832*t133-768*t106-192*t194*l1+128*t118+320*t5
     #+128*t118*t106
      t237 = z*t17
      t258 = 64*t9+384*t112*t187+256*t136*t17*s1l1-128*t118*t180-128*t15
     #0*s2l1*l1+384*t118*t237-384*t118*t215-192*t208*z-128*t187*t136+128
     #*t9*t17+256*t124*s2l1+384*z+896*t237+128*s2l1+128*s1l1+128*t115+64
     #*t4+320*t132-256*t39
      taux7 = t153+t190+t225+t258
      t1 = -1+l1
      t2 = l1**2
      t6 = -4*t1
      t7 = 2*t1
      t9 = -t6
      t13 = s2l1*l1
      t14 = s1l1*l1
      t17 = 1-z
      t18 = t17**2
      t22 = s2l1*t2
      t24 = 4*l1
      t25 = s1l1*t2
      t27 = 2*t13
      t30 = 2*s2l1
      t31 = 2*t22
      t32 = 6*l1
      t34 = 6*s1l1
      t35 = 2*t25
      t36 = 2*t14
      taux8 = -128*(3*t1*t2+t6*t7-t9*t7)*(t13+t14-s2l1+1)*t18*t17+128*(t
     #1*(3*t2+5*t22+t24+11*t25-t27)+t6*(t30-6+t31+t32-4*t13-t34+t35+t36)
     #+t9*(8+t34-4*t14+t27-t35-t31-8*l1))*t18-128*(t1*(-s2l1+t13+1+t24+1
     #1*t14)+t6*(4*s2l1-4+t24-t27-8*s1l1+t36)+t9*(-6*t14+8+10*s1l1-t32-t
     #27+t30))*t17+128*t1*(-s2l1+3*s1l1+1)+128*t6*(4*s2l1-4*s1l1)+128*t9
     #*(2+4*s1l1)

CCC Coefficients

      tcoff13 = taux7
      tcoff14 = taux7
      tcoff15 = taux8
      tcoff16 = taux8

CCC Jet functions

      tFjet13 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff16)

CCC Sum

      qgRVnum4e2 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFjet
     #16*tcoff16
      RETURN
      END

      double precision function qgRVnum4e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z*l1
      t3 = z**2
      t4 = t3*z
      t5 = l1**2
      t8 = t3*t5
      t10 = t5*l1
      t13 = s2l1**2
      t16 = s1l1**2
      t17 = t16*t10
      t22 = t13*z
      t25 = 64+448*t1-320*t4*t5+960*t8+320*t4*t10+64*t13*t10+64*t17+128*
     #s1l1+128*s2l1-64*z+192*t22*t5
      t26 = z*t10
      t29 = t3*t10
      t32 = s1l1*s2l1
      t37 = t13*t3
      t40 = s2l1*t3
      t43 = s2l1*z
      t46 = s1l1*z
      t55 = -128*t26*t13+64*t29*t16-256*t32*t5+128*t32*l1-128*t37*t5-128
     #*t40*t5+128*t43*t5-128*t46*t5+64*t29*t13+128*s1l1*t3*t5+64*t22
      t75 = 128*t43-640*t29+320*t26-256*t3*l1-192*t16*t5+128*t32-64*t13*
     #t5-128*t22*l1+64*t37*l1-128*t43*l1-128*t32*t8
      t76 = z*t5
      t79 = t10*s1l1
      t89 = l5**2
      t90 = dsqrt(z)
      t91 = dsqrt(l5)
      t92 = t90*t91
      t94 = dsqrt(1-l5)
      t95 = t92*t94
      t96 = dsqrt(l1)
      t97 = t96*l1
      t99 = dsqrt(1-l1)
      t101 = t97*t99*l6
      t106 = t91*l5
      t107 = t90*z*t106
      t112 = t96*t99*l6
      t115 = t90*t106
      t117 = t94*t97*t99
      t120 = l5*z
      t121 = l6**2
      t122 = t5*t121
      t125 = t89*z
      t128 = l1*t121
      t132 = t94*t96*t99
      t135 = t115*t94
      t138 = l1*l6
      t143 = -1+2*l5-t89-8*t95*t101+2*l1+8*t107*t94*t101+8*t95*t112-t5-4
     #*t115*t117+16*t120*t122-16*t125*t122-16*t120*t128+4*t115*t132+8*t1
     #35*t101-16*t125*t138+16*t120*t138
      t144 = t5*l6
      t176 = 16*t125*t144-16*t120*t144+16*t125*t128+2*l5*t5+2*t89*l1-4*l
     #5*l1-t89*t5-4*t92*t132+4*t92*t117-8*t135*t112+6*t125*l1-6*t125*t5+
     #6*t120*t5-t89*t3*t5-6*t120*l1-4*t107*t117
      t180 = (1-l1+t1)**2
      t196 = 384*t32*t76+128*t40*t79+64*t13+256*t46*l1+128*t79*s2l1-64*(
     #-1+z)*(t143+t176)/t180+128*t40*l1+192*t16*z*t5-128*z*t17-640*t76+1
     #92*t16*l1-256*t43*t79
      taux7 = t25+t55+t75+t196
      t1 = -1+l1
      t2 = l1**2
      t5 = -4*t1
      t8 = -t5
      t12 = s2l1*l1
      t13 = s1l1*l1
      t16 = 1-z
      t17 = t16**2
      t20 = s2l1*t2
      t21 = 2*t20
      t22 = s1l1*t2
      t28 = 4*l1
      t30 = 2*t12
      t38 = 2*t13
      taux8 = -128*(t1*t2+3*t5*t2+2*t8*t1)*(t12+t13-s2l1+1)*t17*t16+128*
     #(t1*(l1+t21+4*t22-t12+t2)+t5*(3*t2+5*t20+t28+11*t22-t30)+t8*(2*s2l
     #1-6+t21+6*l1-4*t12-6*s1l1+2*t22+t38))*t17-128*(t1*(-t12+3*t13+l1)+
     #t5*(-s2l1+t12+1+t28+11*t13)+t8*(4*s2l1-4+t28-t30-8*s1l1+t38))*t16+
     #128*t5*(-s2l1+3*s1l1+1)+128*t8*(4*s2l1-4*s1l1)

CCC Coefficients

      tcoff13 = taux7
      tcoff14 = taux7
      tcoff15 = taux8
      tcoff16 = taux8

CCC Jet functions

      tFjet13 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff13)
      tFjet14 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff15)
      tFjet16 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff16)

CCC Sum

      qgRVnum4e3 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFjet
     #16*tcoff16
      RETURN
      END

      double precision function qgRVnum5e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = t1*z
      t3 = t2*l1
      t4 = 5120*t3
      t5 = l1**2
      t6 = t5*l1
      t7 = t1*t6
      t9 = z*t6
      t10 = 128*t9
      t11 = t1*l1
      t13 = t1**2
      t14 = t13*z
      t15 = t14*t5
      t17 = t5**2
      t18 = t1*t17
      t19 = 256*t18
      t20 = t2*t17
      t21 = 384*t20
      t22 = z*t17
      t23 = 64*t22
      t24 = t14*t6
      t26 = t13*t17
      t27 = 256*t26
      t28 = t14*t17
      t29 = 64*t28
      t30 = t13*t5
      t32 = t13*l1
      t35 = t13*t6
      t37 = t5*t2
      t39 = t1*t5
      t41 = t2*t6
      t43 = z*t5
      t45 = -t4+512*t7-t10+4992*t11+64*t15-t19+t21+t23-128*t24-t27+t29-3
     #84*t30+128*t32+2496*t2+512*t35+3072*t37-5248*t39-768*t41+2496*t43
      t46 = s1l1**2
      t51 = 384*t15
      t58 = t14*l1
      t68 = -t4+1152*t7-256*t9+4864*t11+t51-512*t18+768*t20+128*t22-384*
     #t24-512*t26+128*t28-128*t58-1280*t30+384*t32+2432*t2+1408*t35+3840
     #*t37-5376*t39-1920*t41+2432*t43
      t74 = 6*t18
      t75 = 4*t20
      t76 = 4*t22
      t79 = z*l1
      t86 = 15*t41
      t88 = -t30+38*t79+t32+39*t2+19*t1-t35+19*t5+77*t37-132*t39-t86+37*
     #t43
      t98 = -2944*t3+640*t7-t10+2560*t11+t51-t19+t21+t23-256*t24-t27+t29
      t110 = -256*t58-128*t13+64*t14-1152*t30+640*t32+1280*t2+896*t35+24
     #32*t37-2944*t39-1152*t41+1280*t43
      t112 = s2l1**2
      t119 = -27*t3-21*t7+13*t9-9*t11+t74-t75-t76+t26+t13-3*t6+t17
      t129 = 6*t30+40*t79-4*t32+17*t2+20*t1-4*t35+20*t5-t37+4*t39+t86-29
     #*t43
      t134 = t17*l1
      t147 = t17*t5
      t153 = t13*t1
      t170 = t2-10*t1+22*t6+43*t153*t134-20*t153*t5-73*t153*t17-165*t1*t
     #147-11*t153*t147-165*t13*t147+220*t2*t147-163*z*t134
      t183 = 24*t134+435*t7-139*t9+44*t11+163*t15-587*t18+938*t20+191*t2
     #2-333*t24-827*t26+383*t28
      t198 = -357*t30+103*t32+550*t13*t134+66*z*t147+61*t153*t6+350*t37-
     #191*t39-752*t41+65*t43-10*t5-11*t147
      t203 = (1-l1+t79)**2
      t209 = (-l1-z+t79)**2
      taux9 = (t45*t46+(t68*s2l1+64*z*(-97*t3+51*t7-53*t9+58*t11+t74-t75
     #-t76+t26+18*t6+t17+t88))*s1l1+(t98+t110)*t112+64*z*(t119+t129)*s2l
     #1-64*z*(-239*t14*t134-670*t2*t134-20*t79-25*t17-87*t3+455*t1*t134+
     #706*t35-20*t13-40*t58+66*t14*t147+t170+t183+t198)/t203)/t209
      t1 = l1**2
      t2 = t1*l1
      t3 = z*t2
      t6 = z**2
      t7 = z*t6
      t8 = s2l1*t7
      t11 = t6**2
      t12 = t11*z
      t13 = s2l1*t12
      t16 = s2l1*z
      t17 = t1**2
      t20 = t6*t17
      t23 = s1l1*t12
      t29 = t12*t17
      t32 = t11*t17
      t35 = t7*t17
      t46 = t2*t6
      t51 = -768*t3*s2l1-21760*t8*l1-1024*t13*t2+256*t16*t17-1024*t20*s2
     #l1-768*t23*t2-1024*s1l1*t17*t6+256*t29*s1l1-1024*t32*s1l1+1536*t35
     #*s1l1+1536*t35*s2l1+256*t29*s2l1-1024*t32*s2l1+768*t23*t1+3328*t46
     #*s2l1-256*t23*l1
      t53 = s1l1*t7
      t58 = t11*t2
      t61 = t11*s1l1
      t64 = s1l1*z
      t67 = s2l1*t11
      t74 = t11*t1
      t79 = s2l1*t6
      t90 = -256*t11-35328*t53*l1-1024*t13*l1+7936*t58*s1l1+25856*t61*l1
     #+256*t64*t17+3840*t67*t2+48128*t53*t1+1536*t13*t1-22784*t74*s1l1-9
     #984*t67*t1-17152*t79*t1+20480*t8*t1-5376*t8*t2+5120*t16*t1+12544*t
     #67*l1
      t99 = s1l1*t6
      t102 = t7*t2
      t121 = 17664*s1l1*t2*t6+4864*t64*t1-5632*t3*s1l1-30976*t99*t1-1920
     #0*t102*s1l1+10240*t79*l1+9728*t99*l1+2560*t7-4864*t12-6144*t7*l1-7
     #936*t46-512*t3+5120*t6*l1-16896*t12*t1+9216*t20-13824*t35
      t146 = -2304*z*t17+9472*t12*t2+9216*t32-2304*t29+5120*t8+14592*t12
     #*l1+256*t13+32512*t74-13568*t11*l1-27904*t58-11776*t7*t1-6400*t1*t
     #6+26880*t102-9984*t61+2560*t1*z+4864*t53-5376*t67
      t151 = (-l1-z+z*l1)**2
      taux10 = (t51+t90+t121+t146)/t151

CCC Coefficients

      tcoff17 = taux9
      tcoff18 = taux9
      tcoff19 = taux10
      tcoff20 = taux10

CCC Jet functions

      tFjet17 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff20)

CCC Sum

      qgRVnum5e0 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet
     #20*tcoff20
      RETURN
      END

      double precision function qgRVnum5e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = l1**2
      t3 = t1*t2
      t5 = t1**2
      t6 = t5*z
      t7 = t6*t2
      t9 = t2**2
      t10 = t1*t9
      t11 = 256*t10
      t12 = t1*z
      t13 = t12*t9
      t14 = 256*t13
      t15 = z*t9
      t16 = 384*t15
      t17 = t2*l1
      t18 = t6*t17
      t20 = 128*t9
      t21 = t5*t9
      t22 = 384*t21
      t23 = t6*t9
      t24 = 128*t23
      t25 = z*l1
      t28 = 256*t17
      t29 = 7488*t3-128*t7-t11-t14+t16+256*t18-t20+t22-t24+832*t25-4992*
     #t12+t28
      t30 = t12*t17
      t32 = t5*t17
      t34 = t5*l1
      t36 = z*t2
      t38 = t12*l1
      t42 = t1*t17
      t44 = z*t17
      t46 = t1*l1
      t48 = t5*t2
      t50 = t12*t2
      t51 = 4224*t50
      t53 = 320*t30-704*t32-128*t34-3456*t36+8832*t38-128*t2+960*t1+704*
     #t42-832*t44-9536*t46+448*t48-t51-256*z
      t55 = s1l1**2
      t59 = 3072*t36
      t69 = 1792*t48+6400*t3-t59+1280*t30-2048*t32-4352*t50-8320*t46+768
     #*t42+7552*t38-3968*t12-1024*t44+512*t1
      t79 = 768*t7
      t80 = t6*l1
      t82 = 1024*t25-512*t34-256*t9+t28+768*t18+768*t15-512*t10+768*t21-
     #512*t13-256*t23-t79+256*t80
      t101 = 3*t9
      t103 = 12*t15
      t104 = 18*t10
      t105 = 3*t21
      t106 = 12*t13
      t107 = 90*t44-16*t1-37*t25+4*l1-21*t2-t101-32*t17+t103-t104-t105+t
     #106
      t122 = 1728*t48+6336*t3-t59+960*t30-1344*t32-t51-6144*t46+64*t42+6
     #592*t38-3072*t12-192*t44
      t128 = 192*t5-960*t34-t20-128*t6+512*t18+t16-t11+t22-t14-t24-t79+5
     #12*t80
      t130 = s2l1**2
      t143 = -15*t48-91*t3+98*t36-56*t30+11*t32+48*t50+51*t46+102*t42+24
     #*t38-28*t12-80*t44
      t150 = -2*t5-44*t1-84*t25+9*t34-40*t2-t101+23*t17+t103-t104-t105+t
     #106
      t161 = t5*t1
      t166 = t9*t2
      t183 = t9*l1
      t191 = 42*t161*t166+42*t166+120*t80+120*t9+630*t1*t166-840*t12*t16
     #6+630*t5*t166-2015*t5*t183+844*t6*t183+2560*t12*t183-109*t183
      t208 = 2944*t21-1276*t23+60*t25-1825*t1*t183+692*z*t183-252*t6*t16
     #6-147*t161*t183-29*t12-81*t17+2895*t30-2486*t32
      t220 = -361*t34-246*t36+387*t38+60*t5+28*t2+36*t1-1814*t42+590*t44
     #-206*t46+1262*t48-1409*t50
      t225 = (1-l1+t25)**2
      t231 = (-l1-z+t25)**2
      taux9 = ((t29+t53)*t55+((t69+t82)*s2l1+64*z*(-3*t48+246*t3-93*t36+
     #14*t30+6*t32-129*t50-148*t46-78*t42+181*t38-78*t12+t107))*s1l1+(t1
     #22+t128)*t130+64*z*(t143+t150)*s2l1-64*z*(820*t3-515*t7+2404*t10-3
     #576*t13-844*t15+1079*t18-183*t161*t17+60*t161*t2-252*z*t166+228*t1
     #61*t9+t191+t208+t220)/t225)/t231
      t1 = z**2
      t2 = s1l1*t1
      t5 = s2l1*z
      t6 = l1**2
      t7 = t6**2
      t10 = s2l1*t1
      t13 = t1**2
      t14 = t6*l1
      t15 = t13*t14
      t18 = t1*z
      t19 = t18*t14
      t22 = t1*t14
      t26 = t13*z
      t27 = t26*s2l1
      t30 = s1l1*z
      t35 = t13*t6
      t41 = z*t14
      t44 = t13*s1l1
      t47 = s1l1*t18
      t52 = s2l1*t13
      t55 = -55808*t2*l1-10496*t5*t7-51712*t10*l1-44288*t15*s1l1+114688*
     #t19*s1l1-31232*t22*s2l1-155136*t19+5120*t27*l1-10496*t30*t7-7680*t
     #27*t6+113920*t35*s1l1+24320*t26+5120*t27*t14+47360*t41*s1l1-126976
     #*t44*l1-245248*t47*t6-28160*t5*t6-21504*t52*t14
      t56 = s2l1*t18
      t59 = s1l1*t14
      t87 = 36352*t56*t14-116736*t59*t1+166656*t6*t2-34048*t30*t6+6144*t
     #30*l1-95744*t56*t6+86272*t10*t6-1024*t13+4608*t2+155648*t15-10752*
     #t18+72704*t13*l1-9728*z*t6+20736*t18*l1+51200*t22+256*t41-20480*t1
     #*l1-172032*t35
      t100 = t1*t7
      t102 = t18*t7
      t108 = t13*t7
      t110 = t26*t7
      t114 = s1l1*t7
      t120 = 73984*t18*t6+20992*t1*t6-26624*t47+24576*t52-4864*t59+2560*
     #s1l1*t6+86784*t26*t6-55296*t100+82944*t102+13824*z*t7-51968*t26*t1
     #4-55296*t108+13824*t110-72960*t26*l1+2304*t114+49920*t44-1280*t27+
     #2304*s2l1*t7
      t126 = s1l1*t26
      t157 = -2560*s2l1*t14-26112*t56+18944*t100*s2l1+3840*t126*t14+1753
     #60*t47*l1+102400*t56*l1+13824*t41*s2l1+45312*t52*t6+18944*t114*t1-
     #55808*t52*l1-16896*t102*s1l1-1280*t110*s2l1+7424*t108*s2l1-16896*t
     #102*s2l1-3840*t126*t6+1280*t126*l1-1280*t110*s1l1+7424*t108*s1l1
      t162 = (-l1-z+z*l1)**2
      taux10 = (t55+t87+t120+t157)/t162

CCC Coefficients

      tcoff17 = taux9
      tcoff18 = taux9
      tcoff19 = taux10
      tcoff20 = taux10

CCC Jet functions

      tFjet17 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff20)

CCC Sum

      qgRVnum5e1 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet
     #20*tcoff20
      RETURN
      END

      double precision function qgRVnum5e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
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
      t13 = 18432*t3*t9*t10
      t16 = 18432*t3*t9*l6
      t17 = t1*z
      t18 = t17*l1
      t22 = t1*l1
      t23 = l1**2
      t24 = z*t23
      t26 = t1*t23
      t28 = t17*t23
      t30 = t2*t23
      t32 = t23*l1
      t33 = z*t32
      t35 = t1*t32
      t37 = t17*t32
      t39 = t2*t32
      t40 = 6*t18-1-2*z+3*l1-t22+6*t24-4*t26+2*t28-t30-t1-3*t23-4*t33+6*
     #t35-4*t37+t39+t32
      t46 = l5**2
      t48 = dsqrt(z)
      t51 = dsqrt(1-l5)
      t52 = t48*t1*t51
      t53 = dsqrt(l1)
      t55 = t22-2*t5-z+l1-1
      t59 = dsqrt(1-l1)
      t60 = t59*t4
      t62 = t60*t8*l6
      t71 = dsqrt(l5)
      t86 = t48*t17*t51*t53
      t105 = t23**2
      t106 = 64*t105
      t108 = 256*l1+384*t32+3520*t24+2944*t22-512*t18-2560*t5+64*t30-115
     #2*t33+1280*t35-768*t37+384*t39-t106+512*z
      t112 = 128*t3
      t113 = t2*z
      t114 = t113*t23
      t116 = t1*t105
      t117 = 640*t116
      t118 = t17*t105
      t119 = 640*t118
      t120 = z*t105
      t121 = 320*t120
      t122 = t113*t32
      t124 = t2*t105
      t125 = 320*t124
      t126 = t113*t105
      t127 = 64*t126
      t130 = 2048*t28+64*t17-5120*t26-t112+64*t114-t117+t119+t121-128*t1
     #22-t125+t127-1216*t1-576*t23
      t132 = s1l1**2
      t144 = t113*l1
      t148 = 640*t32+4352*t24+4224*t22-2944*t18-1024*t5-768*t30-2176*t33
     #+3072*t35-2560*t37+1408*t39-128*t144-128*t105+128*t2
      t152 = 384*t114
      t161 = 3712*t28+512*t17-7168*t26-t112+t152-1280*t116+1280*t118+640
     #*t120-384*t122-640*t124+128*t126+2304*t1-512*t23
      t175 = 10*t28+44*t24-49*t26+6*t30+3*t37-6*t39-7*t18-39*t33+57*t22+
     #27*t35-52*t5
      t183 = 2*t105
      t184 = 2*t124
      t185 = 8*t118
      t186 = 8*t120
      t187 = 12*t116
      t188 = 8*z+3*t1+2*t17-2*t3-11*t23+4*l1+15*t32+t183+t184-t185-t186+
     #t187
      t206 = 64*t113+256*t32+3712*t24+5504*t22-6912*t18+1152*t5-1152*t30
     #-1024*t33+1792*t35-1792*t37+1024*t39-256*t144-t106
      t215 = 1152*z-64*t2+4864*t28+3200*t17-7808*t26+512*t3+t152-t117+t1
     #19+t121-256*t122-t125+t127+1152*t1
      t217 = s2l1**2
      t230 = -83*t28-115*t24+158*t26+12*t30+55*t37-8*t39+33*t18+101*t33-
     #121*t22-117*t35+96*t5
      t238 = 36*z+68*t1+2*t2+3*t17-8*t3+28*t23-31*t32+t183+t184-t185-t18
     #6+t187
      t243 = t105*l1
      t246 = t105*t23
      t253 = t2*t1
      t283 = 312*t113*t246-921*z*t243+2280*t1*t243-3010*t17*t243-885*t11
     #3*t243-52*t246+155*t243+111*t32+351*t24+323*t22-453*t18
      t296 = -76*t5-1122*t30-833*t33+2316*t35-3184*t37+2321*t39-82*t144-
     #184*t105-41*t2+1586*t28+23*t17
      t308 = -1115*t26+288*t3+371*t114-3005*t116+4030*t118+1167*t120-855
     #*t122-2970*t124+1139*t126-50*t1-30*t23
      t316 = (-l1-z+t5)**2
      taux9 = ((-t13+t16-1152*t40*t4*z*t8)*t46+(-9216*t52*t53*t55*t62+46
     #08*t52*t53*t55*t59*t9)*t71*l5+(t13-t16-2304*(-3*t22+t5+z-2*l1+1+t2
     #6-2*t24+t23)*t4*t1*t8)*l5+(-9216*t86*t62+4608*t86*t60*t8)*t71+(t10
     #8+t130)*t132+((t148+t161)*s2l1+64*z*(t175+t188))*s1l1+(t206+t215)*
     #t217+64*z*(t230+t238)*s2l1-64*z*(2235*t2*t243-780*t2*t246+1040*t17
     #*t246-780*t1*t246-52*t253*t246-177*t253*t105+312*z*t246-41*t253*t2
     #3+124*t253*t32+146*t253*t243+t283+t296+t308)*t8)/t316
      t1 = z**2
      t2 = t1*z
      t3 = l1**2
      t4 = t3**2
      t5 = t2*t4
      t7 = t1**2
      t8 = t7*z
      t10 = t3*l1
      t11 = s1l1*t10
      t13 = s1l1*z
      t16 = s2l1*t1
      t19 = s1l1*t1
      t22 = t1*t10
      t25 = s2l1*t8
      t28 = s2l1*z
      t31 = t1*t4
      t38 = z*t10
      t41 = s2l1*t7
      t44 = s2l1*t2
      t51 = t2*t10
      t54 = t7*t10
      t61 = 102912*t13*t3-351360*t19*t3-159744*t38*s1l1+57856*t41*t10-12
     #9664*t44*t10+322944*t11*t1-31488*t13*l1-281728*t51*s1l1+97920*t54*
     #s1l1+512*t1*t3-171520*t2*t3
      t67 = t4*l1
      t72 = s1l1*t7
      t77 = s1l1*t8
      t80 = s1l1*t4
      t85 = t4*t8
      t89 = t7*t4
      t103 = s1l1*t2
      t109 = -27776*t89*s2l1+78336*t5*s1l1+2816*t85*s1l1-26624*t89*s1l1-
     #147072*t8*t3+12160*s2l1*t10-22400*t19+47360*t103+323712*t51-307072
     #*t54-129664*t7*l1
      t138 = t7*t3
      t146 = -79872*t72-14976*s2l1*t4+384*s1l1*t67-15360*t80-37120*t41+4
     #4160*t44-12160*s1l1*t3+310272*t138-2048*t77*l1+6144*t77*t3-8192*t2
     #5*l1
      t189 = -165376*t44*l1-295680*t103*l1+12672*t25*t3-194304*t138*s1l1
     #+82304*t41*l1+1920*t41*t67-1920*t13*t67-1920*t28*t67+59264*t28*t3+
     #10880*t2+5120*t7
      t195 = (-l1-z+z*l1)**2
      taux10 = (1920*t72*t67-9344*t25*t10-103424*t80*t1-170368*t16*t3+44
     #8768*t103*t3+63104*t28*t4+79104*t5*s2l1-77184*t41*t3+3840*t19*t67-
     #384*t25*t67+135296*t22*s2l1-66304*t38*s2l1+384*s2l1*t67+64256*t13*
     #t4-6528*t77*t10+200960*t72*l1-3840*t103*t67+128256*t19*l1+18688*t1
     #*l1+91264*t16*l1+116736*t8*l1+7808*z*t3+3200*t85*s2l1-182016*t5-38
     #912*t8-124544*t22+27136*t11+t146+t61+t189+2048*t25-384*t77*t67+995
     #84*t8*t10+3840*t16*t67-102656*t31*s2l1-3840*t44*t67-30336*z*t4+175
     #616*t44*t3+121344*t31+8320*t38-5760*t2*l1-30336*t85+121344*t89+t10
     #9)/t195

CCC Coefficients

      tcoff17 = taux9
      tcoff18 = taux9
      tcoff19 = taux10
      tcoff20 = taux10

CCC Jet functions

      tFjet17 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff20)

CCC Sum

      qgRVnum5e2 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet
     #20*tcoff20
      RETURN
      END

      double precision function qgRVnum5e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = l1**2
      t3 = t1*t2
      t5 = t2*l1
      t6 = t1*t5
      t7 = 192*t6
      t8 = t1**2
      t9 = t2**2
      t10 = t8*t9
      t11 = 64*t10
      t12 = t1*z
      t13 = t12*t9
      t14 = 256*t13
      t15 = t1*t9
      t16 = 384*t15
      t17 = t8*t5
      t19 = z*t2
      t21 = t1*l1
      t23 = z*l1
      t25 = t12*t2
      t26 = 128*t25
      t27 = z*t5
      t30 = z*t9
      t31 = 256*t30
      t34 = 64*t9
      t35 = t12*t5
      t38 = 704*t3-t7-t11+t14-t16+64*t17-1024*t19-448*t21+704*t23-t26+32
     #0*t27+448*t2+t31-128*t5-256*l1-t34-64*t35-256*z
      t39 = s1l1**2
      t48 = t8*t2
      t57 = 1152*t3-384*t6-128*t10+512*t13-768*t15+256*t17-1408*t19-t26-
     #128*t48+896*t27+512*t2+512*t30-384*t5-128*t9-512*t1-384*t35
      t62 = 11*t21
      t71 = 17*t3-9*t6-25*t19-t62+19*t23-3*t25+9*t27+11*t2-3*t5-8*l1+3*t
     #35-8*z
      t80 = t12*l1
      t86 = t8*l1
      t90 = 64*t3-t7-t11+t14-t16+192*t17-128*t19+512*t21-576*t80+256*t25
     #-192*t48+576*t27+384*t12+64*t86+t31-256*t5-t34-320*t35
      t91 = s2l1**2
      t104 = 30*t3-21*t6-18*t19-t62+4*t23+7*t80-14*t25+21*t27+2*t2-7*t5+
     #2*t1+7*t35
      t109 = t8*z
      t116 = t9*t2
      t119 = t8*t1
      t130 = t9*l1
      t146 = -10*t86-5*t109*t2+104*t109*t5-270*t109*t9-138*z*t116+23*t11
     #9*t9+23*t119*t116+345*t1*t116-460*t12*t116+345*t8*t116-855*t8*t130
     #+309*t109*t130+1250*t12*t130-1020*t1*t130+441*z*t130-138*t109*t116
     #+106*t80+481*t3-1049*t6
      t168 = 954*t10-1576*t13+1359*t15-564*t17-184*t19-130*t21+34*t23-47
     #5*t25+441*t27-594*t30+1131*t35+168*t48+15*t2-63*t5+104*t9-46*t119*
     #t130-79*t130+19*t1+23*t116-5*t12
      t172 = (1-l1+t23)**2
      t178 = (-l1-z+t23)**2
      taux9 = (t38*t39+(t57*s2l1+64*z*t71)*s1l1+t90*t91+64*z*t104*s2l1-6
     #4*z*(t146+t168)/t172)/t178
      t1 = l1**2
      t2 = t1**2
      t3 = s1l1*t2
      t5 = t1*l1
      t6 = s1l1*t5
      t8 = z**2
      t9 = t8**2
      t10 = z*t9
      t13 = s1l1*t9
      t15 = t10*s2l1
      t17 = t2*l1
      t24 = t8*z
      t25 = s1l1*t24
      t28 = s2l1*z
      t32 = z*t5
      t35 = s2l1*t24
      t38 = s2l1*t9
      t43 = s1l1*z
      t48 = t8*t2
      t51 = s1l1*t10
      t54 = t24*t2
      t57 = t10*t2
      t60 = t24*t5
      t63 = 116992*t32*s2l1+103424*t35*l1+66688*t38*t1-44416*t38*l1+7040
     #*t43*t17+201984*t3*t8+199168*t48*s2l1+4480*t51*t5-150272*t54*s2l1-
     #5248*t57*s2l1+340352*t60*s1l1
      t65 = s1l1*t8
      t71 = t8*t5
      t86 = t9*t2
      t106 = 50688*t86*s2l1-428416*t6*t8-147456*t54*s1l1+378880*t65*t1-3
     #840*t57*s1l1+46464*t86*s1l1+1024*t51*l1+29056*t57-116224*t86+35712
     #*t65+6784*t1*z
      t116 = t9*t5
      t127 = s2l1*t8
      t136 = t1*t9
      t149 = -14080*t127*t17+7040*t28*t17-156544*t43*t1+4096*t15*l1-2272
     #00*t136+13056*t8*l1-65920*t8*t1+19584*t38+198912*t24*t1-77568*t10*
     #t5-18816*s2l1*t5
      t189 = -106240*t116*s1l1-63104*t127*l1-116224*t48-7040*t38*t17-128
     #256*t43*t2+238336*t32*s1l1-7552*t15*t1+128384*t136*s1l1+19456*t10+
     #87424*t10*t1-23168*t35
      t195 = (-l1-z+z*l1)**2
      taux10 = (174336*t54+t189-58368*t10*l1-153984*t65*l1+53504*t43*l1-
     #58752*t28*t1-85504*t38*t5-366464*t25*t1-14080*t65*t17+1408*t15*t17
     #+225408*t35*t5-304896*t60+t63+200960*t25*l1+18816*s1l1*t1+169088*t
     #127*t1-169472*t35*t1+91392*t9*l1+29056*z*t2+14080*t25*t17-46080*t2
     #4*l1-7040*t13*t17+149760*t71-124032*t28*t2+t106+31104*t3-48512*t6-
     #4736*t9-101504*t13*l1-1408*s2l1*t17+39936*t13-1024*t15+256768*t116
     #-1408*s1l1*t17+29696*s2l1*t2-3072*t51*t1+14080*t35*t17+1408*t51*t1
     #7+t149-246400*t71*s2l1+6272*t24-29056*t25+8320*t15*t5-24064*t32)/t
     #195

CCC Coefficients

      tcoff17 = taux9
      tcoff18 = taux9
      tcoff19 = taux10
      tcoff20 = taux10

CCC Jet functions

      tFjet17 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff17)
      tFjet18 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff19)
      tFjet20 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff20)

CCC Sum

      qgRVnum5e3 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet
     #20*tcoff20
      RETURN
      END

      double precision function qgRVnum6e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = t1*l1
      t3 = z*l1
      t4 = -z-t2+t3+t1
      t5 = 768*t4
      t6 = s1l1**2
      t11 = l1**2
      t13 = z*t11
      t15 = t1*t11
      t25 = s2l1**2
      t32 = t1**2
      t33 = t32*z
      t34 = t11**2
      t35 = t34*l1
      t45 = t1*z
      t46 = t11*l1
      t57 = -1-t33*t35+t35+3*l1+z-2*t1-3*t34+3*t33*t34-4*t11+16*t13-30*t
     #15-40*t45*t46+20*t32*t46+4*t32*l1-12*t45*l1+40*t1*t46
      t88 = -20*z*t46-14*t32*t11+30*t45*t11+2*t33*t11-30*t1*t34+30*t45*t
     #34+15*z*t34-4*t33*t46-15*t32*t34-5*z*t35+10*t35*t1+12*t2-7*t3+5*t3
     #2*t35+2*t45+4*t46-10*t45*t35
      t92 = (1-l1+t3)**2
      t100 = s13e0((1-l1)*(1-z))
      taux11 = -(t5*t6+(s2l1*t5+192*z*(-2+6*t3-6*t2+2*t11-4*t13+2*t15-2*
     #z+4*t1))*s1l1+384*t4*t25-192*z*(2+2*t2-2*l1-2*t1)*s2l1+192*z*(t57+
     #t88)/t92)*t100
      t1 = s2l1*z
      t3 = z**2
      t4 = t3*z
      t5 = l1**2
      t6 = t5*l1
      t9 = t3**2
      t30 = s2l1*t4
      t32 = s2l1*t3
      t35 = s1l1*t3
      t38 = -4*t1-6*t4*t6+2*t9*t6+8*t9*l1+2*z*t5-12*t4*l1+6*t3*t6-2*z*t6
     #+2*t3*l1-6*t9*t5+14*t4*t5-10*t3*t5-4*t30-8*t32*l1-20*t35*l1
      t41 = s1l1*t4
      t48 = s1l1*z
      t65 = 4*t1*l1+12*t41*l1+4*t30*l1+8*t35*t5-4*t48*t5+8*t48*l1-4*t41*
     #t5-4*t48+2*z*l1+8*t32+4*t4+12*t35-2*z-4*t9+2*t3-8*t41
      t70 = s13e0((1-l1)*(1-z))
      taux12 = -384*(t38+t65)*t70

CCC Coefficients

      tcoff21 = taux11
      tcoff22 = taux11
      tcoff23 = taux12
      tcoff24 = taux12

CCC Jet functions

      tFjet21 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff24)

CCC Sum

      qgRVnum6e0 = tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet
     #24*tcoff24
      RETURN
      END

      double precision function qgRVnum6e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = t1*l1
      t8 = z*l1
      t9 = -z+t8-t2+t1
      t11 = s2l1**2
      t13 = l1**2
      t14 = t13**2
      t15 = t14*l1
      t16 = z*t15
      t19 = t1**2
      t20 = t19*t15
      t22 = t19*z
      t23 = t22*t15
      t24 = t1*z
      t25 = t24*t15
      t27 = t1*t15
      t29 = 12*t2
      t31 = t13*l1
      t35 = t24*t31
      t37 = t19*t31
      t39 = t19*l1
      t41 = z*t13
      t43 = -1-5*t16-2*t1+5*t20-t23-10*t25+10*t27+t29-7*t8+4*t31-3*t14+2
     #*t24-40*t35+20*t37+4*t39+16*t41
      t44 = t24*l1
      t46 = t1*t31
      t48 = z*t31
      t50 = t19*t13
      t53 = t24*t13
      t55 = t1*t13
      t57 = z*t14
      t59 = t22*t31
      t61 = t22*t14
      t63 = t22*t13
      t65 = t1*t14
      t67 = t24*t14
      t69 = t19*t14
      t71 = 3*l1
      t72 = -12*t44+40*t46-20*t48-14*t50-4*t13+30*t53+t15-30*t55+15*t57-
     #4*t59+3*t61+2*t63-30*t65+30*t67-15*t69+z+t71
      t76 = (1-l1+t8)**2
      t77 = 1/t76
      t80 = 768*t9
      t81 = s1l1**2
      t85 = 2*z
      t86 = 4*t1
      t97 = 1-l1
      t98 = 1-z
      t99 = t97*t98
      t100 = s13e1(t99)
      t102 = -t98
      t103 = t102*t77
      t104 = l6**2
      t107 = 6144*t39*t103*t104
      t110 = 6144*t39*t103*l6
      t119 = -6*t44+1+t85-t71+t2-6*t41+4*t55-2*t53+t50+t1+3*t13+4*t48-6*
     #t46+4*t35-t37-t31
      t125 = l5**2
      t127 = dsqrt(z)
      t130 = dsqrt(1-l5)
      t131 = t127*t1*t130
      t132 = dsqrt(l1)
      t134 = t2-2*t8-z+l1-1
      t137 = dsqrt(t97)
      t138 = t137*t102
      t140 = t138*t77*l6
      t149 = dsqrt(l5)
      t152 = 3*t2
      t153 = 2*l1
      t154 = 2*t41
      t164 = t127*t24*t130*t132
      t172 = 384*l1
      t173 = 192*t13
      t174 = 384*t41
      t175 = 192*t55
      t227 = 4+20*t16+10*t1-20*t20+4*t23+40*t25-40*t27-46*t2+27*t8-13*t3
     #1+11*t14-8*t24+133*t35-64*t37-12*t39-56*t41
      t245 = 42*t44-135*t46+67*t48+44*t50+13*t13-100*t53-4*t15+105*t55-5
     #4*t57+12*t59-10*t61-6*t63+106*t65-104*t67+51*t69-4*z-11*l1
      t251 = s13e0(t99)
      taux11 = (-192*z*(2-2*l1-2*t1+2*t2)*s2l1+384*t9*t11+192*z*(t43+t72
     #)*t77+t80*t81+(t80*s2l1+192*z*(-2+6*t8-t85+t86+2*t13-4*t41+2*t55-6
     #*t2))*s1l1)*t100-((t107-t110-384*t119*t102*z*t77)*t125+(3072*t131*
     #t132*t134*t140-1536*t131*t132*t134*t137*t103)*t149*l5+(-t107+t110-
     #768*(t152-t8-z+t153-1-t55+t154-t13)*t102*t1*t77)*l5+(3072*t164*t14
     #0-1536*t164*t138*t77)*t149+(192-t172+t173-t174+t175+1920*t2-1152*t
     #1+1920*z-1536*t8)*t81+((-t172+384*t13-768*t41+384*t55+1536*t2-1536
     #*t1+2304*z-1152*t8)*s2l1+192*z*(-5*t13-5*t55+10*t41-10*t8+t29+8*z+
     #7-t153-8*t1))*s1l1+(t173-t174+t175+768*t2-960*t1+768*z-768*t8)*t11
     #-192*z*(-t13-t55+t154-t152+z-6+t71+t86)*s2l1+192*z*(t227+t245)*t77
     #)*t251
      t2 = z**2
      t3 = t2*z
      t4 = s1l1*t3
      t7 = t2**2
      t9 = s2l1*t2
      t10 = t9*l1
      t12 = s1l1*t2
      t13 = t12*l1
      t15 = s2l1*z
      t16 = t15*l1
      t18 = t4*l1
      t20 = s2l1*t3
      t21 = t20*l1
      t23 = l1**2
      t24 = t12*t23
      t28 = t23*l1
      t29 = t3*t28
      t31 = t7*t28
      t33 = t7*l1
      t35 = -2*z-8*t4+4*t3-4*t7-8*t10-20*t13+4*t16+12*t18+4*t21+8*t24+2*
     #t2-4*t15-6*t29+2*t31+8*t33
      t36 = z*t23
      t38 = t3*l1
      t40 = t2*t28
      t42 = z*t28
      t44 = t2*l1
      t46 = t7*t23
      t48 = t3*t23
      t50 = t2*t23
      t54 = s1l1*z
      t56 = z*l1
      t59 = t54*t23
      t61 = t54*l1
      t63 = t4*t23
      t65 = 2*t36-12*t38+6*t40-2*t42+2*t44-6*t46+14*t48-10*t50-4*t20+12*
     #t12-4*t54+2*t56+8*t9-4*t59+8*t61-4*t63
      t69 = (1-l1)*(1-z)
      t70 = s13e1(t69)
      t84 = 10*z+40*t4-20*t3+20*t7+40*t10+98*t13-21*t16-60*t18-21*t21-50
     #*t24+2*s2l1*l1
      t97 = -6*t42*s2l1-7*t2+22*t15+36*t29-12*t31-40*t33-8*t36+59*t38-36
     #*t40+12*t42-11*t44
      t109 = s1l1*t28
      t113 = 32*t46-72*t48+48*t50+20*t20-60*t12+25*t54-8*t56-39*t9+2*s2l
     #1*t28+2*t109+2*s2l1*t23
      t133 = 25*t59-38*t61+25*t63-2*s1l1-4*t9*t23+3*t20*t23-2*t20*t28-t1
     #5*t23+6*t109*t2-6*t42*s1l1+6*t40*s2l1-2*t29*s1l1
      t136 = s13e0(t69)
      taux12 = 384*(t35+t65)*t70-384*(t84+t97+t113+t133)*t136

CCC Coefficients

      tcoff21 = taux11
      tcoff22 = taux11
      tcoff23 = taux12
      tcoff24 = taux12

CCC Jet functions

      tFjet21 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff24)

CCC Sum

      qgRVnum6e1 = tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet
     #24*tcoff24
      RETURN
      END

      double precision function qgRVnum6e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = t1*l1
      t8 = z*l1
      t9 = -z+t8-t2+t1
      t11 = s2l1**2
      t13 = l1**2
      t14 = t13**2
      t15 = t14*l1
      t16 = 12*t2
      t19 = t1**2
      t20 = t13*l1
      t21 = t19*t20
      t23 = t1*z
      t27 = 3*l1
      t28 = t23*t20
      t30 = t23*t13
      t32 = t19*t13
      t34 = t1*t13
      t36 = z*t20
      t38 = z*t13
      t40 = -1+t15+t16-7*t8-2*t1+20*t21+2*t23+4*t20-4*t13+t27-40*t28+30*
     #t30-14*t32-30*t34-20*t36+16*t38
      t41 = t19*z
      t42 = t41*t14
      t44 = t19*l1
      t46 = t1*t20
      t48 = t23*l1
      t50 = z*t15
      t52 = t1*t15
      t54 = t23*t15
      t56 = t41*t15
      t57 = t19*t14
      t59 = t19*t15
      t61 = t41*t20
      t63 = z*t14
      t66 = t23*t14
      t68 = t1*t14
      t70 = t41*t13
      t72 = 3*t42+4*t44+40*t46-12*t48-5*t50+10*t52-10*t54-t56-15*t57+5*t
     #59-4*t61+15*t63+z-3*t14+30*t66-30*t68+2*t70
      t76 = (1-l1+t8)**2
      t77 = 1/t76
      t80 = 768*t9
      t81 = s1l1**2
      t84 = 6*t8
      t85 = 2*z
      t86 = 4*t1
      t87 = 2*t13
      t88 = 4*t38
      t89 = 2*t34
      t97 = 1-l1
      t98 = 1-z
      t99 = t97*t98
      t100 = s13e2(t99)
      t102 = dsqrt(z)
      t105 = dsqrt(1-l5)
      t106 = t102*t1*t105
      t107 = dsqrt(l1)
      t109 = t2-2*t8-z+l1-1
      t112 = dsqrt(t97)
      t113 = -t98
      t114 = t112*t113
      t116 = t114*t77*l6
      t121 = t113*t77
      t125 = 3072*t106*t107*t109*t116-1536*t106*t107*t109*t112*t121
      t126 = dsqrt(l5)
      t127 = t126*l5
      t129 = 3*t2
      t130 = 2*t38
      t137 = t102*t23*t105*t107
      t143 = 3072*t137*t116-1536*t137*t114*t77
      t145 = l6**2
      t148 = 6144*t44*t121*t145
      t151 = 6144*t44*t121*l6
      t160 = -6*t48+1+t85-t27+t2-6*t38+4*t34-2*t30+t32+t1+3*t13+4*t36-6*
     #t46+4*t28-t21-t20
      t162 = z*t77
      t166 = l5**2
      t168 = 384*t38
      t169 = 192*t13
      t171 = 1536*t8
      t173 = 192*t34
      t174 = 384*l1
      t179 = 960*t1
      t189 = 1536*t2
      t193 = 10*t8
      t194 = 2*l1
      t197 = 8*z
      t205 = t129-t8-z+t194-1-t34+t130-t13
      t207 = t1*t77
      t215 = 10*t1
      t227 = 4-4*t15-46*t2+27*t8+t215-64*t21-8*t23-13*t20+13*t13-11*l1+1
     #33*t28-100*t30+44*t32+105*t34+67*t36-56*t38
      t240 = 4*z
      t245 = -10*t42-12*t44-135*t46+42*t48+20*t50-40*t52+40*t54+4*t56+51
     #*t57-20*t59+12*t61-54*t63-t240+11*t14-104*t66+106*t68-6*t70
      t251 = s13e1(t99)
      t269 = 1152*l1
      t270 = 576*t13
      t271 = 576*t34
      t272 = 1152*t38
      t275 = 384*z
      t313 = -5+5*t15+57*t2-47*t8-t215+53*t21+6*t23+23*t20-23*t13+16*l1-
     #134*t28+93*t30-32*t32-129*t34-98*t36+87*t38
      t331 = 9*t42+8*t44+164*t46-34*t48-25*t50+50*t52-50*t54-5*t56-52*t5
     #7+25*t59-8*t61+73*t63+10*z-16*t14+118*t66-132*t68+4*t70
      t337 = s13e0(t99)
      s1 = (-192*z*(2-2*l1-2*t1+2*t2)*s2l1+384*t9*t11+192*z*(t40+t72)*t7
     #7+t80*t81+(t80*s2l1+192*z*(-2+t84-t85+t86+t87-t88+t89-6*t2))*s1l1)
     #*t100
      s3 = (t125*t127-192*z*(-t13-6-t34+z-t129+t130+t86+t27)*s2l1+t143*t
     #126+(t148-t151-384*t160*t113*t162)*t166+(-t168+192+t169+1920*t2-t1
     #71-1152*t1+t173-t174+1920*z)*t81+(-t168-768*t8-t179+t169+768*t2+t1
     #73+768*z)*t11+((384*t13+384*t34-768*t38-1536*t1+2304*z+t189-t174-1
     #152*t8)*s2l1+192*z*(-t193-t194-5*t13+7+10*t38+t16+t197-5*t34-8*t1)
     #)*s1l1+(-t148+t151-768*t205*t113*t207)*l5+192*z*(t227+t245)*t77)*t
     #251
      s4 = -((-t148+t151+384*t160*t113*t162)*t166-t125*t127+(t148-t151+7
     #68*t205*t113*t207)*l5-t143*t126+(-t269+t270+t271-t272-384*t2-384*t
     #1-t275+t171+576)*t81+((-t269+1152*t13+1152*t34-2304*t38-1152*t2+t2
     #75+2304*t8)*s2l1+192*z*(t89-t88+t87+t84-t240+4-6*l1))*s1l1+(t270+t
     #271-t272-t189+t179-t275+t171)*t11-192*z*(t89-t88+t87+t193-2*t2-t19
     #7+4-8*l1)*s2l1+192*z*(t313+t331)*t77)*t337
      s2 = s3+s4
      taux11 = s1+s2
      t2 = s2l1*z
      t4 = z**2
      t5 = t4**2
      t7 = s2l1*t4
      t9 = t4*z
      t10 = l1**2
      t11 = t9*t10
      t13 = t5*t10
      t15 = t4*t10
      t17 = t10*l1
      t18 = z*t17
      t20 = z*t10
      t22 = t5*l1
      t24 = t4*t17
      t26 = t9*l1
      t28 = t4*l1
      t30 = z*l1
      t32 = s1l1*t4
      t33 = t32*t10
      t35 = -2*z-4*t2-4*t5+8*t7+14*t11-6*t13-10*t15-2*t18+2*t20+8*t22+6*
     #t24-12*t26+2*t28+2*t30+8*t33
      t36 = s1l1*z
      t37 = t36*t10
      t39 = s1l1*t9
      t40 = t39*t10
      t44 = 2*t4
      t45 = s2l1*t9
      t49 = t5*t17
      t51 = t9*t17
      t53 = t7*l1
      t55 = t36*l1
      t57 = t32*l1
      t59 = t2*l1
      t61 = t39*l1
      t63 = t45*l1
      t65 = -4*t37-4*t40-4*t36+4*t9+t44-4*t45-8*t39+12*t32+2*t49-6*t51-8
     #*t53+8*t55-20*t57+4*t59+12*t61+4*t63
      t69 = (1-l1)*(1-z)
      t70 = s13e2(t69)
      t75 = s1l1*t17
      t76 = t75*t4
      t78 = t45*t17
      t86 = 10*z+22*t2+20*t5+6*t76-2*t78-39*t7-72*t11+32*t13+48*t15+12*t
     #18-8*t20
      t92 = t18*s2l1
      t99 = -40*t22-36*t24+59*t26-11*t28-8*t30-6*t92-50*t33+25*t37+25*t4
     #0+25*t36-20*t9
      t101 = s2l1*l1
      t103 = t7*t10
      t105 = t2*t10
      t106 = t45*t10
      t108 = t18*s1l1
      t110 = t51*s1l1
      t112 = t24*s2l1
      t115 = s2l1*t17
      t119 = 2*t101-4*t103-t105+3*t106-6*t108-2*t110+6*t112-2*s1l1+2*t11
     #5-7*t4+20*t45
      t133 = 40*t39-60*t32-12*t49+2*t75+2*s2l1*t10+36*t51+40*t53-38*t55+
     #98*t57-21*t59-60*t61-21*t63
      t136 = s13e1(t69)
      t149 = -15*z-41*t2-32*t5-33*t76+11*t78+60*t7+131*t11-58*t13-88*t15
     #-26*t18+15*t20
      t162 = 64*t22+78*t24-98*t26+30*t28+4*t30+33*t92+123*t33-78*t37-56*
     #t40-56*t36+32*t9-11*t101
      t174 = 30*t103-15*t105-15*t106+33*t108+11*t110-33*t112+11*s1l1-11*
     #t115+t44-32*t45-64*t39
      t189 = 96*t32+26*t49-11*t75-78*t51-11*s1l1*l1+11*s1l1*t10-70*t53+7
     #9*t55-164*t57+45*t59+96*t61+36*t63
      t192 = s13e0(t69)
      taux12 = 384*(t35+t65)*t70+384*(t86+t99+t119+t133)*t136-384*(t149+
     #t162+t174+t189)*t192

CCC Coefficients

      tcoff21 = taux11
      tcoff22 = taux11
      tcoff23 = taux12
      tcoff24 = taux12

CCC Jet functions

      tFjet21 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff24)

CCC Sum

      qgRVnum6e2 = tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet
     #24*tcoff24
      RETURN
      END

      double precision function qgRVnum6e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = t1*l1
      t8 = z*l1
      t9 = -z+t8-t2+t1
      t11 = s2l1**2
      t13 = t1*z
      t14 = t13*l1
      t16 = l1**2
      t17 = t16*l1
      t18 = t1*t17
      t20 = z*t17
      t22 = t1**2
      t23 = t22*t16
      t25 = t13*t16
      t27 = t16**2
      t28 = t22*t27
      t30 = 3*l1
      t31 = t22*z
      t32 = t31*t27
      t34 = t31*t16
      t36 = t1*t27
      t38 = t13*t27
      t40 = z*t27
      t42 = t31*t17
      t44 = t27*l1
      t45 = t22*t44
      t47 = -1+z-12*t14+40*t18-20*t20-14*t23+30*t25-15*t28+t30+3*t32+2*t
     #34-30*t36+30*t38+15*t40-4*t42+5*t45
      t48 = z*t44
      t50 = t1*t44
      t52 = t13*t44
      t54 = t31*t44
      t55 = 12*t2
      t61 = t1*t16
      t63 = t13*t17
      t65 = t22*t17
      t67 = t22*l1
      t69 = z*t16
      t72 = -5*t48+10*t50-10*t52-t54+t55-2*t1-7*t8+t44-3*t27+4*t17-4*t16
     #-30*t61-40*t63+20*t65+4*t67+16*t69+2*t13
      t76 = (1-l1+t8)**2
      t77 = 1/t76
      t80 = 768*t9
      t81 = s1l1**2
      t84 = 6*t8
      t85 = 2*z
      t86 = 4*t1
      t87 = 2*t16
      t88 = 4*t69
      t89 = 2*t61
      t97 = 1-l1
      t98 = 1-z
      t99 = t97*t98
      t100 = s13e3(t99)
      t102 = dsqrt(z)
      t105 = dsqrt(1-l5)
      t106 = t102*t1*t105
      t107 = dsqrt(l1)
      t109 = t2-2*t8-z+l1-1
      t112 = dsqrt(t97)
      t113 = -t98
      t114 = t112*t113
      t116 = t114*t77*l6
      t121 = t113*t77
      t125 = 3072*t106*t107*t109*t116-1536*t106*t107*t109*t112*t121
      t126 = dsqrt(l5)
      t127 = t126*l5
      t129 = 3*t2
      t130 = 2*t69
      t137 = t102*t13*t105*t107
      t143 = 3072*t137*t116-1536*t137*t114*t77
      t145 = l6**2
      t148 = 6144*t67*t121*t145
      t151 = 6144*t67*t121*l6
      t160 = -6*t14+1+t85-t30+t2-6*t69+4*t61-2*t25+t23+t1+3*t16+4*t20-6*
     #t18+4*t63-t65-t17
      t162 = z*t77
      t166 = l5**2
      t168 = 384*t69
      t169 = 192*t16
      t171 = 1536*t8
      t173 = 192*t61
      t174 = 384*l1
      t179 = 960*t1
      t189 = 1536*t2
      t193 = 10*t8
      t194 = 2*l1
      t197 = 8*z
      t205 = t129-t8-z+t194-1-t61+t130-t16
      t207 = t1*t77
      t212 = 4*z
      t227 = 4-t212+42*t14-135*t18+67*t20+44*t23-100*t25+51*t28-11*l1-10
     #*t32-6*t34+106*t36-104*t38-54*t40+12*t42-20*t45
      t233 = 10*t1
      t245 = 20*t48-40*t50+40*t52+4*t54-46*t2+t233+27*t8-4*t44+11*t27-13
     #*t17+13*t16+105*t61+133*t63-64*t65-12*t67-56*t69-8*t13
      t251 = s13e2(t99)
      t253 = 576*t16
      t254 = 576*t61
      t257 = 1152*t69
      t258 = 1152*l1
      t259 = 384*z
      t314 = -5+10*z-34*t14+164*t18-98*t20-32*t23+93*t25-52*t28+16*l1+9*
     #t32+4*t34-132*t36+118*t38+73*t40-8*t42+25*t45
      t331 = -25*t48+50*t50-50*t52-5*t54+57*t2-t233-47*t8+5*t44-16*t27+2
     #3*t17-23*t16-129*t61-134*t63+53*t65+8*t67+87*t69+6*t13
      t337 = s13e1(t99)
      t350 = 2-t85-108*t18+88*t20-16*t25+20*t28-12*l1-2*t32+80*t36-60*t3
     #8-50*t40-10*t45+10*t48
      t364 = -20*t50+20*t52+2*t54-10*t2+22*t8-2*t44+12*t27-26*t17+26*t16
     #+58*t61+56*t63-10*t65-68*t69
      t367 = s13e0(t99)
      s1 = (-192*z*(2-2*l1-2*t1+2*t2)*s2l1+384*t9*t11+192*z*(t47+t72)*t7
     #7+t80*t81+(t80*s2l1+192*z*(-2+t84-t85+t86+t87-t88+t89-6*t2))*s1l1)
     #*t100+(t125*t127-192*z*(-t16-6-t61+z-t129+t130+t86+t30)*s2l1+t143*
     #t126+(t148-t151-384*t160*t113*t162)*t166+(-t168+192+t169+1920*t2-t
     #171-1152*t1+t173-t174+1920*z)*t81+(-t168-768*t8-t179+t169+768*t2+t
     #173+768*z)*t11+((384*t16+384*t61-768*t69-1536*t1+2304*z+t189-t174-
     #1152*t8)*s2l1+192*z*(-t193-t194-5*t16+7+10*t69+t55+t197-5*t61-8*t1
     #))*s1l1+(-t148+t151-768*t205*t113*t207)*l5+192*z*(t227+t245)*t77)*
     #t251
      taux11 = s1+((t253+t254-384*t1-384*t2-t257-t258+t171-t259+576)*t81
     #+(-t148+t151+384*t160*t113*t162)*t166-t125*t127+(t254-t257+t179-t1
     #89+t253+t171-t259)*t11+(t148-t151+768*t205*t113*t207)*l5-192*z*(t8
     #7-2*t2-t197+t193-t88-8*l1+4+t89)*s2l1-t143*t126+((1152*t61-1152*t2
     #-t258+1152*t16-2304*t69+2304*t8+t259)*s2l1+192*z*(t87-t88-t212+4+t
     #89-6*l1+t84))*s1l1+192*z*(t314+t331)*t77)*t337-192*z*(t350+t364)*t
     #77*t367
      t1 = s1l1*z
      t2 = l1**2
      t3 = t1*t2
      t5 = z**2
      t6 = s1l1*t5
      t7 = t6*t2
      t9 = t5*z
      t10 = s1l1*t9
      t11 = t10*t2
      t13 = t1*l1
      t15 = t6*l1
      t17 = s2l1*t5
      t18 = t17*l1
      t20 = s2l1*z
      t21 = t20*l1
      t23 = t10*l1
      t25 = s2l1*t9
      t26 = t25*l1
      t29 = t2*l1
      t30 = t9*t29
      t32 = t5**2
      t33 = t32*l1
      t35 = z*t2
      t37 = t9*l1
      t39 = t5*t29
      t41 = -4*t3+8*t7-4*t11+8*t13-20*t15-8*t18+4*t21+12*t23+4*t26-2*z-6
     #*t30+8*t33+2*t35-12*t37+6*t39
      t42 = z*t29
      t44 = t5*l1
      t49 = t32*t2
      t51 = t9*t2
      t55 = z*l1
      t58 = t5*t2
      t60 = 2*t5
      t63 = t32*t29
      t65 = -2*t42+2*t44-4*t1+8*t17-4*t20-6*t49+14*t51-4*t25-8*t10+2*t55
     #+12*t6-10*t58+t60+4*t9-4*t32+2*t63
      t69 = (1-l1)*(1-z)
      t70 = s13e3(t69)
      t72 = s2l1*t29
      t74 = s1l1*t29
      t85 = 2*t72+2*t74+25*t3-50*t7+25*t11-38*t13+98*t15+40*t18-21*t21-6
     #0*t23-21*t26
      t88 = s2l1*t2
      t98 = 10*z+36*t30+2*t88-40*t33-8*t35+59*t37-36*t39+12*t42-11*t44+2
     #5*t1-39*t17
      t106 = s2l1*l1
      t110 = t42*s1l1
      t112 = t25*t29
      t114 = 22*t20+32*t49-72*t51+20*t25+40*t10-8*t55+2*t106-60*t6+48*t5
     #8-6*t110-2*t112
      t115 = t74*t5
      t117 = t20*t2
      t118 = t17*t2
      t120 = t25*t2
      t122 = t30*s1l1
      t124 = t39*s2l1
      t126 = t42*s2l1
      t133 = 6*t115-t117-4*t118+3*t120-2*t122+6*t124-6*t126-7*t5-20*t9-2
     #*s1l1+20*t32-12*t63
      t136 = s13e2(t69)
      t149 = 20*t72+20*t74+140*t3-148*t7+52*t11-116*t13+120*t15+56*t18-5
     #6*t21-48*t23-20*t26
      t162 = 4*z+72*t30-24*t88-32*t33-28*t35+68*t37-72*t39+24*t42-52*t44
     #+52*t1-28*t17+28*t20
      t170 = s1l1*l1
      t172 = s1l1*t2
      t178 = 40*t49-108*t51+16*t25+32*t10+16*t55+20*t106+44*t170-44*t172
     #-48*t6+96*t58-60*t110-20*t112
      t191 = 60*t115+72*t117-72*t118+24*t120-20*t122+60*t124-60*t126+12*
     #t5-16*t9-20*s1l1+16*t32-24*t63
      t194 = s13e0(t69)
      t207 = -11*t72-11*t74-78*t3+123*t7-56*t11+79*t13-164*t15-70*t18+45
     #*t21+96*t23+36*t26
      t220 = -15*z-78*t30+64*t33+15*t35-98*t37+78*t39-26*t42+30*t44-56*t
     #1+60*t17-41*t20-58*t49
      t233 = 131*t51-32*t25-64*t10+4*t55-11*t106-11*t170+11*t172+96*t6-8
     #8*t58+33*t110+11*t112
      t245 = -33*t115-15*t117+30*t118-15*t120+11*t122-33*t124+33*t126+t6
     #0+32*t9+11*s1l1-32*t32+26*t63
      t248 = s13e1(t69)
      taux12 = 384*(t41+t65)*t70+384*(t85+t98+t114+t133)*t136-384*(t149+
     #t162+t178+t191)*t194+384*(t207+t220+t233+t245)*t248

CCC Coefficients

      tcoff21 = taux11
      tcoff22 = taux11
      tcoff23 = taux12
      tcoff24 = taux12

CCC Jet functions

      tFjet21 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff21)
      tFjet22 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff23)
      tFjet24 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff24)

CCC Sum

      qgRVnum6e3 = tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet
     #24*tcoff24
      RETURN
      END

      double precision function qgRVnum7e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = l1**2
      t3 = t1*t2
      t4 = 512*t3
      t5 = t1*z
      t6 = t2*l1
      t7 = t5*t6
      t9 = z*t2
      t10 = 256*t9
      t11 = t5*l1
      t13 = t1*t6
      t15 = z*t6
      t16 = 512*t15
      t17 = t1*l1
      t18 = 128*t17
      t19 = t1**2
      t20 = t19*t2
      t22 = t5*t2
      t26 = t19*t6
      t27 = 1152*t26
      t34 = t19*z
      t35 = t34*t6
      t37 = -1408*t3-2944*t7+t27+512*t9-128*t11+2944*t13-1024*t15+t18-25
     #6*t20+1152*t22-128*t35
      t38 = t2**2
      t39 = t1*t38
      t41 = t5*t38
      t43 = z*t38
      t45 = t34*t38
      t47 = t19*t38
      t49 = z*l1
      t50 = 128*t49
      t51 = t38*l1
      t52 = t34*t51
      t54 = t5*t51
      t56 = t1*t51
      t58 = z*t51
      t60 = t19*t51
      t63 = -2176*t39+2688*t41+640*t43+256*t45-1408*t47+t50-128*t52-768*
     #t54+512*t56-128*t58+512*t60-128*z
      t66 = 64*t35
      t67 = -t4-704*t7+t10-64*t11+1152*t13-t16+t18+64*t20+192*t22+(t37+t
     #63)*s2l1+t66
      t70 = 320*t43
      t72 = 64*t49
      t73 = 64*t1
      t74 = 64*t52
      t75 = 384*t54
      t76 = 256*t56
      t77 = 64*t58
      t78 = 256*t60
      t79 = 64*z
      t80 = -960*t39+960*t41+t70-320*t47+t72-t73-t74-t75+t76-t77+t78-t79
      t90 = 192*t45
      t94 = -2752*t7-2048*t3-1344*t39+1728*t41+384*t43+t90-960*t47-192*t
     #49-t74-t75+t76-192*t1+t78
      t95 = 64*t5
      t98 = 128*t19*l1
      t106 = 192*t35
      t108 = 64*t34*t2
      t109 = -t77+t95+1280*t26+t98+640*t9-704*t11+2432*t13-768*t15+896*t
     #17-704*t20+2048*t22-t106+t108
      t116 = (t67+t80)*s1l1-1664*t7-1152*t3+(t94+t109)*s2l1-1088*t39+134
     #4*t41+t70+128*t45-704*t47-t50-t74-t75+t76
      t123 = -2304*t7-1536*t3-1152*t39+1536*t41+t70+t90-896*t47-t50-t74-
     #t75+t76-128*t1+t78
      t124 = 448*t9
      t127 = 576*t15
      t131 = -t77+t95+t27+t98+t124-576*t11+1920*t13-t127+640*t17-640*t20
     #+1664*t22-t106+t108
      t133 = s2l1**2
      t143 = -t4-1024*t7+256*t26+t10+1280*t13-t16+256*t22-1024*t39+1152*
     #t41+t70+64*t45-512*t47+t72-t74-t75+t76-t77+t78-t79
      t144 = s1l1**2
      t152 = (t123+t131)*t133-t73+t78-t77+t143*t144+640*t26+t124-192*t11
     #+1664*t13-t127+384*t17-192*t20+896*t22-t66
      t155 = (1-l1+t49)**2
      taux13 = (t116+t152)/t155
      t2 = z**2
      t3 = t2**2
      t6 = l1**2
      t7 = t6*l1
      t8 = t3*t7
      t13 = z*t2
      t16 = t2*t7
      t18 = z*t7
      t22 = s1l1*z
      t24 = s2l1*t2
      t26 = t3*t6
      t30 = s2l1*t13
      t32 = t3*z
      t37 = -4*t3*l1-14*t8+2*t2-10*z*t6+10*t13*l1-26*t16+8*t18-10*t2*l1+
     #2*t22+2*t24+10*t26-24*t13*t6-2*t30+2*t32*t7-2*t32*t6
      t38 = t6**2
      t39 = t2*t38
      t41 = t13*t38
      t45 = t32*t38
      t47 = t3*t38
      t49 = z*l1
      t51 = t13*t7
      t53 = s1l1*t2
      t68 = 8*t39-12*t41-2*z*t38-2*t45+8*t47+4*t49+30*t51-2*t53-14*t8*s1
     #l1-2*t13+26*t2*t6-2*t22*t38+8*t18*s1l1-8*t22*t6+18*t53*t6
      t75 = s2l1*t32
      t78 = s2l1*z
      t87 = s2l1*t3
      t92 = s1l1*t13
      t105 = 42*t30*t7-26*s1l1*t7*t2-2*t75*t6-14*t78*t6+2*t26*s1l1+38*t2
     #4*t6-36*t30*t6+14*t87*t6-4*t87*l1-12*t92*t6+30*t51*s1l1-34*t16*s2l
     #1+2*t53*l1-14*t24*l1+4*t75*t7
      t138 = -2*t78*t38+6*t78*l1+10*t18*s2l1-2*t92*l1+12*t30*l1+8*t39*s2
     #l1+2*s1l1*t32*t7+8*s1l1*t38*t2-12*t41*s2l1-22*t87*t7-2*t45*s1l1+8*
     #t47*s1l1-2*t45*s2l1+8*t47*s2l1-12*t41*s1l1
      t143 = (1-l1+t49)**2
      taux14 = 128*(-1+l1)*(t37+t68+t105+t138)/t143

CCC Coefficients

      tcoff25 = taux13
      tcoff26 = taux13
      tcoff27 = taux14
      tcoff28 = taux14

CCC Jet functions

      tFjet25 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff25)
      tFjet26 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff26)
      tFjet27 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff27)
      tFjet28 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff28)

CCC Sum

      qgRVnum7e0 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFjet
     #28*tcoff28
      RETURN
      END

      double precision function qgRVnum7e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = t1*z
      t3 = t2*l1
      t4 = -1+z
      t5 = z*l1
      t7 = (1-l1+t5)**2
      t8 = 1/t7
      t9 = t4*t8
      t12 = -1+l1
      t13 = t12**2
      t17 = l1**2
      t18 = t17*l1
      t19 = t1*t18
      t22 = z*t17
      t24 = t1**2
      t25 = t24*t18
      t27 = t2*t18
      t28 = 3328*t27
      t29 = t1*t17
      t32 = t17**2
      t33 = t24*t32
      t34 = 1472*t33
      t35 = t24*z
      t36 = t35*t32
      t38 = z*t32
      t40 = t2*t32
      t41 = 3264*t40
      t42 = -128-4352*t19+384*t3-1920*t22-768*t25+t28+3456*t29+1088*t5+t
     #34-192*t36-896*t38-t41
      t43 = t1*t32
      t45 = t2*t17
      t47 = t24*t17
      t49 = t1*l1
      t51 = z*t18
      t53 = t32*l1
      t54 = t24*t53
      t55 = 768*t54
      t56 = z*t53
      t57 = 192*t56
      t58 = t1*t53
      t59 = 768*t58
      t60 = t2*t53
      t61 = 1152*t60
      t62 = t35*t53
      t63 = 192*t62
      t71 = 128*t5
      t73 = 512*t36
      t77 = t35*t18
      t80 = -128+640*t19+1664*t22-1664*t25+2432*t27-128*t29-t71+2048*t33
     #-t73+2560*t38-2048*t40-1024*t43+256*t77-896*t45
      t82 = 256*t49
      t85 = 512*t58
      t88 = 256*z
      t95 = 384*t47-t82-3200*t51-t55-768*t56+t85+512*t60+256*t62-t88+256
     #*l1+256*t1-1024*t17+1536*t18-1024*t32+256*t53
      t98 = 2880*t43-1600*t45+64*t47-1664*t49+1792*t51-t55+t57-t59+t61+t
     #63-384*z+448*t1+(t80+t95)*s2l1
      t115 = -64+64*t3*t9-64*t13*t4*t8+(t42+t98)*s1l1-2368*t19-256*t3-51
     #2*t22-1216*t25+2624*t27+1088*t29+t71+t34-320*t36-512*t38-2496*t40+
     #1856*t43+192*t77-832*t45+256*t47+64*t49
      t132 = -64-192*t19+64*t3+448*t22-1728*t25+2112*t27+64*t29-64*t5+13
     #44*t33-384*t36+960*t38-1344*t40-192*t43+384*t77-1216*t45+960*t47
      t136 = 128*t35*t17
      t137 = t24*l1
      t139 = 384*t54
      t140 = 384*t56
      t141 = 256*t58
      t142 = 256*t60
      t143 = 128*t62
      t144 = 64*z
      t145 = 192*t1
      t150 = 128*t53
      t151 = 192*t49-960*t51-t136-192*t137-t139-t140+t141+t142+t143-t144
     #-t145-128*t17+128*t2+384*t18-384*t32+t150
      t153 = s2l1**2
      t158 = dsqrt(z)
      t161 = 2*t5
      t163 = t158*z*(t49+3*z-t161+l1-1)
      t164 = dsqrt(l1)
      t166 = dsqrt(1-l5)
      t167 = t164*t166
      t169 = -t12
      t170 = dsqrt(t169)
      t171 = t170*t4
      t172 = t8*l6
      t173 = t171*t172
      t182 = dsqrt(l5)
      t183 = t182*l5
      t185 = t169*t4
      t186 = t1*t8
      t187 = l6**2
      t190 = 1024*t185*t186*t187
      t193 = 1024*t185*t186*l6
      t203 = 4*z+t161-l1+2*t17+14*t29-6*t45+4*t51-6*t19+4*t27-t25-10*t22
     #-9*t49-t18
      t209 = l5**2
      t223 = -192+1088*t19+2176*t22-320*t25+576*t27-768*t29-1024*t5+704*
     #t33-128*t36+1600*t38-704*t40-832*t43-128*t45
      t229 = t82-2624*t51-t139-t140+t141+t142+t143+t88+640*l1-1280*t17+1
     #280*t18-640*t32+t150
      t231 = s1l1**2
      t233 = t4*t158
      t234 = t233*t166
      t236 = t164*t170*t169
      t248 = t158*t1*t164*t166
      t256 = 3*t5
      t273 = -128-2304*t19-512*t3-192*t22-2048*t25+t28+576*t29+2112*t33-
     #t73-576*t38-t41+2240*t43+448*t77
      t282 = -1024*t45+768*t47+512*t49+576*t51-t136-64*t137-t55+t57-t59+
     #t61+t63-128*z-384*t1+320*t2
      t286 = t4*z*l1
      t287 = t12*t8
      t290 = 1024*t286*t287*t187
      t293 = 1024*t286*t287*l6
      t304 = t5-1+l1
      t322 = 768*t51+(t132+t151)*t153-512*t54+128*t56-t85+768*t60+t143+(
     #-512*t163*t167*t173+256*t163*t164*t166*t170*t9)*t183-t144+(-t190+t
     #193-64*t203*t4*z*t8)*t209+(t223+t229)*t231-t145+(512*t234*t236*t17
     #2-256*t234*t236*t8)*t182+192*t2+(512*t248*t173-256*t248*t171*t8)*t
     #182+(t190-t193+128*(-t29-t256+2*t22-t17-l1+2)*t4*t186)*l5+(t273+t2
     #82)*s2l1+(-t290+t293+64*(6*t5-6*t22-t29-1+2*l1-t17)*t4*t8)*t209+(5
     #12*t233*t167*t304*t170*t172-256*t234*t164*t304*t170*t8)*t183+(t290
     #-t293+128*t12*(t256-1+l1)*t9)*l5
      taux13 = (t115+t322)*t8
      t2 = z**2
      t3 = t2**2
      t4 = l1**2
      t5 = t4*l1
      t6 = t3*t5
      t9 = z*t5
      t12 = s1l1*z
      t15 = s1l1*t2
      t18 = t2*z
      t19 = s2l1*t18
      t22 = s1l1*t5
      t25 = t3*z
      t26 = s2l1*t25
      t29 = s2l1*z
      t32 = t3*t4
      t35 = s2l1*t2
      t40 = s1l1*t18
      t43 = t18*t5
      t46 = t2*t5
      t52 = 68*t6*s1l1-20*t9*s1l1+26*t12*t4-82*t15*t4-206*t19*t5+106*t22
     #*t2+10*t26*t4+66*t29*t4-12*t32*s1l1-188*t35*t4+180*t19*t4+62*t40*t
     #4-138*t43*s1l1+158*t46*s2l1-4*t15*l1+8*t18
      t53 = t4**2
      t57 = s2l1*t3
      t74 = t2*t53
      t80 = s1l1*t53
      t83 = t18*t53
      t88 = t25*t53
      t91 = t3*t53
      t94 = 2*t12*t53+74*t6+22*t57*l1+68*t35*l1-20*t26*t5+2*t29*t53-30*t
     #29*l1-38*t9*s2l1+6*t40*l1-60*t19*l1-28*t74*s2l1-10*s1l1*t25*t5-28*
     #t80*t2+52*t83*s2l1+110*t57*t5+10*t88*s1l1-38*t91*s1l1
      t123 = 10*t88*s2l1-38*t91*s2l1+52*t83*s1l1+52*z*t4-44*t18*l1+142*t
     #46-44*t9+46*t2*l1-8*t12-8*t35-48*t32+118*t18*t4+8*t19-2*s1l1*l1+2*
     #s2l1*t4+6*s1l1*t4+2*t80
      t134 = z*l1
      t150 = -10*t25*t5+10*t25*t4-48*t74+72*t83+12*z*t53+12*t88-48*t91-2
     #0*t134-6*t22-4*s2l1*t5+18*t3*l1+2*s2l1*t53-8*t2+8*t15-162*t43-132*
     #t2*t4-70*t57*t4
      t155 = (1-l1+t134)**2
      taux14 = 128*(-1+l1)*(t52+t94+t123+t150)/t155

CCC Coefficients

      tcoff25 = taux13
      tcoff26 = taux13
      tcoff27 = taux14
      tcoff28 = taux14

CCC Jet functions

      tFjet25 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff25)
      tFjet26 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff26)
      tFjet27 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff27)
      tFjet28 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff28)

CCC Sum

      qgRVnum7e1 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFjet
     #28*tcoff28
      RETURN
      END

      double precision function qgRVnum7e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = dsqrt(z)
      t4 = dsqrt(l1)
      t7 = dsqrt(1-l5)
      t8 = t2*t1*t4*t7
      t9 = 1-l1
      t10 = dsqrt(t9)
      t11 = -1+z
      t12 = t10*t11
      t13 = z*l1
      t15 = (1-l1+t13)**2
      t16 = 1/t15
      t17 = t16*l6
      t18 = t12*t17
      t25 = dsqrt(l5)
      t27 = l1**2
      t28 = t1*t27
      t30 = t1*z
      t31 = t27*l1
      t32 = t30*t31
      t34 = t1**2
      t35 = t34*t31
      t37 = t34*l1
      t39 = z*t27
      t41 = t30*l1
      t42 = 256*t41
      t43 = t1*t31
      t45 = z*t31
      t47 = t1*l1
      t48 = 128*t47
      t49 = t34*t27
      t51 = t30*t27
      t53 = t34*z
      t54 = t53*t31
      t56 = t53*t27
      t58 = 256-384*t28-1984*t32+1408*t35+128*t37+192*t39+t42+1280*t43-3
     #20*t45-t48-576*t49+640*t51-384*t54+128*t56
      t59 = t27**2
      t60 = t1*t59
      t62 = t30*t59
      t64 = z*t59
      t66 = t53*t59
      t68 = t34*t59
      t71 = t59*l1
      t72 = t53*t71
      t73 = 128*t72
      t74 = t30*t71
      t75 = 768*t74
      t76 = t1*t71
      t77 = 512*t76
      t78 = z*t71
      t79 = 128*t78
      t80 = t34*t71
      t81 = 512*t80
      t85 = -1344*t60+2112*t62+320*t64+384*t66-1472*t68-64*t13-t73-t75+t
     #77-t79+t81+320*t1+256*z-128*t30
      t92 = 64*t37
      t99 = 192*t54
      t100 = 64*t56
      t102 = 128-192*t27-896*t28-1600*t32+768*t35+t92+832*t39+128*t41+21
     #12*t43-1536*t45-320*t49+512*t51-t99+t100-1920*t60
      t107 = 192*t13
      t108 = 64*t72
      t109 = 640*t74
      t110 = 640*t76
      t111 = 320*t78
      t112 = 320*t80
      t114 = 128*z
      t115 = 64*t30
      t118 = 64*t71
      t119 = 1664*t62+1216*t64+192*t66-832*t68-t107-t108-t109+t110-t111+
     #t112+192*t1+t114-t115+448*t31-320*t59+t118
      t121 = s2l1**2
      t136 = 256*t66
      t137 = -192*t43+192*t45-320*t47-128*t49-320*t51-t99+t100-832*t60+1
     #344*t62+192*t64+t136
      t140 = -t9
      t141 = t140**2
      t147 = t140*t11
      t148 = t1*t16
      t149 = l6**2
      t152 = 1024*t147*t148*t149
      t155 = 1024*t147*t148*l6
      t157 = 2*t13
      t158 = 2*t27
      t166 = -4*z-t157+l1-t158-14*t28+6*t51-4*t45+6*t43-4*t32+t35+10*t39
     #+9*t47+t31
      t172 = l5**2
      t177 = t2*z*(t47+3*z-t157+l1-1)
      t178 = t4*t7
      t184 = t11*t16
      t189 = t25*l5
      t191 = 3*t13
      t203 = 1792*t39
      t210 = -64-1408*t27-256*t28+192*t32-64*t35+t203+960*t43-2240*t45+t
     #48-128*t51-1408*t60+640*t62+1344*t64
      t214 = 64*z
      t218 = 64*t66-192*t68-576*t13-t108-t109+t110-t111+t112+t214+832*l1
     #+1152*t31-448*t59+t118
      t220 = s1l1**2
      t223 = t11*z*l1
      t224 = t140*t16
      t227 = 2048*t223*t224*t149
      t230 = 2048*t223*t224*l6
      t247 = t11*t2
      t248 = t247*t7
      t250 = t4*t10*t9
      t262 = t13-1+l1
      t284 = 128*t54
      t287 = 256-1280*t27-768*t32+384*t35+1664*t39+t42+2432*t43-3456*t45
     #-256*t47-384*t51-t284-3328*t60+2304*t62
      t298 = 2560*t64+t136-1024*t68-t73-1280*t74+1280*t76-640*t78+640*t8
     #0+t114+384*l1+1536*t31-768*t59+128*t71
      t302 = t42-2240*t43+1024*t13+t214+(t287+t298)*s2l1-t73+t81-t75+t77
     #-t79-t284-64*t64
      t312 = 256+576*t62-704*t68+t136-64*t60-t203-832*t51+2624*t28+896*t
     #32+320*t35-704*t47+1152*t45
      t315 = -t115+(-t152+t155-64*t166*t11*z*t16)*t172+(512*t177*t178*t1
     #8-256*t177*t4*t7*t10*t184)*t189+(t152-t155+128*(t28+t191-2*t39+t27
     #+l1-2)*t11*t148)*l5+(t210+t218)*t220+(t227-t230+64*(-12*t13+12*t39
     #+2*t28+2-4*l1+t158)*t11*t16)*t172+(-t227+t230-256*t140*(t191-1+l1)
     #*t184)*l5+(-1024*t248*t250*t17+512*t248*t250*t16)*t25-64*t41*t184+
     #(-1024*t247*t178*t262*t10*t17+512*t248*t4*t262*t10*t16)*t189+(t302
     #+t312)*s1l1
      taux13 = (128+(-512*t8*t18+256*t8*t12*t16)*t25+(t58+t85)*s2l1+(t10
     #2+t119)*t121+832*t28-384*t32+576*t35+t92-448*t39+t42+t137-960*t68+
     #t107-t73-t75+t77-t79+t81+128*t141*t11*t16+128*t1+t114+t315)*t16
      t2 = z**2
      t3 = t2*z
      t4 = s2l1*t3
      t6 = l1**2
      t7 = t6**2
      t8 = t7*l1
      t13 = t2**2
      t14 = t13*t7
      t16 = s2l1*z
      t19 = s2l1*t13
      t20 = t6*l1
      t25 = s1l1*t20
      t28 = s1l1*t3
      t31 = t13*z
      t32 = t31*s2l1
      t36 = s1l1*t2
      t39 = s1l1*z
      t44 = z*t20
      t49 = t3*t20
      t51 = t13*t20
      t59 = t2*t20
      t61 = -30*t36*t8-14*t39*t6+102*t36*t6+30*t44*s1l1+224*t2*t6+285*t4
     #9-127*t51-28*t13*l1-91*z*t6+64*t3*l1-253*t59
      t66 = t13*t6
      t74 = t2*t7
      t76 = t3*t7
      t84 = t31*t7
      t88 = s2l1*t2
      t103 = -3*s1l1*t8-23*t84-18*s1l1*t6+10*t88-12*t36+49*t16*l1+20*t66
     #*s1l1-16*t32*t6-15*t19*t8+32*t32*t20+106*t19*t6
      t123 = z*l1
      t140 = s1l1*t31
      t148 = -2*t16-12*t3-9*s2l1*t6+77*t14*s2l1-112*t76*s1l1-16*t84*s1l1
     #+68*t14*s1l1-90*t28*t6+3*t140*t8-15*s1l1*t13*t8+295*t88*t6
      t178 = s1l1*t7
      t189 = -103*t88*l1+6*t36*l1-260*t59*s2l1+10*t39+16*t140*t20+88*t17
     #8*t2-118*t76*s2l1-19*t84*s2l1+4*t178+8*t25+9*s2l1*t20
      t195 = (1-l1+t123)**2
      taux14 = 128*(-1+l1)*(37*t123+10*t2-2*z-3*s2l1*t8+t148+t189-12*t4+
     #74*t66+92*t14+92*t74-138*t76-73*t2*l1-23*t16*t7-36*t19*l1+s2l1*t7+
     #15*t16*t8-2*s1l1+t103-30*t88*t8-154*t25*t2+30*t28*t8+88*t4*l1-191*
     #t3*t6+342*t4*t20-32*t39*t7+30*t4*t8+11*s1l1*l1+2*s2l1*l1+3*t32*t8-
     #181*t19*t20+16*t31*t20-279*t4*t6-9*t39*l1+58*t44*s2l1-102*t51*s1l1
     #-97*t16*t6+82*t74*s2l1+202*t49*s1l1-8*t28*l1-16*t31*t6+15*t39*t8+7
     #9*t44-23*z*t7+t61)/t195

CCC Coefficients

      tcoff25 = taux13
      tcoff26 = taux13
      tcoff27 = taux14
      tcoff28 = taux14

CCC Jet functions

      tFjet25 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff25)
      tFjet26 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff26)
      tFjet27 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff27)
      tFjet28 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff28)

CCC Sum

      qgRVnum7e2 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFjet
     #28*tcoff28
      RETURN
      END

      double precision function qgRVnum7e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = t1*z
      t3 = l1**2
      t4 = t2*t3
      t6 = t1*t3
      t8 = t3*l1
      t9 = t2*t8
      t11 = t1**2
      t12 = t11*t8
      t14 = z*t3
      t15 = 384*t14
      t16 = t1*t8
      t18 = z*t8
      t19 = 768*t18
      t20 = t1*l1
      t23 = t3**2
      t24 = t23*l1
      t25 = 64*t24
      t29 = z*t24
      t30 = 256*t29
      t31 = t11*t24
      t32 = 64*t31
      t34 = t2*t24
      t35 = 256*t34
      t36 = t1*t24
      t37 = 384*t36
      t38 = z*t23
      t42 = t2*t23
      t44 = t1*t23
      t45 = 960*t44
      t46 = -1280*t18+t25+640*t16+640*t14+640*t8-t30+t32-640*t3-t35+t37+
     #960*t38-320*t23+320*l1+320*t42-t45
      t47 = s1l1**2
      t52 = 768*t36
      t54 = 64-128*t4+512*t6+1536*t9-384*t12-t15-1920*t16+t19+128*t20+t4
     #6*t47+192*t29-768*t31+1152*t34-t52-640*t38
      t57 = -1+z
      t59 = t57*z*l1
      t60 = -1+l1
      t61 = z*l1
      t63 = (1-l1+t61)**2
      t64 = 1/t63
      t65 = t60*t64
      t66 = l6**2
      t69 = 1024*t59*t65*t66
      t72 = 1024*t59*t65*l6
      t81 = l5**2
      t84 = dsqrt(z)
      t85 = t57*t84
      t87 = dsqrt(1-l5)
      t88 = dsqrt(l1)
      t91 = t61-1+l1
      t92 = -t60
      t93 = dsqrt(t92)
      t95 = t64*l6
      t99 = t85*t87
      t106 = dsqrt(l5)
      t118 = t88*t93*t92
      t141 = t11*t23
      t142 = 128*t141
      t146 = -1792*t18+128*t24-256*t9+1280*t16+512*t14+768*t8-512*t29+12
     #8*t31-512*t3-512*t34+t52+1664*t38-512*t23+128*l1-t142+896*t42-1920
     #*t44+128+256*t61
      t148 = 512*t9
      t149 = 768*t6
      t150 = 1024*t16
      t152 = 128*t38
      t153 = 384*t42
      t154 = 384*t44
      t156 = 512*t18+t146*s2l1+t148+t149-t150-768*t14-t152+t142-t153+t15
     #4+128+512*t61
      t159 = t11*z
      t170 = 128*t61
      t172 = -384*t18-640*t9-t149+896*t16+t15+128*t12+t152+384*t20-t142+
     #t153-t154+128*z+128-t170+384*t4
      t174 = t60**2
      t178 = 64*z
      t183 = -t19+t25-t148-576*t6+t150+448*t14+192*t8-t30+t32-64*t3-t35
      t190 = 64+t37+64*t12+704*t38+192*t20-192*t23-t142+576*t42-t45+t178
     #-t170+192*t4
      t192 = s2l1**2
      t194 = -2496*t42+2112*t44+(t69-t72+64*(-6*t61+6*t14+t6+1-2*l1+t3)*
     #t57*t64)*t81+64*t61+(-512*t85*t87*t88*t91*t93*t95+256*t99*t88*t91*
     #t93*t64)*t106*l5+(-t69+t72-128*t60*(3*t61-1+l1)*t57*t64)*l5+(-512*
     #t99*t118*t95+256*t99*t118*t64)*t106+t156*s1l1+1216*t141-192*t159*t
     #23+192*t159*t24+t172*s2l1+64*t174*t57*t64+t178+(t183+t190)*t192
      taux13 = (t54+t194)*t64
      t2 = s2l1*z
      t3 = l1**2
      t4 = t3**2
      t7 = z**2
      t8 = t7*t4
      t11 = t7**2
      t12 = t11*z
      t13 = s1l1*t12
      t14 = t3*l1
      t17 = t12*t4
      t20 = t11*t4
      t23 = t7*z
      t24 = t23*t4
      t31 = s2l1*t23
      t32 = t4*l1
      t42 = z*l1
      t48 = s1l1*t23
      t51 = s1l1*t7
      t54 = t12*s2l1
      t60 = s1l1*t4
      t64 = -24*s2l1*t4+4*s1l1*t3-34*t42+11*s2l1*t32-8*t12*t14-110*t48*t
     #32+110*t51*t32-11*t54*t32+55*s1l1*t11*t32-35*t60+11*s1l1*t32
      t66 = t11*t3
      t70 = s1l1*z
      t74 = s2l1*t7
      t83 = s2l1*t11
      t86 = z*t14
      t90 = s1l1*t14
      t110 = 164*t90*t7+4*t70*t3-34*t23*l1+55*t83*t32+110*t74*t32+13*t17
     #+13*z*t4+78*t24-52*t8+188*t24*s1l1-242*t60*t7
      t127 = t23*t14
      t130 = t11*t14
      t134 = t7*t14
      t151 = 70*t130+60*z*t3+146*t134-46*t86+50*t7*l1+108*t23*t3+18*t11*
     #l1-7*s2l1*l1-26*t51*t3+23*t70*l1+52*t74*l1
      t190 = -19*s1l1*l1+32*t90+7*s1l1+10*t31+12*t48*l1+148*t70*t4+14*t8
     #3*l1+8*t54*t3+10*t51+8*s2l1*t14+12*s2l1*t3
      t196 = (1-l1+t42)**2
      taux14 = 128*(-1+l1)*(-52*t20+t110+10*t23+110*t83*t14+210*t24*s2l1
     #-116*t127*s1l1+226*t134*s2l1-140*t7*t3-55*t70*t32+50*t2*t3-116*t86
     #*s1l1-16*t54*t14+8*t17*s1l1-16*t51*l1+44*t130*s1l1+115*t2*t4-67*t2
     #0*s1l1-42*t83*t3-28*t31*l1-162*t127+t64+t151-8*t13*t14+t190-11*t13
     #*t32-31*t2*l1-86*t86*s2l1+24*t48*t3+8*t12*t3-55*t2*t32+19*t17*s2l1
     #-242*t31*t14-110*t31*t32-100*t20*s2l1-6*t66*s1l1+7*t2-220*t8*s2l1-
     #162*t74*t3+134*t31*t3+7*z-36*t66-4*t70-4*t7-4*t74)/t196

CCC Coefficients

      tcoff25 = taux13
      tcoff26 = taux13
      tcoff27 = taux14
      tcoff28 = taux14

CCC Jet functions

      tFjet25 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff25)
      tFjet26 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff26)
      tFjet27 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff27)
      tFjet28 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff28)

CCC Sum

      qgRVnum7e3 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFjet
     #28*tcoff28
      RETURN
      END

      double precision function qgRVnum8e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = 64*z
      t2 = z**2
      t3 = t2*z
      t4 = l1**2
      t5 = t4*l1
      t6 = t3*t5
      t8 = -1+z
      t10 = t8*z*l1
      t11 = -1+l1
      t12 = z*l1
      t14 = (1-l1+t12)**2
      t15 = 1/t14
      t16 = t11*t15
      t17 = l6**2
      t20 = 1024*t10*t16*t17
      t22 = z*t4
      t24 = t2*t4
      t32 = 1024*t10*t16*l6
      t34 = l5**2
      t36 = t11*t2
      t37 = t8*t15
      t40 = 1024*t36*t37*t17
      t43 = 1024*t36*t37*l6
      t44 = 3*t12
      t53 = 64*t3
      t54 = dsqrt(z)
      t55 = t8*t54
      t57 = dsqrt(1-l5)
      t58 = t55*t57
      t59 = dsqrt(l1)
      t60 = -t11
      t61 = dsqrt(t60)
      t63 = t59*t61*t60
      t64 = t15*l6
      t72 = dsqrt(l5)
      t75 = t2*l1
      t77 = 2*t12
      t79 = t54*z*(t75+3*z-t77+l1-1)
      t80 = t59*t57
      t82 = t61*t8
      t83 = t82*t64
      t92 = t72*l5
      t94 = t2**2
      t95 = t94*t5
      t97 = 256*t22
      t98 = t2*t5
      t100 = z*t5
      t102 = 384*t75
      t103 = t3*t4
      t107 = t3*l1
      t109 = -64-t1-576*t6+(-t20+64*(6*t12-6*t22-t24-1+2*l1-t4)*t8*t15+t
     #32)*t34+(-t40+t43-128*(t24+t44-2*t22+t4+l1-2)*t2*t37)*l5-128*t2+t5
     #3+(512*t58*t63*t64-256*t58*t63*t15)*t72+(-512*t79*t80*t83+256*t79*
     #t59*t57*t61*t37)*t92+192*t95+t97+576*t98-192*t100+t102+704*t103+12
     #8*t94*l1-384*t107
      t116 = 128*z
      t119 = 320*t24
      t120 = 64*t6
      t122 = 128*t98
      t123 = 64*t100
      t125 = 192*t2
      t129 = -128+(-128-256*t24+t97+256*t75-256*t12-t116)*s2l1+t119+t120
     #-64*t22-t122+t123-t102-256*t103+t125+192*t107+64*t12-192*z
      t131 = 192*t24
      t132 = 192*t22
      t134 = 192*t12
      t136 = s2l1**2
      t138 = t4**2
      t142 = t12-1+l1
      t163 = s1l1**2
      t168 = t54*t2*t59*t57
      t176 = t11**2
      t195 = -4*z-t77+l1-2*t4-14*t24+6*t103-4*t100+6*t98-4*t6+t95+10*t22
     #+9*t75+t5
      t205 = -128-t119-t120+t97+t122-t123+t102+64*t103-320*t2-64*t107+t5
     #3-t134-t116
      t209 = -256*t94*t4-128*t12+t129*s1l1+(-64-t131+t132+320*t75-t125-t
     #134-t1)*t136-192*t2*t138+(512*t55*t80*t142*t61*t64-256*t58*t59*t14
     #2*t61*t15)*t92+(t20-t32+128*t11*(t44-1+l1)*t37)*l5+(-64-t131+t132+
     #192*t75-320*t12+t116)*t163-704*t24+(512*t168*t83-256*t168*t82*t15)
     #*t72-64*t176*t8*t15+192*t3*t138+64*z*t138-64*t94*t138+(t40-t43+64*
     #t195*t8*z*t15)*t34+t205*s2l1+64*t107*t37
      t213 = s13e0(-l1*t8)
      taux15 = -(t109+t209)*t213
      t2 = z**2
      t3 = t2*z
      t4 = s1l1*t3
      t7 = l1**2
      t8 = t7*l1
      t15 = s1l1*t2
      t17 = t2**2
      t24 = s1l1*z
      t26 = s2l1*t3
      t28 = s2l1*t2
      t35 = -6*t4*l1+2*z*t8-6*t2*l1-12*t3*t7-2*t15-4*t17*l1+10*t3*l1+4*t
     #17*t7+2*t24-2*t26+2*t28+12*t2*t7-6*t2*t8+2*t2
      t53 = s2l1*z
      t64 = -4*t15*t7+2*t24*t7+2*t4*t7-8*t24*l1-4*z*t7-2*t3-2*t26*t7+4*t
     #28*t7-2*t17*t8-2*t53*l1+14*t15*l1+2*t28*l1-2*t53*t7+6*t3*t8
      t69 = s13e0(l1*(1-z))
      taux16 = -128*(1-l1)*(t35+t64)*t69

CCC Coefficients

      tcoff29 = taux15
      tcoff30 = taux15
      tcoff31 = taux16
      tcoff32 = taux16

CCC Jet functions

      tFjet29 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff29)
      tFjet30 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff30)
      tFjet31 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff32)

CCC Sum

      qgRVnum8e0 = tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFjet
     #32*tcoff32
      RETURN
      END

      double precision function qgRVnum8e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = l1**2
      t3 = t1*t2
      t6 = t1*z
      t7 = 64*t6
      t8 = 64*z
      t9 = 320*t3
      t10 = t2*l1
      t11 = t6*t10
      t12 = 64*t11
      t13 = t2*z
      t14 = 256*t13
      t15 = t1*t10
      t16 = 128*t15
      t17 = z*t10
      t18 = 64*t17
      t19 = t1*l1
      t20 = 384*t19
      t21 = t6*t2
      t24 = t6*l1
      t26 = z*l1
      t27 = 192*t26
      t28 = 128*z
      t29 = -128-t9-t12+t14+t16-t18+t20+64*t21-320*t1-64*t24+t7-t27-t28
      t31 = -1+z
      t32 = dsqrt(z)
      t33 = t31*t32
      t35 = dsqrt(1-l5)
      t36 = dsqrt(l1)
      t37 = t35*t36
      t39 = t26-1+l1
      t40 = 1-l1
      t41 = dsqrt(t40)
      t44 = (1-l1+t26)**2
      t45 = 1/t44
      t46 = t45*l6
      t48 = t33*t37*t39*t41*t46
      t50 = t33*t35
      t54 = t50*t36*t39*t41*t45
      t57 = dsqrt(l5)
      t58 = t57*l5
      t60 = t2**2
      t61 = z*t60
      t63 = t1**2
      t64 = t63*t60
      t66 = -t40
      t67 = t66*t1
      t68 = t31*t45
      t69 = l6**2
      t70 = t68*t69
      t72 = 1024*t67*t70
      t73 = t68*l6
      t75 = 1024*t67*t73
      t77 = 2*t26
      t78 = 2*t2
      t84 = t63*t10
      t87 = -4*z-t77+l1-t78-14*t3+6*t21-4*t17+6*t15-4*t11+t84+10*t13+9*t
     #19+t10
      t89 = z*t45
      t93 = l5**2
      t96 = t31*z*l1
      t97 = t66*t45
      t99 = t96*t97*t69
      t100 = 1024*t99
      t101 = 6*t26
      t103 = 2*l1
      t109 = t96*t97*l6
      t110 = 1024*t109
      t114 = 192*t3
      t115 = 192*t13
      t117 = 320*t26
      t119 = s1l1**2
      t122 = t36*t41*t40
      t124 = t50*t122*t46
      t127 = t50*t122*t45
      t133 = -64-704*t3-128*t1+t7-t8+t29*s2l1+(512*t48-256*t54)*t58+64*t
     #61-64*t64+(t72-t75+64*t87*t31*t89)*t93+(-t100+64*(t101-6*t13-t3-1+
     #t103-t2)*t31*t45+t110)*t93+192*t84+(-64-t114+t115+192*t19-t117+t28
     #)*t119+(512*t124-256*t127)*t57+t14+576*t15-192*t17
      t135 = t63*l1
      t137 = 384*t24
      t138 = t63*t2
      t140 = 128*t26
      t141 = 256*t3
      t142 = 256*t19
      t143 = 256*t26
      t147 = 256*t21
      t148 = 192*t1
      t151 = 192*z
      t152 = -128+(-128-t141+t14+t142-t143-t28)*s2l1+t9+t12-64*t13-t16+t
     #18-t20-t147+t148+192*t24+64*t26-t151
      t154 = 320*t19
      t156 = s2l1**2
      t158 = t1*t60
      t162 = t32*t1*t36*t35
      t163 = t41*t31
      t164 = t163*t46
      t165 = t162*t164
      t168 = t162*t163*t45
      t175 = t32*z*(t19+3*z-t77+l1-1)
      t177 = t175*t37*t164
      t182 = t175*t36*t35*t41*t68
      t186 = 3*t26
      t189 = t66*(t186-1+l1)*t68
      t193 = t6*t60
      t195 = t24*t68
      t205 = t66**2
      t207 = t205*t31*t45
      t209 = t20+704*t21+128*t135-t137-256*t138-t140+t152*s1l1+(-64-t114
     #+t115+t154-t148-t27-t8)*t156-192*t158+(512*t165-256*t168)*t57+(-51
     #2*t177+256*t182)*t58+(t100-t110+128*t189)*l5+192*t193+64*t195+(-t7
     #2+t75-128*(t3+t186-2*t13+t2+l1-2)*t1*t68)*l5-576*t11-64*t207
      t212 = -l1*t31
      t213 = s13e1(t212)
      t216 = 128*t6
      t224 = 2*t40*t1
      t226 = 1024*t224*t70
      t228 = 1024*t224*t73
      t229 = 2*t3
      t240 = 256*t17
      t247 = -t147+832*t3-512*t19+t240+256-576*t13-512*t15+t117+128*t24+
     #256*t11+448*t1+256*z-t216
      t255 = 128+2048*t3+t148-t216+t28+(-1024*t124+512*t127)*t57-256*t61
     #+256*t64+(-t226+t228-128*(-t229-t101+4*t13-t78-t103+4)*t1*t68)*l5-
     #640*t84+t247*s2l1-768*t13-1920*t15+640*t17-960*t19-2048*t21-384*t1
     #35
      t258 = 384*t26
      t272 = 8*z+4*t26-t103+4*t2+28*t3-12*t21+8*t17-12*t15+8*t11-2*t84-2
     #0*t13-18*t19-2*t10
      t279 = 256*t1
      t281 = 128*t10
      t286 = 256*t15
      t295 = 512*t21-t9+t154-128*t17+256-t115+t286+t258-t137-128*t11-t27
     #9+(-t142-512*t17+256+512*t13+t286+t140+256*l1-512*t2+256*t10+t28)*
     #s2l1+t151
      t297 = 2048*t99
      t298 = 2048*t109
      t332 = 1024*t24+768*t138+t258+768*t158+(t226-t228+64*t272*t31*t89)
     #*t93+(128*t3-t20-t240+128+t16+t143+t279-128*t2+t281+t28)*t156+t295
     #*s1l1+(-t297+t298-256*t189)*l5+(-1024*t48+512*t54)*t58-768*t193+(1
     #024*t177-512*t182)*t58-128*t195+1920*t11+128*t207+(-1024*t165+512*
     #t168)*t57+(t141-t20-t240+128*t13+t16+t143+384*l1-384*t2+t281-t28)*
     #t119+(t297-t298+64*(-12*t26+12*t13+t229+2-4*l1+t78)*t31*t45)*t93
      t334 = s13e0(t212)
      taux15 = (t133+t209)*t213-(t255+t332)*t334
      t1 = 1-l1
      t2 = z**2
      t3 = s2l1*t2
      t4 = l1**2
      t5 = t3*t4
      t7 = t2*z
      t8 = s2l1*t7
      t9 = t8*t4
      t11 = z*t4
      t13 = t2**2
      t14 = t4*l1
      t15 = t13*t14
      t17 = t7*t14
      t19 = s1l1*z
      t20 = t19*l1
      t22 = s1l1*t7
      t23 = t22*t4
      t25 = t2*t4
      t29 = t19*t4
      t31 = s1l1*t2
      t32 = t31*t4
      t34 = t22*l1
      t36 = s2l1*z
      t37 = t36*t4
      t39 = 4*t5-2*t9-4*t11-2*t15+6*t17-8*t20+2*t23+12*t25+2*t3-2*t8+2*t
     #29-4*t32-6*t34-2*t37
      t40 = t3*l1
      t43 = t13*t4
      t45 = t31*l1
      t47 = t7*l1
      t49 = t13*l1
      t52 = t36*l1
      t54 = t7*t4
      t56 = t2*l1
      t58 = z*t14
      t60 = t2*t14
      t64 = 2*t40+2*t19+4*t43+14*t45+10*t47-4*t49-2*t31-2*t52-12*t54-6*t
     #56+2*t58-6*t60+2*t2-2*t7
      t65 = t39+t64
      t68 = l1*(1-z)
      t69 = s13e1(t68)
      t90 = -12*t5+6*t9+12*t11+8*t15-24*t17+22*t20-8*t23-36*t25-6*t3+6*t
     #8-12*t29+18*t32+18*t34+6*t37-8*t40-6*t19-12*t43-40*t45-28*t47
      t100 = s1l1*t14
      t120 = 12*t49+6*t31+6*t52+36*t54+16*t56-8*t58+24*t60-2*s2l1*t14-2*
     #t100+2*s1l1*t4+6*t58*s1l1-6*t2+2*t17*s1l1-6*t100*t2+2*t8*t14-6*t60
     #*s2l1+6*t7+6*t58*s2l1+2*t8*l1
      t127 = s13e0(t68)
      taux16 = 128*t1*t65*t69-128*(t1*(t90+t120)-2*t1*t65)*t127

CCC Coefficients

      tcoff29 = taux15
      tcoff30 = taux15
      tcoff31 = taux16
      tcoff32 = taux16

CCC Jet functions

      tFjet29 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff29)
      tFjet30 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff30)
      tFjet31 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff32)

CCC Sum

      qgRVnum8e1 = tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFjet
     #32*tcoff32
      RETURN
      END

      double precision function qgRVnum8e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = t1*z
      t3 = l1**2
      t4 = t2*t3
      t6 = z*l1
      t7 = 128*t6
      t8 = t1*t3
      t10 = t3*l1
      t11 = t2*t10
      t13 = t1**2
      t14 = t13*t10
      t16 = z*t3
      t17 = 256*t16
      t18 = t1*t10
      t20 = t2*l1
      t21 = 384*t20
      t22 = t13*t3
      t24 = dsqrt(z)
      t26 = t1*l1
      t28 = 2*t6
      t30 = t24*z*(t26+3*z-t28+l1-1)
      t31 = dsqrt(l1)
      t33 = dsqrt(1-l5)
      t34 = t31*t33
      t36 = 1-l1
      t37 = dsqrt(t36)
      t38 = -1+z
      t39 = t37*t38
      t41 = (1-l1+t6)**2
      t42 = 1/t41
      t43 = t42*l6
      t44 = t39*t43
      t45 = t30*t34*t44
      t49 = t38*t42
      t51 = t30*t31*t33*t37*t49
      t54 = dsqrt(l5)
      t55 = t54*l5
      t56 = (-512*t45+256*t51)*t55
      t57 = 128*t18
      t58 = z*t10
      t59 = 64*t58
      t60 = 384*t26
      t63 = 192*t6
      t65 = 64*t11
      t66 = 64*t2
      t67 = 320*t8
      t68 = 128*z
      t69 = -128+t17+t57-t59+t60+64*t4-320*t1-t63-64*t20-t65+t66-t67-t68
      t71 = 64*z
      t72 = t38*t24
      t73 = t72*t33
      t75 = t31*t37*t36
      t77 = t73*t75*t43
      t80 = t73*t75*t42
      t86 = t24*t1*t31*t33
      t87 = t86*t44
      t90 = t86*t39*t42
      t93 = (512*t87-256*t90)*t54
      t94 = 128*t1
      t95 = -64+704*t4-t7-704*t8-576*t11+192*t14+t17+576*t18-t21-256*t22
     #+t56+t69*s2l1-t71+(512*t77-256*t80)*t54+t66+t93-t94
      t96 = 192*t16
      t97 = 320*t26
      t98 = 192*t8
      t99 = 192*t1
      t101 = s2l1**2
      t104 = 320*t6
      t106 = s1l1**2
      t109 = t38*z*l1
      t110 = -t36
      t111 = t110*t42
      t112 = l6**2
      t114 = t109*t111*t112
      t115 = 1024*t114
      t116 = 6*t6
      t118 = 2*l1
      t124 = t109*t111*l6
      t125 = 1024*t124
      t127 = l5**2
      t130 = 256*t4
      t132 = 192*z
      t133 = 256*t26
      t134 = 256*t6
      t135 = 256*t8
      t138 = 64*t6
      t139 = -128+t65-64*t16-t57+t67-t60-t130+t99+t59+192*t20-t132+(-128
     #+t133-t134-t135+t17-t68)*s2l1+t138
      t142 = t3**2
      t143 = t1*t142
      t145 = 3*t6
      t148 = t110*(t145-1+l1)*t49
      t152 = t13*l1
      t154 = t20*t49
      t155 = 64*t154
      t156 = t110*t1
      t157 = t49*t112
      t159 = 1024*t156*t157
      t160 = t49*l6
      t162 = 1024*t156*t160
      t164 = 2*t3
      t172 = -4*z-t28+l1-t164-14*t8+6*t4-4*t58+6*t18-4*t11+t14+10*t16+9*
     #t26+t10
      t174 = z*t42
      t178 = (t159-t162+64*t172*t38*t174)*t127
      t179 = t110**2
      t181 = t179*t38*t42
      t189 = (-t159+t162-128*(t8+t145-2*t16+t3+l1-2)*t1*t49)*l5
      t190 = z*t142
      t192 = t2*t142
      t194 = t13*t142
      t197 = t6-1+l1
      t200 = t72*t34*t197*t37*t43
      t205 = t73*t31*t197*t37*t42
      t209 = (-64+t96+t97-t63-t98-t99-t71)*t101+(-64+192*t26+t68-t98+t96
     #-t104)*t106+t60+(-t115+64*(t116-6*t16-t8-1+t118-t3)*t38*t42+t125)*
     #t127+t139*s1l1-192*t58-192*t143+(t115+128*t148-t125)*l5+128*t152+t
     #155+t178-64*t181+t189+64*t190+192*t192-64*t194+(512*t200-256*t205)
     #*t55
      t212 = -l1*t38
      t213 = s13e2(t212)
      t217 = 128*t10
      t218 = 256*t1
      t219 = 256*t58
      t223 = 384*t6
      t225 = 1920*t11
      t227 = 768*t16
      t231 = 2048*t114
      t232 = 2048*t124
      t235 = 2*t8
      t246 = 128*t2
      t247 = 256*z
      t249 = 512*t26
      t253 = 256*t11
      t254 = 512*t18
      t255 = t247+t219-t130+448*t1-t249+128*t20+832*t8-576*t16+t253+256-
     #t254+t104-t246
      t257 = 128+(t134+128*t8-128*t3+t217+t218-t219+t68-t60+t57+128)*t10
     #1-2048*t4+t223+2048*t8+t225-640*t14-t227-1920*t18+1024*t20+768*t22
     #+(t231-t232+64*(-12*t6+12*t16+t235+2-4*l1+t164)*t38*t42)*t127+(-t2
     #31+t232-256*t148)*l5+t68-t246+t99+t255*s2l1
      t264 = 2*t36*t1
      t266 = 1024*t264*t157
      t268 = 1024*t264*t160
      t281 = 8*z+4*t6-t118+4*t3+28*t8-12*t4+8*t58-12*t18+8*t11-2*t14-20*
     #t16-18*t26-2*t10
      t291 = 128*t58
      t292 = 256*t18
      t293 = 512*t16
      t300 = -t67+512*t4-128*t11+256-t21-t291+t223-t218-t96+(t292+t68+t2
     #93-t133+256-512*t3+256*t10+t7+256*l1-512*t58)*s2l1+t292+t132+t97
      t332 = (-1024*t77+512*t80)*t54-960*t26+(t266-t268+64*t281*t38*t174
     #)*t127+640*t58+768*t143+t300*s1l1-384*t152+(-1024*t87+512*t90)*t54
     #-128*t154+128*t181+(-t266+t268-128*(-t235-t116+4*t16-t164-t118+4)*
     #t1*t49)*l5-256*t190+(-t68-384*t3+128*t16+t57+t134-t219+t217-t60+38
     #4*l1+t135)*t106-768*t192+(-1024*t200+512*t205)*t55+256*t194+(1024*
     #t45-512*t51)*t55
      t334 = s13e1(t212)
      t337 = 448*t6
      t346 = 64*t18
      t349 = 64*t10
      t350 = t17-64*t8+64*t26+t138+t346-t291-192*t3+192*l1+t349-t68-64
      t352 = 1472*t4-t337-1856*t8-t225+576*t14+896*t16+2112*t18-512*t20-
     #512*t22+t56+t350*t106+t93-t99
      t376 = (t133+t134+(t293-t135+t133+t57-t219-256*t3+128*l1+t217-t94-
     #t247)*s2l1-t94-t247)*s1l1+(-t253+t227-1024*t8+768*t26-512*t6+t254-
     #t219+t130-384*t1)*s2l1+896*t26-768*t58-960*t143+256*t152+t155+t178
     #+t189+320*t190+960*t192+(t293-448*t8+t249-t337+t346-t291-64*t3+t34
     #9-t99)*t101-320*t194
      t378 = s13e0(t212)
      taux15 = (t95+t209)*t213+(t257+t332)*t334-(t352+t376)*t378
      t1 = 1-l1
      t2 = s2l1*z
      t3 = t2*l1
      t5 = z**2
      t6 = s2l1*t5
      t7 = l1**2
      t8 = t6*t7
      t10 = t5*z
      t11 = s2l1*t10
      t12 = t11*t7
      t14 = t7*l1
      t15 = t5*t14
      t17 = z*t7
      t19 = t5**2
      t20 = t19*t14
      t22 = s1l1*z
      t23 = t22*l1
      t25 = t10*t14
      t27 = t5*t7
      t31 = s1l1*t10
      t32 = t31*t7
      t34 = t22*t7
      t36 = s1l1*t5
      t37 = t36*t7
      t39 = -2*t3+4*t8-2*t12-6*t15-4*t17-2*t20-8*t23+6*t25+12*t27+2*t6-2
     #*t11+2*t32+2*t34-4*t37
      t40 = t31*l1
      t42 = t2*t7
      t45 = t6*l1
      t47 = t19*t7
      t49 = t10*l1
      t51 = t19*l1
      t54 = t36*l1
      t56 = t10*t7
      t58 = t5*l1
      t60 = z*t14
      t64 = -6*t40-2*t42+2*t22+2*t45+4*t47+10*t49-4*t51-2*t36+14*t54-12*
     #t56-6*t58+2*t60+2*t5-2*t10
      t65 = t39+t64
      t68 = l1*(1-z)
      t69 = s13e2(t68)
      t82 = 8*t32
      t90 = 6*t3-12*t8+6*t12+24*t15+12*t17+8*t20+22*t23-24*t25-36*t27-6*
     #t6+6*t11-t82-12*t34+18*t37+18*t40+6*t42-6*t22-8*t45-12*t47
      t96 = 16*t58
      t98 = s2l1*t14
      t100 = s1l1*t14
      t102 = s1l1*t7
      t104 = t60*s1l1
      t107 = t25*s1l1
      t109 = t100*t5
      t111 = t11*t14
      t113 = t15*s2l1
      t116 = t60*s2l1
      t118 = t11*l1
      t120 = -28*t49+12*t51+6*t36-40*t54+36*t56+t96-8*t60-2*t98-2*t100+2
     #*t102+6*t104-6*t5+2*t107-6*t109+2*t111-6*t113+6*t10+6*t116+2*t118
      t121 = t90+t120
      t124 = -2*t1
      t127 = s13e1(t68)
      t140 = -5*t3+2*t8-5*t12-30*t15-15*t17-10*t20-23*t23+30*t25+38*t27+
     #4*t6-4*t11
      t151 = t82+38*t34-31*t37-12*t40+11*t42+4*t22+6*t45+8*t47+16*t49-8*
     #t51-4*t36
      t165 = 24*t54-31*t56-t96+10*t60+7*t98+7*t100-15*t102-3*t2-3*s1l1-8
     #*s2l1*t7+3*s2l1*l1
      t180 = 11*s1l1*l1-3*z+8*z*l1-21*t104+4*t5-7*t107+21*t109-7*t111+21
     #*t113-4*t10-21*t116-4*t118
      t186 = s13e0(t68)
      taux16 = 128*t1*t65*t69+128*(t1*t121+t124*t65)*t127-128*(t1*(t140+
     #t151+t165+t180)+t124*t121)*t186

CCC Coefficients

      tcoff29 = taux15
      tcoff30 = taux15
      tcoff31 = taux16
      tcoff32 = taux16

CCC Jet functions

      tFjet29 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff29)
      tFjet30 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff30)
      tFjet31 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff32)

CCC Sum

      qgRVnum8e2 = tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFjet
     #32*tcoff32
      RETURN
      END

      double precision function qgRVnum8e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -1+z
      t2 = dsqrt(z)
      t3 = t1*t2
      t5 = dsqrt(1-l5)
      t6 = dsqrt(l1)
      t7 = t5*t6
      t9 = z*l1
      t10 = t9-1+l1
      t11 = 1-l1
      t12 = dsqrt(t11)
      t15 = (1-l1+t9)**2
      t16 = 1/t15
      t17 = t16*l6
      t19 = t3*t7*t10*t12*t17
      t21 = t3*t5
      t25 = t21*t6*t10*t12*t16
      t27 = 512*t19-256*t25
      t28 = dsqrt(l5)
      t29 = t28*l5
      t32 = t6*t12*t11
      t34 = t21*t32*t17
      t37 = t21*t32*t16
      t39 = 512*t34-256*t37
      t41 = z**2
      t42 = 128*t41
      t43 = t41*z
      t44 = l1**2
      t45 = t44*l1
      t46 = t43*t45
      t47 = 64*t46
      t48 = z*t44
      t50 = t41*t45
      t51 = 128*t50
      t52 = t41*t44
      t53 = 320*t52
      t54 = t41*l1
      t55 = 384*t54
      t56 = t43*t44
      t57 = 256*t56
      t58 = 192*t41
      t59 = z*t45
      t60 = 64*t59
      t61 = t43*l1
      t63 = 192*z
      t64 = 256*t54
      t65 = 256*t9
      t66 = 256*t52
      t67 = 256*t48
      t68 = 128*z
      t71 = 64*t9
      t72 = -128+t47-64*t48-t51+t53-t55-t57+t58+t60+192*t61-t63+(-128+t6
     #4-t65-t66+t67-t68)*s2l1+t71
      t75 = 192*t52
      t76 = 192*t48
      t77 = 320*t9
      t79 = s1l1**2
      t82 = t1*z*l1
      t83 = -t11
      t84 = t83*t16
      t85 = l6**2
      t87 = t82*t84*t85
      t88 = 1024*t87
      t89 = 3*t9
      t92 = t1*t16
      t93 = t83*(t89-1+l1)*t92
      t96 = t82*t84*l6
      t97 = 1024*t96
      t98 = t88+128*t93-t97
      t100 = t41**2
      t101 = t100*l1
      t103 = 384*t61
      t104 = t100*t44
      t108 = t2*t41*t6*t5
      t109 = t12*t1
      t110 = t109*t17
      t111 = t108*t110
      t114 = t108*t109*t16
      t117 = (512*t111-256*t114)*t28
      t118 = t44**2
      t119 = t41*t118
      t121 = t43*t118
      t123 = z*t118
      t125 = t100*t118
      t127 = 320*t54
      t128 = 192*t9
      t129 = 64*z
      t131 = s2l1**2
      t133 = t83**2
      t135 = t133*t1*t16
      t136 = 64*t135
      t137 = -64+t27*t29+t39*t28-t42+t72*s1l1+(-64+192*t54+t68-t75+t76-t
     #77)*t79+t98*l5+128*t101-t103-256*t104+t117-192*t119+192*t121+64*t1
     #23-64*t125+(-64+t76+t127-t128-t75-t58-t129)*t131-t136
      t140 = t100*t45
      t145 = t61*t92
      t146 = 64*t145
      t150 = 64*t43
      t151 = -128+t67+t51-t60+t55+64*t56-320*t41-t128-64*t61-t47+t150-t5
     #3-t68
      t153 = 6*t9
      t155 = 2*l1
      t156 = t153-6*t48-t52-1+t155-t44
      t161 = l5**2
      t163 = t83*t41
      t164 = t92*t85
      t166 = 1024*t163*t164
      t168 = 2*t9
      t169 = 2*t44
      t177 = -4*z-t168+l1-t169-14*t52+6*t56-4*t59+6*t50-4*t46+t140+10*t4
     #8+9*t54+t45
      t179 = z*t16
      t182 = t92*l6
      t184 = 1024*t163*t182
      t186 = (t166+64*t177*t1*t179-t184)*t161
      t187 = 128*t9
      t191 = t2*z*(t54+3*z-t168+l1-1)
      t193 = t191*t7*t110
      t198 = t191*t6*t5*t12*t92
      t201 = (-512*t193+256*t198)*t29
      t208 = (-t166-128*(t52+t89-2*t48+t44+l1-2)*t41*t92+t184)*l5
      t209 = -704*t52-576*t46+192*t140+t67+576*t50-192*t59+704*t56+t146+
     #t151*s2l1+(-t88+64*t156*t1*t16+t97)*t161+t186+t150-t187+t201+t55-t
     #129+t208
      t212 = -l1*t1
      t213 = s13e3(t212)
      t215 = 2048*t87
      t216 = 2048*t96
      t224 = 128*t52
      t226 = 128*t45
      t227 = 256*t41
      t228 = 256*t59
      t233 = 128*t59
      t234 = 384*t9
      t235 = 256*t50
      t236 = 512*t48
      t243 = -t53+512*t56-128*t46+256-t103-t233+t234-t227-t76+(t235+t68+
     #t236-t64+256-512*t44+256*t45+t187+256*l1-512*t59)*s2l1+t235+t63+t1
     #27
      t246 = 2*t11*t41
      t248 = 1024*t246*t164
      t250 = 1024*t246*t182
      t251 = 2*t52
      t282 = 8*z+4*t9-t155+4*t44+28*t52-12*t56+8*t59-12*t50+8*t46-2*t140
     #-20*t48-18*t54-2*t45
      t289 = 128+(-t215+t216-256*t93)*l5+(1024*t193-512*t198)*t29+(t65+t
     #224-128*t44+t226+t227-t228+t68-t55+t51+128)*t131+t243*s1l1+(-t248+
     #t250-128*(-t251-t153+4*t48-t169-t155+4)*t41*t92)*l5+(-1024*t34+512
     #*t37)*t28+t58-384*t101+1024*t61+768*t104+768*t119-768*t121-256*t12
     #3+256*t125+(t248-t250+64*t282*t1*t179)*t161+128*t135
      t291 = 1920*t46
      t293 = 768*t48
      t302 = 256*z
      t304 = 512*t54
      t308 = 256*t46
      t309 = 512*t50
      t310 = 128*t43
      t311 = t302+t228-t57+448*t41-t304+128*t61+832*t52-576*t48+t308+256
     #-t309+t77-t310
      t328 = 128*t48
      t332 = 2048*t52+t291-640*t140-t293-1920*t50+640*t59-2048*t56+(-102
     #4*t19+512*t25)*t29-128*t145+t311*s2l1+(t215-t216+64*(-12*t9+12*t48
     #+t251+2-4*l1+t169)*t1*t16)*t161-t310+t234-960*t54+t68+(-1024*t111+
     #512*t114)*t28+(-t68-384*t44+t328+t51+t65-t228+t226-t55+384*l1+t66)
     #*t79
      t334 = s13e2(t212)
      t338 = 64*t50
      t339 = 192*t44
      t340 = 192*l1
      t341 = 64*t45
      t349 = 64*t54
      t350 = 64*t44
      t353 = -128*t140-640*t48+(t76+t338-t233-t339+t340+t341)*t79-t64+70
     #4*t46+t77-1024*t50+960*t52-320*t56+t129+448*t59+(t76+t349-t187+t33
     #8-t224+t129-t233-t350+t341+64)*t131
      t355 = 256*t44
      t356 = 128*l1
      t380 = 64+(-t328+t65+t224+(384*t48+t51-t224-t228-t355+t356+t226+12
     #8)*s2l1+128)*s1l1+(128*t48+128*t54-128*t9-128*t52+128*z+128)*s2l1+
     #(t88-t97-64*t156*t1*t16)*t161-t27*t29-t98*l5+t136-t39*t28-128*t123
     #+384*t119-384*t121+128*t125
      t382 = s13e0(t212)
      t385 = t67+t349-t233-64-64*t52+t338+t341-t339-t68+t340+t71
      t401 = t385*t79+(t293+t57+t309-t228-512*t9-t308+768*t54-1024*t52-3
     #84*t41)*s2l1-t58+256*t101-512*t61-512*t104+t117-960*t119+960*t121+
     #320*t123-320*t125-1856*t52-t291
      t407 = 448*t9
      t416 = 576*t140+896*t48+2112*t50-768*t59+1472*t56+t146+t186-t407+t
     #201+896*t54+(-t350+t236-448*t52-t58+t341-t407-t233+t338+t304)*t131
     #+(-t302-t42+t65+(t236-t66-t228+t226+t64-t42+t356+t51-t355-t302)*s2
     #l1+t64)*s1l1+t208
      t418 = s13e1(t212)
      taux15 = (t137+t209)*t213+(t289+t332)*t334-(t353+t380)*t382+(t401+
     #t416)*t418
      t1 = 1-l1
      t2 = z**2
      t3 = s1l1*t2
      t5 = l1**2
      t6 = t2*t5
      t8 = t2*z
      t9 = t5*l1
      t10 = t8*t9
      t12 = t2**2
      t13 = t12*t9
      t15 = z*t5
      t17 = t2*t9
      t19 = s1l1*z
      t20 = t19*l1
      t22 = t3*t5
      t24 = t19*t5
      t27 = s2l1*t8
      t28 = t27*t5
      t29 = 2*t28
      t30 = s2l1*t2
      t31 = t30*t5
      t33 = s2l1*z
      t34 = t33*t5
      t36 = t3*l1
      t38 = -2*t3+12*t6+6*t10-2*t13-4*t15-6*t17-8*t20-4*t22+2*t24+2*t2-t
     #29+4*t31-2*t34+14*t36
      t39 = t30*l1
      t41 = t33*l1
      t42 = 2*t41
      t44 = t2*l1
      t46 = z*t9
      t48 = s1l1*t8
      t49 = t48*l1
      t51 = t8*t5
      t52 = 12*t51
      t53 = t12*l1
      t55 = t8*l1
      t57 = t12*t5
      t61 = t48*t5
      t64 = 2*t39-t42-2*t8-6*t44+2*t46-6*t49-t52-4*t53+10*t55+4*t57-2*t2
     #7+2*t30+2*t61+2*t19
      t65 = t38+t64
      t68 = l1*(1-z)
      t69 = s13e3(t68)
      t71 = s1l1*t5
      t73 = t46*s2l1
      t75 = t27*l1
      t81 = 12*t15
      t93 = 2*t71+6*t73+2*t75+6*t3-36*t6-24*t10+8*t13+t81+24*t17+22*t20+
     #18*t22-12*t24-6*t2+6*t28-12*t31+6*t34-40*t36-8*t39+6*t41
      t95 = 16*t44
      t96 = s2l1*t9
      t98 = s1l1*t9
      t108 = 8*t61
      t110 = t10*s1l1
      t112 = t98*t2
      t114 = t46*s1l1
      t116 = t27*t9
      t118 = t17*s2l1
      t120 = 6*t8+t95-2*t96-2*t98-8*t46+18*t49+36*t51+12*t53-28*t55-12*t
     #57+6*t27-6*t30-t108-6*t19+2*t110-6*t112+6*t114+2*t116-6*t118
      t121 = t93+t120
      t124 = -2*t1
      t127 = s13e2(t68)
      t140 = -15*t71-21*t73-4*t75-4*t3+38*t6+30*t10-10*t13-15*t15-30*t17
     #-23*t20-31*t22
      t141 = 38*t24
      t151 = t141+4*t2-5*t28+2*t31+11*t34+24*t36+6*t39-5*t41-4*t8-t95+7*
     #t96
      t156 = s2l1*l1
      t159 = z*l1
      t166 = 7*t98-3*s1l1+10*t46+3*t156-12*t49+8*t159-31*t51-8*t53+16*t5
     #5+8*t57-4*t27
      t168 = s1l1*l1
      t173 = s2l1*t5
      t180 = 4*t30+11*t168+t108+4*t19-7*t110+21*t112-8*t173-21*t114-7*t1
     #16-3*t33+21*t118-3*z
      t182 = t140+t151+t166+t180
      t186 = s13e1(t68)
      t200 = 19*t71+18*t73-24*t6-12*t10+4*t13+t81+12*t17+18*t20+19*t22-t
     #141+t29+9*t31-24*t34+3*t39+t42+13*t44
      t216 = -6*t96-6*t98+5*s1l1-4*t46-5*t156-13*t159+t52-18*t168+6*t110
     #-18*t112+13*t173+18*t114+6*t116+5*t33-18*t118+5*z
      t221 = s13e0(t68)
      taux16 = 128*t1*t65*t69+128*(t1*t121+t124*t65)*t127+128*(t1*t182+t
     #124*t121)*t186-128*(t1*(t200+t216)+t124*t182)*t221

CCC Coefficients

      tcoff29 = taux15
      tcoff30 = taux15
      tcoff31 = taux16
      tcoff32 = taux16

CCC Jet functions

      tFjet29 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff29)
      tFjet30 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff30)
      tFjet31 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff31)
      tFjet32 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff32)

CCC Sum

      qgRVnum8e3 = tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFjet
     #32*tcoff32
      RETURN
      END

      double precision function qgRVnum9e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = 1-z
      t3 = t2**2
      t5 = t3**2
      t6 = t5*t3
      t7 = 2048.D0/3.D0*t6
      t9 = t5*t2
      t11 = t3*t2
      t17 = s2l1**2
      t18 = 1024*t17
      t19 = s1l1**2
      t20 = 320*t19
      t21 = 64*s1l1
      t22 = 1472*s2l1
      t23 = s1l1*s2l1
      t24 = 512*t23
      t25 = 320*t17
      t27 = 1024*s2l1
      t28 = 256*t23
      t31 = 128*t17
      t32 = 320*s2l1
      t35 = 64*t17
      t36 = 256*t19
      t37 = 1344*s1l1
      t38 = 768*s2l1
      t39 = 640*t23
      t43 = 576*t19
      t48 = l1**2
      t51 = 256*s2l1
      t54 = t5*t11
      t56 = 192*t17
      t59 = 768*t23
      t62 = 1152*t17
      t63 = 192*t19
      t66 = 384*t23
      t69 = 256*s1l1
      t72 = 896*s2l1
      t76 = t48**2
      t77 = t76*t48
      t84 = t76**2
      t86 = 1024.D0/3.D0*t54
      t88 = 128*z
      t95 = -t54+t6
      t96 = 64*t95
      t97 = t84*t48
      t98 = t96*t97
      t100 = 64.D0/3.D0*t54
      t103 = 320*t6
      t108 = 640.D0/3.D0*t54
      t109 = 256.D0/3.D0*t6
      t110 = 128*t9
      t112 = t84*l1
      t114 = 1792.D0/3.D0*t9
      t117 = 2816.D0/3.D0*t5
      t130 = 320.D0/3.D0*t6
      t136 = 256.D0/3.D0*t54
      t137 = 128*t11
      t141 = t48*l1
      t142 = t76*t141
      t151 = 128.D0/3.D0*t54
      t157 = t76*l1
      t159 = 4096.D0/3.D0*t5
      t165 = t98+1408-2112*z+(-t100-24448.D0/3.D0*t11-2496*t5-t103+3968*
     #t3+21056.D0/3.D0*t9)*t77+(t108-t109-t110)*t112+(t114+256-576*z-889
     #6.D0/3.D0*t11+t117+3520.D0/3.D0*t3)*t48+704*t11+(-704.D0/3.D0*t54-
     #t103+7808.D0/3.D0*t9-2048*t5)*t84-2112*t3+(-t100-8640*t11+2624.D0/
     #3.D0*t9-t130+64*t3+12608.D0/3.D0*t5+10880.D0/3.D0-17024.D0/3.D0*z)
     #*t76+(t136+t137-21824.D0/3.D0*t9+t7+19136.D0/3.D0*t5)*t142+(t109+7
     #680*t3+1280*t5-4352.D0/3.D0*t9-5248+8512*z-6656.D0/3.D0*t11)*t141+
     #(t151+55360.D0/3.D0*t11-2240*t9-28480.D0/3.D0*t3-20672.D0/3.D0*t5+
     #128-t88)*t157+(-t159-512.D0/3.D0-1280*t3+2688*t11-64.D0/3.D0*z)*l1
      t172 = 1088*z
      t186 = 128*t141
      t188 = -64*t48-t186+192*l1
      t190 = 320*t48
      t191 = 320*l1
      t193 = t191-448*t48+t186
      t197 = 2*l6-1
      t198 = t197**2
      t201 = l5-1
      t203 = -1+l1
      t204 = z*l1
      t206 = (1-l1+t204)**2
      t207 = 1/t206
      t212 = dsqrt(z)
      t214 = dsqrt(l5)
      t216 = dsqrt(-t201)
      t219 = dsqrt(l1)
      t223 = l5*z*l1-l1+1-l5+l5*l1
      t226 = dsqrt(-t203)
      t231 = t223**2
      t233 = 4*l1*t198*z*l5*t201*t203*t207-4*t197*t212*t214*t216*t219*t2
     #23*t207*t226+t231*t207
      t234 = dsqrt(t233)
      t236 = 640*t17
      t237 = 768*s1l1
      t238 = 1280*s2l1
      t240 = 640*s2l1
      t241 = 128*t23
      t244 = 64*t19
      t248 = 128*s2l1
      t265 = 2176*t11
      s1 = 5312.D0/3.D0-2304*z-12160.D0/3.D0*t3+(-t7+64.D0/3.D0*t3-1088.
     #D0/3.D0+7360.D0/3.D0*t9+4672.D0/3.D0*t11-832.D0/3.D0*z-9728.D0/3.D
     #0*t5)*l1+(-448-t18+t20+t21-t22-t24+(-704-t25-448*s1l1-t27-t28)*t11
     #+(192+t31+t32)*t5+(832-t35+t36+t37+t38+t39)*t3+(-128+896*t17-t43-1
     #216*s1l1+t38-t28)*t2)*t48+((-640*s1l1+t51+384)*t6+512*t54+(-t56-10
     #88*t19-1536*s1l1+t51-t59)*t11+(896-t62+t63+1920*s1l1-1664*s2l1-t66
     #)*t5+(t25+t20+t69+t51+t66)*t3+(-1792+t18+t43+t72+t59)*t9)*t77+(-30
     #208.D0/3.D0*t6+17024.D0/3.D0*t9+15488.D0/3.D0*t54-768*t5)*t84+(t86
     #-9344.D0/3.D0*t11+576-t88+832.D0/3.D0*t6+16448.D0/3.D0*t5+256.D0/3
     #.D0*t3-3392*t9)*t48
      t275 = s1+t165*t19+(-96896.D0/3.D0*t5-1856*t54+57088.D0/3.D0*t9+61
     #568.D0/3.D0*t11-1408*t6-1088+t172-8704.D0/3.D0*t3)*t157+(8896*t5-4
     #0000.D0/3.D0*t3+20864.D0/3.D0+3392.D0/3.D0*t54-7360*t9-25856.D0/3.
     #D0*z+10688.D0/3.D0*t11+448.D0/3.D0*t6)*t76+(t188*s1l1-t190+t193*s2
     #l1+t191)*t234+(640+t236+t63+t237+t238+t59+(320+t25+128*s1l1+t240+t
     #241)*t3+(-640-t236-t244-512*s1l1-t238-t24)*t2+(-64-t35-t248)*t11)*
     #l1+(-43072.D0/3.D0*t11-7040*t6+4672*t54+2752*t3-22336.D0/3.D0*t9+6
     #4256.D0/3.D0*t5)*t77+(-20096.D0/3.D0*t54-7936.D0/3.D0*t5+42304.D0/
     #3.D0*t6-20800.D0/3.D0*t9+t265)*t142+(-t188*s1l1-t193*s2l1-t191+t19
     #0)*t234+1024.D0/3.D0*t9
      t281 = 5248.D0/3.D0*t9
      t289 = 896*t6
      t301 = 256*z
      t308 = 640.D0/3.D0*t6
      t313 = 128.D0/3.D0*t9
      t320 = 128*t54
      t322 = 128.D0/3.D0*t6
      t342 = 128*t95*t97+2176.D0/3.D0*t11+(t151-2816*t11-t281+512.D0/3.D
     #0*t6+9344*t3+3968.D0/3.D0*t5-18944.D0/3.D0+9472*z)*t141+(-t108-267
     #52.D0/3.D0*t11-t289+4096*t3+25600.D0/3.D0*t9-7808.D0/3.D0*t5)*t77+
     #(1664.D0/3.D0*t54-1280.D0/3.D0*t6-t110)*t112-2176*t3+(-t109+t117+7
     #68*t9+768-t301+2560.D0/3.D0*t3-8576.D0/3.D0*t11)*t48+4352.D0/3.D0-
     #2176*z+(-896*t54+t308+8576.D0/3.D0*t9-2176*t5)*t84+(t313+8704.D0/3
     #.D0*t11-1536*t5-4480.D0/3.D0*t3-896.D0/3.D0-1408.D0/3.D0*z)*l1+(-t
     #320-27520.D0/3.D0*t11+13568.D0/3.D0-t322-3584.D0/3.D0*t3+2432.D0/3
     #.D0*t9-20096.D0/3.D0*z+15616.D0/3.D0*t5)*t76+(640*t54+256*t11+2048
     #0.D0/3.D0*t5+t7-25216.D0/3.D0*t9)*t142+(t320+59648.D0/3.D0*t11+256
     #*t6-28288.D0/3.D0*t3-8192.D0/3.D0*t9-128+t88-23936.D0/3.D0*t5)*t15
     #7
      t353 = 448.D0/3.D0*t54
      t359 = 5504.D0/3.D0*t6
      t369 = 1152*t9
      t376 = 896*z
      t396 = 1664.D0/3.D0*t9
      t398 = 2432.D0/3.D0*t3
      t410 = t342*s2l1+8512.D0/3.D0-3904*z+(t100-30848.D0/3.D0*t11-4736.
     #D0/3.D0*t6+28480.D0/3.D0-9984*t3-2304*t9-38272.D0/3.D0*z+43904.D0/
     #3.D0*t5)*t76+(-t353+t265+33152.D0/3.D0*t6-57152.D0/3.D0*t9+17920.D
     #0/3.D0*t5)*t142+t98+(-t100-35776.D0/3.D0+t359+21312*t3-3136*t9+170
     #88*z+2624*t5-10560*t11)*t141+9536.D0/3.D0*t11+(t353+3008.D0/3.D0*t
     #6-t369)*t112-5312*t3+(t100-36032.D0/3.D0*t3+115648.D0/3.D0*t11+851
     #2.D0/3.D0*t6-896+t376+14656.D0/3.D0*t9-100160.D0/3.D0*t5)*t157+(-t
     #100-5888*t6+27136.D0/3.D0*t9-3136*t5)*t84+(4160.D0/3.D0*t9-1984.D0
     #/3.D0-320.D0/3.D0*z+13184.D0/3.D0*t11-4480*t5-1024*t3)*l1+(-1984.D
     #0/3.D0*t6-18304.D0/3.D0*t11+16832.D0/3.D0*t5-t396+1152-1216*z+t398
     #)*t48+(64*t54-64064.D0/3.D0*t11-25984.D0/3.D0*t6+6208*t3+10880*t9+
     #12864*t5)*t77-704*t5
      t417 = 64*s2l1
      t422 = 576*s1l1
      t426 = 192*s1l1
      t427 = 192*s2l1
      t435 = 1344*t17
      t436 = 704*s1l1
      t437 = 1216*s2l1
      t477 = 512.D0/3.D0*t54
      t488 = 2176.D0/3.D0*t54
      t500 = 768*t54
      t507 = 64.D0/3.D0*t9
      t509 = 1408.D0/3.D0*t6
      t520 = 1280.D0/3.D0*t54
      t530 = t98+5440.D0/3.D0-7616.D0/3.D0*z+(t136-26816.D0/3.D0*t11+236
     #8.D0/3.D0*t6-25792.D0/3.D0+47552.D0/3.D0*t3-6208.D0/3.D0*t9+34112.
     #D0/3.D0*z+2752*t5)*t141-896.D0/3.D0*t5+(t86-t308-t110)*t112+5248.D
     #0/3.D0*t11+(t477+76672.D0/3.D0*t11+6400.D0/3.D0*t6-8960*t3+640*t9-
     #384+384*z-57472.D0/3.D0*t5)*t157+(11200.D0/3.D0*t9-2112*t5-t289-t4
     #88)*t84-9728.D0/3.D0*t3+(-t320-24704.D0/3.D0*t11-1664.D0/3.D0*t6-1
     #9328.D0/3.D0*t3-6656.D0/3.D0*t9+19456.D0/3.D0-26176.D0/3.D0*z+3328
     #0.D0/3.D0*t5)*t76+(t500+384*t11+3904*t6-11584*t9+6528*t5)*t142-t50
     #7+(-t100-8128.D0/3.D0*t11+704-t509+704*t3+t396+768*z+5632.D0/3.D0*
     #t5)*t48+(4160*t3-4800*t6-34880.D0/3.D0*t11+31232.D0/3.D0*t9+6848.D
     #0/3.D0*t5-t520)*t77+(t322-6208.D0/3.D0*t3-8896.D0/3.D0*t5+2048.D0/
     #3.D0*t9-64.D0/3.D0-3776.D0/3.D0*z+11456.D0/3.D0*t11)*l1
      t533 = 64*l1-64*t48
      t549 = -t95
      t585 = t98-t507+(t136-9856.D0/3.D0-7424.D0/3.D0*t11+t322+16256.D0/
     #3.D0*t3-832*t9+4288*z+2752.D0/3.D0*t5)*t141+(-t520-4352*t11-t509+1
     #920*t3+14080.D0/3.D0*t9-t159)*t77+1024.D0/3.D0*t549*t112+2176.D0/3
     #.D0-t172+64*t5+(t477+28352.D0/3.D0*t11-t281+t509-11840.D0/3.D0*t3-
     #12416.D0/3.D0*t5-256+t301)*t157+(-t488+1792.D0/3.D0*t6+t369-1024*t
     #5)*t84-3200.D0/3.D0*t3+(t322-1088*t3-2432.D0/3.D0*t5-t313-640*z+49
     #28.D0/3.D0*t11)*l1+(-t320-10496.D0/3.D0*t11-t109+2688-2048*t3+704.
     #D0/3.D0*t9-3712*z+8512.D0/3.D0*t5)*t76+(t500+t137-t308-4032*t9+100
     #48.D0/3.D0*t5)*t142+896.D0/3.D0*t11+(-t100+128-t130+t398+t114+t376
     #+512.D0/3.D0*t5-3584.D0/3.D0*t11)*t48
      t597 = 448*t19
      t606 = 448*t17
      t620 = 512*t19
      t624 = -256+t606-t597-896*s1l1+t248-t28+(1024-t35+t244+1472*s1l1+t
     #437+t59)*t2+(704-t25+t436-t417+t241)*t5+(-320-t35-t51)*t9+(-1216-t
     #62+768*t19-t21-2880*s2l1-t39)*t3+(64+1280*t17-t620-t37+1984*s2l1)*
     #t11
      t639 = 1920*t17
      t654 = 64-t35-t244+t21+t417+(-256-t606+832*t19+t422-t27)*t2+(1472-
     #t639+2240*s1l1-t238-t28)*t11+(-64+t25-t422+t240)*t9+(-192+t639-128
     #0*t19-2368*s1l1+2432*s2l1+t241)*t3+(-1344+t56+t620+t21-t72+t241)*t
     #5+(t417+320)*t6
      s1 = 11008.D0/3.D0*t11+t410*s1l1-1728*t5+((-64*s1l1+64*s2l1)*t9-32
     #0*t54+(320+t21-t417)*t6)*t112+((-128+t56+t63+t422+t28)*t9+704*t54+
     #(64-t56-t63-t426-t427-t28)*t5+(-640+t427-384*s1l1)*t6)*t84+((256-t
     #435-t244+t436-t437-t39)*t3+(-768+t426-t427)*t6+(256*t17-t36-320*s1
     #l1+t32)*t2-t320+(-768+t435+1472*t19+832*s1l1+448*s2l1+1280*t23)*t1
     #1+(-128+512*t17-960*t19-2496*s1l1+t22-t66)*t5+(1536-768*t17-t63+10
     #88*s1l1-832*s2l1-t28)*t9)*t157+t530*s2l1
      t661 = s1-t96*t97+t585*t17+384*t549*t97+((-t51+t237+448)*t6+(768-7
     #04*t17-t43-1024*s1l1-512*s2l1-t59)*t9+(-384+832*t17+t597+t27+t59)*
     #t5-832*t54+(-t31+128*t19+t69-t51)*t11)*t142+t624*t141+(-3008.D0/3.
     #D0*t54+14656.D0/3.D0*t5-14016*t11-23552.D0/3.D0+10240*z+t359-1280*
     #t9+52288.D0/3.D0*t3)*t141+t654*t76+(-6400.D0/3.D0*t54+9664.D0/3.D0
     #*t6-1088*t9)*t112
      t664 = (l1+z-t204)**2
      taux17 = (t275+t661)/t664*t207
      t1 = 1-z
      t2 = t1**2
      t3 = t2*t1
      t4 = t2**2
      t5 = t4*t3
      t10 = t4*t2
      t12 = t4*t1
      t16 = l1**2
      t17 = t16*l1
      t19 = t10-t5
      t21 = t16**2
      t22 = t21**2
      t23 = t22*t16
      t24 = 256*t19*t23
      t25 = 1536*s1l1
      t26 = 1536*s2l1
      t30 = 1024*s2l1
      t37 = 512*s1l1
      t38 = 512*s2l1
      t45 = t21*t16
      t62 = t22*l1
      t84 = t21*t17
      t89 = 256*s2l1
      t90 = 256*s1l1
      t91 = 768*s1l1
      t103 = s1l1-s2l1
      t104 = 256*t103
      t110 = 512*t5
      t112 = 768*s2l1
      t115 = 3840*s1l1
      t121 = t21*l1
      t125 = 1280*s1l1
      t132 = 1792*s2l1
      s1 = -12544.D0/3.D0+(4352*t5+53248.D0/3.D0+146432.D0/3.D0*t3-14489
     #6.D0/3.D0*t2-17408*t4-20480.D0/3.D0*t10+4864.D0/3.D0*t12-20992*z)*
     #t17+t24+((8960+t25+t26)*t12-2048*t5+(-t30-1536+2560*s1l1)*t10+(-12
     #80+t25-2048*s2l1)*t3+(768+t37+t38)*t2+(-4864-6144*s1l1+t30)*t4)*t4
     #5+(-5120*t5-23552*t4-17920-97792.D0/3.D0*t3+146944.D0/3.D0*t2-4096
     #.D0/3.D0*t10+94720.D0/3.D0*t12+20224*z)*t21+2048*t19*t23
      t139 = s1+(10752*t5-14592*t10+3840*t12)*t62+(-73472.D0/3.D0*t5-179
     #2*t4+124160.D0/3.D0*t10-15104*t12)*t22+(-20480*t5-89600*t4+75008.D
     #0/3.D0*t10+37376*t3-2560*t2+150784.D0/3.D0*t12)*t45+(30208*t5+7475
     #2.D0/3.D0*t4-54528*t10-7680*t3+21248.D0/3.D0*t12)*t84+((-256-256*s
     #2l1)*t4-t89-t90-256+(1024+t91+t30)*t1+(-1536-t91-t26)*t2+(1024+t90
     #+t30)*t3)*t16+((4864+t90+3328*s2l1)*t3+t104*t1+(1280+6400*s1l1-230
     #4*s2l1)*t4+t110+(-2304-2304*s1l1-t112)*t2+(-7424-t115-t112)*t12+(-
     #t91+3072+t112)*t10)*t121+(512+(1024+t38)*t12+(-3072-t125-1280*s2l1
     #)*t4+(2560+2816*s1l1+t89)*t3+t37+t38+(512-t125+t132)*t2+(-1536-t91
     #-t132)*t1)*t17
      t150 = 512.D0/3.D0*t3
      t154 = 1280*z
      t176 = 3072*s2l1
      t191 = 1024*z
      t193 = 60928.D0/3.D0*t4
      t200 = 1280*t5
      t203 = 256*t5
      t238 = 256.D0/3.D0*t5
      t241 = 6656*t10
      t263 = 7936.D0/3.D0*t10
      t272 = (-t203+8704.D0/3.D0*t3+115712.D0/3.D0*t2-21248+16384.D0/3.D
     #0*t10+41984.D0/3.D0*t12+25600*z-39424*t4)*t21+(-16640.D0/3.D0*t12-
     #3584.D0/3.D0*t2+37376.D0/3.D0*t4+1280-t154-20992.D0/3.D0*t3)*l1+(-
     #1280.D0/3.D0*t5-173312.D0/3.D0*t4+96256.D0/3.D0*t10+48896*t3-8448*
     #t2-14336*t12)*t45+(-1792*t5+3840*t4+20992*t10-23040*t12)*t22+t24+(
     #t200+4096.D0/3.D0*t4-116224.D0/3.D0*t10-7168*t3+43264*t12)*t84+(t2
     #38+101632.D0/3.D0*t3+67328.D0/3.D0-21760.D0/3.D0*t4-t241+6400*t12-
     #48896*t2-28416*z)*t17-16640.D0/3.D0+6912*z+(t203+101632*t4-34048.D
     #0/3.D0*t10-224512.D0/3.D0*t3+27904.D0/3.D0*t2+3584-3584*z-85760.D0
     #/3.D0*t12)*t121+12544*t2+(3328.D0/3.D0*t5-14080.D0/3.D0*t10+3584*t
     #12)*t62+2816*t4+(t263-512+12800.D0/3.D0*t12+768*z+39424.D0/3.D0*t3
     #-5632.D0/3.D0*t2-17664*t4)*t16-29440.D0/3.D0*t3
      t352 = (-t110-14336.D0/3.D0+3584.D0/3.D0*t10-t150+38912.D0/3.D0*t2
     #+15104.D0/3.D0*z+34816.D0/3.D0*t12-t193)*t21+(512.D0/3.D0*t10-7936
     #.D0/3.D0*t3-3072*t12+256.D0/3.D0*t2-256.D0/3.D0+256.D0/3.D0*z+1664
     #0.D0/3.D0*t4)*l1+13312.D0/3.D0*t2+(-8704.D0/3.D0*t5-256*t4+t241-10
     #496.D0/3.D0*t12)*t22+t24+(3072*t5+3584.D0/3.D0*t4-46336.D0/3.D0*t1
     #0+33536.D0/3.D0*t12)*t84+(1024.D0/3.D0*t5+21248.D0/3.D0-58624.D0/3
     #.D0*t2-5376*t4-t263+55040.D0/3.D0*t3+1792*t12-22784.D0/3.D0*z)*t17
     #-256.D0/3.D0*t12-4352.D0/3.D0+4352.D0/3.D0*z+5120.D0/3.D0*t4-4096.
     #D0/3.D0*t19*t62+(-5120.D0/3.D0*t5-56576.D0/3.D0*t4+44288.D0/3.D0*t
     #10+26368.D0/3.D0*t3-256*t2-8192.D0/3.D0*t12)*t45+(-t238-768+4864.D
     #0/3.D0*t2-16384.D0/3.D0*t4+1024*t10+4352.D0/3.D0*t3+6656.D0/3.D0*t
     #12+t191)*t16+(2048.D0/3.D0*t5+125440.D0/3.D0*t4-4608*t10-63488.D0/
     #3.D0*t3+2048.D0/3.D0*t2-17408*t12)*t121-4608*t3
      t354 = ((-768-t90-t89)*t4+(1024-t125+t30)*t12+(2560-t112+t25)*t10-
     #2816*t5)*t22+(29696.D0/3.D0*t4-t150-9984*t12+8704.D0/3.D0*t10-1177
     #6.D0/3.D0*t2+1280-t154)*l1+((-4608+t25-t26)*t12-512*t103*t3+3328*t
     #5+(3072+2048*s1l1)*t4+(-3072*s1l1-1792+t30)*t10)*t84+(-256+(2560+t
     #115)*t2+(1024+1792*s1l1-t38)*t12+(4608-t91+t176)*t4+(-7168-4352*s1
     #l1-t176)*t3+(512-t90+t30)*t1+(-t89-1280)*t10-t90-t89)*t21+(-4352.D
     #0/3.D0*t5+46592.D0/3.D0*t12-1536*t10-768+t191+2816.D0/3.D0*t2-t193
     #+22784.D0/3.D0*t3)*t16+(t104*t12+(t89-1280-t90)*t10+t200)*t62+t272
     #*s1l1+4864*z-4352.D0/3.D0*t12+7168*t4+(24832.D0/3.D0*t5+332032.D0/
     #3.D0*t4+22528.D0/3.D0*t10-119552.D0/3.D0*t3-21248.D0/3.D0*t2+3840-
     #3840*z-250112.D0/3.D0*t12)*t121-40448.D0/3.D0*t3+35840.D0/3.D0*t2+
     #t352*s2l1
      t356 = z*l1
      t358 = (l1+z-t356)**2
      t362 = (1-l1+t356)**2
      taux18 = (t139+t354)/t358/t362

CCC Coefficients

      tcoff33 = taux17
      tcoff34 = taux17
      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet33 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff33)
      tFjet34 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff34)
      tFjet35 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff35)
      tFjet36 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff36)

CCC Sum

      qgRVnum9e0 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFjet
     #36*tcoff36
      RETURN
      END

      double precision function qgRVnum9e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      t2 = t1**2
      t5 = z*l1
      t7 = (1-l1+t5)**2
      t8 = 1/t7
      t11 = 1-z
      t12 = t11**2
      t13 = t12**2
      t14 = t13*t12
      t15 = 448.D0/3.D0*t14
      t16 = t13*t11
      t17 = 1152*t16
      t19 = t12*t11
      t20 = t13*t19
      t21 = 512.D0/3.D0*t20
      t23 = 23552.D0/3.D0*t12
      t26 = l1**2
      t27 = t26*l1
      t30 = 1088.D0/3.D0*t14
      t37 = t26**2
      t38 = t37*t26
      t41 = 128.D0/3.D0*t16
      t44 = 4928.D0/3.D0*t19
      t45 = 128.D0/3.D0*t20
      t47 = 320.D0/3.D0*t14
      t51 = 832.D0/3.D0*t14
      t55 = 1024.D0/3.D0*t20
      t59 = t37*l1
      t62 = 5824.D0/3.D0*t13
      t64 = 512*t19
      t67 = t37**2
      t68 = t67*l1
      t71 = 1024.D0/3.D0*t16
      t72 = 4352.D0/3.D0*t19
      t74 = 512.D0/3.D0*t13
      t75 = 64*t14
      t79 = 256*t20
      t81 = 192*t14
      t91 = 256*t12
      t96 = 128*t20
      t98 = 256*t14
      t99 = 256*t13
      t101 = t67*t26
      t102 = (t96+256*t16-t98-t99)*t101
      t105 = 1536*t20
      t107 = 1024*z
      t111 = 8896.D0/3.D0*t19
      t113 = t37*t27
      t115 = (t15+14336.D0/3.D0+t17-8960.D0/3.D0*t13-t21+6336*t19-t23-26
     #24*z)*t27+(2560.D0/3.D0*t20-t30-11008.D0/3.D0-6848*t16+11776.D0/3.
     #D0*z+30016.D0/3.D0*t13-29888.D0/3.D0*t19+29440.D0/3.D0*t12)*t38-23
     #68.D0/3.D0*t19-t41+(-1216*t16+4928.D0/3.D0*t13-t44+t45+6016.D0/3.D
     #0*t12+t47+1472.D0/3.D0+5312.D0/3.D0*z)*t26+(-t51+4224*t16+2560*t19
     #-11072.D0/3.D0*t13+2880-t55-14848.D0/3.D0*t12-5248.D0/3.D0*z)*t59+
     #(-2048.D0/3.D0*t20+t62+4096.D0/3.D0*t14-t64-6400.D0/3.D0*t16)*t68+
     #4864.D0/3.D0*t12+(t71-t72+9088.D0/3.D0*t12+t74-t75+384+15680.D0/3.
     #D0*z)*l1+(t79+3136.D0/3.D0*t19-t81-1856.D0/3.D0*t16-768-640*t13-76
     #8*t12-16064.D0/3.D0*z)*t37+(4352.D0/3.D0*t20+4736.D0/3.D0*t19-8128
     #.D0/3.D0*t14+t91+3520*t16-10816.D0/3.D0*t13)*t67-1216.D0/3.D0+t102
     #+192*t13+5120.D0/3.D0*z+(-t105+2240*t14+1024-t107-10432.D0/3.D0*t1
     #2+4096.D0/3.D0*t16-8384.D0/3.D0*t13+t111)*t113
      t116 = s2l1**2
      t119 = 1792*t37
      t121 = 3264*l1
      t124 = 256*t59
      t125 = 6272.D0/3.D0*t37
      t127 = 9920.D0/3.D0*l1
      t128 = 512*t59
      t132 = 6592.D0/3.D0*t37
      t135 = 12544.D0/3.D0*l1
      t138 = 2*l6-1
      t139 = t138**2
      t142 = l5-1
      t144 = -1+l1
      t149 = dsqrt(z)
      t151 = dsqrt(l5)
      t153 = dsqrt(-t142)
      t156 = dsqrt(l1)
      t160 = l5*z*l1-l1+1-l5+l5*l1
      t162 = -t144
      t163 = dsqrt(t162)
      t164 = t8*t163
      t168 = t160**2
      t170 = 4*l1*t139*z*l5*t142*t144*t8-4*t138*t149*t151*t153*t156*t160
     #*t164+t168*t8
      t171 = dsqrt(t170)
      t182 = 17536.D0/3.D0*t13
      t183 = 512*z
      t185 = 10432*t16
      t191 = 1280.D0/3.D0*t20
      t192 = 2176.D0/3.D0*t14
      t197 = 256*t19
      t198 = 320*t14
      t207 = 8320.D0/3.D0*t13
      t213 = 2304*z
      t214 = 1280*t12
      t215 = 7552.D0/3.D0*t13
      t218 = 1280.D0/3.D0*t19
      t220 = 128*t14
      t224 = 256.D0/3.D0*t20
      t244 = 2432.D0/3.D0*t13
      t248 = 17792.D0/3.D0*t12
      t253 = 768*t13
      t257 = (-t191+t192-6016.D0/3.D0*t16+t62)*t68+(1408.D0/3.D0*t20+t19
     #7-t198+t91+2752.D0/3.D0*t16-3520.D0/3.D0*t13)*t67+(-576*t14-7360.D
     #0/3.D0*t12+12800.D0/3.D0*t16+1472.D0/3.D0*t19-t207-t21)*t113+t102+
     #(1024.D0/3.D0*t14+22592.D0/3.D0*t19-11648.D0/3.D0*t16-2048+t213+t2
     #14-t215+t45)*t38-t218-768+2752.D0/3.D0*z+(t220+640.D0/3.D0*t16-173
     #44.D0/3.D0*t19+5824*t13+26816.D0/3.D0-t224-9216*t12-8832*z)*t59+37
     #76.D0/3.D0*t12+(t75-1472.D0/3.D0*t16-24832.D0/3.D0-2752.D0/3.D0*t1
     #9-2624.D0/3.D0*t13+t45+27328.D0/3.D0*t12+5440*z)*t37+(1024*t16-128
     #*t13+1664+13312.D0/3.D0*t12-2432*t19-t47+29696.D0/3.D0*z)*t27+(-t2
     #44-4864.D0/3.D0+30464.D0/3.D0*z+1024.D0/3.D0*t19-896.D0/3.D0*t16+t
     #248)*t26+(-4096.D0/3.D0*t19+6208.D0/3.D0+896*t12+t253+3072*z)*l1
      t258 = s1l1**2
      t263 = 18688.D0/3.D0*z
      t265 = 1792.D0/3.D0*t16
      t268 = 960*t20
      t282 = 1856*t20
      t292 = 5504.D0/3.D0*t19
      t293 = 1408.D0/3.D0*t14
      t313 = 3136*t16
      t314 = 18560.D0/3.D0*t13
      t320 = 192*t20
      t322 = 768*t16
      t324 = (t320-t99-704*t14+t322)*t101
      t327 = 1088*t16
      t335 = 26752.D0/3.D0*t19
      t344 = 10496.D0/3.D0*t19
      t348 = (-t268-768*t19+12416.D0/3.D0*t14-5888*t16+10432.D0/3.D0*t13
     #)*t68+(38144.D0/3.D0*t16+2752*t19-8192*t14+5696.D0/3.D0*t20-27520.
     #D0/3.D0*t13)*t67+(-t282+17984.D0/3.D0*t14+1536-1536*z-6144*t12-544
     #0*t16+2368.D0/3.D0*t13+5120*t19)*t113+1984.D0/3.D0*t13-t292+(-t293
     #-24064.D0/3.D0*t16+23744*t12+t268-7040+52672.D0/3.D0*t13-80192.D0/
     #3.D0*t19+7808*z)*t38+(-1216.D0/3.D0*t20-t75+10432.D0/3.D0*t16-2214
     #4.D0/3.D0*t13+8896+55936.D0/3.D0*t19-22400*t12-7040*z)*t59-128+723
     #2.D0/3.D0*z+7808.D0/3.D0*t12+(320*t20-3712.D0/3.D0*t14+11392.D0/3.
     #D0+t313-t314+15872.D0/3.D0*t19-8896*t12-57088.D0/3.D0*z)*t37+t324-
     #320.D0/3.D0*t16+(-t320+30592.D0/3.D0*t19+t30+t327+5120-4736*t13-20
     #672.D0/3.D0*t12+21184.D0/3.D0*z)*t27+(t81-1664*t16+4096*t13-6656.D
     #0/3.D0+t45-t335+13248*t12+39296.D0/3.D0*z)*t26+(-64*t16-1600.D0/3.
     #D0+8704*t12+3328.D0/3.D0*t13-64.D0/3.D0*t14-t344+44416.D0/3.D0*z)*
     #l1
      t351 = t162*t2
      t352 = l6**2
      t353 = t8*t352
      t355 = 1024.D0/3.D0*t351*t353
      t356 = t8*l6
      t358 = 1024.D0/3.D0*t351*t356
      t360 = 2*t5
      t362 = t1*t26
      t364 = t1*z
      t374 = z*t26
      t376 = t1*l1
      t378 = 4*z+t360-l1+2*t26+14*t362-6*t364*t26+4*z*t27-6*t1*t27+4*t36
     #4*t27-t2*t27-10*t374-9*t376-t27
      t383 = l5**2
      t387 = t376+3*z-t360+l1-1
      t389 = t149*t364*t387*t156
      t390 = t153*t163
      t394 = t390*t8
      t398 = t151*l5
      t400 = 3136*s2l1
      t406 = s1l1*s2l1
      t429 = 2*s2l1
      t431 = (s1l1+1+s2l1)*(t429+1)
      s1 = 64.D0/3.D0*t2*t1*l1*t8+t115*t116+((-15808.D0/3.D0*t27+t119-10
     #88.D0/3.D0-12992.D0/3.D0*t26-t121)*s1l1-t124+(-2624.D0/3.D0+t125-5
     #440.D0/3.D0*t27-t127-t128-2368*t26)*s2l1+t132-11840.D0/3.D0*t27-27
     #52.D0/3.D0-4288*t26-t135)*t171+(39040.D0/3.D0*t13+34112.D0/3.D0*t1
     #9+9728.D0/3.D0*t14+14912.D0/3.D0-55616.D0/3.D0*t16-50048.D0/3.D0*t
     #12+2880*t20-3328*z)*t59+(t182+512-t183+9024*t20-t185-15040.D0/3.D0
     #*t12-21760.D0/3.D0*t14+21952.D0/3.D0*t19)*t113+t257*t258
      t444 = s1+(9728.D0/3.D0-2944.D0/3.D0*t14-t182+40768.D0/3.D0*t19-22
     #016.D0/3.D0*t12+t263+1472*t20-t265)*t27+t348*s2l1-1664.D0/3.D0*t11
     #6+(t355-t358+64.D0/3.D0*t378*t364*t8)*t383+(512.D0/3.D0*t389*t390*
     #t356-256.D0/3.D0*t389*t394)*t398+((-5120-t400-3328*s1l1)*t14+(1113
     #6+25024.D0/3.D0*t116+15488*s1l1+12160*t406+60736.D0/3.D0*s2l1+1299
     #2.D0/3.D0*t258)*t16-64.D0/3.D0*(-235-278*t406-163*t116-216*s1l1-28
     #8*s2l1-115*t258)*t12+(-5056-25216.D0/3.D0*t116-48128.D0/3.D0*s1l1-
     #35072.D0/3.D0*t406-56704.D0/3.D0*s2l1-10048.D0/3.D0*t258)*t13+2560
     #*t20-512*t431*t11+(-24064.D0/3.D0-6016.D0/3.D0*t116-640.D0/3.D0*s1
     #l1-2688*t406-2816*s2l1-2048.D0/3.D0*t258)*t19)*t113-576*t16
      t445 = dlog(z)
      t450 = 2176*t19
      t453 = 2048.D0/3.D0*t14
      t475 = -t20+t14
      t477 = 64*t475*t101
      t488 = 640.D0/3.D0*t14
      t489 = 128*t16
      t496 = 640*t16
      t503 = 896*t14
      t504 = 2176.D0/3.D0*t20
      t516 = 768*t20
      t523 = 64.D0/3.D0*t16
      t524 = 64.D0/3.D0*t20
      t527 = 1664.D0/3.D0*t16
      t536 = 6848.D0/3.D0*t13
      t539 = 128.D0/3.D0*t14
      t541 = 8896.D0/3.D0*t13
      t547 = t477+5440.D0/3.D0-7616.D0/3.D0*z+(t224-26816.D0/3.D0*t19+23
     #68.D0/3.D0*t14-25792.D0/3.D0+47552.D0/3.D0*t12-6208.D0/3.D0*t16+34
     #112.D0/3.D0*z+2752*t13)*t27-896.D0/3.D0*t13+(t55-t488-t489)*t68+52
     #48.D0/3.D0*t19+(t21+76672.D0/3.D0*t19+6400.D0/3.D0*t14-8960*t12+t4
     #96-384+384*z-57472.D0/3.D0*t13)*t59+(11200.D0/3.D0*t16-2112*t13-t5
     #03-t504)*t67-9728.D0/3.D0*t12+(-t96-24704.D0/3.D0*t19-1664.D0/3.D0
     #*t14-19328.D0/3.D0*t12-6656.D0/3.D0*t16+19456.D0/3.D0-26176.D0/3.D
     #0*z+33280.D0/3.D0*t13)*t37+(t516+384*t19+3904*t14-11584*t16+6528*t
     #13)*t113-t523+(-t524-8128.D0/3.D0*t19+704-t293+704*t12+t527+768*z+
     #5632.D0/3.D0*t13)*t26+(4160*t12-4800*t14-34880.D0/3.D0*t19+31232.D
     #0/3.D0*t16+t536-t191)*t38+(t539-6208.D0/3.D0*t12-t541+2048.D0/3.D0
     #*t16-64.D0/3.D0-3776.D0/3.D0*z+11456.D0/3.D0*t19)*l1
      t553 = 5504.D0/3.D0*t14
      t554 = 1280*t16
      t562 = 1088*z
      t566 = 5312.D0/3.D0+(-20096.D0/3.D0*t20-7936.D0/3.D0*t13+42304.D0/
     #3.D0*t14-20800.D0/3.D0*t16+t450)*t113+(-t453+64.D0/3.D0*t12-1088.D
     #0/3.D0+7360.D0/3.D0*t16+4672.D0/3.D0*t19-832.D0/3.D0*z-9728.D0/3.D
     #0*t13)*l1+11008.D0/3.D0*t19+(-6400.D0/3.D0*t20+9664.D0/3.D0*t14-t3
     #27)*t68+(-43072.D0/3.D0*t19-7040*t14+4672*t20+2752*t12-22336.D0/3.
     #D0*t16+64256.D0/3.D0*t13)*t38-12160.D0/3.D0*t12+t71+t547*s2l1+(-30
     #08.D0/3.D0*t20+14656.D0/3.D0*t13-14016*t19-23552.D0/3.D0+10240*z+t
     #553-t554+52288.D0/3.D0*t12)*t27+(-96896.D0/3.D0*t13-t282+57088.D0/
     #3.D0*t16+61568.D0/3.D0*t19-1408*t14-1088+t562-8704.D0/3.D0*t12)*t5
     #9
      t568 = -t475
      t580 = 128*z
      t596 = 640.D0/3.D0*t20
      t597 = 256.D0/3.D0*t14
      t601 = 2816.D0/3.D0*t13
      t619 = 128*t19
      t637 = 4096.D0/3.D0*t13
      t642 = t477+1408-2112*z+(-t524-24448.D0/3.D0*t19-2496*t13-t198+396
     #8*t12+21056.D0/3.D0*t16)*t38+(t596-t597-t489)*t68+(t265+256-576*z-
     #t111+t601+3520.D0/3.D0*t12)*t26+704*t19+(-704.D0/3.D0*t20-t198+780
     #8.D0/3.D0*t16-2048*t13)*t67-2112*t12+(-t524-8640*t19+2624.D0/3.D0*
     #t16-t47+64*t12+12608.D0/3.D0*t13+10880.D0/3.D0-17024.D0/3.D0*z)*t3
     #7+(t224+t619-21824.D0/3.D0*t16+t453+19136.D0/3.D0*t13)*t113+(t597+
     #7680*t12+1280*t13-4352.D0/3.D0*t16-5248+8512*z-6656.D0/3.D0*t19)*t
     #27+(t45+55360.D0/3.D0*t19-2240*t16-28480.D0/3.D0*t12-20672.D0/3.D0
     #*t13+128-t580)*t59+(-t637-512.D0/3.D0-t214+2688*t19-64.D0/3.D0*z)*
     #l1
      t665 = 5248.D0/3.D0*t16
      t668 = 256*z
      t691 = 2432.D0/3.D0*t12
      t692 = 896*z
      t696 = t477-t523+(t224-9856.D0/3.D0-7424.D0/3.D0*t19+t539+16256.D0
     #/3.D0*t12-832*t16+4288*z+2752.D0/3.D0*t13)*t27+(-t191-4352*t19-t29
     #3+1920*t12+14080.D0/3.D0*t16-t637)*t38+1024.D0/3.D0*t568*t68+2176.
     #D0/3.D0-t562+64*t13+(t21+28352.D0/3.D0*t19-t665+t293-11840.D0/3.D0
     #*t12-12416.D0/3.D0*t13-256+t668)*t59+(-t504+1792.D0/3.D0*t14+t17-1
     #024*t13)*t67-3200.D0/3.D0*t12+(t539-1088*t12-t244-t41-640*z+t44)*l
     #1+(-t96-t344-t597+2688-2048*t12+704.D0/3.D0*t16-3712*z+8512.D0/3.D
     #0*t13)*t37+(t516+t619-t488-4032*t16+10048.D0/3.D0*t13)*t113+896.D0
     #/3.D0*t19+(-t524+128-t47+t691+t265+t692+t74-3584.D0/3.D0*t19)*t26
      t699 = 128*t27
      t704 = 320*l1
      t715 = 9344*t12
      t730 = 2560.D0/3.D0*t12
      t735 = 896*t20
      t746 = 27520.D0/3.D0*t19
      t748 = 2432.D0/3.D0*t16
      t753 = 640*t20
      t764 = 128*t475*t101+2176.D0/3.D0*t19+(t45-2816*t19-t665+512.D0/3.
     #D0*t14+t715+3968.D0/3.D0*t13-18944.D0/3.D0+9472*z)*t27+(-t596-t335
     #-t503+4096*t12+25600.D0/3.D0*t16-7808.D0/3.D0*t13)*t38+(1664.D0/3.
     #D0*t20-1280.D0/3.D0*t14-t489)*t68-2176*t12+(-t597+t601+t322+768-t6
     #68+t730-8576.D0/3.D0*t19)*t26+4352.D0/3.D0-2176*z+(-t735+t488+8576
     #.D0/3.D0*t16-2176*t13)*t67+(t41+8704.D0/3.D0*t19-1536*t13-4480.D0/
     #3.D0*t12-896.D0/3.D0-1408.D0/3.D0*z)*l1+(-t96-t746+13568.D0/3.D0-t
     #539-3584.D0/3.D0*t12+t748-20096.D0/3.D0*z+15616.D0/3.D0*t13)*t37+(
     #t753+t197+20480.D0/3.D0*t13+t453-25216.D0/3.D0*t16)*t113+(t96+5964
     #8.D0/3.D0*t19+t98-28288.D0/3.D0*t12-8192.D0/3.D0*t16-128+t580-2393
     #6.D0/3.D0*t13)*t59
      t775 = 448.D0/3.D0*t20
      t808 = 1024*t12
      t817 = 64*t20
      t826 = t764*s2l1+8512.D0/3.D0-3904*z+(t524-30848.D0/3.D0*t19-4736.
     #D0/3.D0*t14+28480.D0/3.D0-9984*t12-2304*t16-38272.D0/3.D0*z+43904.
     #D0/3.D0*t13)*t37+(-t775+t450+33152.D0/3.D0*t14-57152.D0/3.D0*t16+1
     #7920.D0/3.D0*t13)*t113+t477+(-t524-35776.D0/3.D0+t553+21312*t12-t3
     #13+17088*z+2624*t13-10560*t19)*t27+9536.D0/3.D0*t19+(t775+3008.D0/
     #3.D0*t14-t17)*t68-5312*t12+(t524-36032.D0/3.D0*t12+115648.D0/3.D0*
     #t19+8512.D0/3.D0*t14-896+t692+14656.D0/3.D0*t16-100160.D0/3.D0*t13
     #)*t59+(-t524-5888*t14+27136.D0/3.D0*t16-3136*t13)*t67+(4160.D0/3.D
     #0*t16-1984.D0/3.D0-320.D0/3.D0*z+13184.D0/3.D0*t19-4480*t13-t808)*
     #l1+(-1984.D0/3.D0*t14-18304.D0/3.D0*t19+16832.D0/3.D0*t13-t527+115
     #2-1216*z+t691)*t26+(t817-64064.D0/3.D0*t19-25984.D0/3.D0*t14+6208*
     #t12+10880*t16+12864*t13)*t38-704*t13
      t828 = -1728*t13+384*t568*t101+(8896*t13-40000.D0/3.D0*t12+20864.D
     #0/3.D0+3392.D0/3.D0*t20-7360*t16-25856.D0/3.D0*z+10688.D0/3.D0*t19
     #+t15)*t37+(t55-9344.D0/3.D0*t19+576-t580+t51+16448.D0/3.D0*t13+256
     #.D0/3.D0*t12-3392*t16)*t26+(-64*l1+64*t26)*t170+t642*t258+(-30208.
     #D0/3.D0*t14+17024.D0/3.D0*t16+15488.D0/3.D0*t20-t253)*t67+t696*t11
     #6-t213+((-64*t26-t699+192*l1)*s1l1-320*t26+(t704-448*t26+t699)*s2l
     #1+t704)*t171+t826*s1l1
      t833 = -t362-3*t5+2*t374-t26-l1+2
      t839 = t149*t2
      t840 = t839*t156
      t842 = t164*l6
      t853 = -896*l1-2816.D0/3.D0*t26-896.D0/3.D0-1216.D0/3.D0*t27+256*t
     #37
      t858 = 4608*t12
      t860 = 1792*t20
      t877 = 16832.D0/3.D0*t116
      t894 = 1728*s2l1
      t895 = 896.D0/3.D0*t258
      t913 = -27136.D0/3.D0-t877-41344.D0/3.D0*s1l1-10240*t406+64.D0/3.D
     #0*(-14*s2l1+t116-15)*t16+64.D0/3.D0*(38*s1l1+95*s2l1+29*t116+34*t4
     #06+67)*t13+(-1472.D0/3.D0-2752.D0/3.D0*t116-1472*s1l1-6272.D0/3.D0
     #*t406-t894-t895)*t19+(-16192.D0/3.D0-1792*t116-11264.D0/3.D0*s1l1-
     #1408.D0/3.D0*t406-18304.D0/3.D0*s2l1-t895)*t12+(31616.D0/3.D0+1580
     #8.D0/3.D0*t116+12544*s1l1+7296*t406+44032.D0/3.D0*s2l1+10112.D0/3.
     #D0*t258)*t11-14272*s2l1-15424.D0/3.D0*t258
      t922 = t390*t156
      t923 = t387*t839
      t956 = 1408*t16
      t970 = 6272.D0/3.D0*t19
      t1025 = 512*t13
      t1037 = (t860+t292-8320.D0/3.D0*t14+512*t12+19072.D0/3.D0*t16-2009
     #6.D0/3.D0*t13)*t67+(2432.D0/3.D0*t14-t248+3328*t19+896*t16+1024-t1
     #07-1280*t20-t99)*t113-640*t19+(t191+t192-17152.D0/3.D0-18304.D0/3.
     #D0*t16+t263+5248*t13-t746+44032.D0/3.D0*t12)*t38+(-3328.D0/3.D0*t2
     #0-t64+6272.D0/3.D0*t14-12416.D0/3.D0*t16+11648.D0/3.D0*t13)*t68+t8
     #08+(-t220+t17+9088.D0/3.D0*t19-1664*t13-t224+13568.D0/3.D0-t730+31
     #232.D0/3.D0*z)*t27+(-t79-t539+t956+4864.D0/3.D0*t13+35072.D0/3.D0+
     #17408.D0/3.D0*t19-19584*t12-31360.D0/3.D0*z)*t59+512.D0/3.D0+t183+
     #256.D0/3.D0*t13+(-896.D0/3.D0*t14+5248*t12+t496-15104.D0/3.D0-1753
     #6.D0/3.D0*z-2560*t13+t79-3328.D0/3.D0*t19)*t37+(t79+512*t16-512*t1
     #4-t1025)*t101+(-2048.D0/3.D0+33536.D0/3.D0*z+t220+1024.D0/3.D0*t13
     #-2176.D0/3.D0*t16+t23-t970)*t26+(t1025+9728.D0/3.D0+4352.D0/3.D0*t
     #12-t218+20992.D0/3.D0*z)*l1
      t1059 = (-t753-t197+4544.D0/3.D0*t14-3840*t16+9664.D0/3.D0*t13)*t6
     #8+(704*t20+4288.D0/3.D0*t19+3584.D0/3.D0*t14+2624*t16-5952*t13)*t6
     #7+(-t79-15808.D0/3.D0*t14+512-t183-t858+7040.D0/3.D0*t16+26176.D0/
     #3.D0*t13-t72)*t113-1600.D0/3.D0+1408*z+(t817+2560.D0/3.D0*t14+t956
     #-5184*t13-22144.D0/3.D0-17152.D0/3.D0*t19+12480*t12-4416*z)*t37+(t
     #817+4096*t14-16256.D0/3.D0+t265+18560.D0/3.D0*z-48320.D0/3.D0*t13+
     #2240*t19+14528*t12)*t38-t970+2816*t12+(-2752.D0/3.D0*t14-108608.D0
     #/3.D0*t12-13504.D0/3.D0*t16-t96+51392.D0/3.D0-49664.D0/3.D0*z+3616
     #0.D0/3.D0*t13+38464.D0/3.D0*t19)*t59+t1037*s2l1+(-1728*t19+t536+56
     #32.D0/3.D0+11904*z+4032*t12-t748)*l1+(-3392.D0/3.D0*t14+2240.D0/3.
     #D0*t16+8576.D0/3.D0*t13-64-10048.D0/3.D0*t19+t715+70144.D0/3.D0*z)
     #*t27+t324+1472.D0/3.D0*t13+(-6464*t19+t30+61312.D0/3.D0*t12-6400.D
     #0/3.D0*t13-18368.D0/3.D0+24256*z+1984.D0/3.D0*t16)*t26
      t1062 = 2*t406
      t1068 = 8*s1l1
      t1077 = -t445*(t566+t828)+(-t355+t358-128.D0/3.D0*t833*t2*t8)*l5+(
     #-512.D0/3.D0*t840*t153*t842+256.D0/3.D0*t840*t394)*t151+t853*t170+
     #(t554-t99-t735-t220)*t101+512.D0/3.D0*s1l1+(-512-t314+t858-6720*t1
     #9-t860-5312.D0/3.D0*t14+t185-23104.D0/3.D0*z)*t37+((1088.D0/3.D0+t
     #121-t119+4416*t26+15616.D0/3.D0*t27)*s1l1+(2624.D0/3.D0+t128+t127-
     #t125+7360.D0/3.D0*t26+5248.D0/3.D0*t27)*s2l1+2752.D0/3.D0+t124+t13
     #5-t132+13120.D0/3.D0*t26+11648.D0/3.D0*t27)*t171+t913*l1+(-7936*t2
     #0+5312*t14+t91-2624.D0/3.D0*t19+27712.D0/3.D0*t16-17984.D0/3.D0*t1
     #3)*t67+(-512.D0/3.D0*t922*t923*t356+256.D0/3.D0*t922*t923*t8)*t398
     #+(-3968.D0/3.D0*t14-16256.D0/3.D0*t16+4032*t20+t541-t197)*t68+t105
     #9*s1l1+(-t79-256*(5+3*s2l1+t1062+3*s1l1+t116+t258)*t16+64*(20+8*s2
     #l1+t1068)*t14-256*(-1-t116-s1l1-t1062-s2l1-t258)*t13)*t101
      t1120 = -s2l1-1
      t1133 = -9536-2688*t116-18816*s1l1-33536.D0/3.D0*t406-64.D0/3.D0*(
     #18*s2l1+17*t116+3+t1068+18*t406)*t16+(34816.D0/3.D0+12992.D0/3.D0*
     #t116+39616.D0/3.D0*s1l1+17408.D0/3.D0*t406+14080*s2l1+8192.D0/3.D0
     #*t258)*t19+(-4160-3584.D0/3.D0*t116-8384.D0/3.D0*s1l1-2560.D0/3.D0
     #*t406-13952.D0/3.D0*s2l1-384*t258)*t13+(-44800.D0/3.D0-4224*t116-2
     #3616*s1l1-31744.D0/3.D0*t406-51712.D0/3.D0*s2l1-22912.D0/3.D0*t258
     #)*t12-512.D0/3.D0*t1120*t14+(40448.D0/3.D0+2624*t116+25280*s1l1+36
     #608.D0/3.D0*t406+44416.D0/3.D0*s2l1+31232.D0/3.D0*t258)*t11-35200.
     #D0/3.D0*s2l1-26368.D0/3.D0*t258
      t1180 = -37696.D0/3.D0-14464.D0/3.D0*t116-80960.D0/3.D0*s1l1-53504
     #.D0/3.D0*t406+64.D0/3.D0*(-22-23*s2l1-t1068)*t14+(-46720.D0/3.D0-8
     #96*t116-27264*s1l1-25600.D0/3.D0*t406-11456*s2l1-31808.D0/3.D0*t25
     #8)*t12+(58496.D0/3.D0+17344.D0/3.D0*t116+37120*s1l1+58624.D0/3.D0*
     #t406+21568*s2l1+16704*t258)*t11+(-4352-5696*s1l1-3200.D0/3.D0*t406
     #-2368*s2l1-2752.D0/3.D0*t258)*t13+(8384.D0/3.D0+3520.D0/3.D0*t116+
     #1984*s1l1+3968.D0/3.D0*t406+t400+1088.D0/3.D0*t258)*t16+(21568.D0/
     #3.D0-7744.D0/3.D0*t116+12736*s1l1+2432.D0/3.D0*t406+1024*s2l1+1222
     #4.D0/3.D0*t258)*t19-48704.D0/3.D0*s2l1-41536.D0/3.D0*t258
      t1211 = 576*s1l1
      t1228 = (23104+42368.D0/3.D0*t116+141056.D0/3.D0*s1l1+30720*t406+1
     #18016.D0/3.D0*s2l1+52480.D0/3.D0*t258)*t12+(-2368.D0/3.D0*s2l1+224
     #0.D0/3.D0-704*s1l1)*t14+(-2304-4864.D0/3.D0*t116-43904.D0/3.D0*s1l
     #1-28672.D0/3.D0*t406-7424*s2l1-7936*t258)*t11+384*t20+(-107776.D0/
     #3.D0-56128.D0/3.D0*t116-48704*s1l1-27520*t406-154688.D0/3.D0*s2l1-
     #34496.D0/3.D0*t258)*t19-4928.D0/3.D0-3712.D0/3.D0*t406-3392.D0/3.D
     #0*t116-t1211-1856*s2l1-320.D0/3.D0*t258+(58432.D0/3.D0+11776.D0/3.
     #D0*t116+50816.D0/3.D0*s1l1+12416.D0/3.D0*t406+54976.D0/3.D0*s2l1+3
     #776.D0/3.D0*t258)*t13+(-4096+9280.D0/3.D0*t116+1216.D0/3.D0*s1l1+3
     #072*t406+3200*s2l1+3328.D0/3.D0*t258)*t16
      t1283 = (6016.D0/3.D0*s2l1+5504.D0/3.D0+6400.D0/3.D0*s1l1)*t14-t10
     #5+(25856+46400.D0/3.D0*t116+59072.D0/3.D0*s1l1+18944*t406+109760.D
     #0/3.D0*s2l1+5824.D0/3.D0*t258)*t19-256*t116-512*t406-256-768*s1l1-
     #256*t258-768*s2l1+(-24704.D0/3.D0+1408*t116+14464.D0/3.D0*s1l1+832
     #0.D0/3.D0*t406-t894+9344.D0/3.D0*t258)*t13+(-50176.D0/3.D0-41728.D
     #0/3.D0*t116-22016*s1l1-62464.D0/3.D0*t406-30592*s2l1-6144*t258)*t1
     #2+(-10624.D0/3.D0-19520.D0/3.D0*t116-10048*s1l1-21376.D0/3.D0*t406
     #-13312*s2l1-2560*t258)*t16-128.D0/3.D0*(-61-146*t406-92*t116-145*s
     #1l1-183*s2l1-54*t258)*t11
      t1309 = -4928.D0/3.D0+(-38080.D0/3.D0*t19-7040.D0/3.D0+14016*t16-6
     #272*t20-23872.D0/3.D0*t13+7808.D0/3.D0*z+2624*t14+37888.D0/3.D0*t1
     #2)*t38-t853*t170+t1133*t26+2496*z+1408.D0/3.D0*t258+12928.D0/3.D0*
     #t12+t1180*t27+512.D0/3.D0*t406+t1228*t59+(-3520+37888.D0/3.D0*z-51
     #2*t20-7104*t19+38144.D0/3.D0*t12+10240.D0/3.D0*t16-896*t13-2240.D0
     #/3.D0*t14)*t26+t1283*t38+(256*t431*t19+64.D0/3.D0*(-139-182*t406-9
     #1*t116-151*s1l1-163*s2l1-91*t258)*t13+(17792.D0/3.D0+5632.D0/3.D0*
     #t116+4736*s1l1+12800.D0/3.D0*t406+4608*s2l1+7168.D0/3.D0*t258)*t16
     #+1152*t20+(-13120.D0/3.D0-5696.D0/3.D0*s2l1-5312.D0/3.D0*s1l1)*t14
     #)*t68
      t1331 = 39616.D0/3.D0*s2l1
      t1347 = t2*z
      t1349 = t149*t1347*t156
      t1357 = t1347*t144
      t1359 = 1024.D0/3.D0*t1357*t353
      t1361 = 1024.D0/3.D0*t1357*t356
      t1440 = 32320.D0/3.D0+25088.D0/3.D0*t116+44800.D0/3.D0*s1l1+13440*
     #t406+(11776+8576.D0/3.D0*t116+704.D0/3.D0*s1l1-5248.D0/3.D0*t406+t
     #1331-20480.D0/3.D0*t258)*t12+(-7936.D0/3.D0-3328.D0/3.D0*t116-4544
     #*s1l1-2816.D0/3.D0*t406-8320.D0/3.D0*s2l1-4352.D0/3.D0*t258)*t13+(
     #-320+512*s2l1+320*s1l1)*t14+(-58496.D0/3.D0-36928.D0/3.D0*t116-178
     #56*s1l1-14464*t406-91072.D0/3.D0*s2l1-3008.D0/3.D0*t258)*t11+(1728
     #+t877+32896.D0/3.D0*s1l1+25088.D0/3.D0*t406+7296*s2l1+18880.D0/3.D
     #0*t258)*t19+(128-5312.D0/3.D0*t116-t1211-2048*t406-8576.D0/3.D0*s2
     #l1-1856.D0/3.D0*t258)*t16+18752*s2l1+14080.D0/3.D0*t258
      t1450 = (256*(-t116-t1062-1-t258)*t12+64.D0/3.D0*(41-86*t406-74*t1
     #16-67*s1l1-129*s2l1-12*t258)*t19+(-10688-t877-33152.D0/3.D0*s1l1-1
     #0752*t406-13568*s2l1-14272.D0/3.D0*t258)*t16-2304*t20+(18368.D0/3.
     #D0+19456.D0/3.D0*t116+9472*s1l1+30848.D0/3.D0*t406+t1331+10240.D0/
     #3.D0*t258)*t13+(18752.D0/3.D0+9344.D0/3.D0*s2l1+3008*s1l1)*t14)*t6
     #7-3328.D0/3.D0*s2l1-64.D0/3.D0*t1120*(-t429+6)*t13+(512.D0/3.D0*t1
     #349*t153*t842-256.D0/3.D0*t1349*t394)*t151+(t1359-t1361-64.D0/3.D0
     #*t2*t378*t8)*t383+t215-64.D0/3.D0*t1347*l1*t8+(-t1359+t1361+128.D0
     #/3.D0*t833*t1347*t8)*l5+64.D0/3.D0*(-37-44*s2l1-7*t116-20*s1l1-14*
     #t406)*t19-64.D0/3.D0*(46+39*t258+26*t116+104*s1l1+92*t406+72*s2l1)
     #*t11+64.D0/3.D0*(74+34*t116+84*s1l1+68*t406+108*s2l1+14*t258)*t12+
     #t1440*t37+(-1152+30592.D0/3.D0*z+1024*t14+17152.D0/3.D0*t12+t207-6
     #784.D0/3.D0*t19-2816*t16)*l1-4544*t19
      t1454 = (l1+z-t5)**2
      taux17 = (t444+t1077+t1309+t1450)/t1454*t8
      t1 = 1-z
      t2 = t1**2
      t3 = t2**2
      t4 = t3*t1
      t6 = t3*t2
      t8 = t2*t1
      t9 = t3*t8
      t16 = l1**2
      t17 = t16*l1
      t19 = s1l1+1+s2l1
      t33 = t16**2
      t34 = t33**2
      t35 = t34*t16
      t49 = 139+115*s2l1+151*s1l1
      t62 = t33*t17
      t70 = 7936.D0/3.D0*t8
      t80 = 16+16*s2l1
      t84 = 16*s1l1
      t117 = t34*l1
      t127 = t33*l1
      t191 = t33*t16
      t193 = 5120*t4
      s2 = 44800.D0/3.D0+(-180608.D0/3.D0+2432.D0/3.D0*t4+60800.D0/3.D0*
     #t6-15616*t9-178688*t8+507904.D0/3.D0*t2+64512*t3+216064.D0/3.D0*z)
     #*t17+(-1024*t19*t3-512*(-6*s1l1-10-6*s2l1)*t4+(-5632-2304*s1l1-230
     #4*s2l1)*t6+1536*t9)*t35+((-23936.D0/3.D0-58880.D0/3.D0*s1l1-2048.D
     #0/3.D0*s2l1)*t4+2048*t19*t1+(152704.D0/3.D0+32896*s1l1+19584*s2l1)
     #*t8-512.D0/3.D0*t49*t2-16896*t9+(-30848+44416.D0/3.D0*s1l1+1024.D0
     #/3.D0*s2l1)*t3+(26496+7808*s2l1+22016*s1l1)*t6)*t62
      s1 = s2+(-148096.D0/3.D0*t4+156032.D0/3.D0*t3-16896*t2+15616.D0/3.
     #D0*t9+5376*t6+15488.D0/3.D0-28672.D0/3.D0*z-t70)*t16+((39040.D0/3.
     #D0+15872.D0/3.D0*s1l1+6144*s2l1)*t3+13184.D0/3.D0+13184.D0/3.D0*s1
     #l1+9728.D0/3.D0*s2l1-128.D0/3.D0*t80*t6-128.D0/3.D0*(51+19*s2l1+t8
     #4)*t4+(-28672.D0/3.D0-8960*s1l1-7936*s2l1)*t1+(17152+50432.D0/3.D0
     #*s1l1+39680.D0/3.D0*s2l1)*t2+(-20992-47104.D0/3.D0*s1l1-12032*s2l1
     #)*t8)*t16-51968.D0/3.D0*z+((-71168.D0/3.D0-11648*s1l1-8064*s2l1)*t
     #4+(60160.D0/3.D0+23168.D0/3.D0*s1l1+20864.D0/3.D0*s2l1)*t6-7168*t9
     #+256.D0/3.D0*t49*t3-1024*t19*t8)*t117
      t196 = s1+(-91904.D0/3.D0*t9-71680.D0/3.D0*t6+805888.D0/3.D0*t4+24
     #832.D0/3.D0-963328.D0/3.D0*t3-60416.D0/3.D0*z+269312.D0/3.D0*t8+79
     #36*t2)*t127+2560.D0/3.D0*s2l1+(19200*t9+82048.D0/3.D0-366592.D0/3.
     #D0*t4+8320*t6-41216.D0/3.D0*z+68352*t8-97792*t2+301312.D0/3.D0*t3)
     #*t33+(2048*t19*t2+(-37760.D0/3.D0-27904*s1l1-13056*s2l1)*t3+256.D0
     #/3.D0*(-77-29*s2l1-29*s1l1)*t8+(33536+57728.D0/3.D0*s1l1+16384.D0/
     #3.D0*s2l1)*t4+(-93824.D0/3.D0-26624.D0/3.D0*s2l1-15488*s1l1)*t6+14
     #848*t9)*t34+1408.D0/3.D0*s1l1+((-110720.D0/3.D0-2816.D0/3.D0*s1l1-
     #4480*s2l1)*t4-512.D0/3.D0*(13+25*s2l1+25*s1l1)*t1+10240*t9+(-15104
     #.D0/3.D0-14080.D0/3.D0*s2l1-49408.D0/3.D0*s1l1)*t6+(-85376-224128.
     #D0/3.D0*s1l1-92288.D0/3.D0*s2l1)*t8+(116992.D0/3.D0+53504*s1l1+299
     #52*s2l1)*t2-1024*s1l1-1024-1024*s2l1+(243968.D0/3.D0+34432*s1l1+45
     #824.D0/3.D0*s2l1)*t3)*t191+t193-75776.D0/3.D0*t3-41472*t2
      t219 = 512*s2l1
      t252 = 141440.D0/3.D0*t8
      t267 = 2560*t9
      t293 = 1024*t8
      t299 = 1024*t3
      t320 = 12800.D0/3.D0*z
      t326 = 1280*t9
      t339 = 2048*z
      t345 = 29696.D0/3.D0*t3
      t346 = 2048*t2
      t350 = 7424.D0/3.D0*t8
      t355 = 3072*t4
      t357 = (t326-1280*t6+t299-t355)*t35
      t383 = 1280.D0/3.D0*t9
      t392 = (-16640.D0/3.D0*t9+18560*t6-38656.D0/3.D0*t3+t293+2944*t4)*
     #t117+(6400.D0/3.D0*t9+15872.D0/3.D0-319232.D0/3.D0*t6+74240.D0/3.D
     #0*t4-t320+492928.D0/3.D0*t3-201856.D0/3.D0*t8-27392*t2)*t191+(t326
     #-53120*t4-59008.D0/3.D0*t6+134528.D0/3.D0-98432.D0/3.D0*t8+141440*
     #t3-113792.D0/3.D0*z-234112.D0/3.D0*t2)*t33-44416*t2+(-6400*t9+3934
     #72.D0/3.D0*t6+77312.D0/3.D0*t2-2048+t339-85120.D0/3.D0*t3-306176.D
     #0/3.D0*t4-18560*t8)*t62-t345+(-t346+8960*t9-74880*t6+129920.D0/3.D
     #0*t4+t350+21760*t3)*t34+t357+103424.D0/3.D0*t8+(58112.D0/3.D0*t4-1
     #30304.D0/3.D0*t3+8704.D0/3.D0+89600.D0/3.D0*t8-30592.D0/3.D0*t2-30
     #592.D0/3.D0*z)*l1+19456-76160.D0/3.D0*z+(-t326+115840.D0/3.D0*t6+9
     #8944*t4+9344.D0/3.D0-923008.D0/3.D0*t3-16000*z+199808*t8-87424.D0/
     #3.D0*t2)*t127+(-27136.D0/3.D0*t6+153856.D0/3.D0*t3-75776.D0/3.D0*t
     #8-45824.D0/3.D0*t4-24320.D0/3.D0*t2+13696.D0/3.D0-8960*z)*t16+(-t3
     #83+68864.D0/3.D0*t6-234368.D0/3.D0+302080.D0/3.D0*z-47360.D0/3.D0*
     #t4+173568*t2+67328.D0/3.D0*t3-371200.D0/3.D0*t8)*t17
      t394 = dlog(z)
      t427 = 256*t9
      t439 = 1280*z
      t456 = -t9+t6
      t458 = 256*t456*t35
      t465 = 256.D0/3.D0*t9
      t468 = 6656*t6
      t490 = 7936.D0/3.D0*t6
      t493 = 39424.D0/3.D0*t8
      t499 = (-t427+8704.D0/3.D0*t8+115712.D0/3.D0*t2-21248+16384.D0/3.D
     #0*t6+41984.D0/3.D0*t4+25600*z-39424*t3)*t33+(-16640.D0/3.D0*t4-358
     #4.D0/3.D0*t2+37376.D0/3.D0*t3+1280-t439-20992.D0/3.D0*t8)*l1+(-t38
     #3-173312.D0/3.D0*t3+96256.D0/3.D0*t6+48896*t8-8448*t2-14336*t4)*t1
     #91+(-1792*t9+3840*t3+20992*t6-23040*t4)*t34+t458+(t326+4096.D0/3.D
     #0*t3-116224.D0/3.D0*t6-7168*t8+43264*t4)*t62+(t465+101632.D0/3.D0*
     #t8+67328.D0/3.D0-21760.D0/3.D0*t3-t468+6400*t4-48896*t2-28416*z)*t
     #17-16640.D0/3.D0+6912*z+(t427+101632*t3-34048.D0/3.D0*t6-224512.D0
     #/3.D0*t8+27904.D0/3.D0*t2+3584-3584*z-85760.D0/3.D0*t4)*t127+12544
     #*t2+(3328.D0/3.D0*t9-14080.D0/3.D0*t6+3584*t4)*t117+2816*t3+(t490-
     #512+12800.D0/3.D0*t4+768*z+t493-5632.D0/3.D0*t2-17664*t3)*t16-2944
     #0.D0/3.D0*t8
      t518 = 1024*z
      t520 = 60928.D0/3.D0*t3
      t534 = 512.D0/3.D0*t8
      t572 = 5120.D0/3.D0*t9
      t596 = (-512*t9-14336.D0/3.D0+3584.D0/3.D0*t6-t534+38912.D0/3.D0*t
     #2+15104.D0/3.D0*z+34816.D0/3.D0*t4-t520)*t33+(512.D0/3.D0*t6-t70-t
     #355+256.D0/3.D0*t2-256.D0/3.D0+256.D0/3.D0*z+16640.D0/3.D0*t3)*l1+
     #13312.D0/3.D0*t2+(-8704.D0/3.D0*t9-256*t3+t468-10496.D0/3.D0*t4)*t
     #34+t458+(3072*t9+3584.D0/3.D0*t3-46336.D0/3.D0*t6+33536.D0/3.D0*t4
     #)*t62+(1024.D0/3.D0*t9+21248.D0/3.D0-58624.D0/3.D0*t2-5376*t3-t490
     #+55040.D0/3.D0*t8+1792*t4-22784.D0/3.D0*z)*t17-256.D0/3.D0*t4-4352
     #.D0/3.D0+4352.D0/3.D0*z+5120.D0/3.D0*t3-4096.D0/3.D0*t456*t117+(-t
     #572-56576.D0/3.D0*t3+44288.D0/3.D0*t6+26368.D0/3.D0*t8-256*t2-8192
     #.D0/3.D0*t4)*t191+(-t465-768+4864.D0/3.D0*t2-16384.D0/3.D0*t3+1024
     #*t6+4352.D0/3.D0*t8+6656.D0/3.D0*t4+t518)*t16+(2048.D0/3.D0*t9+125
     #440.D0/3.D0*t3-4608*t6-63488.D0/3.D0*t8+2048.D0/3.D0*t2-17408*t4)*
     #t127-4608*t8
      t602 = 8704.D0/3.D0*t6
      t609 = (-5120*t9-23552*t3-17920-97792.D0/3.D0*t8+146944.D0/3.D0*t2
     #-4096.D0/3.D0*t6+94720.D0/3.D0*t4+20224*z)*t33+(-20480*t9-89600*t3
     #+75008.D0/3.D0*t6+37376*t8-2560*t2+150784.D0/3.D0*t4)*t191+(4352*t
     #9+53248.D0/3.D0+146432.D0/3.D0*t8-144896.D0/3.D0*t2-17408*t3-20480
     #.D0/3.D0*t6+4864.D0/3.D0*t4-20992*z)*t17+(-73472.D0/3.D0*t9-1792*t
     #3+124160.D0/3.D0*t6-15104*t4)*t34-12544.D0/3.D0+t499*s1l1+(10752*t
     #9-14592*t6+3840*t4)*t117+(24832.D0/3.D0*t9+332032.D0/3.D0*t3+22528
     #.D0/3.D0*t6-119552.D0/3.D0*t8-21248.D0/3.D0*t2+3840-3840*z-250112.
     #D0/3.D0*t4)*t127+(-4352.D0/3.D0*t9+46592.D0/3.D0*t4-1536*t6-768+t5
     #18+2816.D0/3.D0*t2-t520+22784.D0/3.D0*t8)*t16+(30208*t9+74752.D0/3
     #.D0*t3-54528*t6-7680*t8+21248.D0/3.D0*t4)*t62+7168*t3+t596*s2l1-40
     #448.D0/3.D0*t8+2048*t456*t35+(t345-t534-9984*t4+t602-11776.D0/3.D0
     #*t2+1280-t439)*l1-4352.D0/3.D0*t4+4864*z+35840.D0/3.D0*t2
      t695 = 9472.D0/3.D0-39424.D0/3.D0*t2+(t267-5120.D0/3.D0*t6-158464.
     #D0/3.D0*t4+81920*t3-26368.D0/3.D0+76288.D0/3.D0*z+4352.D0/3.D0*t2-
     #56320.D0/3.D0*t8)*t33+t357+(-2048.D0/3.D0*t4+8064*t6-18944+58880.D
     #0/3.D0*z-t572+19456*t3-205312.D0/3.D0*t8+189952.D0/3.D0*t2)*t17+(-
     #20480.D0/3.D0*t9+3968*t6-29440.D0/3.D0*t3+t293+13696*t4)*t117+(-27
     #008.D0/3.D0*t4-t602+t383+14080.D0/3.D0-49408.D0/3.D0*t2+10752*t3-7
     #936*z+34048.D0/3.D0*t8)*t16-3328*z+(-2816.D0/3.D0*t6+16000.D0/3.D0
     #+35840.D0/3.D0*t4-60800.D0/3.D0*t3-31360.D0/3.D0*z+t493-32128.D0/3
     #.D0*t2)*l1+(58880.D0/3.D0*t2-15360*t9+135040.D0/3.D0*t6-2048+t339+
     #9728.D0/3.D0*t4-33408*t8-47872.D0/3.D0*t3)*t62+1024.D0/3.D0*t4+(-t
     #346+43520.D0/3.D0*t9-51712.D0/3.D0*t6-64256.D0/3.D0*t4+t350+25344*
     #t3)*t34-6144*t3+(25600.D0/3.D0*t9+15872.D0/3.D0-45824*t6-40064.D0/
     #3.D0*t4-t320+103168.D0/3.D0*t3+t252-39680*t2)*t191+(-10240.D0/3.D0
     #*t9+33920.D0/3.D0+38528.D0/3.D0*t6+213376.D0/3.D0*t4-21120*z-12070
     #4*t3+30080*t8-7040.D0/3.D0*t2)*t127+46208.D0/3.D0*t8
      s3 = (-62336.D0/3.D0-23296*s1l1-14592*s2l1+(-199168.D0/3.D0-148864
     #.D0/3.D0*s1l1-77056.D0/3.D0*s2l1)*t2+(-16640.D0/3.D0-22144.D0/3.D0
     #*s1l1-6656.D0/3.D0*s2l1)*t4+(-13824+9856.D0/3.D0*s1l1-14336.D0/3.D
     #0*s2l1)*t3+(149248.D0/3.D0+63616.D0/3.D0*s1l1+17408*s2l1)*t8+(4940
     #8+145280.D0/3.D0*s1l1+73216.D0/3.D0*s2l1)*t1+(-t219-896*s1l1+3712)
     #*t6)*t33+(-6784.D0/3.D0-2176.D0/3.D0*s1l1-10240.D0/3.D0*s2l1+128.D
     #0/3.D0*(53+37*s2l1+t84)*t6+(-24704.D0/3.D0-6400.D0/3.D0*s1l1-3584*
     #s2l1)*t4+(39680.D0/3.D0+29696.D0/3.D0*s1l1+24320.D0/3.D0*s2l1)*t3+
     #(-39680.D0/3.D0-36352.D0/3.D0*s1l1-7168*s2l1)*t8+(20480.D0/3.D0+10
     #752*s1l1-t219)*t2+(1024-13312.D0/3.D0*s1l1+4608*s2l1)*t1)*t17
      s2 = s3+t252+(35584.D0/3.D0+29440.D0/3.D0*s2l1+38656.D0/3.D0*s1l1+
     #(82688.D0/3.D0+51840*s1l1+40832.D0/3.D0*s2l1)*t8+(-82688.D0/3.D0-2
     #1632*s1l1-11904*s2l1)*t1+(-174080.D0/3.D0-152960.D0/3.D0*s1l1-1702
     #4*s2l1)*t3-t267+(15104-40832.D0/3.D0*s1l1-5504.D0/3.D0*s2l1)*t2+(4
     #4544+20096*s1l1+29824.D0/3.D0*s2l1)*t4+(-29696.D0/3.D0+2176.D0/3.D
     #0*s2l1+4736*s1l1)*t6)*t127
      s1 = s2+(-37248*t3+12416.D0/3.D0-33280.D0/3.D0*z-13568.D0/3.D0*t2+
     #32896.D0/3.D0*t8+106496.D0/3.D0*t4-30976.D0/3.D0*t6)*l1+(-44288*t9
     #+126464.D0/3.D0*t6-35584.D0/3.D0*t3+t293+38912.D0/3.D0*t4)*t117+(9
     #216*t9-5120*t6+t299-t193)*t35+(71424*t9+9728.D0/3.D0-223744.D0/3.D
     #0*t6-509824.D0/3.D0*t4-6656.D0/3.D0*z+217856*t3-29312.D0/3.D0*t8-1
     #14688.D0/3.D0*t2)*t191+t392*s1l1
      t740 = s1-t394*t609+(-106752*t9+474496.D0/3.D0*t6+71168.D0/3.D0*t2
     #-2048+t339-210304.D0/3.D0*t3+89216.D0/3.D0*t4-98176.D0/3.D0*t8)*t6
     #2-128.D0/3.D0*(60+60*s2l1+39*s1l1)*t2+t695*s2l1+(-128.D0/3.D0*(79-
     #12*s2l1+14*s1l1)*t8+20864.D0/3.D0+5120*s2l1+7296*s1l1+256.D0/3.D0*
     #t80*t4-128.D0/3.D0*(57+73*s2l1+t84)*t3+128.D0/3.D0*(294+117*s2l1+1
     #77*s1l1)*t2-128.D0/3.D0*(318+173*s2l1+283*s1l1)*t1)*l1+128.D0/3.D0
     #*(55+55*s2l1+t84)*t8+128.D0/3.D0*(10+10*s2l1+21*s1l1)*t1+(-t346+27
     #6736.D0/3.D0*t9-361088.D0/3.D0*t6-6144*t4+19712.D0/3.D0*t8+89216.D
     #0/3.D0*t3)*t34-128.D0/3.D0*t80*t3
      t742 = z*l1
      t744 = (l1+z-t742)**2
      t748 = (1-l1+t742)**2
      taux18 = (t196+t740)/t744/t748

CCC Coefficients

      tcoff33 = taux17
      tcoff34 = taux17
      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet33 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff33)
      tFjet34 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff34)
      tFjet35 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff35)
      tFjet36 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff36)

CCC Sum

      qgRVnum9e1 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFjet
     #36*tcoff36
      RETURN
      END

      double precision function qgRVnum9e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -1+l1
      t2 = z**2
      t3 = t2**2
      t4 = t1*t3
      t5 = z*l1
      t7 = (1-l1+t5)**2
      t8 = 1/t7
      t9 = l6**2
      t10 = t8*t9
      t12 = 1024.D0/3.D0*t4*t10
      t13 = t8*l6
      t15 = 1024.D0/3.D0*t4*t13
      t16 = l1**2
      t17 = t2*t16
      t19 = z*t16
      t21 = t17+3*t5-2*t19+t16+l1-2
      t28 = t16**2
      t30 = t28*l1
      t31 = 256*t30
      t32 = t16*l1
      t37 = 512*t30
      t50 = 2*l6-1
      t51 = t50**2
      t54 = l5-1
      t60 = dsqrt(z)
      t62 = dsqrt(l5)
      t64 = dsqrt(-t54)
      t67 = dsqrt(l1)
      t71 = l5*z*l1-l1+1-l5+l5*l1
      t74 = dsqrt(-t1)
      t75 = t8*t74
      t79 = t71**2
      t81 = 4*l1*t51*z*l5*t54*t1*t8-4*t50*t60*t62*t64*t67*t71*t75+t79*t8
      t82 = dsqrt(t81)
      t84 = t3*z
      t88 = t84*t1
      t90 = 1024.D0/3.D0*t88*t10
      t92 = 1024.D0/3.D0*t88*t13
      t98 = 1-z
      t99 = t98**2
      t100 = t99**2
      t101 = t100*t99
      t102 = 512.D0/3.D0*t101
      t103 = t100*t98
      t104 = 128.D0/3.D0*t103
      t106 = t99*t98
      t109 = t100*t106
      t110 = 256*t109
      t115 = 38144.D0/3.D0*t103
      t118 = 2560.D0/3.D0*t109
      t122 = t28*t16
      t124 = 256*t101
      t127 = 15872.D0/3.D0*t106
      t129 = 512.D0/3.D0*t109
      t130 = 7808*z
      t133 = 2560.D0/3.D0*t101
      t134 = 1664.D0/3.D0*t103
      t136 = 1024.D0/3.D0*t109
      t144 = 7168.D0/3.D0*t106
      t145 = 1024.D0/3.D0*t103
      t147 = 128.D0/3.D0*t101
      t152 = 384*t101
      t153 = 128*t109
      t154 = 512*t103
      t155 = 256*t100
      t157 = t28**2
      t158 = t157*t16
      t159 = (t152-t153-t154+t155)*t158
      t161 = 1408.D0/3.D0*t100
      t163 = 2048.D0/3.D0*t109
      t165 = 896.D0/3.D0*t100
      t166 = 256*t106
      t168 = t157*l1
      t173 = 896.D0/3.D0*t101
      t175 = 17024.D0/3.D0*z
      t176 = 128.D0/3.D0*t109
      t179 = 6400.D0/3.D0*t101
      t181 = 4352.D0/3.D0*t109
      t186 = 512*z
      t189 = 1536*t109
      t193 = t28*t32
      t195 = (t102+t104+57088.D0/3.D0+13952.D0/3.D0*t100+10240.D0/3.D0*t
     #106-25728*t99-t110-71680.D0/3.D0*z)*t28+(-2048*t101+t115+45568.D0/
     #3.D0*t99+22400.D0/3.D0*t100-t118+1024.D0/3.D0-3328.D0/3.D0*z-98432
     #.D0/3.D0*t106)*t122+2304+(t124-3328*t103-32384.D0/3.D0+3968*t100-t
     #127+5376*t99+t129-t130)*t32+(t133-t134-61312.D0/3.D0*t100+t136+267
     #52*t106+27392.D0/3.D0*t99-47488.D0/3.D0+72448.D0/3.D0*z)*t30+6016.
     #D0/3.D0*t106+(-11648.D0/3.D0-7040.D0/3.D0*t100+t144+t145+1792.D0/3
     #.D0*t99+t147-2944.D0/3.D0*z)*l1-3456*t99+t159+t104-6784.D0/3.D0*z-
     #t161+(-1664*t101+t163+2816.D0/3.D0*t103+t165-t166)*t168+(2944.D0/3
     #.D0*t103+9088.D0/3.D0*t100-18688.D0/3.D0*t106-t173+1280.D0/3.D0*t9
     #9-4096.D0/3.D0-t175-t176)*t16+(t179+3712.D0/3.D0*t106-t181+16000.D
     #0/3.D0*t103-21760.D0/3.D0*t100)*t157+(t102+512-t186+10880*t100-160
     #00*t103+t189+25600.D0/3.D0*t106-5632*t99)*t193
      t197 = dlog(z)
      t199 = 2*t5
      t202 = t2*z
      t213 = t2*l1
      t215 = -4*z-t199+l1-2*t16-14*t17+6*t202*t16-4*z*t32+6*t2*t32-4*t20
     #2*t32+t3*t32+10*t19+9*t213+t32
      t219 = t90-t92+64.D0/3.D0*t3*t215*t8
      t220 = l5**2
      t235 = 10240.D0/3.D0*t103
      t240 = 640*t109
      t254 = 4608*t99
      t257 = 4352.D0/3.D0*t106
      t261 = 64*t109
      t262 = 1408*t103
      t270 = 1792.D0/3.D0*t103
      t277 = 6272.D0/3.D0*t106
      t287 = 1792*t109
      t288 = 5504.D0/3.D0*t106
      t290 = 512*t99
      t291 = 19072.D0/3.D0*t103
      t296 = 17792.D0/3.D0*t99
      t299 = 1024*z
      t303 = 640*t106
      t304 = 1280.D0/3.D0*t109
      t305 = 2176.D0/3.D0*t101
      t307 = 18688.D0/3.D0*z
      t309 = 27520.D0/3.D0*t106
      t314 = 512*t106
      t317 = 11648.D0/3.D0*t100
      t320 = 1024*t99
      t321 = 128*t101
      t322 = 1152*t103
      t323 = 9088.D0/3.D0*t106
      t325 = 256.D0/3.D0*t109
      t326 = 2560.D0/3.D0*t99
      t338 = 640*t103
      t344 = 512*t101
      t345 = 512*t100
      t350 = 2176.D0/3.D0*t103
      t351 = 23552.D0/3.D0*t99
      t355 = 1280.D0/3.D0*t106
      t359 = (t287+t288-8320.D0/3.D0*t101+t290+t291-20096.D0/3.D0*t100)*
     #t157+(2432.D0/3.D0*t101-t296+3328*t106+896*t103+1024-t299-1280*t10
     #9-t155)*t193-t303+(t304+t305-17152.D0/3.D0-18304.D0/3.D0*t103+t307
     #+5248*t100-t309+44032.D0/3.D0*t99)*t122+(-3328.D0/3.D0*t109-t314+6
     #272.D0/3.D0*t101-12416.D0/3.D0*t103+t317)*t168+t320+(-t321+t322+t3
     #23-1664*t100-t325+13568.D0/3.D0-t326+31232.D0/3.D0*z)*t32+(-t110-t
     #147+t262+4864.D0/3.D0*t100+35072.D0/3.D0+17408.D0/3.D0*t106-19584*
     #t99-31360.D0/3.D0*z)*t30+512.D0/3.D0+t186+256.D0/3.D0*t100+(-t173+
     #5248*t99+t338-15104.D0/3.D0-17536.D0/3.D0*z-2560*t100+t110-3328.D0
     #/3.D0*t106)*t28+(t110+t154-t344-t345)*t158+(-2048.D0/3.D0+33536.D0
     #/3.D0*z+t321+1024.D0/3.D0*t100-t350+t351-t277)*t16+(t345+9728.D0/3
     #.D0+4352.D0/3.D0*t99-t355+20992.D0/3.D0*z)*l1
      t362 = 6848.D0/3.D0*t100
      t365 = 2432.D0/3.D0*t103
      t372 = 9344*t99
      t376 = 192*t109
      t378 = 768*t103
      t380 = (t376-t155-704*t101+t378)*t158
      t383 = 1088.D0/3.D0*t101
      t387 = 1984.D0/3.D0*t103
      t390 = (-t240-t166+4544.D0/3.D0*t101-3840*t103+9664.D0/3.D0*t100)*
     #t168+(704*t109+4288.D0/3.D0*t106+3584.D0/3.D0*t101+2624*t103-5952*
     #t100)*t157+(-t110-15808.D0/3.D0*t101+512-t186-t254+7040.D0/3.D0*t1
     #03+26176.D0/3.D0*t100-t257)*t193-1600.D0/3.D0+1408*z+(t261+t133+t2
     #62-5184*t100-22144.D0/3.D0-17152.D0/3.D0*t106+12480*t99-4416*z)*t2
     #8+(t261+4096*t101-16256.D0/3.D0+t270+18560.D0/3.D0*z-48320.D0/3.D0
     #*t100+2240*t106+14528*t99)*t122-t277+2816*t99+(-2752.D0/3.D0*t101-
     #108608.D0/3.D0*t99-13504.D0/3.D0*t103-t153+51392.D0/3.D0-49664.D0/
     #3.D0*z+36160.D0/3.D0*t100+38464.D0/3.D0*t106)*t30+t359*s2l1+(-1728
     #*t106+t362+5632.D0/3.D0+11904*z+4032*t99-t365)*l1+(-3392.D0/3.D0*t
     #101+2240.D0/3.D0*t103+8576.D0/3.D0*t100-64-10048.D0/3.D0*t106+t372
     #+70144.D0/3.D0*z)*t32+t380+1472.D0/3.D0*t100+(-6464*t106+t383+6131
     #2.D0/3.D0*t99-6400.D0/3.D0*t100-18368.D0/3.D0+24256*z+t387)*t16
      t394 = 256*t99
      t401 = 5824.D0/3.D0*t100
      t405 = 320*t101
      t414 = 8320.D0/3.D0*t100
      t417 = 256*t103
      t419 = (t153+t417-t124-t155)*t158
      t423 = 2304*z
      t424 = 1280*t99
      t425 = 7552.D0/3.D0*t100
      t437 = 64*t101
      t440 = 2624.D0/3.D0*t100
      t449 = 320.D0/3.D0*t101
      t453 = 2432.D0/3.D0*t100
      t460 = 896*t99
      t461 = 768*t100
      t462 = 3072*z
      t465 = (-t304+t305-6016.D0/3.D0*t103+t401)*t168+(1408.D0/3.D0*t109
     #+t166-t405+t394+2752.D0/3.D0*t103-3520.D0/3.D0*t100)*t157+(-576*t1
     #01-7360.D0/3.D0*t99+12800.D0/3.D0*t103+1472.D0/3.D0*t106-t414-t129
     #)*t193+t419+(1024.D0/3.D0*t101+22592.D0/3.D0*t106-11648.D0/3.D0*t1
     #03-2048+t423+t424-t425+t176)*t122-t355-768+2752.D0/3.D0*z+(t321+64
     #0.D0/3.D0*t103-17344.D0/3.D0*t106+5824*t100+26816.D0/3.D0-t325-921
     #6*t99-8832*z)*t30+3776.D0/3.D0*t99+(t437-1472.D0/3.D0*t103-24832.D
     #0/3.D0-2752.D0/3.D0*t106-t440+t176+27328.D0/3.D0*t99+5440*z)*t28+(
     #1024*t103-128*t100+1664+13312.D0/3.D0*t99-2432*t106-t449+29696.D0/
     #3.D0*z)*t32+(-t453-4864.D0/3.D0+30464.D0/3.D0*z+1024.D0/3.D0*t106-
     #896.D0/3.D0*t103+t296)*t16+(-4096.D0/3.D0*t106+6208.D0/3.D0+t460+t
     #461+t462)*l1
      t466 = s1l1**2
      t469 = 18560.D0/3.D0*t100
      t472 = 10432*t103
      t482 = 17536.D0/3.D0*t100
      t490 = 1280*t103
      t491 = 896*t109
      t494 = -1088+t219*t220+(39040.D0/3.D0*t100+34112.D0/3.D0*t106+9728
     #.D0/3.D0*t101+14912.D0/3.D0-55616.D0/3.D0*t103-50048.D0/3.D0*t99+2
     #880*t109-3328*z)*t30+(-3520+37888.D0/3.D0*z-512*t109-7104*t106+381
     #44.D0/3.D0*t99+t235-896*t100-2240.D0/3.D0*t101)*t16+t390*s1l1+(-79
     #36*t109+5312*t101+t394-2624.D0/3.D0*t106+27712.D0/3.D0*t103-17984.
     #D0/3.D0*t100)*t157+t465*t466+12928.D0/3.D0*t99+(-512-t469+t254-672
     #0*t106-t287-5312.D0/3.D0*t101+t472-23104.D0/3.D0*z)*t28+(-t90+t92-
     #128.D0/3.D0*t21*t84*t8)*l5+(t482+512-t186+9024*t109-t472-15040.D0/
     #3.D0*t99-21760.D0/3.D0*t101+21952.D0/3.D0*t106)*t193-576*t103+(t49
     #0-t155-t491-t321)*t158
      t495 = 4544*t106
      t497 = 16256.D0/3.D0*t103
      t499 = 8896.D0/3.D0*t100
      t502 = t74*t64
      t503 = t502*t67
      t505 = t213+3*z-t199+l1-1
      t506 = t60*t3
      t507 = t505*t506
      t514 = -512.D0/3.D0*t503*t507*t13+256.D0/3.D0*t503*t507*t8
      t515 = t62*l5
      t517 = 960*t109
      t518 = 768*t106
      t526 = 5696.D0/3.D0*t109
      t530 = 1856*t109
      t532 = 1536*z
      t540 = 1408.D0/3.D0*t101
      t557 = 320*t109
      t559 = 3136*t103
      t564 = 320.D0/3.D0*t103
      t566 = 1088*t103
      t572 = 192*t101
      t575 = 26752.D0/3.D0*t106
      t582 = 3328.D0/3.D0*t100
      t583 = 64.D0/3.D0*t101
      t584 = 10496.D0/3.D0*t106
      t588 = (-t517-t518+12416.D0/3.D0*t101-5888*t103+10432.D0/3.D0*t100
     #)*t168+(t115+2752*t106-8192*t101+t526-27520.D0/3.D0*t100)*t157+(-t
     #530+17984.D0/3.D0*t101+1536-t532-6144*t99-5440*t103+2368.D0/3.D0*t
     #100+5120*t106)*t193+1984.D0/3.D0*t100-t288+(-t540-24064.D0/3.D0*t1
     #03+23744*t99+t517-7040+52672.D0/3.D0*t100-80192.D0/3.D0*t106+t130)
     #*t122+(-1216.D0/3.D0*t109-t437+10432.D0/3.D0*t103-22144.D0/3.D0*t1
     #00+8896+55936.D0/3.D0*t106-22400*t99-7040*z)*t30-128+7232.D0/3.D0*
     #z+7808.D0/3.D0*t99+(t557-3712.D0/3.D0*t101+11392.D0/3.D0+t559-t469
     #+t127-8896*t99-57088.D0/3.D0*z)*t28+t380-t564+(-t376+30592.D0/3.D0
     #*t106+t383+t566+5120-4736*t100-20672.D0/3.D0*t99+21184.D0/3.D0*z)*
     #t32+(t572-1664*t103+4096*t100-6656.D0/3.D0+t176-t575+13248*t99+392
     #96.D0/3.D0*z)*t16+(-64*t103-1600.D0/3.D0+8704*t99+t582-t583-t584+4
     #4416.D0/3.D0*z)*l1
      t601 = t60*t84*t67
      t603 = t75*l6
      t606 = t502*t8
      t609 = 512.D0/3.D0*t601*t64*t603-256.D0/3.D0*t601*t606
      t612 = 256*t28
      t620 = 6784.D0/3.D0*t106
      t633 = 64.D0/3.D0*t3*t2*l1*t8
      t634 = 448.D0/3.D0*t101
      t650 = 4928.D0/3.D0*t106
      t651 = 6016.D0/3.D0*t99
      t652 = 5312.D0/3.D0*z
      t655 = 832.D0/3.D0*t101
      t664 = 6400.D0/3.D0*t103
      t667 = 4864.D0/3.D0*t99
      t669 = 512.D0/3.D0*t100
      t675 = 640*t100
      t692 = 8896.D0/3.D0*t106
      t695 = (t634+14336.D0/3.D0+t322-8960.D0/3.D0*t100-t129+6336*t106-t
     #351-2624*z)*t32+(t118-t383-11008.D0/3.D0-6848*t103+11776.D0/3.D0*z
     #+30016.D0/3.D0*t100-29888.D0/3.D0*t106+29440.D0/3.D0*t99)*t122-236
     #8.D0/3.D0*t106-t104+(-1216*t103+4928.D0/3.D0*t100-t650+t176+t651+t
     #449+1472.D0/3.D0+t652)*t16+(-t655+4224*t103+2560*t106-11072.D0/3.D
     #0*t100+2880-t136-14848.D0/3.D0*t99-5248.D0/3.D0*z)*t30+(-t163+t401
     #+4096.D0/3.D0*t101-t314-t664)*t168+t667+(t145-t257+9088.D0/3.D0*t9
     #9+t669-t437+384+15680.D0/3.D0*z)*l1+(t110+3136.D0/3.D0*t106-t572-1
     #856.D0/3.D0*t103-768-t675-768*t99-16064.D0/3.D0*z)*t28+(t181+4736.
     #D0/3.D0*t106-8128.D0/3.D0*t101+t394+3520*t103-10816.D0/3.D0*t100)*
     #t157-1216.D0/3.D0+t419+192*t100+5120.D0/3.D0*z+(-t189+2240*t101+10
     #24-t299-10432.D0/3.D0*t99+4096.D0/3.D0*t103-8384.D0/3.D0*t100+t692
     #)*t193
      t696 = s2l1**2
      t716 = t425-t495+(-3968.D0/3.D0*t101-t497+4032*t109+t499-t166)*t16
     #8+t514*t515+t588*s2l1+2496*z+(-38080.D0/3.D0*t106-7040.D0/3.D0+140
     #16*t103-6272*t109-23872.D0/3.D0*t100+7808.D0/3.D0*z+2624*t101+3788
     #8.D0/3.D0*t99)*t122+t609*t62+(1216.D0/3.D0*t32-t612+896*l1+2816.D0
     #/3.D0*t16+896.D0/3.D0)*t81+(-1152+30592.D0/3.D0*z+1024*t101+17152.
     #D0/3.D0*t99+t414-t620-2816*t103)*l1+(9728.D0/3.D0-2944.D0/3.D0*t10
     #1-t482+40768.D0/3.D0*t106-22016.D0/3.D0*t99+t307+1472*t109-t270)*t
     #32+t633+t695*t696+((-15808.D0/3.D0*t32+1792*t28-1088.D0/3.D0-12992
     #.D0/3.D0*t16-3264*l1)*s1l1-t31+(-2624.D0/3.D0+6272.D0/3.D0*t28-544
     #0.D0/3.D0*t32-9920.D0/3.D0*l1-t37-2368*t16)*s2l1+6592.D0/3.D0*t28-
     #11840.D0/3.D0*t32-2752.D0/3.D0-4288*t16-12544.D0/3.D0*l1)*t82
      t722 = (s1l1+1+s2l1)*(2*s1l1+1)
      t726 = s1l1*s2l1
      t776 = 14656.D0/3.D0*t100
      s1 = 12544.D0/3.D0+(-t12+t15-128.D0/3.D0*t21*t3*t8)*l5+(9152.D0/3.
     #D0*t16+256.D0/3.D0*t28-512.D0/3.D0-t31+20800.D0/3.D0*t32+(7808.D0/
     #3.D0*t28+512.D0/3.D0+8512.D0/3.D0*t16+5184*t32-t37+1472*l1)*s1l1+(
     #1856*t16+5696*t32-1920*t28-384+1728*l1)*s2l1+6208.D0/3.D0*l1)*t82+
     #64.D0/3.D0*t84*l1*t8+(t90-t92+128.D0/3.D0*t21*t84*t8)*l5+t195*s2l1
      t785 = s1-t197*(t494+t716)+(256*t722*t106+64.D0/3.D0*(22*s2l1+140*
     #t726+94*t696+34-38*s1l1+46*t466)*t100+(-3584.D0/3.D0*t696-7168.D0/
     #3.D0*t726-13120.D0/3.D0-832*s1l1-3584.D0/3.D0*t466-1856*s2l1)*t103
     #-t189+(4928+2624.D0/3.D0*s1l1+3392.D0/3.D0*s2l1)*t101)*t168+((-140
     #8.D0/3.D0+t37-1856*l1-7040.D0/3.D0*t28-9920.D0/3.D0*t16-12928.D0/3
     #.D0*t32)*s1l1+(256.D0/3.D0-2112*l1+2176*t28-6976.D0/3.D0*t16-14464
     #.D0/3.D0*t32)*s2l1-128+t31-7360.D0/3.D0*l1+512.D0/3.D0*t28-3520*t1
     #6-18176.D0/3.D0*t32)*t82-64.D0/3.D0*(44*t696+138*t726+168*s2l1+124
     #+183*s1l1+48*t466)*t98+(-6784.D0/3.D0*t99+1472.D0/3.D0*t109-3392*t
     #106+t305+t776+1664.D0/3.D0-19840.D0/3.D0*z-10048.D0/3.D0*t103)*t16
     #-t514*t515-t609*t62
      t788 = 2176*t726
      t798 = 704.D0/3.D0*t103
      t804 = 2880*t696
      t806 = 512*t726
      t808 = 2432*t466
      t812 = 448*t696
      t813 = 256*t726
      t815 = 640.D0/3.D0*t466
      t827 = 12*s2l1
      t828 = 6*t696
      t848 = t804+(1600.D0/3.D0*t696-t806-12352.D0/3.D0-6784*s1l1-t808-2
     #176.D0/3.D0*s2l1)*t106+(-t812-t813+1024+9152.D0/3.D0*s1l1+t815-576
     #*s2l1)*t100+(3392.D0/3.D0*t696+4480*t726+9088+13696*s1l1+13696.D0/
     #3.D0*t466+7168*s2l1)*t99-64.D0/3.D0*(8*t726+24+t827+t828+23*s1l1)*
     #t103-64.D0/3.D0*(-5*s2l1-15)*t101+(-7040.D0/3.D0*t696-9600*t726-11
     #392-19968*s1l1-22720.D0/3.D0*t466-36736.D0/3.D0*s2l1)*t98+30848.D0
     #/3.D0*t726+25216.D0/3.D0+47744.D0/3.D0*s1l1+22016.D0/3.D0*t466+328
     #96.D0/3.D0*s2l1
      t850 = 2*t726
      t861 = 64.D0/3.D0*t696
      t871 = 9408*s1l1
      t882 = 1344*t466
      t902 = 768*s1l1
      t916 = 14080.D0/3.D0*s2l1
      t933 = (-1792.D0/3.D0*s1l1+2816+5056.D0/3.D0*s2l1)*t101+t181+(1088
     #*t696+5888.D0/3.D0*t726+21760.D0/3.D0+53248.D0/3.D0*s1l1+7232.D0/3
     #.D0*t466+35072.D0/3.D0*s2l1)*t106+256*t696+t806+256+t902+256*t466+
     #768*s2l1+(-7360.D0/3.D0*t696-3584.D0/3.D0*t726-1408-3904.D0/3.D0*s
     #1l1-12544.D0/3.D0*t466-14144.D0/3.D0*s2l1)*t100+(4672*t696-6784*t7
     #26-8320-26880*s1l1-12480*t466-t916)*t99+(-1792.D0/3.D0*t696-t788-1
     #5104.D0/3.D0+6208.D0/3.D0*s1l1+448*t466-8576.D0/3.D0*s2l1)*t103-12
     #8.D0/3.D0*(-104*t726-70-192*s1l1+44*s2l1+66*t696-170*t466)*t98
      t949 = 1024.D0/3.D0*t696+4928.D0/3.D0*s2l1+t788+64.D0/3.D0*(41*t69
     #6+72*t726+111+97*s1l1+152*s2l1+10*t466)*t99+(-t394+9344*t109+t798-
     #22400.D0/3.D0*t101-19648.D0/3.D0*t100+14080.D0/3.D0*t106)*t157+t84
     #8*t16+(256*(t696+t850+1+t466)*t99+64.D0/3.D0*(-58*s2l1-248*t726+t8
     #28-232-318*s1l1-254*t466)*t106+(t861-256.D0/3.D0*t726-3584.D0/3.D0
     #-2368*s1l1-1088.D0/3.D0*t466-3328.D0/3.D0*s2l1)*t103+9152.D0/3.D0*
     #t109+(-9920.D0/3.D0*t696-2560*t726+6272+t871+3776.D0/3.D0*t466+620
     #8.D0/3.D0*s2l1)*t100+(-10304.D0/3.D0+832.D0/3.D0*s2l1-256*s1l1)*t1
     #01)*t157+t882-17984.D0/3.D0*t99+64.D0/3.D0*(1+s2l1)*(2*s2l1+7)*t10
     #0+t933*t122+(128.D0/3.D0*t99+23680*t106+64256.D0/3.D0*t103-95744.D
     #0/3.D0*t100-17728.D0/3.D0+29888.D0/3.D0*z-11840.D0/3.D0*t101-9344.
     #D0/3.D0*t109)*t30+8768.D0/3.D0*s1l1+(-1536*t103+t155+t572+1088*t10
     #9)*t158
      t961 = -4480.D0/3.D0*t32-704*t16+64-1408.D0/3.D0*l1+t612
      t963 = 64.D0/3.D0*t103
      t976 = 7360.D0/3.D0*t103
      t984 = 1600.D0/3.D0*t100
      t988 = 256*z
      t991 = 256.D0/3.D0*t101
      t993 = 1088*t100
      t994 = 2176.D0/3.D0*t106
      t1001 = 2176.D0/3.D0*t109
      t1002 = 384*t106
      t1006 = (t155-t261-t417)*t158
      t1007 = 2048.D0/3.D0*t103
      t1008 = 2752.D0/3.D0*t100
      t1010 = 1088*t99
      t1011 = 64.D0/3.D0*t109
      t1017 = 3200.D0/3.D0*t101
      t1019 = 768*t109
      t1025 = t963+(t344+t497+2304-t304-31616.D0/3.D0*t106+3008*t100-108
     #8.D0/3.D0*t99-2560*z)*t122+960-1856.D0/3.D0*z+(-t321-3904.D0/3.D0*
     #t103+t325+t675-1216*t106+3136*t99-16256.D0/3.D0-t652)*t32+(-t976-5
     #440*t100+27712.D0/3.D0*t106+5696*t99+t129-7744+30784.D0/3.D0*z)*t3
     #0-3136.D0/3.D0*t99+(-6464.D0/3.D0-t984+2240.D0/3.D0*t106-512.D0/3.
     #D0*t103+832*t99+t583+t988)*l1+(t991+4288.D0/3.D0*t103+8256+t993+t9
     #94-27712.D0/3.D0*t99-t153-21760.D0/3.D0*z)*t28+(-t394+9344.D0/3.D0
     #*t103+768*t101+t993-t1001-t1002)*t157+t355+t1006+(t1007-3136.D0/3.
     #D0+t1008-8320.D0/3.D0*t106+t1010-t1011-t988)*t16+(t136-t564-t572-4
     #480.D0/3.D0*t100)*t168+(-t1017-19136.D0/3.D0*t103+t984+t1019+3648*
     #t106+1408.D0/3.D0*t99)*t193-128.D0/3.D0*t100
      t1028 = t197**2
      t1032 = 1408*t101
      t1033 = 1088*z
      t1038 = 128*t32
      t1043 = 320*l1
      t1049 = t109-t101
      t1053 = -t1049
      t1055 = 64*t1053*t158
      t1056 = 7424.D0/3.D0*t106
      t1057 = 16256.D0/3.D0*t99
      t1065 = 4096.D0/3.D0*t100
      t1072 = 5248.D0/3.D0*t103
      t1078 = 1024*t100
      t1082 = 640*z
      t1090 = 128*t106
      t1091 = 640.D0/3.D0*t101
      t1097 = 2432.D0/3.D0*t99
      t1098 = 896*z
      t1102 = t1055-t963+(t325-9856.D0/3.D0-t1056+t147+t1057-832*t103+42
     #88*z+t1008)*t32+(-t304-4352*t106-t540+1920*t99+14080.D0/3.D0*t103-
     #t1065)*t122+1024.D0/3.D0*t1049*t168+2176.D0/3.D0-t1033+64*t100+(t1
     #29+28352.D0/3.D0*t106-t1072+t540-11840.D0/3.D0*t99-12416.D0/3.D0*t
     #100-256+t988)*t30+(-t1001+1792.D0/3.D0*t101+t322-t1078)*t157-3200.
     #D0/3.D0*t99+(t147-t1010-t453-t104-t1082+t650)*l1+(-t153-t584-t991+
     #2688-2048*t99+t798-3712*z+8512.D0/3.D0*t100)*t28+(t1019+t1090-t109
     #1-4032*t103+10048.D0/3.D0*t100)*t193+896.D0/3.D0*t106+(-t1011+128-
     #t449+t1097+t270+t1098+t669-3584.D0/3.D0*t106)*t16
      t1115 = 3968.D0/3.D0*t100
      t1119 = 640.D0/3.D0*t109
      t1120 = 896*t101
      t1126 = 1664.D0/3.D0*t109
      t1128 = 128*t103
      t1132 = 2816.D0/3.D0*t100
      t1153 = 2048.D0/3.D0*t101
      t1160 = 128*z
      t1164 = 128*t1053*t158+t994+(t176-2816*t106-t1072+t102+t372+t1115-
     #18944.D0/3.D0+9472*z)*t32+(-t1119-t575-t1120+4096*t99+25600.D0/3.D
     #0*t103-7808.D0/3.D0*t100)*t122+(t1126-1280.D0/3.D0*t101-t1128)*t16
     #8-2176*t99+(-t991+t1132+t378+768-t988+t326-8576.D0/3.D0*t106)*t16+
     #4352.D0/3.D0-2176*z+(-t491+t1091+8576.D0/3.D0*t103-2176*t100)*t157
     #+(t104+8704.D0/3.D0*t106-1536*t100-4480.D0/3.D0*t99-896.D0/3.D0-14
     #08.D0/3.D0*z)*l1+(-t153-t309+13568.D0/3.D0-t147-3584.D0/3.D0*t99+t
     #365-20096.D0/3.D0*z+15616.D0/3.D0*t100)*t28+(t240+t166+20480.D0/3.
     #D0*t100+t1153-25216.D0/3.D0*t103)*t193+(t153+59648.D0/3.D0*t106+t1
     #24-28288.D0/3.D0*t99-8192.D0/3.D0*t103-128+t1160-23936.D0/3.D0*t10
     #0)*t30
      t1167 = 30848.D0/3.D0*t106
      t1175 = 448.D0/3.D0*t109
      t1176 = 2176*t106
      t1182 = 5504.D0/3.D0*t101
      t1190 = 3008.D0/3.D0*t101
      t1206 = 4160.D0/3.D0*t103
      t1226 = t1164*s2l1+8512.D0/3.D0-3904*z+(t1011-t1167-4736.D0/3.D0*t
     #101+28480.D0/3.D0-9984*t99-2304*t103-38272.D0/3.D0*z+43904.D0/3.D0
     #*t100)*t28+(-t1175+t1176+33152.D0/3.D0*t101-57152.D0/3.D0*t103+179
     #20.D0/3.D0*t100)*t193+t1055+(-t1011-35776.D0/3.D0+t1182+21312*t99-
     #t559+17088*z+2624*t100-10560*t106)*t32+9536.D0/3.D0*t106+(t1175+t1
     #190-t322)*t168-5312*t99+(t1011-36032.D0/3.D0*t99+115648.D0/3.D0*t1
     #06+8512.D0/3.D0*t101-896+t1098+14656.D0/3.D0*t103-100160.D0/3.D0*t
     #100)*t30+(-t1011-5888*t101+27136.D0/3.D0*t103-3136*t100)*t157+(t12
     #06-1984.D0/3.D0-320.D0/3.D0*z+13184.D0/3.D0*t106-4480*t100-t320)*l
     #1+(-1984.D0/3.D0*t101-18304.D0/3.D0*t106+16832.D0/3.D0*t100-t134+1
     #152-1216*z+t1097)*t16+(t261-64064.D0/3.D0*t106-25984.D0/3.D0*t101+
     #6208*t99+10880*t103+12864*t100)*t122-704*t100
      t1281 = t1055+5440.D0/3.D0-7616.D0/3.D0*z+(t325-26816.D0/3.D0*t106
     #+2368.D0/3.D0*t101-25792.D0/3.D0+47552.D0/3.D0*t99-6208.D0/3.D0*t1
     #03+34112.D0/3.D0*z+2752*t100)*t32-t165+(t136-t1091-t1128)*t168+524
     #8.D0/3.D0*t106+(t129+76672.D0/3.D0*t106+t179-8960*t99+t338-384+384
     #*z-57472.D0/3.D0*t100)*t30+(11200.D0/3.D0*t103-2112*t100-t1120-t10
     #01)*t157-9728.D0/3.D0*t99+(-t153-24704.D0/3.D0*t106-1664.D0/3.D0*t
     #101-19328.D0/3.D0*t99-6656.D0/3.D0*t103+19456.D0/3.D0-26176.D0/3.D
     #0*z+33280.D0/3.D0*t100)*t28+(t1019+t1002+3904*t101-11584*t103+6528
     #*t100)*t193-t963+(-t1011-8128.D0/3.D0*t106+704-t540+704*t99+t134+7
     #68*z+5632.D0/3.D0*t100)*t16+(4160*t99-4800*t101-34880.D0/3.D0*t106
     #+31232.D0/3.D0*t103+t362-t304)*t122+(t147-6208.D0/3.D0*t99-t499+t1
     #007-64.D0/3.D0-3776.D0/3.D0*z+11456.D0/3.D0*t106)*l1
      t1290 = (-96896.D0/3.D0*t100-t530+57088.D0/3.D0*t103+61568.D0/3.D0
     #*t106-t1032-1088+t1033-8704.D0/3.D0*t99)*t30+((-64*t16-t1038+192*l
     #1)*s1l1-320*t16+(t1043-448*t16+t1038)*s2l1+t1043)*t82+384*t1049*t1
     #58+11008.D0/3.D0*t106+t1102*t696+(8896*t100-40000.D0/3.D0*t99+2086
     #4.D0/3.D0+3392.D0/3.D0*t109-7360*t103-25856.D0/3.D0*z+10688.D0/3.D
     #0*t106+t634)*t28+t1226*s1l1+t1281*s2l1+t145+(-64*l1+64*t16)*t81+(-
     #6400.D0/3.D0*t109+9664.D0/3.D0*t101-t566)*t168
      t1298 = 832.D0/3.D0*z
      t1344 = 704.D0/3.D0*t109
      t1377 = t1055+1408-2112*z+(-t1011-24448.D0/3.D0*t106-2496*t100-t40
     #5+3968*t99+21056.D0/3.D0*t103)*t122+(t1119-t991-t1128)*t168+(t270+
     #256-576*z-t692+t1132+3520.D0/3.D0*t99)*t16+704*t106+(-t1344-t405+7
     #808.D0/3.D0*t103-2048*t100)*t157-2112*t99+(-t1011-8640*t106+2624.D
     #0/3.D0*t103-t449+64*t99+12608.D0/3.D0*t100+10880.D0/3.D0-t175)*t28
     #+(t325+t1090-21824.D0/3.D0*t103+t1153+19136.D0/3.D0*t100)*t193+(t9
     #91+7680*t99+1280*t100-4352.D0/3.D0*t103-5248+8512*z-6656.D0/3.D0*t
     #106)*t32+(t176+55360.D0/3.D0*t106-2240*t103-28480.D0/3.D0*t99-2067
     #2.D0/3.D0*t100+128-t1160)*t30+(-t1065-512.D0/3.D0-t424+2688*t106-6
     #4.D0/3.D0*z)*l1
      t1379 = 5312.D0/3.D0+(-30208.D0/3.D0*t101+17024.D0/3.D0*t103+15488
     #.D0/3.D0*t109-t461)*t157+(-t1153+64.D0/3.D0*t99-1088.D0/3.D0+t976+
     #4672.D0/3.D0*t106-t1298-9728.D0/3.D0*t100)*l1+(t136-9344.D0/3.D0*t
     #106+576-t1160+t655+16448.D0/3.D0*t100+256.D0/3.D0*t99-3392*t103)*t
     #16-1728*t100+(-20096.D0/3.D0*t109-7936.D0/3.D0*t100+42304.D0/3.D0*
     #t101-20800.D0/3.D0*t103+t1176)*t193+(-43072.D0/3.D0*t106-7040*t101
     #+4672*t109+2752*t99-22336.D0/3.D0*t103+64256.D0/3.D0*t100)*t122-12
     #160.D0/3.D0*t99+(-3008.D0/3.D0*t109+t776-14016*t106-23552.D0/3.D0+
     #10240*z+t1182-t490+52288.D0/3.D0*t99)*t32+t1377*t466-t423
      t1393 = 7424.D0/3.D0*t100
      t1397 = 6272.D0/3.D0*t99
      t1402 = 15488.D0/3.D0*t106
      t1409 = 33664.D0/3.D0*t99
      t1425 = 5504.D0/3.D0*t726
      t1432 = 576*s1l1
      t1433 = 384*t466
      t1453 = 13504*s1l1
      t1456 = -9280.D0/3.D0*t696+t1119+(-1216.D0/3.D0*t696-19456.D0/3.D0
     #*t726-10112.D0/3.D0-50240.D0/3.D0*s1l1-7360*t466-19328.D0/3.D0*s2l
     #1)*t98+(-896*t696+t1425+1984+2752.D0/3.D0*s1l1-4864*t466+8704.D0/3
     #.D0*s2l1)*t106+(t861-1408.D0/3.D0*t726-320-t1432-t1433+896*s2l1)*t
     #103+(t804+t806+3712.D0/3.D0+11584.D0/3.D0*s1l1+4160*t466+15616.D0/
     #3.D0*s2l1)*t99+(-2752.D0/3.D0*t696-4096.D0/3.D0*t726-1856+1728*s1l
     #1+4288.D0/3.D0*t466-4224*s2l1)*t100+(448.D0/3.D0*s2l1+1600.D0/3.D0
     #)*t101+4864*t726+1664+t1453+25408.D0/3.D0*t466-320*s2l1
      t1458 = 64.D0/3.D0*(-46-62*s2l1-16*t696-12*t726-17*s1l1)*t106+t102
     #5*t696+1600.D0/3.D0*t103+t1028*(t1290+t1379)/2+(-26112*t99+10432*t
     #106+1920*t101+6208*t100+44416.D0/3.D0-63296.D0/3.D0*z+t526-28480.D
     #0/3.D0*t103)*t28-t1393+(7552*t100+512-t186-6464.D0/3.D0*t103+26816
     #.D0/3.D0*t101-t1397-30272.D0/3.D0*t109-8000.D0/3.D0*t106)*t193+t14
     #02+(t12-t15+64.D0/3.D0*t215*t202*t8)*t220+(24640.D0/3.D0*t100+t140
     #9-25024.D0/3.D0-6592.D0/3.D0*z-t181-t1206-14656*t106+960*t101)*t32
     #-11072.D0/3.D0*z+t1456*t28
      t1469 = 192*s2l1
      t1470 = 4*t726
      t1577 = 19072.D0/3.D0*t106
      t1582 = (t317-3840*t106+t124-t1409+t490+10752-63872.D0/3.D0*z-t153
     #)*t28+(-t290+3712*t103+1280*t101+5248.D0/3.D0*t100-t491+14336.D0/3
     #.D0*t106)*t157+(-t262+t582-t323-t321+12800.D0/3.D0*t99-9856-t462+t
     #176)*t32+(t344+21376.D0/3.D0*t103-16384.D0/3.D0-t1119-t1167-6272.D
     #0/3.D0*t100+42880.D0/3.D0*t99+14848.D0/3.D0*z)*t122+(-t152+t1126+t
     #104-t1393-t314)*t168+(-t323-8192.D0/3.D0+2560.D0/3.D0*t103+3968.D0
     #/3.D0*t99+t165-11392.D0/3.D0*z)*t16+(-t321-t235-5504*t100+t153+192
     #00*t106-8576*t99-6272.D0/3.D0+25984.D0/3.D0*z)*t30+1536-t423+t303+
     #(-t1078+t144-256.D0/3.D0*t103-t460-7936.D0/3.D0-t1082)*l1+(10624.D
     #0/3.D0*t100-22784.D0/3.D0*t103-3328.D0/3.D0*t99+t240-t1577-t1032+1
     #024-t299)*t193+(t345-t153-t154)*t158-t1397
      t1601 = (-4736+t153+143360.D0/3.D0*t106+2816*t101-16128*t99+4096*t
     #103+14848*z-99584.D0/3.D0*t100)*t30+(-t1017-2176*t103+t620+29696.D
     #0/3.D0*t100-35328*t99+70016.D0/3.D0-45056*z-t153)*t28+(-9728.D0/3.
     #D0*t103+3456*t100-4864*t106-29824.D0/3.D0+1536*t101+11648.D0/3.D0*
     #t99-22784.D0/3.D0*z+t176)*t32+(-8064*t101+12032*t103-9728-t1119-13
     #9648.D0/3.D0*t106+12416*t100+120320.D0/3.D0*t99+8960*z)*t122+3072*
     #t106+(-4224*t100+13696.D0/3.D0*t106-3712.D0/3.D0*t99-8960.D0/3.D0-
     #2816.D0/3.D0*z+t490)*l1+(4608*t100-20480.D0/3.D0*t106-t1153-t104-1
     #152*t99-384-11008*z)*t16-t1057+8704.D0/3.D0+t1582*s2l1-2048.D0/3.D
     #0*t100-4608*z+t159+(-t124+t1126-t270+3200.D0/3.D0*t100-t518)*t168+
     #(-10624.D0/3.D0*t101+32768.D0/3.D0*t103-39808.D0/3.D0*t100-t491+67
     #84*t106)*t157+(26624.D0/3.D0*t101+1536-t532+59008.D0/3.D0*t100-217
     #60*t103+t240-t518-8192*t99)*t193
      t1614 = t60*t202*t505*t67
      t1669 = (-23296.D0/3.D0-t1011-t1577-14464.D0/3.D0*t100+t991+22528.
     #D0/3.D0*z+19264*t99+21568.D0/3.D0*t103)*t122+(t1119-t572-t314+448.
     #D0/3.D0*t103-2944.D0/3.D0*t100)*t168+(t387+13120.D0/3.D0*t100-t495
     #+t991-9024*t99-t1011+21440.D0/3.D0-20736*z)*t28+(-t291+1024-t299+2
     #3488.D0/3.D0*t100+t325-33856.D0/3.D0*t106-1216.D0/3.D0*t101-4736.D
     #0/3.D0*t99)*t193+(4352.D0/3.D0*t100-t1056+5504.D0/3.D0*t99-3968.D0
     #/3.D0*t103-t583-4800-4736.D0/3.D0*z)*t32+3328.D0/3.D0-2240*z+1984.
     #D0/3.D0*t106+(384*t100-7616.D0/3.D0*t106+t667-5312.D0/3.D0-3392*z+
     #t350)*t16+(-t437-8768.D0/3.D0*t103-18880.D0/3.D0*t100+64576.D0/3.D
     #0*t106-60224.D0/3.D0*t99+t176+17728.D0/3.D0-1856*z)*t30+(-1600*t99
     #+8192.D0/3.D0*t106-2432.D0/3.D0-t1298-t1115)*l1-t651+t1006+(-t394+
     #t664+t344-t440-t1344+t1402)*t157
      t1706 = 1984.D0/3.D0*t696
      t1712 = (-18304.D0/3.D0*t696+t813-6720-5440.D0/3.D0*s1l1+6016.D0/3
     #.D0*t466-37184.D0/3.D0*s2l1)*t99+(-t1469+64*s1l1-5056.D0/3.D0)*t10
     #1+(22016.D0/3.D0*t696+7168*t726+8000+t1453-512.D0/3.D0*t466+58432.
     #D0/3.D0*s2l1)*t98-1088.D0/3.D0*t109+(-1728*t696-14464.D0/3.D0*t726
     #+3968.D0/3.D0+8896.D0/3.D0*s1l1+10432.D0/3.D0*t466-20672.D0/3.D0*s
     #2l1)*t106-7552*s2l1-18176.D0/3.D0*t726-6016.D0/3.D0*t696-11392.D0/
     #3.D0-9856*s1l1-12160.D0/3.D0*t466+(10688.D0/3.D0*t696+3584*t726+89
     #6-21440.D0/3.D0*s1l1+320*t466+28736.D0/3.D0*s2l1)*t100+(-t1706-128
     #*t726+6272.D0/3.D0+4544.D0/3.D0*s1l1-t815-6848.D0/3.D0*s2l1)*t103
      t1714 = t506*t67
      t1770 = 9664.D0/3.D0*t696+64.D0/3.D0*(-t827-16-t696)*t103+64.D0/3.
     #D0*(76+18*t726+33*s1l1+79*s2l1+17*t696)*t100+(-t1706-t1425-11200.D
     #0/3.D0-4288*s1l1-t1433-3648*s2l1)*t106+(2176.D0/3.D0*t696+11776.D0
     #/3.D0*t726+17792.D0/3.D0+8512*s1l1+t808+15424.D0/3.D0*s2l1)*t99+(-
     #1984*t696-17920.D0/3.D0*t726-20864.D0/3.D0-t871-11008.D0/3.D0*t466
     #-7488*s2l1)*t98+6528*t726+16384.D0/3.D0+23552.D0/3.D0*s1l1+3008*t4
     #66+8192*s2l1
      t1827 = 15872.D0/3.D0*t696+64.D0/3.D0*(-21-15*s2l1+5*s1l1)*t101+(t
     #812+2176.D0/3.D0*t726+384-1984.D0/3.D0*s1l1-64.D0/3.D0*t466+3776.D
     #0/3.D0*s2l1)*t103+(6592.D0/3.D0*t696+2944.D0/3.D0*t726-10048.D0/3.
     #D0-29056.D0/3.D0*s1l1-3392*t466+3712*s2l1)*t106-320.D0/3.D0*t109+(
     #-9152.D0/3.D0*t696+512.D0/3.D0*t726+18944.D0/3.D0+16512*s1l1+27328
     #.D0/3.D0*t466-9728.D0/3.D0*s2l1)*t99+(-3008.D0/3.D0*t696-4352*t726
     #-28736.D0/3.D0-51968.D0/3.D0*s1l1-30272.D0/3.D0*t466-17536.D0/3.D0
     #*s2l1)*t98+(-1664.D0/3.D0*t696-t813+3520.D0/3.D0+2432*s1l1+t882-58
     #88.D0/3.D0*s2l1)*t100+9856*t726+10176+16000*s1l1+21760.D0/3.D0*t46
     #6+42880.D0/3.D0*s2l1
      s1 = ((320*s1l1-6976.D0/3.D0*s2l1-3520.D0/3.D0)*t101+(1600*t696+11
     #008.D0/3.D0*t726+22528.D0/3.D0+t1432+1216*t466+5696*s2l1)*t103-64.
     #D0/3.D0*(-t1469-t1470+70*t696-74-360*s1l1-74*t466)*t99+(5888.D0/3.
     #D0*t696-2432.D0/3.D0*t726-25024.D0/3.D0-7808*s1l1-3328.D0/3.D0*t46
     #6-6784.D0/3.D0*s2l1)*t100-9088.D0/3.D0*t109-512*t722*t98+(-896.D0/
     #3.D0*t696+7424*t726+3968+t902+23168.D0/3.D0*t466-t916)*t106)*t193-
     #t219*t220+t1601*s1l1-t633+(-72064.D0/3.D0*t106-32960.D0/3.D0*t103-
     #10496.D0/3.D0+9728.D0/3.D0*z-6848.D0/3.D0*t101+46208.D0/3.D0*t99+6
     #720*t109+56000.D0/3.D0*t100)*t122+(-512.D0/3.D0*t1614*t502*t13+256
     #.D0/3.D0*t1614*t606)*t515+t1669*t466
      t1829 = s1+t1712*t30+(512.D0/3.D0*t1714*t64*t603-256.D0/3.D0*t1714
     #*t606)*t62+(t557-256*(-t1470-2*t466-7-2*t696-3*s2l1-3*s1l1)*t103+6
     #4*(-29-8*s2l1-8*s1l1)*t101-256*(1+t696+s1l1+t850+s2l1+t466)*t100)*
     #t158+t1770*l1+(-1472*t99-3328.D0/3.D0-7040.D0/3.D0*z-t1190+3264*t1
     #06+3072*t103-4352*t100)*l1+(5888.D0/3.D0*t101-t161-14656.D0/3.D0*t
     #109-t166+3648*t103)*t168+t1827*t32
      t1833 = (l1+z-t5)**2
      taux17 = (t785+t949+t1458+t1829)/t1833*t8
      t1 = 1-z
      t2 = t1**2
      t3 = t2**2
      t5 = t3*t1
      t7 = t3*t2
      t9 = t2*t1
      t12 = t3*t9
      t17 = l1**2
      t18 = t17*l1
      t28 = 56576.D0/3.D0*s2l1
      t29 = 38*s2l1
      t51 = 6144*t2
      t57 = t17**2
      t58 = t57**2
      t64 = -s1l1-1-s2l1
      t65 = 2*t64
      t87 = t57*t18
      t95 = 4096.D0/3.D0*t3
      t102 = t57*t17
      t104 = 3*t64
      t124 = 91904.D0/3.D0*t12
      t128 = 2048*t9
      t134 = t58*l1
      t178 = t57*l1
      t181 = 4096*z
      t219 = 384*t12*t58*t18
      t224 = 24064.D0/3.D0*z
      t237 = 3072*t3
      t238 = 3328*t12
      t239 = 12288*t5
      t242 = t58*t17
      t245 = 17664*t3
      t246 = 6656.D0/3.D0*t5
      t249 = 2048.D0/3.D0*t12
      t257 = 46336.D0/3.D0*t9
      t267 = (10368*t12+4096-t181-382336.D0/3.D0*t7-64640.D0/3.D0*t5+182
     #272.D0/3.D0*t3-70656*t2+342656.D0/3.D0*t9)*t87+143744.D0/3.D0*t2+(
     #-2048*t12+199040.D0/3.D0*t5+66688.D0/3.D0*t7+26240+268160.D0/3.D0*
     #t9-145536*t3-239872.D0/3.D0*z-186880.D0/3.D0*t2)*t57-105472.D0/3.D
     #0*t9+(-19456*t5-141184.D0/3.D0*t9+131072.D0/3.D0*t3-22144+151936.D
     #0/3.D0*t2+146176.D0/3.D0*z)*l1-61952.D0/3.D0+98432.D0/3.D0*z+(-108
     #416.D0/3.D0-115328.D0/3.D0*t7-99840*t5+2176*t12-186880.D0/3.D0*t9+
     #189312*t3+211328.D0/3.D0*z+84736.D0/3.D0*t2)*t178+t219+(288256.D0/
     #3.D0*t7-33280.D0/3.D0-212608.D0/3.D0*t3+130688.D0/3.D0*t5-11008.D0
     #/3.D0*t12+t224-147200*t9+380800.D0/3.D0*t2)*t102+(8704*t7-40960.D0
     #/3.D0+95872.D0/3.D0*z+48128.D0/3.D0*t5-23040*t3-73216.D0/3.D0*t9+1
     #36448.D0/3.D0*t2)*t17+(-t237-t238+t239+2688*t7)*t242+(-22528*t7+73
     #344-t245-t246-172416*t2+390016.D0/3.D0*t9+t249-107904*z)*t18+30208
     #.D0/3.D0*t3+(-14848*t12+t51+146944.D0/3.D0*t5+258560.D0/3.D0*t7-t2
     #57-80000*t3)*t58+(-t128-82432.D0/3.D0*t7+30848.D0/3.D0*t12+36352*t
     #3-44160*t5)*t134
      t269 = 6656*s2l1
      t270 = dlog(z)
      t283 = 1024*t9
      t287 = -1+l1
      t295 = 2*l6-1
      t296 = t295**2
      t299 = l5-1
      t301 = z*l1
      t303 = (1-l1+t301)**2
      t304 = 1/t303
      t309 = dsqrt(z)
      t311 = dsqrt(l5)
      t313 = dsqrt(-t299)
      t316 = dsqrt(l1)
      t320 = l5*z*l1-l1+1-l5+l5*l1
      t323 = dsqrt(-t287)
      t328 = t320**2
      t331 = dsqrt(4*l1*t296*z*l5*t299*t287*t304-4*t295*t309*t311*t313*t
     #316*t320*t304*t323+t328*t304)
      t346 = 2048*z
      t361 = 12800.D0/3.D0*z
      t367 = 1280*t12
      t385 = 29696.D0/3.D0*t3
      t386 = 2048*t2
      t390 = 7424.D0/3.D0*t9
      t395 = 1024*t3
      t396 = 3072*t5
      t398 = (t367-1280*t7+t395-t396)*t242
      t424 = 1280.D0/3.D0*t12
      t429 = 67328.D0/3.D0*t3
      t433 = (-16640.D0/3.D0*t12+18560*t7-38656.D0/3.D0*t3+t283+2944*t5)
     #*t134+(6400.D0/3.D0*t12+15872.D0/3.D0-319232.D0/3.D0*t7+74240.D0/3
     #.D0*t5-t361+492928.D0/3.D0*t3-201856.D0/3.D0*t9-27392*t2)*t102+(t3
     #67-53120*t5-59008.D0/3.D0*t7+134528.D0/3.D0-98432.D0/3.D0*t9+14144
     #0*t3-113792.D0/3.D0*z-234112.D0/3.D0*t2)*t57-44416*t2+(-6400*t12+3
     #93472.D0/3.D0*t7+77312.D0/3.D0*t2-2048+t346-85120.D0/3.D0*t3-30617
     #6.D0/3.D0*t5-18560*t9)*t87-t385+(-t386+8960*t12-74880*t7+129920.D0
     #/3.D0*t5+t390+21760*t3)*t58+t398+103424.D0/3.D0*t9+(58112.D0/3.D0*
     #t5-130304.D0/3.D0*t3+8704.D0/3.D0+89600.D0/3.D0*t9-30592.D0/3.D0*t
     #2-30592.D0/3.D0*z)*l1+19456-76160.D0/3.D0*z+(-t367+115840.D0/3.D0*
     #t7+98944*t5+9344.D0/3.D0-923008.D0/3.D0*t3-16000*z+199808*t9-87424
     #.D0/3.D0*t2)*t178+(-27136.D0/3.D0*t7+153856.D0/3.D0*t3-75776.D0/3.
     #D0*t9-45824.D0/3.D0*t5-24320.D0/3.D0*t2+13696.D0/3.D0-8960*z)*t17+
     #(-t424+68864.D0/3.D0*t7-234368.D0/3.D0+302080.D0/3.D0*z-47360.D0/3
     #.D0*t5+173568*t2+t429-371200.D0/3.D0*t9)*t18
      t435 = 141440.D0/3.D0*t9
      t442 = 7936.D0/3.D0*t9
      t473 = 5120*t5
      t498 = 5120.D0/3.D0*t12
      t511 = 8704.D0/3.D0*t7
      t523 = 39424.D0/3.D0*t9
      t528 = 15360*t12
      t560 = 9472.D0/3.D0-39424.D0/3.D0*t2+(2560*t12-5120.D0/3.D0*t7-158
     #464.D0/3.D0*t5+81920*t3-26368.D0/3.D0+76288.D0/3.D0*z+4352.D0/3.D0
     #*t2-56320.D0/3.D0*t9)*t57+t398+(-2048.D0/3.D0*t5+8064*t7-18944+588
     #80.D0/3.D0*z-t498+19456*t3-205312.D0/3.D0*t9+189952.D0/3.D0*t2)*t1
     #8+(-20480.D0/3.D0*t12+3968*t7-29440.D0/3.D0*t3+t283+13696*t5)*t134
     #+(-27008.D0/3.D0*t5-t511+t424+14080.D0/3.D0-49408.D0/3.D0*t2+10752
     #*t3-7936*z+34048.D0/3.D0*t9)*t17-3328*z+(-2816.D0/3.D0*t7+16000.D0
     #/3.D0+35840.D0/3.D0*t5-60800.D0/3.D0*t3-31360.D0/3.D0*z+t523-32128
     #.D0/3.D0*t2)*l1+(58880.D0/3.D0*t2-t528+135040.D0/3.D0*t7-2048+t346
     #+9728.D0/3.D0*t5-33408*t9-47872.D0/3.D0*t3)*t87+1024.D0/3.D0*t5+(-
     #t386+43520.D0/3.D0*t12-51712.D0/3.D0*t7-64256.D0/3.D0*t5+t390+2534
     #4*t3)*t58-6144*t3+(25600.D0/3.D0*t12+15872.D0/3.D0-45824*t7-40064.
     #D0/3.D0*t5-t361+103168.D0/3.D0*t3+t435-39680*t2)*t102+(-10240.D0/3
     #.D0*t12+33920.D0/3.D0+38528.D0/3.D0*t7+213376.D0/3.D0*t5-21120*z-1
     #20704*t3+30080*t9-7040.D0/3.D0*t2)*t178+46208.D0/3.D0*t9
      s1 = (-t124-71680.D0/3.D0*t7+805888.D0/3.D0*t5+24832.D0/3.D0-96332
     #8.D0/3.D0*t3-60416.D0/3.D0*z+269312.D0/3.D0*t9+7936*t2)*t178-75776
     #.D0/3.D0*t3+(-44288*t12+126464.D0/3.D0*t7-35584.D0/3.D0*t3+t283+38
     #912.D0/3.D0*t5)*t134+128.D0/3.D0*t287*(-62*t18+24*t57-9-44*l1-71*t
     #17)*t331+(71424*t12+9728.D0/3.D0-223744.D0/3.D0*t7-509824.D0/3.D0*
     #t5-6656.D0/3.D0*z+217856*t3-29312.D0/3.D0*t9-114688.D0/3.D0*t2)*t1
     #02+(-106752*t12+474496.D0/3.D0*t7+71168.D0/3.D0*t2-2048+t346-21030
     #4.D0/3.D0*t3+89216.D0/3.D0*t5-98176.D0/3.D0*t9)*t87+t433*s1l1+t435
     #+(-148096.D0/3.D0*t5+156032.D0/3.D0*t3-16896*t2+15616.D0/3.D0*t12+
     #5376*t7+15488.D0/3.D0-28672.D0/3.D0*z-t442)*t17
      t562 = s1+14080+(19200*t12+82048.D0/3.D0-366592.D0/3.D0*t5+8320*t7
     #-41216.D0/3.D0*z+68352*t9-97792*t2+301312.D0/3.D0*t3)*t57-41472*t2
     #+(-180608.D0/3.D0+2432.D0/3.D0*t5+60800.D0/3.D0*t7-15616*t12-17868
     #8*t9+507904.D0/3.D0*t2+64512*t3+216064.D0/3.D0*z)*t18+(-t386+27673
     #6.D0/3.D0*t12-361088.D0/3.D0*t7-6144*t5+19712.D0/3.D0*t9+89216.D0/
     #3.D0*t3)*t58+(9216*t12-5120*t7+t395-t473)*t242+t473-51968.D0/3.D0*
     #z+(-37248*t3+12416.D0/3.D0-33280.D0/3.D0*z-13568.D0/3.D0*t2+32896.
     #D0/3.D0*t9+106496.D0/3.D0*t5-30976.D0/3.D0*t7)*l1+t560*s2l1
      s1 = 26112*t3+(60800-12544*t5-33280.D0/3.D0*t7+584960.D0/3.D0*t9-2
     #27840.D0/3.D0*t3+16384*t12-530048.D0/3.D0*t2-257152.D0/3.D0*z)*t18
     #-148864.D0/3.D0*t9+141824.D0/3.D0*t2+(-128.D0/3.D0*(-143*s1l1-488-
     #89*s2l1)*t9-80768.D0/3.D0*s1l1-25728-t28+256.D0/3.D0*(-33-t29)*t5-
     #128.D0/3.D0*(-18*s1l1-44-140*s2l1)*t3+128.D0/3.D0*(-725*s1l1-1181-
     #470*s2l1)*t2-128.D0/3.D0*(-1062*s1l1-1185-626*s2l1)*t1)*l1+(-32396
     #8.D0/3.D0*t12+t51+381056.D0/3.D0*t5+184192.D0/3.D0*t7-83200.D0/3.D
     #0*t9-176512.D0/3.D0*t3)*t58+((35456+99328.D0/3.D0*s1l1-23552*s2l1)
     #*t5+2048*t65*t1+(-484864.D0/3.D0-292096.D0/3.D0*s1l1-196864.D0/3.D
     #0*s2l1)*t9-512.D0/3.D0*(-426*s1l1-384-318*s2l1)*t2+97792.D0/3.D0*t
     #12+(97024-187648.D0/3.D0*s1l1+70400.D0/3.D0*s2l1)*t3+(-142208.D0/3
     #.D0*s1l1-64896-78080.D0/3.D0*s2l1)*t7)*t87-128.D0/3.D0*(-73*s1l1-5
     #4-54*s2l1)*t2+(123904.D0/3.D0*t7-14848.D0/3.D0-t95+449408.D0/3.D0*
     #t5-71296*t12+5632.D0/3.D0*z-684160.D0/3.D0*t9+114560*t2)*t102
      t564 = s1+(2048*t104*t2+(72320.D0/3.D0+225664.D0/3.D0*s1l1+66688.D
     #0/3.D0*s2l1)*t3+256.D0/3.D0*(205*s1l1+349+205*s2l1)*t9+(-306560.D0
     #/3.D0-110848.D0/3.D0*s1l1-18560.D0/3.D0*s2l1)*t5+(119680.D0/3.D0*s
     #1l1+255232.D0/3.D0+109952.D0/3.D0*s2l1)*t7-t124)*t58-4480*s1l1+(-t
     #128-17664*t7+59904*t12+32768*t3-72960*t5)*t134+128.D0/3.D0*(-38*s1
     #l1-111-111*s2l1)*t9+((192256.D0/3.D0-13184.D0/3.D0*s1l1+t28)*t5-51
     #2.D0/3.D0*(-41*s1l1-5-41*s2l1)*t1-54272.D0/3.D0*t12+(55936.D0/3.D0
     #+16384*s2l1+34176*s1l1)*t7+(779008.D0/3.D0+207360*s1l1+214784.D0/3
     #.D0*s2l1)*t9+(-324608.D0/3.D0-161024*s1l1-73472*s2l1)*t2+3072*s1l1
     #+3072+3072*s2l1+(-220032-124160.D0/3.D0*s1l1-40192*s2l1)*t3)*t102+
     #(-167168.D0/3.D0+21120*t7-639616.D0/3.D0*t5+104192.D0/3.D0*t12-704
     #0.D0/3.D0*t9+173184*t3+259328.D0/3.D0*z+44288*t2)*t178+t267*s1l1-t
     #269-t270*t562
      t601 = t270**2
      t634 = 256*t12
      t637 = 16384.D0/3.D0*t7
      t638 = 41984.D0/3.D0*t5
      t646 = 1280*z
      t663 = -t12+t7
      t665 = 256*t663*t242
      t671 = 256.D0/3.D0*t12
      t674 = 6656*t7
      t696 = 7936.D0/3.D0*t7
      t703 = (-t634+8704.D0/3.D0*t9+115712.D0/3.D0*t2-21248+t637+t638+25
     #600*z-39424*t3)*t57+(-16640.D0/3.D0*t5-3584.D0/3.D0*t2+37376.D0/3.
     #D0*t3+1280-t646-20992.D0/3.D0*t9)*l1+(-t424-173312.D0/3.D0*t3+9625
     #6.D0/3.D0*t7+48896*t9-8448*t2-14336*t5)*t102+(-1792*t12+3840*t3+20
     #992*t7-23040*t5)*t58+t665+(t367+t95-116224.D0/3.D0*t7-7168*t9+4326
     #4*t5)*t87+(t671+101632.D0/3.D0*t9+67328.D0/3.D0-21760.D0/3.D0*t3-t
     #674+6400*t5-48896*t2-28416*z)*t18-16640.D0/3.D0+6912*z+(t634+10163
     #2*t3-34048.D0/3.D0*t7-224512.D0/3.D0*t9+27904.D0/3.D0*t2+3584-3584
     #*z-85760.D0/3.D0*t5)*t178+12544*t2+(3328.D0/3.D0*t12-14080.D0/3.D0
     #*t7+3584*t5)*t134+2816*t3+(t696-512+12800.D0/3.D0*t5+768*z+t523-56
     #32.D0/3.D0*t2-t245)*t17-29440.D0/3.D0*t9
      t721 = 1536*t7
      t722 = 1024*z
      t724 = 60928.D0/3.D0*t3
      t738 = 512.D0/3.D0*t9
      t759 = 33536.D0/3.D0*t5
      t797 = (-512*t12-14336.D0/3.D0+3584.D0/3.D0*t7-t738+38912.D0/3.D0*
     #t2+15104.D0/3.D0*z+34816.D0/3.D0*t5-t724)*t57+(512.D0/3.D0*t7-t442
     #-t396+256.D0/3.D0*t2-256.D0/3.D0+256.D0/3.D0*z+16640.D0/3.D0*t3)*l
     #1+13312.D0/3.D0*t2+(-8704.D0/3.D0*t12-256*t3+t674-10496.D0/3.D0*t5
     #)*t58+t665+(3072*t12+3584.D0/3.D0*t3-46336.D0/3.D0*t7+t759)*t87+(1
     #024.D0/3.D0*t12+21248.D0/3.D0-58624.D0/3.D0*t2-5376*t3-t696+55040.
     #D0/3.D0*t9+1792*t5-22784.D0/3.D0*z)*t18-256.D0/3.D0*t5-4352.D0/3.D
     #0+4352.D0/3.D0*z+5120.D0/3.D0*t3-4096.D0/3.D0*t663*t134+(-t498-565
     #76.D0/3.D0*t3+44288.D0/3.D0*t7+26368.D0/3.D0*t9-256*t2-8192.D0/3.D
     #0*t5)*t102+(-t671-768+4864.D0/3.D0*t2-16384.D0/3.D0*t3+1024*t7+435
     #2.D0/3.D0*t9+t246+t722)*t17+(t249+125440.D0/3.D0*t3-4608*t7-63488.
     #D0/3.D0*t9+2048.D0/3.D0*t2-17408*t5)*t178-4608*t9
      t809 = (-5120*t12-23552*t3-17920-97792.D0/3.D0*t9+146944.D0/3.D0*t
     #2-4096.D0/3.D0*t7+94720.D0/3.D0*t5+20224*z)*t57+(-20480*t12-89600*
     #t3+75008.D0/3.D0*t7+37376*t9-2560*t2+150784.D0/3.D0*t5)*t102+(4352
     #*t12+53248.D0/3.D0+146432.D0/3.D0*t9-144896.D0/3.D0*t2-17408*t3-20
     #480.D0/3.D0*t7+4864.D0/3.D0*t5-20992*z)*t18+(-73472.D0/3.D0*t12-17
     #92*t3+124160.D0/3.D0*t7-15104*t5)*t58-12544.D0/3.D0+t703*s1l1+(107
     #52*t12-14592*t7+3840*t5)*t134+(24832.D0/3.D0*t12+332032.D0/3.D0*t3
     #+22528.D0/3.D0*t7-119552.D0/3.D0*t9-21248.D0/3.D0*t2+3840-3840*z-2
     #50112.D0/3.D0*t5)*t178+(-4352.D0/3.D0*t12+46592.D0/3.D0*t5-t721-76
     #8+t722+2816.D0/3.D0*t2-t724+22784.D0/3.D0*t9)*t17+(30208*t12+74752
     #.D0/3.D0*t3-54528*t7-7680*t9+21248.D0/3.D0*t5)*t87+7168*t3+t797*s2
     #l1-40448.D0/3.D0*t9+2048*t663*t242+(t385-t738-9984*t5+t511-11776.D
     #0/3.D0*t2+1280-t646)*l1-4352.D0/3.D0*t5+4864*z+35840.D0/3.D0*t2
      s1 = -5888.D0/3.D0+(-t638+t721-74368.D0/3.D0*t9+t429-51712.D0/3.D0
     #+113152.D0/3.D0*t2+107264.D0/3.D0*z)*l1-512*t5-48512.D0/3.D0*t9+t2
     #19+20992.D0/3.D0*t3+(-t237-3712*t12+t239+3456*t7)*t242+(-4224*t12-
     #768*t7+227584.D0/3.D0*t5-296192.D0/3.D0*t3+179200.D0/3.D0+180992.D
     #0/3.D0*t9-301312.D0/3.D0*t2-322688.D0/3.D0*z)*t57
      t901 = s1+(40832.D0/3.D0*t12+27136*t3-t128-27904.D0/3.D0*t7-47488*
     #t5)*t134+13824*t2+(t759-28672.D0/3.D0+5632.D0/3.D0*t7-43648*t9+202
     #24.D0/3.D0*t3-t249+118400.D0/3.D0*t2+69760.D0/3.D0*z)*t17+20224.D0
     #/3.D0*z+(25216*t12-268160.D0/3.D0*t5-52224*t2+4096-t181+27392*t3-3
     #7632*t7+315776.D0/3.D0*t9)*t87+(-20224.D0/3.D0*t7-10752*t5-70400.D
     #0/3.D0*t3+13184+81536*t9-193024.D0/3.D0*t2+8192.D0/3.D0*t12-73088.
     #D0/3.D0*z)*t18+(-75392.D0/3.D0*t12+t51+275072.D0/3.D0*t5+53248.D0/
     #3.D0*t7-t257-62080*t3)*t58+(-28928.D0/3.D0*t7-111488.D0/3.D0+19712
     #0.D0/3.D0*t3-271360.D0/3.D0*t5+17536.D0/3.D0*t12+186752.D0/3.D0*z+
     #143104.D0/3.D0*t9+52480.D0/3.D0*t2)*t178+(118400.D0/3.D0*t7-33280.
     #D0/3.D0+92800.D0/3.D0*t3+61568*t5-42112.D0/3.D0*t12+t224-578176.D0
     #/3.D0*t9+307072.D0/3.D0*t2)*t102
      s1 = -68864.D0/3.D0+((-121984.D0/3.D0-18816*s1l1-11648*s2l1)*t3-18
     #304*s1l1-18688-13696*s2l1-128.D0/3.D0*(-t29-18)*t7-128.D0/3.D0*(-7
     #8*s1l1-225-53*s2l1)*t5+(32384+31616*s1l1+22912*s2l1)*t1+(-146816.D
     #0/3.D0-55936*s1l1-28928*s2l1)*t2+(60160+159488.D0/3.D0*s1l1+66304.
     #D0/3.D0*s2l1)*t9)*t17+(47104*t3+10752*t7-79360.D0/3.D0+51840*z+166
     #528.D0/3.D0*t2-111104.D0/3.D0*t5-132608.D0/3.D0*t9)*l1+t601*t809/2
     #+(-1024*t104*t3-512*(26*s2l1+26*s1l1+38)*t5+(19712+9728*s1l1+9728*
     #s2l1)*t7-t238)*t242+t901*s2l1+(107264.D0/3.D0*t5-21248.D0/3.D0*t3+
     #51968*t2-t637-42112.D0/3.D0-16384.D0/3.D0*t12+98176.D0/3.D0*z-5030
     #4*t9)*t17+25856*z+(-t237-13952*t12+18432*t5-1408*t7)*t242
      s3 = s1+(-35328*s1l1-31744-26112*s2l1+(-321664.D0/3.D0-175616*s1l1
     #-73472.D0/3.D0*s2l1)*t9+(83968+61568*s1l1+36480*s2l1)*t1+(493184.D
     #0/3.D0+333440.D0/3.D0*s1l1+116480.D0/3.D0*s2l1)*t3+13568.D0/3.D0*t
     #12+(-123136.D0/3.D0+72704*s1l1-t269)*t2+(-257920.D0/3.D0-115712.D0
     #/3.D0*s1l1-62720.D0/3.D0*s2l1)*t5+(-27520.D0/3.D0*s1l1-7168*s2l1+3
     #3152.D0/3.D0)*t7)*t178
      s2 = s3+(13312+3840*s1l1+18688*s2l1+128.D0/3.D0*(-82*s2l1-138-58*s
     #1l1)*t7+(-46336.D0/3.D0-77824.D0/3.D0*s1l1-37120.D0/3.D0*s2l1)*t3+
     #(65792.D0/3.D0+70912.D0/3.D0*s1l1+48640.D0/3.D0*s2l1)*t9+(42496.D0
     #/3.D0+20608.D0/3.D0*s1l1+2688*s2l1)*t5+t424+(-62848.D0/3.D0-93056.
     #D0/3.D0*s1l1+512*s2l1)*t2+(-3456+20096*s1l1-18560*s2l1)*t1)*t18+(6
     #6432+234368.D0/3.D0*s1l1+130688.D0/3.D0*s2l1+(623488.D0/3.D0+40960
     #0.D0/3.D0*s1l1+234496.D0/3.D0*s2l1)*t2+(7424+36224.D0/3.D0*s1l1+36
     #608.D0/3.D0*s2l1)*t5+(53504.D0/3.D0-12416.D0/3.D0*s1l1-2816*s2l1)*
     #t3+(-120832-35200*s1l1-38400*s2l1)*t9-2560.D0/3.D0*t12+(-165504-16
     #2304*s1l1-227200.D0/3.D0*s2l1)*t1+(8192.D0/3.D0*s1l1-1152+8576.D0/
     #3.D0*s2l1)*t7)*t57-128.D0/3.D0*(-38-38*s2l1)*t3
      t1047 = s2+128.D0/3.D0*(19*s1l1+124+124*s2l1)*t1+(109824*t12+4096-
     #t181-80768*t7-410240.D0/3.D0*t5+9728*t3-63488*t2+472064.D0/3.D0*t9
     #)*t87+((239104.D0/3.D0+113536.D0/3.D0*s1l1+93568.D0/3.D0*s2l1)*t5+
     #(-82432.D0/3.D0*s1l1-91648.D0/3.D0*s2l1-189952.D0/3.D0)*t7+t528+25
     #6.D0/3.D0*(-438*s1l1-396-330*s2l1)*t3-1024*t65*t9)*t134-15872.D0/3
     #.D0*t5+(-18176*t7+157696.D0/3.D0+440320.D0/3.D0*t5-22144*t12+51968
     #*t9-142592*t3-108928*z-238720.D0/3.D0*t2)*t57
      t1050 = (l1+z-t301)**2
      taux18 = (t564+t1047)/t1050*t304

CCC Coefficients

      tcoff33 = taux17
      tcoff34 = taux17
      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet33 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff33)
      tFjet34 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff34)
      tFjet35 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff35)
      tFjet36 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff36)

CCC Sum

      qgRVnum9e2 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFjet
     #36*tcoff36
      RETURN
      END

      double precision function qgRVnum9e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = s1l1*s2l1
      t2 = 128*t1
      t3 = 1-z
      t4 = t3**2
      t5 = t4*t3
      t6 = t4**2
      t7 = t6*t5
      t8 = 192*t7
      t9 = 2*t1
      t10 = s1l1**2
      t11 = s2l1**2
      t13 = t6*t3
      t16 = t6*t4
      t17 = 192*t16
      t19 = l1**2
      t20 = t19**2
      t21 = t20**2
      t22 = t21*t19
      t26 = 6*t1
      t36 = 1472*s2l1
      t44 = t21*l1
      t46 = 896.D0/3.D0*t1
      t48 = 448.D0/3.D0*s1l1
      t56 = 448.D0/3.D0*t11
      t64 = 448*s1l1
      t65 = 320.D0/3.D0*t10
      t69 = 1792.D0/3.D0*t1
      t70 = 1088.D0/3.D0*t11
      t78 = 640*t10
      t79 = 896*s2l1
      t80 = 1408*t1
      t83 = 960*t10
      t87 = 320.D0/3.D0*t7
      t88 = (-t46+3200.D0/3.D0*t11-t48-1600.D0/3.D0*t10+1920+2432*s2l1)*
     #t6+(1280.D0/3.D0*s2l1+3520.D0/3.D0)*t16+(-t56+2112*s1l1+5440.D0/3.
     #D0*t10+4864.D0/3.D0-2816.D0/3.D0*s2l1)*t3+(-1408.D0/3.D0*t1-7744.D
     #0/3.D0*t11+t64+t65-11008.D0/3.D0-17536.D0/3.D0*s2l1)*t5+(-t69-t70+
     #512.D0/3.D0*s1l1-t65-1216-2624.D0/3.D0*s2l1)*t13-2048.D0/3.D0*t1-2
     #048.D0/3.D0*t11-1984.D0/3.D0*s1l1-t78-3520.D0/3.D0-t79+(t80+2368*t
     #11-1856*s1l1-t83+4160.D0/3.D0+5376*s2l1)*t4-t87
      t92 = 320.D0/3.D0*t16
      t93 = 7360.D0/3.D0*t4
      t94 = 256.D0/3.D0*t13
      t95 = 6784.D0/3.D0*z
      t96 = 320*t6
      t98 = t20*l1
      t102 = (-256*t13+256*t16)*t22
      t103 = 64.D0/3.D0*t16
      t104 = 1984.D0/3.D0*t5
      t106 = 1472*t4
      t110 = 256.D0/3.D0*t16
      t112 = 640.D0/3.D0*t13
      t113 = 640*z
      t114 = 2560.D0/3.D0*t4
      t115 = 4480.D0/3.D0*t6
      t117 = t20*t19
      t119 = 1984.D0/3.D0*t16
      t120 = 576*t6
      t123 = 128.D0/3.D0*t5
      t124 = 1280.D0/3.D0*t4
      t125 = 768*z
      t128 = 64.D0/3.D0*t4
      t130 = 3904.D0/3.D0*t13
      t135 = 128*z
      t138 = 64.D0/3.D0*z
      t139 = 128.D0/3.D0*t16
      t140 = 128*t13
      t145 = t16*t22
      t149 = 128.D0/3.D0*t13
      t153 = t19*l1
      t155 = 9344.D0/3.D0*t5
      t156 = 896.D0/3.D0*t13
      t158 = 256*t4
      t162 = 1280.D0/3.D0*t16
      t163 = 640*t6
      t167 = 512*t4
      t169 = 128.D0/3.D0*t6
      t172 = 128*t5
      t173 = 9088.D0/3.D0*t4
      t174 = 512.D0/3.D0*t13
      t175 = 9728.D0/3.D0*z
      t176 = 896.D0/3.D0*t6
      t179 = 2048.D0/3.D0*t13
      t180 = 640*t5
      t181 = 2944.D0/3.D0*t6
      t182 = 1408.D0/3.D0*t16
      t185 = 1024.D0/3.D0*t4
      t186 = 896.D0/3.D0*z
      t189 = 768*t6
      t190 = 512.D0/3.D0*t16
      t191 = 1152*t4
      t194 = t20*t153
      t196 = (-t139-t123+t140-1408*t4+5632.D0/3.D0-2560.D0/3.D0*z)*t20-1
     #28*t145+(2432.D0/3.D0*t5-2048.D0/3.D0*t4+t149-896.D0/3.D0-3968.D0/
     #3.D0*z-1280.D0/3.D0*t6)*t153+(-t139-t155-t156+1280-1280*z+t158+358
     #4.D0/3.D0*t6)*t117+(t162+t163-384*t13)*t44+(-t167+256.D0/3.D0-2176
     #.D0/3.D0*z+t169)*t19+(-t172+t110+t173-t174-8576.D0/3.D0+t175+t176)
     #*t98+256.D0/3.D0+t135+(t179-t180-t181-t182)*t21+(-t123+t185-1280.D
     #0/3.D0+t186)*l1+(-t189+t190-t191+9472.D0/3.D0*t5)*t194
      t198 = 1024*t4
      t205 = 4864.D0/3.D0*z
      t208 = (-t92+t93-t94-5824.D0/3.D0+t95-t96)*t98+t102+(-t96+t103+t10
     #4+832.D0/3.D0*t13-t106+1280-1088.D0/3.D0*z)*t20+(t110-9856.D0/3.D0
     #*t5+t112+640-t113+t114+t115)*t117+(-t119+t120+t94)*t44+(t123-t124-
     #128.D0/3.D0-t125)*t19-t128+(1472.D0/3.D0*t16+t130-1856*t6+64*t5)*t
     #21+(t123-t128-192-t135)*l1-t138+t196*s2l1+(t163-t110-t198+5888.D0/
     #3.D0*t5-4480.D0/3.D0*t13)*t194+(t180-t149-2944.D0/3.D0*t4-640.D0/3
     #.D0*t6+256-t205)*t153
      t216 = t2+(t8-256*(t9+t10+t11)*t13-t17)*t22+(64.D0/3.D0*(-27*t11+9
     #*s1l1-t26+33-39*s2l1+21*t10)*t6+(6656.D0/3.D0*t1+4096.D0/3.D0*t11-
     #192*s1l1+2560.D0/3.D0*t10+256+t36)*t13-1088.D0/3.D0*t7+(-640*s2l1-
     #1792.D0/3.D0)*t16)*t44+t88*t20-512.D0/3.D0*s2l1+t64+1792*z+t208*s1
     #l1-896.D0/3.D0*t11+704.D0/3.D0*t10+64.D0/3.D0*(5*s2l1+5*t11)*t5
      t230 = 1024.D0/3.D0*t13
      t231 = dlog(z)
      t232 = z**2
      t233 = t232**2
      t236 = z*l1
      t238 = (1-l1+t236)**2
      t239 = 1/t238
      t241 = 64.D0/3.D0*t233*t232*l1*t239
      t244 = 256*t98
      t249 = 512*t98
      t262 = 2*l6-1
      t263 = t262**2
      t266 = l5-1
      t268 = -1+l1
      t273 = dsqrt(z)
      t275 = dsqrt(l5)
      t277 = dsqrt(-t266)
      t280 = dsqrt(l1)
      t284 = l5*z*l1-l1+1-l5+l5*l1
      t287 = dsqrt(-t268)
      t288 = t239*t287
      t292 = t284**2
      t294 = 4*l1*t263*z*l5*t266*t268*t239-4*t262*t273*t275*t277*t280*t2
     #84*t288+t292*t239
      t295 = dsqrt(t294)
      t299 = 256*t6
      t304 = 7424.D0/3.D0*t6
      t305 = 128*t7
      t314 = 3200.D0/3.D0*t16
      t316 = 6784.D0/3.D0*t5
      t325 = 1536*t16
      t328 = 128.D0/3.D0*t7
      t333 = 640.D0/3.D0*t7
      t345 = 1280*t13
      t350 = 2048.D0/3.D0*t16
      t354 = 16256.D0/3.D0*t4
      t355 = 11648.D0/3.D0*t6
      t357 = 256*t16
      t358 = 33664.D0/3.D0*t4
      t365 = 896*t7
      t369 = 1408*t13
      t370 = 3328.D0/3.D0*t6
      t371 = 9088.D0/3.D0*t5
      t372 = 128*t16
      t374 = 3072*z
      t377 = 512*t16
      t379 = 30848.D0/3.D0*t5
      t385 = 384*t16
      t386 = 1664.D0/3.D0*t7
      t387 = 512*t5
      t391 = 3968.D0/3.D0*t4
      t395 = 10240.D0/3.D0*t13
      t402 = 2304*z
      t403 = 1024*t6
      t404 = 7168.D0/3.D0*t5
      t405 = 896*t4
      t411 = 640*t7
      t412 = 19072.D0/3.D0*t5
      t413 = 1408*t16
      t414 = 1024*z
      t417 = 512*t6
      t418 = 512*t13
      t421 = 6272.D0/3.D0*t4
      t422 = (t355-3840*t5+t357-t358+t345+10752-63872.D0/3.D0*z-t305)*t2
     #0+(-t167+3712*t13+1280*t16+5248.D0/3.D0*t6-t365+14336.D0/3.D0*t5)*
     #t21+(-t369+t370-t371-t372+12800.D0/3.D0*t4-9856-t374+t328)*t153+(t
     #377+21376.D0/3.D0*t13-16384.D0/3.D0-t333-t379-6272.D0/3.D0*t6+4288
     #0.D0/3.D0*t4+14848.D0/3.D0*z)*t117+(-t385+t386+t149-t304-t387)*t44
     #+(-t371-8192.D0/3.D0+2560.D0/3.D0*t13+t391+t176-11392.D0/3.D0*z)*t
     #19+(-t372-t395-5504*t6+t305+19200*t5-8576*t4-6272.D0/3.D0+25984.D0
     #/3.D0*z)*t98+1536-t402+t180+(-t403+t404-t94-t405-7936.D0/3.D0-t113
     #)*l1+(10624.D0/3.D0*t6-22784.D0/3.D0*t13-3328.D0/3.D0*t4+t411-t412
     #-t413+1024-t414)*t194+(t417-t305-t418)*t22-t421
      t427 = (t385-t305-t418+t299)*t22
      t428 = 1792.D0/3.D0*t13
      t430 = 768*t5
      t435 = 39808.D0/3.D0*t6
      t440 = 1536*z
      t446 = (-4736+t305+143360.D0/3.D0*t5+2816*t16-16128*t4+4096*t13+14
     #848*z-99584.D0/3.D0*t6)*t98+(-t314-2176*t13+t316+29696.D0/3.D0*t6-
     #35328*t4+70016.D0/3.D0-45056*z-t305)*t20+(-9728.D0/3.D0*t13+3456*t
     #6-4864*t5-29824.D0/3.D0+t325+11648.D0/3.D0*t4-22784.D0/3.D0*z+t328
     #)*t153+(-8064*t16+12032*t13-9728-t333-139648.D0/3.D0*t5+12416*t6+1
     #20320.D0/3.D0*t4+8960*z)*t117+3072*t5+(-4224*t6+13696.D0/3.D0*t5-3
     #712.D0/3.D0*t4-8960.D0/3.D0-2816.D0/3.D0*z+t345)*l1+(4608*t6-20480
     #.D0/3.D0*t5-t350-t149-t191-384-11008*z)*t19-t354+8704.D0/3.D0+t422
     #*s2l1-2048.D0/3.D0*t6-4608*z+t427+(-t357+t386-t428+3200.D0/3.D0*t6
     #-t430)*t44+(-10624.D0/3.D0*t16+32768.D0/3.D0*t13-t435-t365+6784*t5
     #)*t21+(26624.D0/3.D0*t16+1536-t440+59008.D0/3.D0*t6-21760*t13+t411
     #-t430-8192*t4)*t194
      t449 = 3008.D0/3.D0*t16
      t457 = 1408.D0/3.D0*t6
      t459 = 256*t5
      t463 = t233*z
      t464 = t463*t268
      t465 = l6**2
      t468 = 1024.D0/3.D0*t464*t239*t465
      t469 = t239*l6
      t471 = 1024.D0/3.D0*t464*t469
      t473 = 2*t236
      t475 = t232*t19
      t477 = t232*z
      t487 = z*t19
      t489 = t232*l1
      t491 = -4*z-t473+l1-2*t19-14*t475+6*t477*t19-4*z*t153+6*t232*t153-
     #4*t477*t153+t233*t153+10*t487+9*t489+t153
      t495 = -t468+t471-64.D0/3.D0*t233*t491*t239
      t496 = l5**2
      t498 = 15488.D0/3.D0*t5
      t499 = 2880-t241+(9152.D0/3.D0*t19+256.D0/3.D0*t20-512.D0/3.D0-t24
     #4+20800.D0/3.D0*t153+(7808.D0/3.D0*t20+512.D0/3.D0+8512.D0/3.D0*t1
     #9+5184*t153-t249+1472*l1)*s1l1+(1856*t19+5696*t153-1920*t20-384+17
     #28*l1)*s2l1+6208.D0/3.D0*l1)*t295-11072.D0/3.D0*z+(-1536*t13+t299+
     #t17+1088*t7)*t22+1600.D0/3.D0*t13-t304+t446*s1l1+(-t106-3328.D0/3.
     #D0-7040.D0/3.D0*z-t449+3264*t5+3072*t13-4352*t6)*l1-17984.D0/3.D0*
     #t4+(5888.D0/3.D0*t16-t457-14656.D0/3.D0*t7-t459+3648*t13)*t44+t495
     #*t496+t498
      t501 = 512*z
      t510 = 4352.D0/3.D0*t7
      t511 = 4160.D0/3.D0*t13
      t519 = 2176.D0/3.D0*t16
      t520 = 14656.D0/3.D0*t6
      t526 = t273*t463*t280
      t531 = t277*t287
      t535 = -512.D0/3.D0*t526*t277*t288*l6+256.D0/3.D0*t526*t531*t239
      t545 = t531*t280
      t549 = (t489+3*z-t473+l1-1)*t273*t233
      t556 = 512.D0/3.D0*t545*t549*t469-256.D0/3.D0*t545*t549*t239
      t557 = t275*l5
      t561 = t475+3*t236-2*t487+t19+l1-2
      t572 = 5696.D0/3.D0*t7
      t585 = 64.D0/3.D0*t7
      t586 = 14464.D0/3.D0*t6
      t595 = 1984.D0/3.D0*t13
      t597 = 4544*t5
      t602 = 19072.D0/3.D0*t13
      t604 = 256.D0/3.D0*t7
      t606 = 1216.D0/3.D0*t16
      t611 = 7424.D0/3.D0*t5
      t614 = 4736.D0/3.D0*z
      t620 = 4864.D0/3.D0*t4
      t622 = 2176.D0/3.D0*t13
      t625 = 64*t16
      t635 = 832.D0/3.D0*z
      t636 = 3968.D0/3.D0*t6
      t639 = 6016.D0/3.D0*t4
      t640 = 64*t7
      t641 = 256*t13
      t643 = (t299-t640-t641)*t22
      t644 = 6400.D0/3.D0*t13
      t645 = 2624.D0/3.D0*t6
      t646 = 704.D0/3.D0*t7
      t649 = (-23296.D0/3.D0-t585-t412-t586+t110+22528.D0/3.D0*z+19264*t
     #4+21568.D0/3.D0*t13)*t117+(t333-t17-t387+448.D0/3.D0*t13-t181)*t44
     #+(t595+13120.D0/3.D0*t6-t597+t110-9024*t4-t585+21440.D0/3.D0-20736
     #*z)*t20+(-t602+1024-t414+23488.D0/3.D0*t6+t604-33856.D0/3.D0*t5-t6
     #06-4736.D0/3.D0*t4)*t194+(4352.D0/3.D0*t6-t611+5504.D0/3.D0*t4-396
     #8.D0/3.D0*t13-t103-4800-t614)*t153+3328.D0/3.D0-2240*z+t104+(384*t
     #6-7616.D0/3.D0*t5+t620-5312.D0/3.D0-3392*z+t622)*t19+(-t625-8768.D
     #0/3.D0*t13-18880.D0/3.D0*t6+64576.D0/3.D0*t5-60224.D0/3.D0*t4+t328
     #+17728.D0/3.D0-1856*z)*t98+(-1600*t4+8192.D0/3.D0*t5-2432.D0/3.D0-
     #t635-t636)*l1-t639+t643+(-t158+t644+t377-t645-t646+t498)*t21
      t652 = 704.D0/3.D0*t13
      t661 = 256*t7
      t666 = 38144.D0/3.D0*t13
      t669 = 2560.D0/3.D0*t7
      t676 = 15872.D0/3.D0*t5
      t678 = 512.D0/3.D0*t7
      t679 = 7808*z
      t682 = 2560.D0/3.D0*t16
      t683 = 1664.D0/3.D0*t13
      t685 = 1024.D0/3.D0*t7
      t699 = 2048.D0/3.D0*t7
      t706 = 896.D0/3.D0*t16
      t707 = 17024.D0/3.D0*z
      t710 = 6400.D0/3.D0*t16
      t718 = 1536*t7
      t723 = (t190+t149+57088.D0/3.D0+13952.D0/3.D0*t6+10240.D0/3.D0*t5-
     #25728*t4-t661-71680.D0/3.D0*z)*t20+(-2048*t16+t666+45568.D0/3.D0*t
     #4+22400.D0/3.D0*t6-t669+1024.D0/3.D0-3328.D0/3.D0*z-98432.D0/3.D0*
     #t5)*t117+2304+(t357-3328*t13-32384.D0/3.D0+3968*t6-t676+5376*t4+t6
     #78-t679)*t153+(t682-t683-61312.D0/3.D0*t6+t685+26752*t5+27392.D0/3
     #.D0*t4-47488.D0/3.D0+72448.D0/3.D0*z)*t98+6016.D0/3.D0*t5+(-11648.
     #D0/3.D0-7040.D0/3.D0*t6+t404+t230+1792.D0/3.D0*t4+t139-2944.D0/3.D
     #0*z)*l1-3456*t4+t427+t149-t95-t457+(-1664*t16+t699+2816.D0/3.D0*t1
     #3+t176-t459)*t44+(2944.D0/3.D0*t13+9088.D0/3.D0*t6-18688.D0/3.D0*t
     #5-t706+t124-4096.D0/3.D0-t707-t328)*t19+(t710+3712.D0/3.D0*t5-t510
     #+16000.D0/3.D0*t13-21760.D0/3.D0*t6)*t21+(t190+512-t501+10880*t6-1
     #6000*t13+t718+25600.D0/3.D0*t5-5632*t4)*t194
      t725 = 64.D0/3.D0*t13
      t726 = 16256.D0/3.D0*t13
      t727 = 1280.D0/3.D0*t7
      t734 = 1856.D0/3.D0*z
      t737 = 5312.D0/3.D0*z
      t740 = 7360.D0/3.D0*t13
      t748 = 1600.D0/3.D0*t6
      t751 = 256*z
      t755 = 1088*t6
      t756 = 2176.D0/3.D0*t5
      t758 = 21760.D0/3.D0*z
      t762 = 768*t16
      t763 = 2176.D0/3.D0*t7
      t764 = 384*t5
      t767 = 1280.D0/3.D0*t5
      t768 = 2752.D0/3.D0*t6
      t769 = 8320.D0/3.D0*t5
      t770 = 1088*t4
      t773 = 320.D0/3.D0*t13
      t776 = 19136.D0/3.D0*t13
      t777 = 768*t7
      t782 = t725+(t377+t726+2304-t727-31616.D0/3.D0*t5+3008*t6-1088.D0/
     #3.D0*t4-2560*z)*t117+960-t734+(-t372-t130+t604+t163-1216*t5+3136*t
     #4-16256.D0/3.D0-t737)*t153+(-t740-5440*t6+27712.D0/3.D0*t5+5696*t4
     #+t678-7744+30784.D0/3.D0*z)*t98-3136.D0/3.D0*t4+(-6464.D0/3.D0-t74
     #8+2240.D0/3.D0*t5-t174+832*t4+t103+t751)*l1+(t110+4288.D0/3.D0*t13
     #+8256+t755+t756-27712.D0/3.D0*t4-t305-t758)*t20+(-t158+9344.D0/3.D
     #0*t13+t762+t755-t763-t764)*t21+t767+t643+(t179-3136.D0/3.D0+t768-t
     #769+t770-t585-t751)*t19+(t685-t773-t17-t115)*t44+(-t314-t776+t748+
     #t777+3648*t5+1408.D0/3.D0*t4)*t194-t169
      t787 = 256*t20
      t790 = (7552*t6+512-t501-6464.D0/3.D0*t13+26816.D0/3.D0*t16-t421-3
     #0272.D0/3.D0*t7-8000.D0/3.D0*t5)*t194+(24640.D0/3.D0*t6+t358-25024
     #.D0/3.D0-6592.D0/3.D0*z-t510-t511-14656*t5+960*t16)*t153+(-6784.D0
     #/3.D0*t4+1472.D0/3.D0*t7-3392*t5+t519+t520+1664.D0/3.D0-19840.D0/3
     #.D0*z-10048.D0/3.D0*t13)*t19+t535*t275+(-72064.D0/3.D0*t5-32960.D0
     #/3.D0*t13-10496.D0/3.D0+t175-6848.D0/3.D0*t16+46208.D0/3.D0*t4+672
     #0*t7+56000.D0/3.D0*t6)*t117+t556*t557+(t468-t471+128.D0/3.D0*t561*
     #t463*t239)*l5+(-26112*t4+10432*t5+1920*t16+6208*t6+44416.D0/3.D0-6
     #3296.D0/3.D0*z+t572-28480.D0/3.D0*t13)*t20+(128.D0/3.D0*t4+23680*t
     #5+64256.D0/3.D0*t13-95744.D0/3.D0*t6-17728.D0/3.D0+29888.D0/3.D0*z
     #-11840.D0/3.D0*t16-9344.D0/3.D0*t7)*t98+t649*t10+(-t158+9344*t7+t6
     #52-22400.D0/3.D0*t16-19648.D0/3.D0*t6+14080.D0/3.D0*t5)*t21+t723*s
     #2l1+t782*t11+(-4480.D0/3.D0*t153-704*t19+64-1408.D0/3.D0*l1+t787)*
     #t294
      t793 = 4096.D0/3.D0*t6
      t796 = 384*s1l1
      t800 = 1664.D0/3.D0*t1
      t803 = 640.D0/3.D0*t10
      t806 = 640*t1
      t807 = 448*t11
      t811 = 640.D0/3.D0*t1
      t812 = 576*s1l1
      t817 = 320*s1l1
      t826 = 640.D0/3.D0*t16
      t827 = (512*t1+1280.D0/3.D0*t11+t796+t78+2048.D0/3.D0+896.D0/3.D0*
     #s2l1)*t4+(-t800-1600.D0/3.D0*t11-512*s1l1-t803-512-t79)*t3+t806+t8
     #07+2432.D0/3.D0*s1l1+832.D0/3.D0*t10+1792.D0/3.D0+1024*s2l1+(t811+
     #t56-t812-t78+128.D0/3.D0+2752.D0/3.D0*s2l1)*t5+(-1280.D0/3.D0*t1-t
     #70+t817-5120.D0/3.D0-4096.D0/3.D0*s2l1)*t6-64.D0/3.D0*(-65-5*t11-2
     #5*s2l1)*t13-t826
      t832 = 320*t7
      t842 = 2752.D0/3.D0*s1l1
      t863 = 3*s1l1
      t864 = 15*s2l1
      t870 = 1024.D0/3.D0*t1
      t883 = 1280.D0/3.D0*t10
      t905 = (t870+2496*t11+128*s1l1-t803-2368.D0/3.D0+9728.D0/3.D0*s2l1
     #)*t4+(-t806-2624*t11-256*s1l1-t65+320-3200*s2l1)*t3+(768*t1+2240.D
     #0/3.D0*t11-t48+t883+4480.D0/3.D0+4544.D0/3.D0*s2l1)*t6+(t811-64.D0
     #/3.D0*t11-320.D0/3.D0*s1l1-2048.D0/3.D0-448.D0/3.D0*s2l1)*t13+256*
     #t1+3776.D0/3.D0*t11-1024.D0/3.D0*s1l1-t803-704.D0/3.D0+3968.D0/3.D
     #0*s2l1+(-2560.D0/3.D0*t1-5120.D0/3.D0*t11+t796+t803+640.D0/3.D0-74
     #24.D0/3.D0*s2l1)*t5+64.D0/3.D0*(-10*s2l1-25)*t16+t87
      t930 = 29888.D0/3.D0*t5
      s1 = -64.D0/3.D0*(t9-26*t11+17*s1l1+5*t10+4-22*s2l1)*t3+64.D0/3.D0
     #*(-19*t11+5*s1l1+1-18*s2l1)*t4-t230-t231*(t499+t790)+t793
      t945 = s1+t827*t19+((2432.D0/3.D0+2560.D0/3.D0*s2l1)*t16-t832+(117
     #76.D0/3.D0*t1+4032*t11-7232.D0/3.D0*s1l1-6464.D0/3.D0*t10-1856.D0/
     #3.D0+8640*s2l1)*t5+(-3712*t1-11072.D0/3.D0*t11+t842-4160.D0/3.D0*t
     #10-24832.D0/3.D0-30464.D0/3.D0*s2l1)*t6+(3328.D0/3.D0*t1-512.D0/3.
     #D0*t11+3200.D0/3.D0*s1l1+2560*t10+9152.D0/3.D0-1280*s2l1)*t4+(-896
     #*t1+576*t11+896.D0/3.D0*s1l1-256.D0/3.D0*t10+16448.D0/3.D0+7744.D0
     #/3.D0*s2l1)*t13-128.D0/3.D0*(18*t11+18*t1-t863+3+t864)*t3)*t117+t9
     #05*t153+(64.D0/3.D0*(36*t11-t863+54*t1+3+t864+18*t10)*t5+(-10880.D
     #0/3.D0*t1-6848.D0/3.D0*t11+640*s1l1-t83-1984-13376.D0/3.D0*s2l1)*t
     #13-832.D0/3.D0*t7+(-t69+768*t11-t812-6016.D0/3.D0*t10-5312.D0/3.D0
     #+6016.D0/3.D0*s2l1)*t6+(6400.D0/3.D0*s2l1+3968)*t16)*t21+(t930-t62
     #5-19520.D0/3.D0*t6-10304*t4+14080.D0/3.D0-21568.D0/3.D0*z+1024*t7+
     #960*t13)*t153+(7616.D0/3.D0*t13-24128.D0/3.D0*t7+13376.D0/3.D0*t16
     #+2560.D0/3.D0*t6+192*t5)*t21
      t952 = 320*t5
      t955 = 4928.D0/3.D0*t6
      t958 = 6272.D0/3.D0*t5
      t965 = 4224*t13
      t974 = 13184.D0/3.D0*t5
      t980 = 1088.D0/3.D0*t5
      t981 = 7808.D0/3.D0*t6
      t986 = 10496.D0/3.D0*t5
      t993 = 2048*t4
      t995 = 8576.D0/3.D0*t13
      t998 = 5248.D0/3.D0*z
      t1000 = 1408*z
      t1004 = 1088.D0/3.D0*t16
      t1005 = 4352.D0/3.D0*t13
      t1008 = 704*t5
      t1009 = 3200.D0/3.D0*t4
      t1015 = (-9216*t5+704.D0/3.D0*t16+16000*t4+t965-19264.D0/3.D0+t758
     #-t586)*t98+(13184*t5-t519-32000.D0/3.D0*t4+3776.D0/3.D0*t13+9344.D
     #0/3.D0-t751-16000.D0/3.D0*t6)*t20+(t974+t606-7168*t4+t428+11776.D0
     #/3.D0-8320*z-2880*t6)*t153-t980+(-t981-9344.D0/3.D0*t16+18112.D0/3
     #.D0*t13-t952)*t21+(-t986+3840*t4-1600+t205+3392.D0/3.D0*t6)*l1+t39
     #1+(-21952.D0/3.D0*t6+12352.D0/3.D0*t16-t993+24320.D0/3.D0*t5-t995)
     #*t194+t102-2432.D0/3.D0+t998+(-t185-11584*t5-t776-t325+1408-t1000+
     #18432*t6)*t117+(t1004+t755-t1005)*t44+(-t1008-t1009-3520.D0/3.D0*t
     #13+128-1152*z+2304*t6)*t19
      t1017 = t231**2
      t1023 = 8896.D0/3.D0*t6
      t1026 = 576*t13
      t1036 = 7552.D0/3.D0*t6
      t1040 = 8320.D0/3.D0*t6
      t1072 = 2624.D0/3.D0*t5
      t1077 = -1088+t241+2496*z-t535*t275+(-3968.D0/3.D0*t16-t726+4032*t
     #7+t1023-t459)*t44-t1026+(-38080.D0/3.D0*t5-7040.D0/3.D0+14016*t13-
     #6272*t7-23872.D0/3.D0*t6+7808.D0/3.D0*z+2624*t16+37888.D0/3.D0*t4)
     #*t117+t1036+(-1152+30592.D0/3.D0*z+1024*t16+17152.D0/3.D0*t4+t1040
     #-t316-2816*t13)*l1+((-15808.D0/3.D0*t153+1792*t20-1088.D0/3.D0-129
     #92.D0/3.D0*t19-3264*l1)*s1l1-t244+(-2624.D0/3.D0+6272.D0/3.D0*t20-
     #5440.D0/3.D0*t153-9920.D0/3.D0*l1-t249-2368*t19)*s2l1+6592.D0/3.D0
     #*t20-11840.D0/3.D0*t153-2752.D0/3.D0-4288*t19-12544.D0/3.D0*l1)*t2
     #95+(-t468+t471-128.D0/3.D0*t561*t463*t239)*l5-t495*t496+(-7936*t7+
     #5312*t16+t158-t1072+27712.D0/3.D0*t13-17984.D0/3.D0*t6)*t21
      t1091 = 4608*t4
      t1094 = 4352.D0/3.D0*t5
      t1103 = 4096*t16
      t1119 = 1792*t7
      t1120 = 5504.D0/3.D0*t5
      t1126 = 17792.D0/3.D0*t4
      t1133 = 18688.D0/3.D0*z
      t1135 = 27520.D0/3.D0*t5
      t1144 = 1152*t13
      t1157 = 640*t13
      t1166 = 1024.D0/3.D0*t6
      t1167 = 23552.D0/3.D0*t4
      t1174 = (t1119+t1120-8320.D0/3.D0*t16+t167+t602-20096.D0/3.D0*t6)*
     #t21+(2432.D0/3.D0*t16-t1126+3328*t5+896*t13+1024-t414-1280*t7-t299
     #)*t194-t180+(t727+t519-17152.D0/3.D0-18304.D0/3.D0*t13+t1133+5248*
     #t6-t1135+44032.D0/3.D0*t4)*t117+(-3328.D0/3.D0*t7-t387+6272.D0/3.D
     #0*t16-12416.D0/3.D0*t13+t355)*t44+t198+(-t372+t1144+t371-1664*t6-t
     #604+13568.D0/3.D0-t114+31232.D0/3.D0*z)*t153+(-t661-t139+t369+4864
     #.D0/3.D0*t6+35072.D0/3.D0+17408.D0/3.D0*t5-19584*t4-31360.D0/3.D0*
     #z)*t98+512.D0/3.D0+t501+256.D0/3.D0*t6+(-t706+5248*t4+t1157-15104.
     #D0/3.D0-17536.D0/3.D0*z-2560*t6+t661-3328.D0/3.D0*t5)*t20+(t661+t4
     #18-t377-t417)*t22+(-2048.D0/3.D0+33536.D0/3.D0*z+t372+t1166-t622+t
     #1167-t958)*t19+(t417+9728.D0/3.D0+4352.D0/3.D0*t4-t767+20992.D0/3.
     #D0*z)*l1
      t1177 = 6848.D0/3.D0*t6
      t1180 = 2432.D0/3.D0*t13
      t1187 = 9344*t4
      t1192 = 768*t13
      t1194 = (t8-t299-704*t16+t1192)*t22
      t1195 = 1472.D0/3.D0*t6
      t1202 = (-t411-t459+4544.D0/3.D0*t16-3840*t13+9664.D0/3.D0*t6)*t44
     #+(704*t7+4288.D0/3.D0*t5+3584.D0/3.D0*t16+2624*t13-5952*t6)*t21+(-
     #t661-15808.D0/3.D0*t16+512-t501-t1091+7040.D0/3.D0*t13+26176.D0/3.
     #D0*t6-t1094)*t194-1600.D0/3.D0+t1000+(t640+t682+t369-5184*t6-22144
     #.D0/3.D0-17152.D0/3.D0*t5+12480*t4-4416*z)*t20+(t640+t1103-16256.D
     #0/3.D0+t428+18560.D0/3.D0*z-48320.D0/3.D0*t6+2240*t5+14528*t4)*t11
     #7-t958+2816*t4+(-2752.D0/3.D0*t16-108608.D0/3.D0*t4-13504.D0/3.D0*
     #t13-t305+51392.D0/3.D0-49664.D0/3.D0*z+36160.D0/3.D0*t6+38464.D0/3
     #.D0*t5)*t98+t1174*s2l1+(-1728*t5+t1177+5632.D0/3.D0+11904*z+4032*t
     #4-t1180)*l1+(-3392.D0/3.D0*t16+2240.D0/3.D0*t13+8576.D0/3.D0*t6-64
     #-10048.D0/3.D0*t5+t1187+70144.D0/3.D0*z)*t153+t1194+t1195+(-6464*t
     #5+t1004+61312.D0/3.D0*t4-6400.D0/3.D0*t6-18368.D0/3.D0+24256*z+t59
     #5)*t19
      t1222 = 448.D0/3.D0*t16
      t1223 = 8960.D0/3.D0*t6
      t1236 = 4928.D0/3.D0*t5
      t1239 = 832.D0/3.D0*t16
      t1245 = 5824.D0/3.D0*t6
      t1249 = 512.D0/3.D0*t6
      t1266 = (t305+t641-t357-t299)*t22
      t1270 = 10432.D0/3.D0*t4
      t1273 = 8896.D0/3.D0*t5
      t1276 = (t1222+14336.D0/3.D0+t1144-t1223-t678+6336*t5-t1167-2624*z
     #)*t153+(t669-t1004-11008.D0/3.D0-6848*t13+11776.D0/3.D0*z+30016.D0
     #/3.D0*t6-t930+29440.D0/3.D0*t4)*t117-2368.D0/3.D0*t5-t149+(-1216*t
     #13+t955-t1236+t328+t639+t92+1472.D0/3.D0+t737)*t19+(-t1239+t965+25
     #60*t5-11072.D0/3.D0*t6+2880-t685-14848.D0/3.D0*t4-t998)*t98+(-t699
     #+t1245+4096.D0/3.D0*t16-t387-t644)*t44+t620+(t230-t1094+t173+t1249
     #-t625+384+15680.D0/3.D0*z)*l1+(t661+3136.D0/3.D0*t5-t17-1856.D0/3.
     #D0*t13-768-t163-768*t4-16064.D0/3.D0*z)*t20+(t510+4736.D0/3.D0*t5-
     #8128.D0/3.D0*t16+t158+3520*t13-10816.D0/3.D0*t6)*t21-1216.D0/3.D0+
     #t1266+192*t6+5120.D0/3.D0*z+(-t718+2240*t16+1024-t414-t1270+4096.D
     #0/3.D0*t13-8384.D0/3.D0*t6+t1273)*t194
      t1278 = 18560.D0/3.D0*t6
      t1281 = 10432*t13
      t1285 = 960*t7
      t1288 = 10432.D0/3.D0*t6
      t1296 = 1856*t7
      t1322 = 3136*t13
      t1328 = 1088*t13
      t1336 = 26752.D0/3.D0*t5
      t1341 = 64*t13
      t1346 = (-t1285-t430+12416.D0/3.D0*t16-5888*t13+t1288)*t44+(t666+2
     #752*t5-8192*t16+t572-27520.D0/3.D0*t6)*t21+(-t1296+17984.D0/3.D0*t
     #16+1536-t440-6144*t4-5440*t13+2368.D0/3.D0*t6+5120*t5)*t194+1984.D
     #0/3.D0*t6-t1120+(-t182-24064.D0/3.D0*t13+23744*t4+t1285-7040+52672
     #.D0/3.D0*t6-80192.D0/3.D0*t5+t679)*t117+(-1216.D0/3.D0*t7-t625+104
     #32.D0/3.D0*t13-22144.D0/3.D0*t6+8896+55936.D0/3.D0*t5-22400*t4-704
     #0*z)*t98-128+7232.D0/3.D0*z+7808.D0/3.D0*t4+(t832-3712.D0/3.D0*t16
     #+11392.D0/3.D0+t1322-t1278+t676-8896*t4-57088.D0/3.D0*z)*t20+t1194
     #-t773+(-t8+30592.D0/3.D0*t5+t1004+t1328+5120-4736*t6-20672.D0/3.D0
     #*t4+21184.D0/3.D0*z)*t153+(t17-1664*t13+4096*t6-6656.D0/3.D0+t328-
     #t1336+13248*t4+39296.D0/3.D0*z)*t19+(-t1341-1600.D0/3.D0+8704*t4+t
     #370-t103-t986+44416.D0/3.D0*z)*l1
      t1348 = 17536.D0/3.D0*t6
      t1364 = 320*t16
      t1371 = 1472.D0/3.D0*t5
      t1377 = 1280*t4
      t1387 = 3776.D0/3.D0*t4
      t1388 = 1472.D0/3.D0*t13
      t1395 = 128*t6
      t1401 = 2432.D0/3.D0*t6
      t1409 = (-t727+t519-6016.D0/3.D0*t13+t1245)*t44+(1408.D0/3.D0*t7+t
     #459-t1364+t158+2752.D0/3.D0*t13-3520.D0/3.D0*t6)*t21+(-576*t16-t93
     #+12800.D0/3.D0*t13+t1371-t1040-t678)*t194+t1266+(1024.D0/3.D0*t16+
     #22592.D0/3.D0*t5-11648.D0/3.D0*t13-2048+t402+t1377-t1036+t328)*t11
     #7-t767-768+2752.D0/3.D0*z+(t372+t112-17344.D0/3.D0*t5+5824*t6+2681
     #6.D0/3.D0-t604-9216*t4-8832*z)*t98+t1387+(t625-t1388-24832.D0/3.D0
     #-2752.D0/3.D0*t5-t645+t328+27328.D0/3.D0*t4+5440*z)*t20+(1024*t13-
     #t1395+1664+13312.D0/3.D0*t4-2432*t5-t92+29696.D0/3.D0*z)*t153+(-t1
     #401-4864.D0/3.D0+30464.D0/3.D0*z+1024.D0/3.D0*t5-t156+t1126)*t19+(
     #-4096.D0/3.D0*t5+6208.D0/3.D0+t405+t189+t374)*l1
      t1421 = t1202*s1l1+(39040.D0/3.D0*t6+34112.D0/3.D0*t5+9728.D0/3.D0
     #*t16+14912.D0/3.D0-55616.D0/3.D0*t13-50048.D0/3.D0*t4+2880*t7-3328
     #*z)*t98+12928.D0/3.D0*t4+(-3520+37888.D0/3.D0*z-512*t7-7104*t5+381
     #44.D0/3.D0*t4+t395-896*t6-2240.D0/3.D0*t16)*t19-t597+t1276*t11+(-5
     #12-t1278+t1091-6720*t5-t1119-5312.D0/3.D0*t16+t1281-23104.D0/3.D0*
     #z)*t20+t1346*s2l1+(t1348+512-t501+9024*t7-t1281-15040.D0/3.D0*t4-2
     #1760.D0/3.D0*t16+21952.D0/3.D0*t5)*t194+(1216.D0/3.D0*t153-t787+89
     #6*l1+2816.D0/3.D0*t19+896.D0/3.D0)*t294+t1409*t10-t556*t557+(9728.
     #D0/3.D0-2944.D0/3.D0*t16-t1348+40768.D0/3.D0*t5-22016.D0/3.D0*t4+t
     #1133+1472*t7-t428)*t153+(t345-t299-t365-t372)*t22
      t1433 = -256.D0/3.D0*t19-64.D0/3.D0-64.D0/3.D0*l1+64*t153
      t1437 = 384*z
      t1439 = 2112*t6
      t1445 = 512.D0/3.D0*t153
      t1446 = 128*t20
      t1447 = 320.D0/3.D0*l1
      t1448 = 704.D0/3.D0*t19
      t1452 = 1088.D0/3.D0*t19
      t1454 = 320*t19
      t1456 = 64*l1
      t1467 = 64*t4
      t1470 = 64*t145
      t1476 = 64*t6
      t1485 = 64.D0/3.D0*t6
      t1498 = (-832*t5-t103+256-t751+t128-t1195)*t117+(320*t13+2176.D0/3
     #.D0*t6-t1467+t767)*t194-t1470+(256.D0/3.D0*t5-128*t4+t169-64.D0/3.
     #D0-448*z)*t153-128.D0/3.D0-t138+(t1222+t1476+192*t13)*t44+(t123-70
     #4.D0/3.D0*t4+64.D0/3.D0-t186)*t19+(t1371+t103+640*t4+t149-2176.D0/
     #3.D0+1984.D0/3.D0*z+t1485)*t98+(-t1388-t172-t1166-t110)*t21+(-t146
     #7-64.D0/3.D0-640.D0/3.D0*z)*l1+(-64.D0/3.D0*t5-t1341-1472.D0/3.D0*
     #t4-t1485+1600.D0/3.D0-64*z)*t20
      t1500 = -3328.D0/3.D0+(-5888.D0/3.D0*t6-13312.D0/3.D0*t13-8192.D0/
     #3.D0*t16+960*t4+25472.D0/3.D0*t7-t952)*t194+(-t955-4544.D0/3.D0*t4
     #-1792.D0/3.D0*z-t685-t350+t958+5440.D0/3.D0*t13)*t19+t1015*s2l1+t1
     #017*(t1077+t1421)/2+(-3776*t16-448*t6+12800.D0/3.D0*t7-t149)*t44+t
     #1433*t294+(2624.D0/3.D0*t4+t350-6080.D0/3.D0*t5-64.D0/3.D0-t1437-5
     #120.D0/3.D0*t13+t1439)*l1+(-t641+1216*t16-t1285)*t22+(t1445+(-t144
     #6+t1447+t1448+64+128.D0/3.D0*t153)*s1l1+t1452+t1447+(-320.D0/3.D0-
     #512*t20+t1454+1280.D0/3.D0*t153+t1456)*s2l1-320*t20)*t295+t1498*t1
     #0
      t1523 = (t811+2368.D0/3.D0*t11-t817+3904.D0/3.D0+t36)*t5+(-t870-10
     #88*t11+t842+t883-448-5440.D0/3.D0*s2l1)*t4+t46-256.D0/3.D0*t11+t81
     #2+1216.D0/3.D0*t10+1472.D0/3.D0+64*s2l1+(-128.D0/3.D0*t1+640*t11-3
     #008.D0/3.D0*s1l1-2240.D0/3.D0*t10-1856.D0/3.D0+832*s2l1)*t3+t773+6
     #4.D0/3.D0*(-35-10*t11-20*s2l1)*t6
      t1612 = (-9088.D0/3.D0*t1-14144.D0/3.D0*t11+1216.D0/3.D0*s1l1-3520
     #.D0/3.D0*t10-2112-6848*s2l1)*t4+(1600.D0/3.D0-640.D0/3.D0*s2l1)*t1
     #6+(1536*t1+2816*t11-4160.D0/3.D0*s1l1-1280*t10-3712.D0/3.D0+11584.
     #D0/3.D0*s2l1)*t3-t87+(-t811+7616.D0/3.D0*t11+3904.D0/3.D0*s1l1+569
     #6.D0/3.D0*t10+1792+6848.D0/3.D0*s2l1)*t5-t807+t64+t2+832-576*s2l1+
     #576*t10+(t80+320.D0/3.D0*t11-704*s1l1+448.D0/3.D0*t10+12928.D0/3.D
     #0+3264*s2l1)*t6+(t800-128*t11-64*s1l1+192*t10-12032.D0/3.D0-5312.D
     #0/3.D0*s2l1)*t13
      t1615 = -t7+t16
      t1617 = 64*t1615*t22
      t1627 = -t1615
      t1630 = 1088*z
      t1632 = 5248.D0/3.D0*t13
      t1651 = 2432.D0/3.D0*t4
      t1652 = 896*z
      t1656 = t1617-t725+(t604-9856.D0/3.D0-t611+t139+t354-832*t13+4288*
     #z+t768)*t153+(-t727-4352*t5-t182+1920*t4+14080.D0/3.D0*t13-t793)*t
     #117+1024.D0/3.D0*t1627*t44+2176.D0/3.D0-t1630+t1476+(t678+28352.D0
     #/3.D0*t5-t1632+t182-11840.D0/3.D0*t4-12416.D0/3.D0*t6-256+t751)*t9
     #8+(-t763+1792.D0/3.D0*t16+t1144-t403)*t21-t1009+(t139-t770-t1401-t
     #149-t113+t1236)*l1+(-t305-t986-t110+2688-t993+t652-3712*z+8512.D0/
     #3.D0*t6)*t20+(t777+t172-t826-4032*t13+10048.D0/3.D0*t6)*t194+896.D
     #0/3.D0*t5+(-t585+128-t92+t1651+t428+t1652+t1249-3584.D0/3.D0*t5)*t
     #19
      t1664 = 896*t16
      t1672 = 2816.D0/3.D0*t6
      t1701 = 128*t1615*t22+t756+(t328-2816*t5-t1632+t190+t1187+t636-189
     #44.D0/3.D0+9472*z)*t153+(-t333-t1336-t1664+4096*t4+25600.D0/3.D0*t
     #13-t981)*t117+(t386-t162-t140)*t44-2176*t4+(-t110+t1672+t1192+768-
     #t751+t114-8576.D0/3.D0*t5)*t19+4352.D0/3.D0-2176*z+(-t365+t826+t99
     #5-2176*t6)*t21+(t149+8704.D0/3.D0*t5-1536*t6-4480.D0/3.D0*t4-896.D
     #0/3.D0-1408.D0/3.D0*z)*l1+(-t305-t1135+13568.D0/3.D0-t139-3584.D0/
     #3.D0*t4+t1180-20096.D0/3.D0*z+15616.D0/3.D0*t6)*t20+(t411+t459+204
     #80.D0/3.D0*t6+t350-25216.D0/3.D0*t13)*t194+(t305+59648.D0/3.D0*t5+
     #t357-28288.D0/3.D0*t4-8192.D0/3.D0*t13-128+t135-23936.D0/3.D0*t6)*
     #t98
      t1711 = 448.D0/3.D0*t7
      t1712 = 2176*t5
      t1718 = 5504.D0/3.D0*t16
      t1758 = t1701*s2l1+8512.D0/3.D0-3904*z+(t585-t379-4736.D0/3.D0*t16
     #+28480.D0/3.D0-9984*t4-2304*t13-38272.D0/3.D0*z+43904.D0/3.D0*t6)*
     #t20+(-t1711+t1712+33152.D0/3.D0*t16-57152.D0/3.D0*t13+17920.D0/3.D
     #0*t6)*t194+t1617+(-t585-35776.D0/3.D0+t1718+21312*t4-t1322+17088*z
     #+2624*t6-10560*t5)*t153+9536.D0/3.D0*t5+(t1711+t449-t1144)*t44-531
     #2*t4+(t585-36032.D0/3.D0*t4+115648.D0/3.D0*t5+8512.D0/3.D0*t16-896
     #+t1652+14656.D0/3.D0*t13-100160.D0/3.D0*t6)*t98+(-t585-5888*t16+27
     #136.D0/3.D0*t13-3136*t6)*t21+(t511-1984.D0/3.D0-320.D0/3.D0*z+t974
     #-4480*t6-t198)*l1+(-t119-18304.D0/3.D0*t5+16832.D0/3.D0*t6-t683+11
     #52-1216*z+t1651)*t19+(t640-64064.D0/3.D0*t5-25984.D0/3.D0*t16+6208
     #*t4+10880*t13+12864*t6)*t117-704*t6
      t1807 = t1617+1408-2112*z+(-t585-24448.D0/3.D0*t5-2496*t6-t1364+39
     #68*t4+21056.D0/3.D0*t13)*t117+(t333-t110-t140)*t44+(t428+256-576*z
     #-t1273+t1672+3520.D0/3.D0*t4)*t19+t1008+(-t646-t1364+7808.D0/3.D0*
     #t13-2048*t6)*t21-2112*t4+(-t585-8640*t5+2624.D0/3.D0*t13-t92+t1467
     #+12608.D0/3.D0*t6+10880.D0/3.D0-t707)*t20+(t604+t172-21824.D0/3.D0
     #*t13+t350+19136.D0/3.D0*t6)*t194+(t110+7680*t4+1280*t6-t1005-5248+
     #8512*z-6656.D0/3.D0*t5)*t153+(t328+55360.D0/3.D0*t5-2240*t13-28480
     #.D0/3.D0*t4-20672.D0/3.D0*t6+128-t135)*t98+(-t793-512.D0/3.D0-t137
     #7+2688*t5-t138)*l1
      t1816 = 128*t153
      t1820 = 320*l1
      t1831 = 5312.D0/3.D0+t1656*t11+t1758*s1l1-t402+11008.D0/3.D0*t5+(-
     #6400.D0/3.D0*t7+9664.D0/3.D0*t16-t1328)*t44+t1807*t10+(-3008.D0/3.
     #D0*t7+t520-14016*t5-23552.D0/3.D0+10240*z+t1718-t345+52288.D0/3.D0
     #*t4)*t153+((-64*t19-t1816+192*l1)*s1l1-t1454+(t1820-448*t19+t1816)
     #*s2l1+t1820)*t295+(-30208.D0/3.D0*t16+17024.D0/3.D0*t13+15488.D0/3
     #.D0*t7-t189)*t21+t230
      t1917 = t1617+5440.D0/3.D0-7616.D0/3.D0*z+(t604-26816.D0/3.D0*t5+2
     #368.D0/3.D0*t16-25792.D0/3.D0+47552.D0/3.D0*t4-6208.D0/3.D0*t13+34
     #112.D0/3.D0*z+2752*t6)*t153-t176+(t685-t826-t140)*t44+5248.D0/3.D0
     #*t5+(t678+76672.D0/3.D0*t5+t710-8960*t4+t1157-384+t1437-57472.D0/3
     #.D0*t6)*t98+(11200.D0/3.D0*t13-t1439-t1664-t763)*t21-9728.D0/3.D0*
     #t4+(-t305-24704.D0/3.D0*t5-1664.D0/3.D0*t16-19328.D0/3.D0*t4-6656.
     #D0/3.D0*t13+19456.D0/3.D0-26176.D0/3.D0*z+33280.D0/3.D0*t6)*t20+(t
     #777+t764+3904*t16-11584*t13+6528*t6)*t194-t725+(-t585-8128.D0/3.D0
     #*t5+704-t182+704*t4+t683+t125+5632.D0/3.D0*t6)*t19+(4160*t4-4800*t
     #16-34880.D0/3.D0*t5+31232.D0/3.D0*t13+t1177-t727)*t117+(t139-6208.
     #D0/3.D0*t4-t1023+t179-64.D0/3.D0-3776.D0/3.D0*z+11456.D0/3.D0*t5)*
     #l1
      t1928 = (t685-t155+576-t135+t1239+16448.D0/3.D0*t6+256.D0/3.D0*t4-
     #3392*t13)*t19+(-t350+t128-1088.D0/3.D0+t740+4672.D0/3.D0*t5-t635-9
     #728.D0/3.D0*t6)*l1-1728*t6+(-20096.D0/3.D0*t7-7936.D0/3.D0*t6+4230
     #4.D0/3.D0*t16-20800.D0/3.D0*t13+t1712)*t194+384*t1627*t22+(8896*t6
     #-40000.D0/3.D0*t4+20864.D0/3.D0+3392.D0/3.D0*t7-7360*t13-25856.D0/
     #3.D0*z+10688.D0/3.D0*t5+t1222)*t20-12160.D0/3.D0*t4+(-43072.D0/3.D
     #0*t5-7040*t16+4672*t7+2752*t4-22336.D0/3.D0*t13+64256.D0/3.D0*t6)*
     #t117+t1917*s2l1+(-64*l1+64*t19)*t294+(-96896.D0/3.D0*t6-t1296+5708
     #8.D0/3.D0*t13+61568.D0/3.D0*t5-t413-1088+t1630-8704.D0/3.D0*t4)*t9
     #8
      t1989 = -t1470+(-31232.D0/3.D0*t5+t625+30848.D0/3.D0*t4+1472*t13-9
     #856.D0/3.D0+11200.D0/3.D0*z+1920*t6)*t98+(t1239+t120-t1026)*t44+t1
     #387+(-2048*t5+t639-2240.D0/3.D0+832*z+2240.D0/3.D0*t6)*l1+(20800.D
     #0/3.D0*t5-t625-t1270+t725+64+6464.D0/3.D0*z-t1223)*t20+(t595-t387-
     #t1395-448*t16)*t21+(6272.D0/3.D0*t13-4928*t6+t1364-t770+11648.D0/3
     #.D0*t5)*t194-t980+(t769+t103-17728.D0/3.D0*t4+t428+10240.D0/3.D0-6
     #720*z-4160.D0/3.D0*t6)*t153-768+(-t92-3392.D0/3.D0*t5-11584.D0/3.D
     #0*t13+1024-t414-5824.D0/3.D0*t4+17792.D0/3.D0*t6)*t117+t614+(t1072
     #-3392.D0/3.D0*t4-1216.D0/3.D0*t13+512.D0/3.D0-t734+t299)*t19
      t1991 = t1523*l1+((-448.D0/3.D0-t1447-t1446-t1448+896.D0/3.D0*t153
     #)*s1l1+(64.D0/3.D0-t1456+t787-t1454-256.D0/3.D0*t153)*s2l1-256.D0/
     #3.D0-t1447+64*t20-t1452+t1445)*t295+(6720*t4-13632*t13-384-192*z-3
     #7504.D0/3.D0*t5+t435+2432*t7+t1103)*t98+8128.D0/3.D0*t4-t611+((-14
     #848.D0/3.D0-7040.D0/3.D0*s2l1)*t16+(7936.D0/3.D0*t1+3008.D0/3.D0*t
     #11-2240.D0/3.D0*s1l1+1088.D0/3.D0*t10+1984.D0/3.D0+8000.D0/3.D0*s2
     #l1)*t13-64.D0/3.D0*(-51*t11+24*s1l1-t26+69-72*s2l1+45*t10)*t4+(896
     #0.D0/3.D0*t1+6464.D0/3.D0*t11+t48+8704.D0/3.D0*t10+12224.D0/3.D0+1
     #0688.D0/3.D0*s2l1)*t6+2624.D0/3.D0*t7+(-3968*t1-4096*t11+3328.D0/3
     #.D0*s1l1+128*t10+2432.D0/3.D0-16256.D0/3.D0*s2l1)*t5)*t194+t1612*t
     #98-t1017*t231*(t1831+t1928)/6-t1433*t294+(17728.D0/3.D0*t5+384-t14
     #37+8384*t13-16448.D0/3.D0*t7-10688.D0/3.D0*t6-4864*t4-t762)*t117+(
     #15296.D0/3.D0*t4-t1288+6720*t13-t430-2432*t16-10304.D0/3.D0+18944.
     #D0/3.D0*z-4160.D0/3.D0*t7)*t20+t1989*t11
      t1995 = (l1+z-t236)**2
      taux17 = (t216+t945+t1500+t1991)/t1995*t239
      t4 = 1-z
      t5 = t4**2
      t6 = t5**2
      t7 = t6*t5
      t10 = t5*t4
      t11 = t6*t10
      t12 = 1280*t11
      t32 = t6*t4
      t37 = l1**2
      t38 = t37*l1
      t41 = s1l1+1+s2l1
      t63 = t37**2
      t64 = t63**2
      t66 = dlog(z)
      t67 = t66**2
      t85 = t63*t37
      t102 = 256*t11
      t105 = 16384.D0/3.D0*t7
      t106 = 41984.D0/3.D0*t32
      t114 = 1280*z
      t118 = 1280.D0/3.D0*t11
      t126 = 1792*t11
      t132 = -t11+t7
      t134 = t64*t37
      t135 = 256*t132*t134
      t136 = 4096.D0/3.D0*t6
      t141 = t63*t38
      t143 = 256.D0/3.D0*t11
      t146 = 6656*t7
      t156 = 27904.D0/3.D0*t5
      t160 = t63*l1
      t167 = t64*l1
      t170 = 7936.D0/3.D0*t7
      t173 = 39424.D0/3.D0*t10
      t175 = 17664*t6
      t179 = (-t102+8704.D0/3.D0*t10+115712.D0/3.D0*t5-21248+t105+t106+2
     #5600*z-39424*t6)*t63+(-16640.D0/3.D0*t32-3584.D0/3.D0*t5+37376.D0/
     #3.D0*t6+1280-t114-20992.D0/3.D0*t10)*l1+(-t118-173312.D0/3.D0*t6+9
     #6256.D0/3.D0*t7+48896*t10-8448*t5-14336*t32)*t85+(-t126+3840*t6+20
     #992*t7-23040*t32)*t64+t135+(t12+t136-116224.D0/3.D0*t7-7168*t10+43
     #264*t32)*t141+(t143+101632.D0/3.D0*t10+67328.D0/3.D0-21760.D0/3.D0
     #*t6-t146+6400*t32-48896*t5-28416*z)*t38-16640.D0/3.D0+6912*z+(t102
     #+101632*t6-34048.D0/3.D0*t7-224512.D0/3.D0*t10+t156+3584-3584*z-85
     #760.D0/3.D0*t32)*t160+12544*t5+(3328.D0/3.D0*t11-14080.D0/3.D0*t7+
     #3584*t32)*t167+2816*t6+(t170-512+12800.D0/3.D0*t32+768*z+t173-5632
     #.D0/3.D0*t5-t175)*t37-29440.D0/3.D0*t10
      t197 = 1536*t7
      t198 = 1024*z
      t200 = 60928.D0/3.D0*t6
      t214 = 512.D0/3.D0*t10
      t221 = 7936.D0/3.D0*t10
      t222 = 3072*t32
      t237 = 33536.D0/3.D0*t32
      t240 = 1024.D0/3.D0*t11
      t254 = 5120.D0/3.D0*t11
      t266 = 6656.D0/3.D0*t32
      t269 = 2048.D0/3.D0*t11
      t278 = (-512*t11-14336.D0/3.D0+3584.D0/3.D0*t7-t214+38912.D0/3.D0*
     #t5+15104.D0/3.D0*z+34816.D0/3.D0*t32-t200)*t63+(512.D0/3.D0*t7-t22
     #1-t222+256.D0/3.D0*t5-256.D0/3.D0+256.D0/3.D0*z+16640.D0/3.D0*t6)*
     #l1+13312.D0/3.D0*t5+(-8704.D0/3.D0*t11-256*t6+t146-10496.D0/3.D0*t
     #32)*t64+t135+(3072*t11+3584.D0/3.D0*t6-46336.D0/3.D0*t7+t237)*t141
     #+(t240+21248.D0/3.D0-58624.D0/3.D0*t5-5376*t6-t170+55040.D0/3.D0*t
     #10+1792*t32-22784.D0/3.D0*z)*t38-256.D0/3.D0*t32-4352.D0/3.D0+4352
     #.D0/3.D0*z+5120.D0/3.D0*t6-4096.D0/3.D0*t132*t167+(-t254-56576.D0/
     #3.D0*t6+44288.D0/3.D0*t7+26368.D0/3.D0*t10-256*t5-8192.D0/3.D0*t32
     #)*t85+(-t143-768+4864.D0/3.D0*t5-16384.D0/3.D0*t6+1024*t7+4352.D0/
     #3.D0*t10+t266+t198)*t37+(t269+125440.D0/3.D0*t6-4608*t7-63488.D0/3
     #.D0*t10+2048.D0/3.D0*t5-17408*t32)*t160-4608*t10
      t283 = 29696.D0/3.D0*t6
      t285 = 8704.D0/3.D0*t7
      t292 = (-5120*t11-23552*t6-17920-97792.D0/3.D0*t10+146944.D0/3.D0*
     #t5-4096.D0/3.D0*t7+94720.D0/3.D0*t32+20224*z)*t63+(-20480*t11-8960
     #0*t6+75008.D0/3.D0*t7+37376*t10-2560*t5+150784.D0/3.D0*t32)*t85+(4
     #352*t11+53248.D0/3.D0+146432.D0/3.D0*t10-144896.D0/3.D0*t5-17408*t
     #6-20480.D0/3.D0*t7+4864.D0/3.D0*t32-20992*z)*t38+(-73472.D0/3.D0*t
     #11-1792*t6+124160.D0/3.D0*t7-15104*t32)*t64-12544.D0/3.D0+t179*s1l
     #1+(10752*t11-14592*t7+3840*t32)*t167+(24832.D0/3.D0*t11+332032.D0/
     #3.D0*t6+22528.D0/3.D0*t7-119552.D0/3.D0*t10-21248.D0/3.D0*t5+3840-
     #3840*z-250112.D0/3.D0*t32)*t160+(-4352.D0/3.D0*t11+46592.D0/3.D0*t
     #32-t197-768+t198+2816.D0/3.D0*t5-t200+22784.D0/3.D0*t10)*t37+(3020
     #8*t11+74752.D0/3.D0*t6-54528*t7-7680*t10+21248.D0/3.D0*t32)*t141+7
     #168*t6+t278*s2l1-40448.D0/3.D0*t10+2048*t132*t134+(t283-t214-9984*
     #t32+t285-11776.D0/3.D0*t5+1280-t114)*l1-4352.D0/3.D0*t32+4864*z+35
     #840.D0/3.D0*t5
      t298 = 1024*t10
      t308 = 7*s2l1
      t334 = 5120.D0/3.D0*z
      t367 = t11*t64*t38
      t368 = 1408*t367
      t378 = 2048*t5
      t380 = 75520.D0/3.D0*t10
      t383 = 29440.D0/3.D0*t6
      t402 = 2048*z
      t411 = 14848*t32
      t412 = 1024*t6
      t413 = 2560*t7
      t416 = (-42240*t7-2048.D0/3.D0+275840.D0/3.D0*t10-46592*t6+t334+79
     #36.D0/3.D0*t11-20224*t5-202880.D0/3.D0*t32)*t85-78080.D0/3.D0*t5+(
     #8704*t32-19968*t6+113536.D0/3.D0*t10+22144-54528*t5-141184.D0/3.D0
     #*z)*l1+(19968*t7+52864*t32+20992.D0/3.D0-79232.D0/3.D0*t10-15104.D
     #0/3.D0*t6-28160.D0/3.D0*z-t156-4480.D0/3.D0*t11)*t160+(-135680.D0/
     #3.D0*t32-36608.D0/3.D0*t7-210688.D0/3.D0*t10+151552.D0/3.D0*t6-446
     #72+256384.D0/3.D0*z+1024*t11+311296.D0/3.D0*t5)*t63+48128.D0/3.D0*
     #t10-t368+(30208.D0/3.D0*t7+44032.D0/3.D0*t32+15104*t6-125824.D0/3.
     #D0*t10+141056.D0/3.D0*t5-17280+131968.D0/3.D0*z-t240)*t38+35072.D0
     #/3.D0-24832*z+(-t378-50432*t7+t380-95872*t32+27136.D0/3.D0*t11-t38
     #3)*t64+(-19456.D0/3.D0*t6+60416.D0/3.D0*t7-9600*t11-t298+176896.D0
     #/3.D0*t32)*t167+(-12800*t6-3840*t7-7168*t32+111616.D0/3.D0*t10-142
     #592.D0/3.D0*t5+59264.D0/3.D0-143488.D0/3.D0*z)*t37+(-16768.D0/3.D0
     #*t11+286592.D0/3.D0*t32+2048-t402+115712.D0/3.D0*t6-68608*t10+2662
     #4.D0/3.D0*t5+183296.D0/3.D0*t7)*t141-13568.D0/3.D0*t6+(17152.D0/3.
     #D0*t11-t411+t412-t413)*t134
      t445 = -1+l1
      t453 = 2*l6-1
      t454 = t453**2
      t457 = l5-1
      t459 = z*l1
      t461 = (1-l1+t459)**2
      t462 = 1/t461
      t467 = dsqrt(z)
      t469 = dsqrt(l5)
      t471 = dsqrt(-t457)
      t474 = dsqrt(l1)
      t478 = l5*z*l1-l1+1-l5+l5*l1
      t481 = dsqrt(-t445)
      t486 = t478**2
      t489 = dsqrt(4*l1*t454*z*l5*t457*t445*t462-4*t453*t467*t469*t471*t
     #474*t478*t462*t481+t486*t462)
      t504 = 67328.D0/3.D0*t6
      t511 = 384*t367
      t513 = 3072*t6
      t515 = 12288*t32
      t530 = 2048*t10
      t547 = 4096*z
      t563 = 6144*t5
      t566 = 46336.D0/3.D0*t10
      t583 = 24064.D0/3.D0*z
      s1 = -5888.D0/3.D0+(-t106+t197-74368.D0/3.D0*t10+t504-51712.D0/3.D
     #0+113152.D0/3.D0*t5+107264.D0/3.D0*z)*l1-512*t32-48512.D0/3.D0*t10
     #+t511+20992.D0/3.D0*t6+(-t513-3712*t11+t515+3456*t7)*t134+(-4224*t
     #11-768*t7+227584.D0/3.D0*t32-296192.D0/3.D0*t6+179200.D0/3.D0+1809
     #92.D0/3.D0*t10-301312.D0/3.D0*t5-322688.D0/3.D0*z)*t63
      t588 = s1+(40832.D0/3.D0*t11+27136*t6-t530-27904.D0/3.D0*t7-47488*
     #t32)*t167+13824*t5+(t237-28672.D0/3.D0+5632.D0/3.D0*t7-43648*t10+2
     #0224.D0/3.D0*t6-t269+118400.D0/3.D0*t5+69760.D0/3.D0*z)*t37+20224.
     #D0/3.D0*z+(25216*t11-268160.D0/3.D0*t32-52224*t5+4096-t547+27392*t
     #6-37632*t7+315776.D0/3.D0*t10)*t141+(-20224.D0/3.D0*t7-10752*t32-7
     #0400.D0/3.D0*t6+13184+81536*t10-193024.D0/3.D0*t5+8192.D0/3.D0*t11
     #-73088.D0/3.D0*z)*t38+(-75392.D0/3.D0*t11+t563+275072.D0/3.D0*t32+
     #53248.D0/3.D0*t7-t566-62080*t6)*t64+(-28928.D0/3.D0*t7-111488.D0/3
     #.D0+197120.D0/3.D0*t6-271360.D0/3.D0*t32+17536.D0/3.D0*t11+186752.
     #D0/3.D0*z+143104.D0/3.D0*t10+52480.D0/3.D0*t5)*t160+(118400.D0/3.D
     #0*t7-33280.D0/3.D0+92800.D0/3.D0*t6+61568*t32-42112.D0/3.D0*t11+t5
     #83-578176.D0/3.D0*t10+307072.D0/3.D0*t5)*t85
      t692 = (10368*t11+4096-t547-382336.D0/3.D0*t7-64640.D0/3.D0*t32+18
     #2272.D0/3.D0*t6-70656*t5+342656.D0/3.D0*t10)*t141+143744.D0/3.D0*t
     #5+(-2048*t11+199040.D0/3.D0*t32+66688.D0/3.D0*t7+26240+268160.D0/3
     #.D0*t10-145536*t6-239872.D0/3.D0*z-186880.D0/3.D0*t5)*t63-105472.D
     #0/3.D0*t10+(-19456*t32-141184.D0/3.D0*t10+131072.D0/3.D0*t6-22144+
     #151936.D0/3.D0*t5+146176.D0/3.D0*z)*l1-61952.D0/3.D0+98432.D0/3.D0
     #*z+(-108416.D0/3.D0-115328.D0/3.D0*t7-99840*t32+2176*t11-186880.D0
     #/3.D0*t10+189312*t6+211328.D0/3.D0*z+84736.D0/3.D0*t5)*t160+t511+(
     #288256.D0/3.D0*t7-33280.D0/3.D0-212608.D0/3.D0*t6+130688.D0/3.D0*t
     #32-11008.D0/3.D0*t11+t583-147200*t10+380800.D0/3.D0*t5)*t85+(8704*
     #t7-40960.D0/3.D0+95872.D0/3.D0*z+48128.D0/3.D0*t32-23040*t6-73216.
     #D0/3.D0*t10+136448.D0/3.D0*t5)*t37+(-t513-3328*t11+t515+2688*t7)*t
     #134+(-22528*t7+73344-t175-t266-172416*t5+390016.D0/3.D0*t10+t269-1
     #07904*z)*t38+30208.D0/3.D0*t6+(-14848*t11+t563+146944.D0/3.D0*t32+
     #258560.D0/3.D0*t7-t566-80000*t6)*t64+(-t530-82432.D0/3.D0*t7+30848
     #.D0/3.D0*t11+36352*t6-44160*t32)*t167
      s1 = -48896.D0/3.D0+(123904.D0/3.D0*t7-14848.D0/3.D0-t136+449408.D
     #0/3.D0*t32-71296*t11+5632.D0/3.D0*z-684160.D0/3.D0*t10+114560*t5)*
     #t85+(-18176*t7+157696.D0/3.D0+440320.D0/3.D0*t32-22144*t11+51968*t
     #10-142592*t6-108928*z-238720.D0/3.D0*t5)*t63+(-167168.D0/3.D0+2112
     #0*t7-639616.D0/3.D0*t32+104192.D0/3.D0*t11-7040.D0/3.D0*t10+173184
     #*t6+259328.D0/3.D0*z+44288*t5)*t160+25856*z+128.D0/3.D0*t445*(310*
     #t37-48*t63+184*l1+51+214*t38)*t489+(60800-12544*t32-33280.D0/3.D0*
     #t7+584960.D0/3.D0*t10-227840.D0/3.D0*t6+16384*t11-530048.D0/3.D0*t
     #5-257152.D0/3.D0*z)*t38-148864.D0/3.D0*t10-15872.D0/3.D0*t32
      t711 = s1+t588*s2l1+(107264.D0/3.D0*t32-21248.D0/3.D0*t6+51968*t5-
     #t105-42112.D0/3.D0-16384.D0/3.D0*t11+98176.D0/3.D0*z-50304*t10)*t3
     #7+26112*t6+(-t513-13952*t11+18432*t32-1408*t7)*t134+(-t530-17664*t
     #7+59904*t11+32768*t6-72960*t32)*t167+(-323968.D0/3.D0*t11+t563+381
     #056.D0/3.D0*t32+184192.D0/3.D0*t7-83200.D0/3.D0*t10-176512.D0/3.D0
     #*t6)*t64+t692*s1l1+(109824*t11+4096-t547-80768*t7-410240.D0/3.D0*t
     #32+9728*t6-63488*t5+472064.D0/3.D0*t10)*t141+(47104*t6+10752*t7-79
     #360.D0/3.D0+51840*z+166528.D0/3.D0*t5-111104.D0/3.D0*t32-132608.D0
     #/3.D0*t10)*l1+141824.D0/3.D0*t5
      t777 = -t41
      t783 = 53632.D0/3.D0*t7
      s2 = 69376.D0/3.D0+(128.D0/3.D0*(67*s1l1+36*s2l1+109)*t7-t12+(-989
     #44.D0/3.D0*s1l1-1024.D0/3.D0+74240.D0/3.D0*s2l1)*t4+(22272*s1l1+20
     #48-5120*s2l1)*t6+(-103168.D0/3.D0*s1l1-27264-7552*s2l1)*t10+(51968
     #*s1l1+101888.D0/3.D0-2432*s2l1)*t5+(-19840.D0/3.D0*s1l1-1408+23168
     #.D0/3.D0*s2l1)*t32-2816*s1l1-16896-24448*s2l1)*t38+7424.D0/3.D0*t3
     #2+(2048*t41*t5+(8960.D0/3.D0*s1l1+108416.D0/3.D0+118528.D0/3.D0*s2
     #l1)*t6+256.D0/3.D0*(-367*s1l1-367*s2l1-403)*t10+(22528.D0/3.D0*s1l
     #1+164864.D0/3.D0-20864.D0/3.D0*s2l1)*t32+(-39552*s1l1-252800.D0/3.
     #D0-51456*s2l1)*t7+76544.D0/3.D0*t11)*t64
      s1 = s2-t67*t66*t292/6+(-14848.D0/3.D0*t6-12416*t7-35584*t11-t298+
     #161920.D0/3.D0*t32)*t167+((68992.D0/3.D0*s1l1+121856.D0/3.D0+13952
     #.D0/3.D0*s2l1)*t6+88448.D0/3.D0+84224.D0/3.D0*s1l1+70016.D0/3.D0*s
     #2l1-128.D0/3.D0*(t308-53)*t7-128.D0/3.D0*(127*s1l1+87*s2l1+383)*t3
     #2+(-140416.D0/3.D0*s1l1-150400.D0/3.D0-96640.D0/3.D0*s2l1)*t4+(613
     #12*s1l1+167296.D0/3.D0+65792.D0/3.D0*s2l1)*t5+(-160768.D0/3.D0*s1l
     #1-162304.D0/3.D0-6144*s2l1)*t10)*t37+t416*s1l1-t66*t711
      s2 = s1+(78208.D0/3.D0-151424.D0/3.D0*z-93824.D0/3.D0*t6+50944.D0/
     #3.D0*t32-62848*t5+150400.D0/3.D0*t10-4864*t7)*l1+(5504*t32-2432*t7
     #-223616.D0/3.D0*t10+123392.D0/3.D0*t6-7168*t11+204800.D0/3.D0*t5-7
     #2704.D0/3.D0+154112.D0/3.D0*z)*t38-128.D0/3.D0*(7+7*s2l1)*t6+(1331
     #2.D0/3.D0*s1l1+7168.D0/3.D0*s2l1+8704.D0/3.D0+(487040.D0/3.D0*s1l1
     #+265472.D0/3.D0-107648.D0/3.D0*s2l1)*t10+(-27392.D0/3.D0*s1l1-8076
     #8.D0/3.D0-23552.D0/3.D0*s2l1)*t4+(-88192*s1l1-406144.D0/3.D0-384*s
     #2l1)*t6-7424.D0/3.D0*t11+(-95616*s1l1+35456.D0/3.D0+75520.D0/3.D0*
     #s2l1)*t5+(96640.D0/3.D0*s1l1+209408.D0/3.D0+38144.D0/3.D0*s2l1)*t3
     #2+(-27776.D0/3.D0+20480.D0/3.D0*s2l1+6528*s1l1)*t7)*t160
      t810 = s2+((-111616.D0/3.D0*s1l1-183424.D0/3.D0-120064.D0/3.D0*s2l
     #1)*t32+(191872.D0/3.D0+100736.D0/3.D0*s1l1+128000.D0/3.D0*s2l1)*t7
     #-32512.D0/3.D0*t11+256.D0/3.D0*(100*s1l1+76*s2l1+82)*t6-1024*t777*
     #t10)*t167+(-213760.D0/3.D0*t32+t783-292864.D0/3.D0*t10+11136*t11-1
     #74592.D0/3.D0+291968.D0/3.D0*z+76928*t6+389248.D0/3.D0*t5)*t63+(88
     #32*t11-16384*t32+t412+6528*t7)*t134+9984*s2l1+(-1024*t41*t6-512*(-
     #35*s1l1-38-35*s2l1)*t32+(-20224-12672*s2l1-12672*s1l1)*t7+t126)*t1
     #34
      t841 = 5376*t7
      t845 = 158464.D0/3.D0*t32
      t850 = 1280*t7
      s1 = (-1792.D0/3.D0*t7+5632*t32-10112*t6+59776.D0/3.D0*t10-106624.
     #D0/3.D0*t5+44672.D0/3.D0-98944.D0/3.D0*z)*l1+1280.D0/3.D0*t32-1254
     #4.D0/3.D0*t6+(7552.D0/3.D0*t11-41216*t32-6400.D0/3.D0*t7-112000.D0
     #/3.D0*t10+31104*t6-135680.D0/3.D0+222208.D0/3.D0*z+260864.D0/3.D0*
     #t5)*t63+32000.D0/3.D0*t10+(29440.D0/3.D0*t32+10624.D0/3.D0*t7-1318
     #40.D0/3.D0*t10+61952.D0/3.D0*t6-6272+81920.D0/3.D0*z-4096.D0/3.D0*
     #t11+92800.D0/3.D0*t5)*t38-t368+(21376.D0/3.D0*t11-t411+t412-t841)*
     #t134
      t890 = s1+12800.D0/3.D0+(-45952.D0/3.D0*t11-t298+t845-13312.D0/3.D
     #0*t6+t783)*t167-41984.D0/3.D0*t5+(-t850-10112.D0/3.D0*t32-43264.D0
     #/3.D0*t6+32640*t10-31360*t5+29696.D0/3.D0-99712.D0/3.D0*z+t240)*t3
     #7-40192.D0/3.D0*z+(2048-t402+14336.D0/3.D0*t5+89088*t32-43904.D0/3
     #.D0*t11+73856.D0/3.D0*t7-71424*t10+23424*t6)*t141+(-t378-24320*t7+
     #t380-275968.D0/3.D0*t32+55936.D0/3.D0*t11+512*t6)*t64+(9856*t7+632
     #32.D0/3.D0-160384.D0/3.D0*t10+54784.D0/3.D0*t6-64256.D0/3.D0*z-124
     #16.D0/3.D0*t11-18816*t5+39936*t32)*t160+(-66560.D0/3.D0*t7-2048.D0
     #/3.D0+118784*t10-61824*t6+t334+24704.D0/3.D0*t11-62720.D0/3.D0*t5-
     #138880.D0/3.D0*t32)*t85
      t960 = 8064*t7
      t1038 = 6400.D0/3.D0*t11
      t1041 = 12800.D0/3.D0*z
      t1067 = 7424.D0/3.D0*t10
      t1072 = (t12-t850+t412-t222)*t134
      t1105 = (-16640.D0/3.D0*t11+18560*t7-38656.D0/3.D0*t6+t298+2944*t3
     #2)*t167+(t1038+15872.D0/3.D0-319232.D0/3.D0*t7+74240.D0/3.D0*t32-t
     #1041+492928.D0/3.D0*t6-201856.D0/3.D0*t10-27392*t5)*t85+(t12-53120
     #*t32-59008.D0/3.D0*t7+134528.D0/3.D0-98432.D0/3.D0*t10+141440*t6-1
     #13792.D0/3.D0*z-234112.D0/3.D0*t5)*t63-44416*t5+(-6400*t11+393472.
     #D0/3.D0*t7+77312.D0/3.D0*t5-2048+t402-85120.D0/3.D0*t6-306176.D0/3
     #.D0*t32-18560*t10)*t141-t283+(-t378+8960*t11-74880*t7+129920.D0/3.
     #D0*t32+t1067+21760*t6)*t64+t1072+103424.D0/3.D0*t10+(58112.D0/3.D0
     #*t32-130304.D0/3.D0*t6+8704.D0/3.D0+89600.D0/3.D0*t10-30592.D0/3.D
     #0*t5-30592.D0/3.D0*z)*l1+19456-76160.D0/3.D0*z+(-t12+115840.D0/3.D
     #0*t7+98944*t32+9344.D0/3.D0-923008.D0/3.D0*t6-16000*z+199808*t10-8
     #7424.D0/3.D0*t5)*t160+(-27136.D0/3.D0*t7+153856.D0/3.D0*t6-75776.D
     #0/3.D0*t10-45824.D0/3.D0*t32-24320.D0/3.D0*t5+13696.D0/3.D0-8960*z
     #)*t37+(-t118+68864.D0/3.D0*t7-234368.D0/3.D0+302080.D0/3.D0*z-4736
     #0.D0/3.D0*t32+173568*t5+t504-371200.D0/3.D0*t10)*t38
      t1107 = 141440.D0/3.D0*t10
      t1143 = 5120*t32
      t1224 = 9472.D0/3.D0-39424.D0/3.D0*t5+(2560*t11-5120.D0/3.D0*t7-t8
     #45+81920*t6-26368.D0/3.D0+76288.D0/3.D0*z+4352.D0/3.D0*t5-56320.D0
     #/3.D0*t10)*t63+t1072+(-2048.D0/3.D0*t32+t960-18944+58880.D0/3.D0*z
     #-t254+19456*t6-205312.D0/3.D0*t10+189952.D0/3.D0*t5)*t38+(-20480.D
     #0/3.D0*t11+3968*t7-t383+t298+13696*t32)*t167+(-27008.D0/3.D0*t32-t
     #285+t118+14080.D0/3.D0-49408.D0/3.D0*t5+10752*t6-7936*z+34048.D0/3
     #.D0*t10)*t37-3328*z+(-2816.D0/3.D0*t7+16000.D0/3.D0+35840.D0/3.D0*
     #t32-60800.D0/3.D0*t6-31360.D0/3.D0*z+t173-32128.D0/3.D0*t5)*l1+(58
     #880.D0/3.D0*t5-15360*t11+135040.D0/3.D0*t7-2048+t402+9728.D0/3.D0*
     #t32-33408*t10-47872.D0/3.D0*t6)*t141+1024.D0/3.D0*t32+(-t378+43520
     #.D0/3.D0*t11-51712.D0/3.D0*t7-64256.D0/3.D0*t32+t1067+25344*t6)*t6
     #4-6144*t6+(25600.D0/3.D0*t11+15872.D0/3.D0-45824*t7-40064.D0/3.D0*
     #t32-t1041+103168.D0/3.D0*t6+t1107-39680*t5)*t85+(-10240.D0/3.D0*t1
     #1+33920.D0/3.D0+38528.D0/3.D0*t7+213376.D0/3.D0*t32-21120*z-120704
     #*t6+30080*t10-7040.D0/3.D0*t5)*t160+46208.D0/3.D0*t10
      s1 = (-91904.D0/3.D0*t11-71680.D0/3.D0*t7+805888.D0/3.D0*t32+24832
     #.D0/3.D0-963328.D0/3.D0*t6-60416.D0/3.D0*z+269312.D0/3.D0*t10+7936
     #*t5)*t160-75776.D0/3.D0*t6+(-44288*t11+126464.D0/3.D0*t7-35584.D0/
     #3.D0*t6+t298+38912.D0/3.D0*t32)*t167+128.D0/3.D0*t445*(-62*t38+24*
     #t63-9-44*l1-71*t37)*t489+(71424*t11+9728.D0/3.D0-223744.D0/3.D0*t7
     #-509824.D0/3.D0*t32-6656.D0/3.D0*z+217856*t6-29312.D0/3.D0*t10-114
     #688.D0/3.D0*t5)*t85+(-106752*t11+474496.D0/3.D0*t7+71168.D0/3.D0*t
     #5-2048+t402-210304.D0/3.D0*t6+89216.D0/3.D0*t32-98176.D0/3.D0*t10)
     #*t141+t1105*s1l1+t1107+(-148096.D0/3.D0*t32+156032.D0/3.D0*t6-1689
     #6*t5+15616.D0/3.D0*t11+t841+15488.D0/3.D0-28672.D0/3.D0*z-t221)*t3
     #7
      t1226 = s1+14080+(19200*t11+82048.D0/3.D0-366592.D0/3.D0*t32+8320*
     #t7-41216.D0/3.D0*z+68352*t10-97792*t5+301312.D0/3.D0*t6)*t63-41472
     #*t5+(-180608.D0/3.D0+2432.D0/3.D0*t32+60800.D0/3.D0*t7-15616*t11-1
     #78688*t10+507904.D0/3.D0*t5+64512*t6+216064.D0/3.D0*z)*t38+(-t378+
     #276736.D0/3.D0*t11-361088.D0/3.D0*t7-6144*t32+19712.D0/3.D0*t10+89
     #216.D0/3.D0*t6)*t64+(9216*t11-5120*t7+t412-t1143)*t134+t1143-51968
     #.D0/3.D0*z+(-37248*t6+12416.D0/3.D0-33280.D0/3.D0*z-13568.D0/3.D0*
     #t5+32896.D0/3.D0*t10+106496.D0/3.D0*t32-30976.D0/3.D0*t7)*l1+t1224
     #*s2l1
      s1 = t890*s2l1+128.D0/3.D0*(-131*s1l1-299*s2l1-299)*t4+128.D0/3.D0
     #*(7*s1l1-15*s2l1-15)*t10+(-128.D0/3.D0*(337*s1l1+337*s2l1+847)*t10
     #+25600*s1l1+56320.D0/3.D0*s2l1+25088+256.D0/3.D0*(-8+t308)*t32-128
     #.D0/3.D0*(-53*s1l1-66*s2l1-239)*t6+128.D0/3.D0*(889*s1l1+588*s2l1+
     #1357)*t5-128.D0/3.D0*(1069*s1l1+635*s2l1+1185)*t4)*l1-128.D0/3.D0*
     #(-22*s1l1-153*s2l1-153)*t5+((512*s1l1+6400.D0/3.D0+70528*s2l1)*t32
     #+2048*t777*t4+(43904*s1l1+85632+45440*s2l1)*t10-512.D0/3.D0*(76*s1
     #l1+52*s2l1+58)*t5-83200.D0/3.D0*t11+(74240.D0/3.D0*s1l1-315392.D0/
     #3.D0-262784.D0/3.D0*s2l1)*t6+(38272*s1l1+171136.D0/3.D0+80384.D0/3
     #.D0*s2l1)*t7)*t141+(t960-6656.D0/3.D0+360704.D0/3.D0*t10-86656*t6+
     #9728.D0/3.D0*z+31616*t11-25728*t5-45312*t32)*t85-77056.D0/3.D0*z+(
     #-40576.D0/3.D0*t7+58624*t32+65152.D0/3.D0-26880*t10+10112.D0/3.D0*
     #t6-67712.D0/3.D0*z-24448*t5-53504.D0/3.D0*t11)*t160
      s2 = s1-37120.D0/3.D0*t6+26368*t10+7168*s1l1+(-50688*t11+80640*t32
     #+2048-t402+155648.D0/3.D0*t6-74112*t10+17408.D0/3.D0*t5-46720.D0/3
     #.D0*t7)*t141
      t1302 = s2-97280.D0/3.D0*t5+t67*t1226/2+((-4096.D0/3.D0*s1l1-21017
     #6.D0/3.D0-152704.D0/3.D0*s2l1)*t32-512.D0/3.D0*(-28*s1l1-28*s2l1-3
     #7)*t4+38656.D0/3.D0*t11+(-78208.D0/3.D0*s1l1-26624.D0/3.D0*s2l1-70
     #40)*t7+(-291712.D0/3.D0*s1l1-139392+6272*s2l1)*t10+(201088.D0/3.D0
     #*s1l1+88576.D0/3.D0-1920*s2l1)*t5-1024*s1l1-1024-1024*s2l1+(128.D0
     #/3.D0*s1l1+506368.D0/3.D0+126208.D0/3.D0*s2l1)*t6)*t85+(-69760.D0/
     #3.D0*t6-56320*t5+t413+65024.D0/3.D0-11392.D0/3.D0*t32-153472.D0/3.
     #D0*z+147584.D0/3.D0*t10+7168.D0/3.D0*t11)*t37+(t1038+(70016.D0/3.D
     #0*s1l1+350336.D0/3.D0+129536.D0/3.D0*s2l1)*t10+(412160.D0/3.D0*s1l
     #1+136448+157312.D0/3.D0*s2l1)*t4+(-2944*s1l1-21376.D0/3.D0-13568.D
     #0/3.D0*s2l1)*t7+(-112640*s1l1-187648-202624.D0/3.D0*s2l1)*t5+(-755
     #2*s1l1+9856.D0/3.D0-23680.D0/3.D0*s2l1)*t32+(5376*s1l1-68864.D0/3.
     #D0-3328.D0/3.D0*s2l1)*t6-61440*s1l1-49536-87040.D0/3.D0*s2l1)*t63+
     #(-t378+41344.D0/3.D0*t7+84736.D0/3.D0*t10-244352.D0/3.D0*t32+17190
     #4.D0/3.D0*t11-47488.D0/3.D0*t6)*t64
      t1305 = (l1+z-t459)**2
      taux18 = (t810+t1302)/t1305*t462

CCC Coefficients

      tcoff33 = taux17
      tcoff34 = taux17
      tcoff35 = taux18
      tcoff36 = taux18

CCC Jet functions

      tFjet33 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff33)
      tFjet34 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff34)
      tFjet35 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff35)
      tFjet36 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff36)

CCC Sum

      qgRVnum9e3 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFjet
     #36*tcoff36
      RETURN
      END

      double precision function qgRVnum10e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = s2l1**2
      t2 = 256.D0/3.D0*t1
      t3 = s1l1**2
      t4 = 256.D0/3.D0*t3
      t6 = 320.D0/3.D0*s2l1
      t8 = 1-z
      t10 = 5*s1l1
      t11 = 3*s2l1
      t13 = t8**2
      t14 = t13**2
      t15 = t14*t13
      t18 = 1024.D0/3.D0*t1
      t20 = 320*s2l1
      t21 = s1l1*s2l1
      t22 = 256.D0/3.D0*t21
      t25 = 512.D0/3.D0*t3
      t26 = 448.D0/3.D0*s1l1
      t27 = 448.D0/3.D0*s2l1
      t28 = 512.D0/3.D0*t21
      t30 = t13*t8
      t32 = 64*t1
      t33 = 320.D0/3.D0*t3
      t34 = 192*s1l1
      t35 = 832.D0/3.D0*s2l1
      t38 = t14*t30
      t39 = 448.D0/3.D0*t38
      t40 = 448.D0/3.D0*t1
      t41 = 64*t3
      t42 = 64.D0/3.D0*s1l1
      t44 = t14*t8
      t47 = l1**2
      t48 = t47**2
      t55 = 512.D0/3.D0*t1
      t56 = 128*s2l1
      t59 = 128*t38
      t60 = 448.D0/3.D0*t3
      t61 = 576*s1l1
      t62 = 128.D0/3.D0*t21
      t65 = 320.D0/3.D0*t1
      t66 = 256.D0/3.D0*s1l1
      t67 = 256.D0/3.D0*s2l1
      t71 = 320*t3
      t73 = 256*t21
      t79 = 64.D0/3.D0*t1
      t80 = 64.D0/3.D0*t3
      t81 = 2*s2l1
      t82 = 7*s1l1
      t86 = 896.D0/3.D0*t1
      t90 = -t86+1024.D0/3.D0*t3+1408.D0/3.D0*s1l1-1472.D0/3.D0*s2l1+t62
     #-64.D0/3.D0
      t93 = 128*t3
      t98 = 192*t1
      t99 = 512.D0/3.D0*s1l1
      t100 = 704.D0/3.D0*s2l1
      t107 = 192*t3
      t108 = t32-t107-t66+t100+256.D0/3.D0
      t110 = 64.D0/3.D0*s2l1
      t111 = t79+t80-t59+64.D0/3.D0*(1-t81+t82)*t15+t90*t13+(256*t1-t93-
     #704.D0/3.D0-512*s1l1+192*s2l1-t22)*t30+(-t98-t93+320.D0/3.D0-t99-t
     #100-t28)*t44+(192+t40+t4+t99+1088.D0/3.D0*s2l1+640.D0/3.D0*t21)*t1
     #4+t108*t8-t42-t110-64.D0/3.D0
      t113 = 640.D0/3.D0*t1
      t119 = 128*t1
      t120 = 832.D0/3.D0*s1l1
      t123 = 704.D0/3.D0*t1
      t127 = 256*s1l1
      t135 = t113+640.D0/3.D0+t41-64.D0/3.D0*(s2l1+5)*t15+(t34-t56+128-t
     #2)*t44+(-640.D0/3.D0-t119-t33-t120-t20-t28)*t8+(192-t123+t4+1088.D
     #0/3.D0*s1l1-t6-t62)*t13+(448.D0/3.D0-t25-t127+640.D0/3.D0*s2l1-t22
     #)*t14+(-832.D0/3.D0+320*t1+t25-t26+t6+t28)*t30+t127+1280.D0/3.D0*s
     #2l1+t73
      t138 = 2*l6-1
      t139 = t138**2
      t142 = l5-1
      t144 = -1+l1
      t145 = z*l1
      t147 = (1-l1+t145)**2
      t148 = 1/t147
      t153 = dsqrt(z)
      t155 = dsqrt(l5)
      t157 = dsqrt(-t142)
      t160 = dsqrt(l1)
      t164 = l5*z*l1-l1+1-l5+l5*l1
      t167 = dsqrt(-t144)
      t168 = t148*t167
      t172 = t164**2
      t174 = 4*l1*t139*z*l5*t142*t144*t148-4*t138*t153*t155*t157*t160*t1
     #64*t168+t172*t148
      t177 = t47*l1
      t178 = 128.D0/3.D0*t177
      t182 = 320.D0/3.D0*t47
      t186 = dsqrt(t174)
      t194 = 64*s1l1
      t208 = z**2
      t209 = t208**2
      t211 = t153*t209*t160
      t216 = t157*t167
      t217 = t216*t148
      t227 = 512.D0/3.D0*s2l1
      t230 = 384*s2l1
      t239 = -256+t2-640.D0/3.D0*t3-t61-t227-t73
      t243 = 64*s2l1
      t244 = 64-t119+t25+128.D0/3.D0*t38+64.D0/3.D0*(t81+7-3*s1l1)*t15+(
     #256-t119+t60+t61-t227+t62)*t14+(-256-t123-t71-t61-t230-1024.D0/3.D
     #0*t21)*t30+(384+t98+t107+t61+1664.D0/3.D0*s2l1+384*t21)*t13+t239*t
     #8+(-384+t55+t41-t99+t27+t22)*t44+t120-t243+t22
      t246 = t208*z
      t248 = t208*l1
      t250 = 2*t145
      t253 = t153*t246*(t248+3*z-t250+l1-1)*t160
      t254 = t148*l6
      s1 = -64+((-t2+t4+320.D0/3.D0*s1l1-t6)*t8-64.D0/3.D0*(1+t10-t11)*t
     #15+(t18-t4-256.D0/3.D0-320*s1l1+t20+t22)*t13+(256-t18-t25+t26-t27-
     #t28)*t30+(-128-t32+t33+t34-t35)*t14+t39+(t40+t41-512.D0/3.D0-t42+t
     #27+t22)*t44)*t48*l1+(64.D0/3.D0*(-3-t11+4*s1l1)*t15+(384-t55-t41+t
     #26-t56-t22)*t44-t59+(-192+t55-t60-t61+t35-t62)*t14+(-t65-t33-t66-t
     #67-128*t21)*t13+(t65+t71+1280.D0/3.D0*s1l1+t73)*t30)*t48*t47+t111*
     #t48+t135*t47+64.D0/3.D0*t47*t174
      t267 = s1+((128.D0/3.D0-t178-64*t47)*s1l1+(-128.D0/3.D0+t178-t182)
     #*s2l1-t182)*t186+(t39+(-640.D0/3.D0+t55+t93+t99+t27+t28)*t44+(128.
     #D0/3.D0*t1-128.D0/3.D0*t3-t66+t67)*t30+(-t113-t4+320.D0/3.D0+t194-
     #t35-t28)*t14-64.D0/3.D0*(t82-t81+2)*t15)*t48*t177-64.D0/3.D0*(1+s2
     #l1)*(t81+3)*t14+(-512.D0/3.D0*t211*t157*t168*l6+256.D0/3.D0*t211*t
     #217)*t155+t244*t177+(512.D0/3.D0*t253*t216*t254-256.D0/3.D0*t253*t
     #217)*t155*l5-64.D0/3.D0*t209*z*l1*t148
      t268 = t144*t209
      t269 = l6**2
      t272 = 1024.D0/3.D0*t268*t148*t269
      t274 = 1024.D0/3.D0*t268*t254
      t277 = t208*t47
      t288 = z*t47
      t291 = -4*z-t250+l1-2*t47-14*t277+6*t246*t47-4*z*t177+6*t208*t177-
     #4*t246*t177+t209*t177+10*t288+9*t248+t177
      t296 = l5**2
      t306 = t48**2
      t352 = (-t272+t274-64.D0/3.D0*t291*t246*t148)*t296+((64.D0/3.D0*s1
     #l1-64.D0/3.D0*s2l1)*t44-64.D0/3.D0*(4+s1l1-s2l1)*t15+256.D0/3.D0*t
     #38)*t306*l1-t239*t8-t22+(t272-t274+128.D0/3.D0*(t277+3*t145-2*t288
     #+t47+l1-2)*t209*t148)*l5+(-384-t2-t4-1472.D0/3.D0*s1l1-1408.D0/3.D
     #0*s2l1-t73)*t13+(-t79-t80+64.D0/3.D0+(704.D0/3.D0-t86+t25+1792.D0/
     #3.D0*s1l1-t35+t22)*t30+(-896.D0/3.D0+t32-t67-t62-704.D0/3.D0*s1l1)
     #*t14-t108*t8-t90*t13+(320.D0/3.D0+t79+t67)*t44+t42+t110)*l1+(64.D0
     #/3.D0*(6+t10-t81)*t15-t39+(-64.D0/3.D0+t32+t41+t194+t243+t22)*t14+
     #(128.D0/3.D0-t32-t41-t99-t110-t22)*t44)*t306+(-64.D0/3.D0*t38+64.D
     #0/3.D0*t15)*t306*t47+(t230+256+t119+t22+t26)*t30+t119+t243-t25-t12
     #0
      t354 = zqe0(l1,z,l3)
      t357 = (l1+z-t145)**2
      taux19 = -(t267+t352)*t354/t357*t148
      t1 = 1-z
      t2 = t1**2
      t3 = t2*t1
      t4 = t2**2
      t5 = t4*t3
      t6 = t4*t2
      t9 = l1**2
      t10 = t9**2
      t11 = t10**2
      t14 = -s1l1+s2l1
      t15 = 2*t14
      t16 = t4*t1
      t19 = 2*s1l1
      t20 = 2*s2l1
      t31 = 10*s1l1
      t32 = 4*s2l1
      t36 = 8*s1l1
      t37 = 6*s2l1
      t41 = 1792.D0/3.D0*t5
      t44 = 14*s1l1
      t48 = 4*s1l1
      t59 = t9*l1
      t65 = 8*s2l1
      t75 = 34*s1l1
      t79 = 512*t5
      t83 = 6*s1l1
      t88 = 128.D0/3.D0*t15*t1
      t90 = 10*s2l1
      t96 = 128.D0/3.D0*(-12-t44-t20)*t2
      t100 = 14*s2l1
      t103 = 128.D0/3.D0*(-28-t31-t100)*t3
      t109 = 128.D0/3.D0*(-4-t37)*t1
      t117 = 256.D0/3.D0*s1l1
      t118 = 256.D0/3.D0*s2l1
      t122 = 16*s1l1
      t125 = 128.D0/3.D0*(8+t122-t20)*t2
      t131 = 5*s1l1
      t140 = 128.D0/3.D0*(8+t83+t65)*t1
      t158 = 12*s2l1
      s2 = (256.D0/3.D0*t5-256.D0/3.D0*t6)*t11*t9+(128.D0/3.D0*t15*t16-1
     #28.D0/3.D0*(-t19-8+t20)*t6-1024.D0/3.D0*t5)*t11*l1+(-128.D0/3.D0*(
     #-6-t19-t20)*t4+128.D0/3.D0*(-12-t31+t32)*t6-128.D0/3.D0*(-t36+t37+
     #8)*t16+t41)*t11+(128.D0/3.D0*(-18-t44+t20)*t4+128.D0/3.D0*(28-t48+
     #t37)*t16-t41-512.D0/3.D0*t14*t3-128.D0/3.D0*(t32-4-t44)*t6)*t10*t5
     #9
      s1 = s2+(128.D0/3.D0*(-6-t48-t32)*t2-128.D0/3.D0*(42+t44+t65)*t16+
     #128.D0/3.D0*(6-t36+t37)*t6-256.D0/3.D0*(-5+t48-t37)*t3-128.D0/3.D0
     #*(-20-t75+t37)*t4+t79)*t10*t9+(-t41+128.D0/3.D0*(t83+t37+24)*t16+t
     #88+128.D0/3.D0*(4-18*s1l1+t90)*t4-t96-128.D0/3.D0*(-2-t31+t37)*t6+
     #t103)*t10*l1+(t109+128.D0/3.D0*(-t44-2+t32)*t6+128.D0/3.D0*(28+20*
     #s1l1+t100)*t3+256.D0/3.D0+t117+t118-128.D0/3.D0*(8-t48+t100)*t4+t7
     #9-t125-128.D0/3.D0*(20-t48+t20)*t16)*t10+(256.D0/3.D0*(-1+t131-t20
     #)*t3-128.D0/3.D0*(6+t19+t65)*t2+t140-256.D0/3.D0-t117-t118-128.D0/
     #3.D0*(-t37-t122-42)*t16-512.D0/3.D0*t5+128.D0/3.D0*(-22-t75+t32)*t
     #4-256.D0/3.D0*(t20+7-3*s1l1)*t6)*t59
      t188 = s1+256.D0/3.D0+(256.D0/3.D0*(-7*s1l1-8)*t16-t88-t103+128.D0
     #/3.D0*(t122-t158-10)*t4-256.D0/3.D0*(-s2l1-5)*t6+t96)*t9+(512.D0/3
     #.D0*(-2-s2l1)*t16-256.D0/3.D0*(-13-t37-t131)*t4-128.D0/3.D0*(24*s1
     #l1+t90+28)*t3-t109+t125-256.D0/3.D0-t117-t118)*l1+128.D0/3.D0*(t83
     #+t158+12)*t2+256.D0/3.D0*(-s1l1-4-t32)*t3-256.D0/3.D0*(-1-s2l1)*t4
     #-t140+t117+t118
      t189 = zqe0(l1,z,l3)
      t191 = z*l1
      t193 = (l1+z-t191)**2
      t196 = (1-l1+t191)**2
      taux20 = -t188*t189/t193/t196

CCC Coefficients

      tcoff37 = taux19
      tcoff38 = taux19
      tcoff39 = taux20
      tcoff40 = taux20

CCC Jet functions

      tFjet37 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff38)
      tFjet39 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff40)

CCC Sum

      qgRVnum10e0 = tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tFje
     #t40*tcoff40
      RETURN
      END

      double precision function qgRVnum10e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = s2l1**2
      t2 = 64.D0/3.D0*t1
      t3 = s1l1**2
      t4 = 64.D0/3.D0*t3
      t5 = 1-z
      t6 = t5**2
      t7 = t6*t5
      t8 = t6**2
      t9 = t8*t7
      t10 = 128*t9
      t11 = 2*s2l1
      t12 = 7*s1l1
      t14 = t8*t6
      t17 = 256*t1
      t18 = 128*t3
      t19 = 512*s1l1
      t20 = 192*s2l1
      t21 = s1l1*s2l1
      t22 = 256.D0/3.D0*t21
      t25 = 192*t1
      t26 = 512.D0/3.D0*t21
      t27 = 512.D0/3.D0*s1l1
      t28 = 704.D0/3.D0*s2l1
      t30 = t8*t5
      t32 = 64*t1
      t33 = 192*t3
      t34 = 256.D0/3.D0*s1l1
      t35 = t28+t32-t33-t34+256.D0/3.D0
      t37 = 64.D0/3.D0*s1l1
      t39 = 896.D0/3.D0*t1
      t40 = 1024.D0/3.D0*t3
      t41 = 1472.D0/3.D0*s2l1
      t42 = 128.D0/3.D0*t21
      t43 = 1408.D0/3.D0*s1l1-t39+t40-t41+t42-64.D0/3.D0
      t45 = 64.D0/3.D0*s2l1
      t46 = 448.D0/3.D0*t1
      t47 = 256.D0/3.D0*t3
      t49 = 640.D0/3.D0*t21
      t52 = t2+t4-t10-64.D0/3.D0+64.D0/3.D0*(1-t11+t12)*t14+(t17-t18-704
     #.D0/3.D0-t19+t20-t22)*t7+(320.D0/3.D0-t25-t18-t26-t27-t28)*t30+t35
     #*t5-t37+t43*t6-t45+(192+t27+t46+t47+1088.D0/3.D0*s2l1+t49)*t8
      t53 = l1**2
      t54 = t53**2
      t56 = -1+l1
      t57 = z**2
      t58 = t57**2
      t59 = t56*t58
      t60 = z*l1
      t62 = (1-l1+t60)**2
      t63 = 1/t62
      t64 = l6**2
      t65 = t63*t64
      t67 = 1024.D0/3.D0*t59*t65
      t68 = t57*t53
      t70 = z*t53
      t76 = t63*l6
      t78 = 1024.D0/3.D0*t59*t76
      t81 = 1+s2l1
      t86 = 5*s1l1
      t90 = 64*t3
      t93 = 448.D0/3.D0*t9
      t94 = 64*s1l1
      t95 = 64*s2l1
      t99 = t54**2
      t101 = dsqrt(z)
      t103 = dsqrt(l1)
      t104 = t101*t58*t103
      t105 = 1-l5
      t106 = dsqrt(t105)
      t108 = -t56
      t109 = dsqrt(t108)
      t110 = t109*t63
      t112 = t104*t106*t110*l6
      t114 = t106*t109
      t115 = t114*t63
      t116 = t104*t115
      t119 = dsqrt(l5)
      t122 = 2*l6-1
      t123 = t122**2
      t138 = l5*z*l1-l1+1-l5+l5*l1
      t143 = t138**2
      t145 = -4*l1*t123*z*l5*t105*t56*t63-4*t122*t101*t119*t106*t103*t13
     #8*t110+t143*t63
      t146 = t53*t145
      t149 = 2*t60
      t150 = 2*t53
      t152 = t57*z
      t153 = t152*t53
      t155 = t53*l1
      t156 = z*t155
      t158 = t57*t155
      t160 = t152*t155
      t162 = t58*t155
      t164 = t57*l1
      t166 = -4*z-t149+l1-t150-14*t68+6*t153-4*t156+6*t158-4*t160+t162+1
     #0*t70+9*t164+t155
      t171 = l5**2
      t181 = t99*l1
      t183 = 640.D0/3.D0*t1
      t184 = 832.D0/3.D0*s2l1
      t190 = 512.D0/3.D0*t1
      t191 = 448.D0/3.D0*s2l1
      t194 = 128.D0/3.D0*t1
      t195 = 128.D0/3.D0*t3
      t196 = 256.D0/3.D0*s2l1
      t200 = t54*t155
      t202 = 128*t1
      t205 = t58*z*l1*t63
      t207 = -64+t52*t54+(t67+128.D0/3.D0*(t68+3*t60-2*t70+t53+l1-2)*t58
     #*t63-t78)*l5-64.D0/3.D0*t81*(t11+3)*t8+(64.D0/3.D0*(6+t86-t11)*t14
     #+(128.D0/3.D0-t45-t32-t90-t27-t22)*t30-t93+(-64.D0/3.D0+t94+t32+t9
     #0+t95+t22)*t8)*t99+(-512.D0/3.D0*t112+256.D0/3.D0*t116)*t119+64.D0
     #/3.D0*t146+(-t67-64.D0/3.D0*t166*t152*t63+t78)*t171+((64.D0/3.D0*s
     #1l1-64.D0/3.D0*s2l1)*t30-64.D0/3.D0*(4+s1l1-s2l1)*t14+256.D0/3.D0*
     #t9)*t181+(t93+(320.D0/3.D0-t183-t47-t26+t94-t184)*t8-64.D0/3.D0*(t
     #12-t11+2)*t14+(-640.D0/3.D0+t27+t190+t18+t191+t26)*t30+(-t34+t194-
     #t195+t196)*t7)*t200+t202+t95-64.D0/3.D0*t205
      t208 = 832.D0/3.D0*s1l1
      t209 = 512.D0/3.D0*t3
      t210 = 576*s1l1
      t211 = 448.D0/3.D0*t3
      t212 = 512.D0/3.D0*s2l1
      t216 = 704.D0/3.D0*t1
      t217 = 320*t3
      t218 = 384*s2l1
      t219 = 1024.D0/3.D0*t21
      t228 = 1664.D0/3.D0*s2l1
      t232 = 256.D0/3.D0*t1
      t233 = 640.D0/3.D0*t3
      t234 = 256*t21
      t235 = -256-t212+t232-t233-t210-t234
      t237 = 64+(256+t210-t202+t211-t212+t42)*t8-t202+t209+128.D0/3.D0*t
     #9-t95+(-256-t210-t216-t217-t218-t219)*t7+(-384-t27+t190+t90+t191+t
     #22)*t30+t208+64.D0/3.D0*(t11+7-3*s1l1)*t14+(384+t210+t25+t33+t228+
     #384*t21)*t6+t235*t5+t22
      t243 = t101*t152*(t164+3*z-t149+l1-1)*t103
      t245 = t243*t114*t76
      t247 = t243*t115
      t250 = t119*l5
      t252 = 704.D0/3.D0*s1l1
      t257 = 1792.D0/3.D0*s1l1
      t266 = 320.D0/3.D0*t3
      t267 = 320*s2l1
      t270 = 256*s1l1
      t272 = 320.D0/3.D0*s2l1
      t275 = 128*s2l1
      t276 = 192*s1l1
      t282 = 1280.D0/3.D0*s2l1
      t283 = 448.D0/3.D0*s1l1
      t290 = (-640.D0/3.D0-t208-t202-t266-t267-t26)*t5+t183+640.D0/3.D0+
     #t90+t270+(192+1088.D0/3.D0*s1l1-t216+t47-t272-t42)*t6+(-t275+128+t
     #276-t232)*t30-64.D0/3.D0*(s2l1+5)*t14+t234+t282+(-832.D0/3.D0-t283
     #+320*t1+t209+t272+t26)*t7+(448.D0/3.D0+640.D0/3.D0*s2l1-t209-t270-
     #t22)*t8
      t292 = 320.D0/3.D0*s1l1
      t295 = 3*s2l1
      t299 = 1024.D0/3.D0*t1
      t310 = t54*l1
      t318 = 320.D0/3.D0*t1
      t319 = 128*t21
      t322 = 1280.D0/3.D0*s1l1
      t328 = t54*t53
      t333 = 1472.D0/3.D0*s1l1
      t336 = -t9+t14
      t338 = t99*t53
      t340 = 128.D0/3.D0*t155
      t344 = 320.D0/3.D0*t53
      t348 = dsqrt(t145)
      s1 = -t208-t209-t22+t237*t155+(512.D0/3.D0*t245-256.D0/3.D0*t247)*
     #t250+(-t2-t4+t45+(-t252-896.D0/3.D0+t32-t42-t196)*t8+(320.D0/3.D0+
     #t2+t196)*t30+64.D0/3.D0+(704.D0/3.D0+t257-t39+t209-t184+t22)*t7-t4
     #3*t6-t35*t5+t37)*l1+t290*t53
      t352 = s1+((-t272-t232+t47+t292)*t5-64.D0/3.D0*(1+t86-t295)*t14+(2
     #56+t283-t299-t209-t191-t26)*t7+(t46+t90-512.D0/3.D0-t37+t191+t22)*
     #t30+t93+(t299-t47-256.D0/3.D0-320*s1l1+t267+t22)*t6+(-128-t32+t266
     #+t276-t184)*t8)*t310+(64.D0/3.D0*(-3-t295+4*s1l1)*t14+(-192-t210+t
     #190-t211+t184-t42)*t8+(-t34-t318-t266-t196-t319)*t6+(t322+t318+t21
     #7+t234)*t7+(384-t275-t190-t90+t283-t22)*t30-t10)*t328+(t218+t202+t
     #22+256+t283)*t7+(-384-1408.D0/3.D0*s2l1-t232-t47-t333-t234)*t6+64.
     #D0/3.D0*t336*t338+((128.D0/3.D0-t340-64*t53)*s1l1-t344+(-t344-128.
     #D0/3.D0+t340)*s2l1)*t348-t235*t5
      t354 = zqe0(l1,z,l3)
      t355 = dlog(t108)
      t357 = sqe1(l1)
      t358 = zqe1(l1,z,l3)
      t359 = sbqe0(l1)
      t360 = dlog(z)
      t365 = 6*s2l1
      t366 = 14*s1l1
      t370 = 1408.D0/3.D0*t1
      t371 = 3776.D0/3.D0*s1l1
      t374 = -384+t370-t371-t22+1984.D0/3.D0*s2l1-2048.D0/3.D0*t3
      t377 = 896.D0/3.D0*t21
      t378 = 256*t3
      t381 = 1216.D0/3.D0*s2l1
      t384 = 1664.D0/3.D0*t21
      t385 = 2560.D0/3.D0*s2l1
      t386 = 1408.D0/3.D0*t3
      t391 = 128.D0/3.D0-t190+1600.D0/3.D0*s1l1-t381+1664.D0/3.D0*t3
      t397 = 2*t108*t58
      t399 = 1024.D0/3.D0*t397*t65
      t401 = 1024.D0/3.D0*t397*t76
      t405 = 2*l1
      t413 = 896.D0/3.D0*t3
      t414 = 8*s2l1
      t415 = 10*s1l1
      t419 = 1408.D0/3.D0*t9
      t444 = 384*s1l1
      s1 = 640.D0/3.D0+(384*t9-t94+64.D0/3.D0*(t365-4-t366)*t14+t374*t6+
     #(2624.D0/3.D0-t17+2944.D0/3.D0*s1l1+t377+t378)*t7+(-1088.D0/3.D0+t
     #183+t257+t22+t381+t233)*t30+(-1408.D0/3.D0-t39-t333-t384-t385-t386
     #)*t8+t391*t5+t95-t18)*t54-128.D0/3.D0*t146+(t399-t401+128.D0/3.D0*
     #(-2*t68-6*t60+4*t70-t150-t405+4)*t58*t63)*l5-t232+t275
      t455 = s1+128.D0/3.D0*t205+t210+t413+(64.D0/3.D0*(-18+t414-t415)*t
     #14+t419+(256.D0/3.D0-t232-t252-t26-t272-t233)*t8+(-512.D0/3.D0+448
     #*s1l1-t22-t95+t195)*t30)*t99+t377+(-t419+(1856.D0/3.D0-t190-t257-t
     #22-t28-t233)*t30+(64.D0/3.D0-t202+t276-t319-t20)*t7+(-1024.D0/3.D0
     #+t299+t292+1408.D0/3.D0*t21+t228+t386)*t8-64.D0/3.D0*(t365-8-t366)
     #*t14)*t200+(t94+(-896+384*t1-3520.D0/3.D0*s1l1-t26+t20-t378)*t7+(8
     #96.D0/3.D0*s2l1-t194+2432.D0/3.D0+t22+t444)*t8-t391*t5-t374*t6+(-t
     #212-t194-256)*t30-t95+t18)*l1
      t468 = 1280.D0/3.D0*t1
      t475 = 2176.D0/3.D0*s2l1
      t483 = 640*t3
      t487 = -1280.D0/3.D0-t468-t19-64.D0/3.D0*(-14-t11)*t14+(t20-1472.D
     #0/3.D0+t202-1024.D0/3.D0*s1l1)*t30+(1408.D0/3.D0+t39+t257+t219+t47
     #5+t209)*t5+(-1408.D0/3.D0+t183-t19+t22-t196+t209)*t6+(-256.D0/3.D0
     #-t232+2368.D0/3.D0*s1l1+t42-t282+t413)*t8+(1600.D0/3.D0-t39-t208-t
     #49+t95-t483)*t7-512*t21-t385-t18
      t489 = 4*s2l1
      t494 = 1216*s1l1
      t500 = 2176.D0/3.D0*t3
      t510 = 640+t232+1088*s1l1+640*t21+t475+t40
      t514 = -640.D0/3.D0-t10+t232-t210+64.D0/3.D0*(-18-t489+6*s1l1)*t14
     #-t377+(-2560.D0/3.D0+t183-t494+t26+t184-t195)*t8+(704+t299+4736.D0
     #/3.D0*s1l1+1280.D0/3.D0*t21+2368.D0/3.D0*s2l1+t500)*t7+(-896-t370-
     #t371-2432.D0/3.D0*t21-3968.D0/3.D0*s2l1-512*t3)*t6+t510*t5+(3328.D
     #0/3.D0-t183+t283-t22-t28-t18)*t30-t275-t413
      t549 = 640.D0/3.D0*t53
      t577 = 8*z+4*t60-t405+4*t53+28*t68-12*t153+8*t156-12*t158+8*t160-2
     #*t162-20*t70-18*t164-2*t155
      s1 = (-1024.D0/3.D0*t245+512.D0/3.D0*t247)*t250+(1024.D0/3.D0*t112
     #-512.D0/3.D0*t116)*t119+(2560.D0/3.D0+t299+2624.D0/3.D0*s1l1+t384+
     #3584.D0/3.D0*s2l1+t18)*t6+t487*t53+t514*t155-256.D0/3.D0*t336*t338
     #+(-1600.D0/3.D0-832*s2l1-t39-t270-t26)*t7
      s2 = s1+((-128.D0/3.D0+t202-t270+t275-t18)*t5-64.D0/3.D0*(-4+10*s2
     #l1-t415)*t14+(256-t468+t322-t26-1024.D0/3.D0*s2l1-t209)*t6+(-1472.
     #D0/3.D0+t299+t208+t49+t45+t483)*t7+(t190-1664.D0/3.D0*s1l1+t319+16
     #00.D0/3.D0*s2l1-t47)*t8-t419+(1984.D0/3.D0-t202-t292+t22-t275-t195
     #)*t30)*t310-64.D0/3.D0*t81*(-t489-6)*t8+((-256.D0/3.D0+256.D0/3.D0
     #*t155+128*t53)*s1l1+(128.D0/3.D0-t340+t549)*s2l1-64.D0/3.D0+64.D0/
     #3.D0*t155+t549)*t348
      t603 = s2+((-128.D0/3.D0*s1l1+t22+t195+t196+t194)*t30-64.D0/3.D0*(
     #-14+t489-2*s1l1)*t14-896.D0/3.D0*t9)*t181+(-t399+t401-64.D0/3.D0*t
     #577*t152*t63)*t171-t510*t5+(64.D0/3.D0*(t414+4-8*s1l1)*t14+(-3328.
     #D0/3.D0+t190-t292+t191+t47)*t30+1280.D0/3.D0*t9+(2176.D0/3.D0-t39+
     #t494-t26-t41+t195)*t8+(128.D0/3.D0+t202+t444+t234+t275+384*t3)*t6+
     #(-512.D0/3.D0-t194-3968.D0/3.D0*s1l1-t234+128.D0/3.D0*s2l1-t500)*t
     #7)*t328
      t608 = (l1+z-t60)**2
      taux19 = ((t207+t352)*(t354*t355+t357-t358+t359*t360)-(t455+t603)*
     #t354)/t608*t63
      t1 = 6*s2l1
      t3 = 1-z
      t5 = 128.D0/3.D0*(-4-t1)*t3
      t6 = 4*s1l1
      t7 = 2*s2l1
      t9 = t3**2
      t10 = t9**2
      t11 = t10*t3
      t14 = 14*s1l1
      t15 = 4*s2l1
      t17 = t10*t9
      t20 = 14*s2l1
      t24 = t9*t3
      t25 = t10*t24
      t26 = 512*t25
      t27 = 256.D0/3.D0*s1l1
      t28 = 256.D0/3.D0*s2l1
      t29 = 16*s1l1
      t32 = 128.D0/3.D0*(8+t29-t7)*t9
      t38 = l1**2
      t39 = t38**2
      t47 = 8*s2l1
      t51 = 34*s1l1
      t55 = 8*s1l1
      t60 = t39*t38
      t62 = 2*s1l1
      t66 = 1792.D0/3.D0*t25
      t67 = 10*s1l1
      t75 = t39**2
      t79 = t75*t38
      t81 = 5*s1l1
      t99 = 6*s1l1
      t102 = 128.D0/3.D0*(8+t99+t47)*t3
      t104 = t38*l1
      t106 = 12*s2l1
      t115 = 128.D0/3.D0*(-12-t14-t7)*t9
      t119 = 18*s1l1
      t120 = 10*s2l1
      t126 = 128.D0/3.D0*(-28-t67-t20)*t24
      t127 = -s1l1+s2l1
      t128 = 2*t127
      t130 = 128.D0/3.D0*t128*t3
      t132 = t39*l1
      t147 = t39*t104
      t156 = t75*l1
      t179 = 24*s1l1
      t185 = -1-s2l1
      s1 = 256.D0/3.D0+(t5-128.D0/3.D0*(20-t6+t7)*t11+128.D0/3.D0*(-t14-
     #2+t15)*t17-128.D0/3.D0*(8-t6+t20)*t10+t26+256.D0/3.D0+t27+t28-t32+
     #128.D0/3.D0*(28+20*s1l1+t20)*t24)*t39+(128.D0/3.D0*(-6-t6-t15)*t9-
     #256.D0/3.D0*(-5+t6-t1)*t24-128.D0/3.D0*(42+t14+t47)*t11-128.D0/3.D
     #0*(-20-t51+t1)*t10+t26+128.D0/3.D0*(6-t55+t1)*t17)*t60+(-128.D0/3.
     #D0*(-6-t62-t7)*t10+t66+128.D0/3.D0*(-12-t67+t15)*t17-128.D0/3.D0*(
     #-t55+t1+8)*t11)*t75+t27+t28+(256.D0/3.D0*t25-256.D0/3.D0*t17)*t79+
     #(256.D0/3.D0*(-1+t81-t7)*t24-256.D0/3.D0*(t7+7-3*s1l1)*t17-128.D0/
     #3.D0*(6+t62+t47)*t9-128.D0/3.D0*(-t1-t29-42)*t11-512.D0/3.D0*t25-2
     #56.D0/3.D0-t27-t28+128.D0/3.D0*(-22-t51+t15)*t10+t102)*t104
      t188 = s1+128.D0/3.D0*(t99+t106+12)*t9+(-t66+128.D0/3.D0*(t99+t1+2
     #4)*t11-t115-128.D0/3.D0*(-2-t67+t1)*t17+128.D0/3.D0*(4-t119+t120)*
     #t10+t126+t130)*t132+(128.D0/3.D0*(-18-t14+t7)*t10+128.D0/3.D0*(28-
     #t6+t1)*t11-512.D0/3.D0*t127*t24-128.D0/3.D0*(t15-4-t14)*t17-t66)*t
     #147+(128.D0/3.D0*t128*t11-128.D0/3.D0*(-t62-8+t7)*t17-1024.D0/3.D0
     #*t25)*t156+(256.D0/3.D0*(-7*s1l1-8)*t11-256.D0/3.D0*(-s2l1-5)*t17-
     #t126+128.D0/3.D0*(t29-t106-10)*t10-t130+t115)*t38+256.D0/3.D0*(-s1
     #l1-4-t15)*t24+(512.D0/3.D0*(-2-s2l1)*t11-t5+t32-256.D0/3.D0*(-13-t
     #1-t81)*t10-256.D0/3.D0-t27-t28-128.D0/3.D0*(t179+t120+28)*t24)*l1-
     #t102-256.D0/3.D0*t185*t10
      t189 = zqe0(l1,z,l3)
      t191 = dlog(1-l1)
      t193 = sqe1(l1)
      t194 = zqe1(l1,z,l3)
      t195 = sbqe0(l1)
      t196 = dlog(z)
      t205 = -t127
      t206 = 10*t205
      t215 = 12*s1l1
      t219 = 56*s1l1
      t220 = 16*s2l1
      t225 = 30*s2l1
      t229 = 3072*t25
      t232 = 74*s1l1
      t236 = 26*s2l1
      t252 = 40*s2l1
      t253 = 72*s1l1
      t257 = 46*s1l1
      t266 = 174*s1l1
      t267 = 22*s2l1
      t280 = 128.D0/3.D0*t206*t3
      t289 = 128.D0/3.D0*(56+68*s1l1+t1)*t9
      t293 = 58*s2l1
      t296 = 128.D0/3.D0*(128+t257+t293)*t24
      t301 = 128.D0/3.D0*(12-t99+t267)*t3
      t310 = 1024.D0/3.D0*s1l1
      t311 = 1024.D0/3.D0*s2l1
      t320 = 128.D0/3.D0*(-52-86*s1l1-t7)*t9
      t334 = 30*s1l1
      t337 = 128.D0/3.D0*(-40-t334-t252)*t3
      t338 = 1280.D0/3.D0*s1l1
      t339 = 1280.D0/3.D0*s2l1
      t340 = 82*s1l1
      t363 = 5*s2l1
      s2 = (-t26-128.D0/3.D0*(-t62-t7-12)*t17)*t79+(128.D0/3.D0*t206*t11
     #-128.D0/3.D0*(44+t29-t1)*t17+5632.D0/3.D0*t25)*t156+(-128.D0/3.D0*
     #(32+t215+t106)*t10+128.D0/3.D0*(60+t219-t220)*t17-128.D0/3.D0*(-44
     #+38*s1l1-t225)*t11-t229)*t75+(128.D0/3.D0*(90+t232-t120)*t10+128.D
     #0/3.D0*(-140-t236+t29)*t11+t229-2560.D0/3.D0*t205*t24-128.D0/3.D0*
     #(22+t232-t20)*t17)*t147
      s1 = s2+(128.D0/3.D0*(28+t119+18*s2l1)*t9-128.D0/3.D0*(-t252-208-t
     #253)*t11+128.D0/3.D0*(t257-24-24*s2l1)*t17-256.D0/3.D0*(22-22*s1l1
     #+t236)*t24-128.D0/3.D0*(104+t266-t267)*t10-8192.D0/3.D0*t25)*t60+(
     #t229+128.D0/3.D0*(-128-34*s2l1-32*s1l1)*t11+t280+128.D0/3.D0*(-6+9
     #4*s1l1-36*s2l1)*t10-t289-128.D0/3.D0*(-t236+10+t219)*t17+t296)*t13
     #2+(t301+128.D0/3.D0*(10+t253-t220)*t17+128.D0/3.D0*(-148-102*s1l1-
     #78*s2l1)*t24-1024.D0/3.D0-t310-t311-128.D0/3.D0*(-42+t179-72*s2l1)
     #*t10-2560*t25-t320-128.D0/3.D0*(-t1-100+t29)*t11)*t39+(256.D0/3.D0
     #*(2-27*s1l1+t1)*t24-128.D0/3.D0*(-32-t215-42*s2l1)*t9+t337+1280.D0
     #/3.D0+t338+t339-128.D0/3.D0*(208+t340+t225)*t11+2560.D0/3.D0*t25+1
     #28.D0/3.D0*(114+t266-t106)*t10-256.D0/3.D0*(15*s1l1-9*s2l1-34)*t17
     #)*t104
      t395 = s1-1280.D0/3.D0+(256.D0/3.D0*(42+t7+35*s1l1)*t11-t280-t296+
     #128.D0/3.D0*(38-t340+48*s2l1)*t10-256.D0/3.D0*(t363+25)*t17+t289)*
     #t38+(512.D0/3.D0*(10+t363)*t11-256.D0/3.D0*(25*s1l1+31*s2l1+66)*t1
     #0-128.D0/3.D0*(-148-122*s1l1-t293)*t24-t301+t320+1024.D0/3.D0+t310
     #+t311)*l1+128.D0/3.D0*(-60*s2l1-t334-60)*t9+256.D0/3.D0*(t81+20+20
     #*s2l1)*t24+1280.D0/3.D0*t185*t10-t337-t338-t339
      t398 = z*l1
      t400 = (l1+z-t398)**2
      t404 = (1-l1+t398)**2
      taux20 = (t188*(t189*t191+t193-t194+t195*t196)-t395*t189)/t400/t40
     #4

CCC Coefficients

      tcoff37 = taux19
      tcoff38 = taux19
      tcoff39 = taux20
      tcoff40 = taux20

CCC Jet functions

      tFjet37 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff38)
      tFjet39 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff40)

CCC Sum

      qgRVnum10e1 = tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tFje
     #t40*tcoff40
      RETURN
      END

      double precision function qgRVnum10e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t3 = 1-z
      t4 = t3**2
      t5 = t4**2
      t6 = t5*t3
      t9 = t5*t4
      t12 = t4*t3
      t13 = t5*t12
      t14 = 256.D0/3.D0*t13
      t16 = l1**2
      t17 = t16**2
      t18 = t17**2
      t19 = t18*l1
      t22 = 2*l6-1
      t23 = t22**2
      t26 = l5-1
      t28 = -1+l1
      t29 = z*l1
      t31 = (1-l1+t29)**2
      t32 = 1/t31
      t37 = dsqrt(z)
      t39 = dsqrt(l5)
      t41 = dsqrt(-t26)
      t44 = dsqrt(l1)
      t48 = l5*z*l1-l1+1-l5+l5*l1
      t50 = -t28
      t51 = dsqrt(t50)
      t52 = t32*t51
      t56 = t48**2
      t58 = 4*l1*t23*z*l5*t26*t28*t32-4*t22*t37*t39*t41*t44*t48*t52+t56*
     #t32
      t59 = t16*t58
      t61 = z**2
      t62 = t61*z
      t64 = t61*l1
      t66 = 2*t29
      t69 = t37*t62*(t64+3*z-t66+l1-1)*t44
      t70 = t41*t51
      t71 = t32*l6
      t73 = t69*t70*t71
      t75 = t70*t32
      t76 = t69*t75
      t79 = t39*l5
      t80 = (512.D0/3.D0*t73-256.D0/3.D0*t76)*t79
      t81 = 448.D0/3.D0*t13
      t82 = s2l1**2
      t83 = 640.D0/3.D0*t82
      t84 = s1l1**2
      t85 = 256.D0/3.D0*t84
      t86 = s1l1*s2l1
      t87 = 512.D0/3.D0*t86
      t88 = 64*s1l1
      t89 = 832.D0/3.D0*s2l1
      t92 = 7*s1l1
      t93 = 2*s2l1
      t97 = 512.D0/3.D0*s1l1
      t98 = 512.D0/3.D0*t82
      t99 = 128*t84
      t100 = 448.D0/3.D0*s2l1
      t103 = 256.D0/3.D0*s1l1
      t104 = 128.D0/3.D0*t82
      t105 = 128.D0/3.D0*t84
      t106 = 256.D0/3.D0*s2l1
      t110 = t16*l1
      t111 = t17*t110
      t113 = 3*s2l1
      t118 = 576*s1l1
      t119 = 448.D0/3.D0*t84
      t120 = 128.D0/3.D0*t86
      t123 = 320.D0/3.D0*t82
      t124 = 320.D0/3.D0*t84
      t125 = 128*t86
      t128 = 1280.D0/3.D0*s1l1
      t129 = 320*t84
      t130 = 256*t86
      t133 = 128*s2l1
      t134 = 64*t84
      t135 = 448.D0/3.D0*s1l1
      t136 = 256.D0/3.D0*t86
      t139 = 128*t13
      t141 = t17*t16
      t143 = 384*s2l1
      t144 = 128*t82
      t147 = t61**2
      t150 = t147*z*l1*t32
      t151 = 64.D0/3.D0*t150
      t152 = 192*t82
      t153 = 192*t84
      t154 = 1664.D0/3.D0*s2l1
      t158 = 512.D0/3.D0*s2l1
      t159 = 256.D0/3.D0*t82
      t160 = 640.D0/3.D0*t84
      t161 = -256-t158+t159-t160-t118-t130
      t165 = 512.D0/3.D0*t84
      t169 = 704.D0/3.D0*t82
      t170 = 1024.D0/3.D0*t86
      t173 = 64*s2l1
      t178 = 832.D0/3.D0*s1l1
      t179 = 64+(384+t118+t152+t153+t154+384*t86)*t4+t161*t3+(256+t118-t
     #144+t119-t158+t120)*t5+t165+128.D0/3.D0*t13+(-384-t97+t98+t134+t10
     #0+t136)*t6+(-256-t118-t169-t129-t143-t170)*t12-t173+64.D0/3.D0*(t9
     #3+7-3*s1l1)*t9+t178+t136-t144
      t183 = 1408.D0/3.D0*s2l1
      t184 = 1472.D0/3.D0*s1l1
      t187 = t28*t147
      t188 = l6**2
      t189 = t32*t188
      t191 = 1024.D0/3.D0*t187*t189
      t193 = 1024.D0/3.D0*t187*t71
      t195 = 2*t16
      t196 = t61*t16
      t198 = t62*t16
      t200 = z*t110
      t202 = t61*t110
      t204 = t62*t110
      t206 = t147*t110
      t207 = z*t16
      t210 = -4*z-t66+l1-t195-14*t196+6*t198-4*t200+6*t202-4*t204+t206+1
     #0*t207+9*t64+t110
      t215 = l5**2
      t216 = (-t191+t193-64.D0/3.D0*t210*t62*t32)*t215
      t217 = 64.D0/3.D0*t82
      t218 = 64.D0/3.D0*t84
      t219 = 448.D0/3.D0*t82
      t221 = 640.D0/3.D0*t86
      t224 = 704.D0/3.D0*s2l1
      t225 = 64*t82
      t226 = t224+t225-t153-t103+256.D0/3.D0
      t228 = 64.D0/3.D0*s1l1
      t229 = 256*t82
      t230 = 512*s1l1
      t231 = 192*s2l1
      t236 = 64.D0/3.D0*s2l1
      t237 = 1408.D0/3.D0*s1l1
      t238 = 896.D0/3.D0*t82
      t239 = 1024.D0/3.D0*t84
      t240 = 1472.D0/3.D0*s2l1
      t241 = t237-t238+t239-t240+t120-64.D0/3.D0
      t246 = t217+t218+(192+t97+t219+t85+1088.D0/3.D0*s2l1+t221)*t5+t226
     #*t3-t228-t139+(t229-t99-704.D0/3.D0-t230+t231-t136)*t12+(320.D0/3.
     #D0-t152-t99-t87-t97-t224)*t6-64.D0/3.D0-t236+t241*t4+64.D0/3.D0*(1
     #-t93+t92)*t9
      t248 = -64+((64.D0/3.D0*s1l1-64.D0/3.D0*s2l1)*t6-64.D0/3.D0*(4+s1l
     #1-s2l1)*t9+t14)*t19+64.D0/3.D0*t59+t80+(t81+(320.D0/3.D0-t83-t85-t
     #87+t88-t89)*t5-64.D0/3.D0*(t92-t93+2)*t9+(-640.D0/3.D0+t97+t98+t99
     #+t100+t87)*t6+(-t103+t104-t105+t106)*t12)*t111+(64.D0/3.D0*(-3-t11
     #3+4*s1l1)*t9+(-192-t118+t98-t119+t89-t120)*t5+(-t103-t123-t124-t10
     #6-t125)*t4+(t128+t123+t129+t130)*t12+(384-t133-t98-t134+t135-t136)
     #*t6-t139)*t141+(t143+t144+t136+256+t135)*t12-t151+t179*t110-t161*t
     #3+(-384-t183-t159-t85-t184-t130)*t4+t216+t246*t17
      t249 = -t13+t9
      t251 = t18*t16
      t260 = (t191-t193+128.D0/3.D0*(t196+3*t29-2*t207+t16+l1-2)*t147*t3
     #2)*l5
      t261 = 5*s1l1
      t271 = 1+s2l1
      t282 = 1792.D0/3.D0*s1l1
      t285 = 704.D0/3.D0*s1l1
      t290 = 128.D0/3.D0*t110
      t294 = 320.D0/3.D0*t16
      t298 = dsqrt(t58)
      t300 = 320.D0/3.D0*s2l1
      t301 = 320.D0/3.D0*s1l1
      t307 = 1024.D0/3.D0*t82
      t313 = 320*s2l1
      t316 = 192*s1l1
      t320 = t17*l1
      t323 = t37*t147*t44
      t326 = t323*t41*t52*l6
      t328 = t323*t75
      t331 = (-512.D0/3.D0*t326+256.D0/3.D0*t328)*t39
      t334 = 640.D0/3.D0*s2l1
      t335 = 256*s1l1
      t338 = 1280.D0/3.D0*s2l1
      t339 = 320*t82
      t350 = (-t133+128+t316-t159)*t6+t83+640.D0/3.D0+t134+(448.D0/3.D0+
     #t334-t165-t335-t136)*t5+t338+t130+(-832.D0/3.D0-t135+t339+t165+t30
     #0+t87)*t12+t335+(-640.D0/3.D0-t178-t144-t124-t313-t87)*t3-64.D0/3.
     #D0*(s2l1+5)*t9+(192+1088.D0/3.D0*s1l1-t169+t85-t300-t120)*t4
      t352 = 64.D0/3.D0*t249*t251+t144+t173+t260+(64.D0/3.D0*(6+t261-t93
     #)*t9+(128.D0/3.D0-t236-t225-t134-t97-t136)*t6-t81+(-64.D0/3.D0+t88
     #+t225+t134+t173+t136)*t5)*t18-t178-64.D0/3.D0*t271*(t93+3)*t5-t165
     #+(-t241*t4+t228+t236-t217-t226*t3+(320.D0/3.D0+t217+t106)*t6+64.D0
     #/3.D0+(704.D0/3.D0+t282-t238+t165-t89+t136)*t12-t218+(-t285-896.D0
     #/3.D0+t225-t120-t106)*t5)*l1+((128.D0/3.D0-t290-64*t16)*s1l1-t294+
     #(-t294-128.D0/3.D0+t290)*s2l1)*t298-t136+((-t300-t159+t85+t301)*t3
     #-64.D0/3.D0*(1+t261-t113)*t9+(256+t135-t307-t165-t100-t87)*t12+(t2
     #19+t134-512.D0/3.D0-t228+t100+t136)*t6+t81+(t307-t85-256.D0/3.D0-3
     #20*s1l1+t313+t136)*t4+(-128-t225+t124+t316-t89)*t5)*t320+t331+t350
     #*t16
      t354 = zqe0(l1,z,l3)
      t355 = dlog(t50)
      t356 = t355**2
      t359 = sqe1(l1)
      t360 = zqe1(l1,z,l3)
      t363 = sqe2(l1)
      t364 = zqe2(l1,z,l3)
      t365 = sbqe0(l1)
      t366 = dlog(z)
      t367 = t366**2
      t370 = sbqe1(l1)
      t374 = 6*s2l1
      t375 = 14*s1l1
      t379 = 1408.D0/3.D0*t86
      t380 = 1408.D0/3.D0*t84
      t387 = 1408.D0/3.D0*t13
      t395 = 2176.D0/3.D0*s2l1
      t396 = -640-640*t86-t159-1088*s1l1-t239-t395
      t399 = 2*t50*t147
      t401 = 1024.D0/3.D0*t399*t189
      t403 = 1024.D0/3.D0*t399*t71
      t406 = 2*l1
      t417 = 8*z+4*t29-t406+4*t16+28*t196-12*t198+8*t200-12*t202+8*t204-
     #2*t206-20*t207-18*t64-2*t110
      t428 = 640.D0/3.D0*t16
      t431 = 64.D0/3.D0*t110
      t436 = 1408.D0/3.D0*t82
      t437 = 3776.D0/3.D0*s1l1
      t440 = t136+384-t436+t437-1984.D0/3.D0*s2l1+2048.D0/3.D0*t84
      t442 = 896.D0/3.D0*s2l1
      t443 = 384*s1l1
      t447 = 1216.D0/3.D0*s2l1
      t449 = -128.D0/3.D0-1600.D0/3.D0*s1l1+t98+t447-1664.D0/3.D0*t84
      t453 = 256*t84
      t458 = 896.D0/3.D0*t13
      t459 = 128.D0/3.D0*s1l1
      t462 = 4*s2l1
      t474 = 1280.D0/3.D0*t82
      t475 = 1024.D0/3.D0*s1l1
      t478 = 640*t84
      t481 = 2560.D0/3.D0*s2l1
      t485 = 896.D0/3.D0*t84
      t493 = -512*t86-t474-t230+(t144-1472.D0/3.D0-t475+t231)*t6+(-t221+
     #1600.D0/3.D0-t238+t173-t178-t478)*t12-t481-1280.D0/3.D0-64.D0/3.D0
     #*(-14-t93)*t9-t99+(t120-256.D0/3.D0+t485-t159-t338+2368.D0/3.D0*s1
     #l1)*t5+(t136-1408.D0/3.D0+t165+t83-t230-t106)*t4+(t170+1408.D0/3.D
     #0+t238+t282+t395+t165)*t3
      t500 = 10*s1l1
      t504 = 1024.D0/3.D0*s2l1
      t507 = 1664.D0/3.D0*s1l1
      s1 = 640.D0/3.D0+(-64.D0/3.D0*(t374-8-t375)*t9+(-1024.D0/3.D0+t379
     #+t307+t301+t154+t380)*t5+(-t136+1856.D0/3.D0-t98-t282-t224-t160)*t
     #6+(-t125+64.D0/3.D0+t316-t144-t231)*t12-t387)*t111-256.D0/3.D0*t24
     #9*t251+t396*t3+(-t401+t403-64.D0/3.D0*t417*t62*t32)*t215-128.D0/3.
     #D0*t59
      t516 = s1+((256.D0/3.D0*t110-256.D0/3.D0+128*t16)*s1l1+t428+(-t290
     #+t428+128.D0/3.D0)*s2l1-64.D0/3.D0+t431)*t298+(t99+(-t158-t104-256
     #)*t6-t173+t88+t440*t4+(-t104+t442+t443+2432.D0/3.D0+t136)*t5+t449*
     #t3+(-t87-896+384*t82+t231-3520.D0/3.D0*s1l1-t453)*t12)*l1+(-t458+(
     #t105+t104+t106+t136-t459)*t6-64.D0/3.D0*(-14+t462-2*s1l1)*t9)*t19+
     #(-1024.D0/3.D0*t73+512.D0/3.D0*t76)*t79+t493*t16+((-128.D0/3.D0-t3
     #35+t144+t133-t99)*t3+(1984.D0/3.D0+t136-t105-t144-t301-t133)*t6-t3
     #87-64.D0/3.D0*(-4+10*s2l1-t500)*t9+(256-t87-t165-t474+t128-t504)*t
     #4+(t125-t507+t98+1600.D0/3.D0*s2l1-t85)*t5+(t221-1472.D0/3.D0+t307
     #+t236+t178+t478)*t12)*t320+128.D0/3.D0*t150
      t522 = 8*s2l1
      t532 = 1664.D0/3.D0*t86
      t537 = 896.D0/3.D0*t86
      t538 = 2944.D0/3.D0*s1l1
      t544 = 384*t13
      t555 = 3584.D0/3.D0*s2l1
      t567 = 1280.D0/3.D0*t86
      t570 = 2176.D0/3.D0*t84
      t573 = 1216*s1l1
      t580 = -t396*t3+(-t136+3328.D0/3.D0-t83-t224+t135-t99)*t6+(-2432.D
     #0/3.D0*t86-896-t436-3968.D0/3.D0*s2l1-t437-512*t84)*t4-640.D0/3.D0
     #-t537+(t567+704+t307+2368.D0/3.D0*s2l1+4736.D0/3.D0*s1l1+t570)*t12
     #+t159-t118-t485+(t87-2560.D0/3.D0+t83-t573+t89-t105)*t5-t139+64.D0
     #/3.D0*(-18-t462+6*s1l1)*t9-t133
      t583 = 128.D0/3.D0*s2l1
      t586 = 384*t84
      s1 = ((-512.D0/3.D0-t136+t105+448*s1l1-t173)*t6+t387+(-t87+256.D0/
     #3.D0-t159-t300-t285-t160)*t5+64.D0/3.D0*(-18+t522-t500)*t9)*t18-64
     #.D0/3.D0*t271*(-6-t462)*t5-t159+t133+((-t532-1408.D0/3.D0-t380-t23
     #8-t184-t481)*t5-t88+t173-t440*t4+(2624.D0/3.D0-t229+t537+t538+t453
     #)*t12-t99+64.D0/3.D0*(t374-4-t375)*t9+t544+(-1088.D0/3.D0+t136+t83
     #+t282+t447+t160)*t6-t449*t3)*t17+(-832*s2l1-t335-t238-1600.D0/3.D0
     #-t87)*t12+(2560.D0/3.D0+t532+t307+2624.D0/3.D0*s1l1+t99+t555)*t4
      t613 = s1+t580*t110+t118+((-t130-512.D0/3.D0-3968.D0/3.D0*s1l1-t10
     #4-t570+t583)*t12+(t130+128.D0/3.D0+t144+t133+t443+t586)*t4+1280.D0
     #/3.D0*t13+(-t87+2176.D0/3.D0-t238+t105+t573-t240)*t5+(-3328.D0/3.D
     #0-t301+t98+t100+t85)*t6+64.D0/3.D0*(t522+4-8*s1l1)*t9)*t141+t485+(
     #t401-t403+128.D0/3.D0*(-2*t196-6*t29+4*t207-t195-t406+4)*t147*t32)
     #*l5+t537+(1024.D0/3.D0*t326-512.D0/3.D0*t328)*t39
      t619 = 1088.D0/3.D0*t82
      t620 = 512*s2l1
      t625 = 128*s1l1
      t628 = 512.D0/3.D0*t9
      t629 = 832.D0/3.D0*t84
      t630 = 256*s2l1
      t634 = 896.D0/3.D0*s1l1
      t635 = 704.D0/3.D0*t84
      t638 = (t619+192-t97-t153+t620)*t4+(t217+t130+640.D0/3.D0+t230+t45
     #3+t158)*t3+(t169+t125-256+t625+t218+t154)*t6-t225+t14+t628-t170-64
     #0.D0/3.D0-t230-t629-t630+(t217-t87+128.D0/3.D0-t128-t134-t620)*t5+
     #(-1792.D0/3.D0*t82+t120-896.D0/3.D0+t634+t635-t154)*t12
      t661 = 832.D0/3.D0*t82
      t664 = (-t217+t567-1600.D0/3.D0+t475+t99+t158)*t4+(t217+320.D0/3.D
     #0+t583)*t6+(-t152-t537+1088.D0/3.D0-t634-t105-t183)*t12+(t225-832.
     #D0/3.D0+t133+t459+t120)*t5+(t661-t120+1856.D0/3.D0+t97+t338)*t3-t1
     #52-t130-1088.D0/3.D0-t128-t165-t338
      t666 = 1088.D0/3.D0*t13
      s1 = t638*t110+((-64+t120+t217-t106+t218)*t6-64.D0/3.D0*(11-t462)*
     #t9+t458)*t19+(t144+t630+128+t120+t459)*t12+(-t123-t120-320.D0/3.D0
     #-t459-t334)*t4+(-t136-t103-t105)*t3+t664*l1
      s2 = s1+(t666+(t339+t87-512.D0/3.D0+t103+t105+t620)*t6+(-t217+t120
     #+192+t97+t134+t630)*t12+(-t661-t136-896.D0/3.D0-t335+t218-t481)*t5
     #-64.D0/3.D0*(4-4*s2l1)*t9)*t111+(64.D0/3.D0*(5-t522)*t9-t666+(t144
     #-64-t105+t158)*t5+(-t98-t125+320-t105)*t6)*t18
      t725 = s2+(64.D0/3.D0*(3-t462)*t9+(-t229-t87+320-t625-t105-t183)*t
     #6-t544+(t217+t87+t128+t134+1792.D0/3.D0*s2l1)*t5+(-t229-192+t85-t5
     #04)*t4+(t436-t136+192-t634-t160+t442)*t12)*t141+(-64.D0/3.D0-64.D0
     #/3.D0*l1)*t58+((-t217-t221-192-t475-t153-t630)*t3-64.D0/3.D0*(-3-t
     #522)*t9+(t619+t379+768+2432.D0/3.D0*s1l1+t629+3200.D0/3.D0*s2l1)*t
     #4+(-576*t82-t130-960-t237-t134-t555)*t12+(t217-t120+1856.D0/3.D0-t
     #459+t105+t583)*t5+t666+(t229+t87-1984.D0/3.D0+t459+t105+t158)*t6)*
     #t320+320.D0/3.D0*t249*t251
      t727 = 128.D0/3.D0*t16
      t728 = 320.D0/3.D0*l1
      t771 = t144-256*t13+448.D0/3.D0-t105+t158+64.D0/3.D0*(-1-t462)*t9+
     #(t217-t379+1280.D0/3.D0-t507-t635-t334)*t4+(t83+t130-1600.D0/3.D0+
     #t97+t334)*t12+(-t307-t87+64-t103-t105-t154)*t6+(t83+t120+576+640.D
     #0/3.D0*s1l1-t218+t395)*t5+(-t229+t221-1216.D0/3.D0+t335+t160-t630)
     #*t3
      t773 = ((320.D0/3.D0*t110+t727+t728)*s1l1+(128.D0/3.D0+t431-t727+6
     #4*l1)*s2l1+128.D0/3.D0+64*t110+t728)*t298+t216+t80+t260+t331-t151+
     #t136+t103+t134-128.D0/3.D0*t271**2*t5+(t219-t97-t119+t158-t628+(-t
     #229-t379-768-t538-t586-896*s2l1)*t4+(-t219+t120-1664.D0/3.D0+t459-
     #t334)*t5+(1664.D0/3.D0*t82+t130+960+t230+t85+1152*s2l1)*t12+(-t120
     #+1024.D0/3.D0-t159-t158-t459)*t6+(-t98+t221+192+t282+t586)*t3)*t16
     #+t771*t17
      t778 = (l1+z-t29)**2
      taux19 = ((t248+t352)*(-t354*t356/2-(t359-t360)*t355+t363-t364-t36
     #5*t367/2+t370*t366)+(t516+t613)*(t354*t355+t359-t360+t365*t366)-(t
     #725+t773)*t354)/t778*t32
      t1 = 6*s2l1
      t3 = 1-z
      t5 = 128.D0/3.D0*(-4-t1)*t3
      t6 = 4*s1l1
      t7 = 2*s2l1
      t9 = t3**2
      t10 = t9**2
      t11 = t10*t3
      t14 = 14*s1l1
      t15 = 4*s2l1
      t17 = t10*t9
      t20 = 14*s2l1
      t24 = t9*t3
      t25 = t10*t24
      t26 = 512*t25
      t27 = 256.D0/3.D0*s1l1
      t28 = 256.D0/3.D0*s2l1
      t29 = 16*s1l1
      t32 = 128.D0/3.D0*(8+t29-t7)*t9
      t38 = l1**2
      t39 = t38**2
      t47 = 8*s2l1
      t51 = 34*s1l1
      t55 = 8*s1l1
      t60 = t39*t38
      t62 = 2*s1l1
      t66 = 1792.D0/3.D0*t25
      t67 = 10*s1l1
      t75 = t39**2
      t79 = t75*t38
      t81 = 5*s1l1
      t99 = 6*s1l1
      t102 = 128.D0/3.D0*(8+t99+t47)*t3
      t104 = t38*l1
      t106 = 12*s2l1
      t115 = 128.D0/3.D0*(-12-t14-t7)*t9
      t119 = 18*s1l1
      t120 = 10*s2l1
      t126 = 128.D0/3.D0*(-28-t67-t20)*t24
      t127 = -s1l1+s2l1
      t128 = 2*t127
      t130 = 128.D0/3.D0*t128*t3
      t132 = t39*l1
      t147 = t39*t104
      t156 = t75*l1
      t179 = 24*s1l1
      t185 = -1-s2l1
      s1 = 256.D0/3.D0+(t5-128.D0/3.D0*(20-t6+t7)*t11+128.D0/3.D0*(-t14-
     #2+t15)*t17-128.D0/3.D0*(8-t6+t20)*t10+t26+256.D0/3.D0+t27+t28-t32+
     #128.D0/3.D0*(28+20*s1l1+t20)*t24)*t39+(128.D0/3.D0*(-6-t6-t15)*t9-
     #256.D0/3.D0*(-5+t6-t1)*t24-128.D0/3.D0*(42+t14+t47)*t11-128.D0/3.D
     #0*(-20-t51+t1)*t10+t26+128.D0/3.D0*(6-t55+t1)*t17)*t60+(-128.D0/3.
     #D0*(-6-t62-t7)*t10+t66+128.D0/3.D0*(-12-t67+t15)*t17-128.D0/3.D0*(
     #-t55+t1+8)*t11)*t75+t27+t28+(256.D0/3.D0*t25-256.D0/3.D0*t17)*t79+
     #(256.D0/3.D0*(-1+t81-t7)*t24-256.D0/3.D0*(t7+7-3*s1l1)*t17-128.D0/
     #3.D0*(6+t62+t47)*t9-128.D0/3.D0*(-t1-t29-42)*t11-512.D0/3.D0*t25-2
     #56.D0/3.D0-t27-t28+128.D0/3.D0*(-22-t51+t15)*t10+t102)*t104
      t188 = s1+128.D0/3.D0*(t99+t106+12)*t9+(-t66+128.D0/3.D0*(t99+t1+2
     #4)*t11-t115-128.D0/3.D0*(-2-t67+t1)*t17+128.D0/3.D0*(4-t119+t120)*
     #t10+t126+t130)*t132+(128.D0/3.D0*(-18-t14+t7)*t10+128.D0/3.D0*(28-
     #t6+t1)*t11-512.D0/3.D0*t127*t24-128.D0/3.D0*(t15-4-t14)*t17-t66)*t
     #147+(128.D0/3.D0*t128*t11-128.D0/3.D0*(-t62-8+t7)*t17-1024.D0/3.D0
     #*t25)*t156+(256.D0/3.D0*(-7*s1l1-8)*t11-256.D0/3.D0*(-s2l1-5)*t17-
     #t126+128.D0/3.D0*(t29-t106-10)*t10-t130+t115)*t38+256.D0/3.D0*(-s1
     #l1-4-t15)*t24+(512.D0/3.D0*(-2-s2l1)*t11-t5+t32-256.D0/3.D0*(-13-t
     #1-t81)*t10-256.D0/3.D0-t27-t28-128.D0/3.D0*(t179+t120+28)*t24)*l1-
     #t102-256.D0/3.D0*t185*t10
      t189 = zqe0(l1,z,l3)
      t190 = 1-l1
      t191 = dlog(t190)
      t192 = t191**2
      t195 = sqe1(l1)
      t196 = zqe1(l1,z,l3)
      t199 = sqe2(l1)
      t200 = zqe2(l1,z,l3)
      t201 = sbqe0(l1)
      t202 = dlog(z)
      t203 = t202**2
      t206 = sbqe1(l1)
      t211 = 9*s2l1
      t215 = 12*s1l1
      t216 = 42*s2l1
      t220 = 174*s1l1
      t228 = 1280.D0/3.D0*s1l1
      t229 = 1280.D0/3.D0*s2l1
      t230 = 30*s1l1
      t231 = 40*s2l1
      t234 = 128.D0/3.D0*(-40-t230-t231)*t3
      t236 = 82*s1l1
      t237 = 30*s2l1
      t243 = 22*s2l1
      t246 = 128.D0/3.D0*(12-t99+t243)*t3
      t247 = 5*s2l1
      t251 = 1024.D0/3.D0*s1l1
      t252 = 1024.D0/3.D0*s2l1
      t254 = 31*s2l1
      t259 = 58*s2l1
      t266 = 128.D0/3.D0*(-52-86*s1l1-t7)*t9
      t269 = 18*s2l1
      t274 = 72*s1l1
      t278 = 46*s1l1
      t284 = 26*s2l1
      t293 = 74*s1l1
      t300 = -t127
      t307 = 3072*t25
      t318 = 128.D0/3.D0*(56+68*s1l1+t1)*t9
      t321 = 128.D0/3.D0*(128+t278+t259)*t24
      t329 = 10*t300
      t331 = 128.D0/3.D0*t329*t3
      t341 = 16*s2l1
      t369 = 34*s2l1
      t374 = 94*s1l1
      t379 = 56*s1l1
      s2 = (-256.D0/3.D0*(15*s1l1-t211-34)*t17-128.D0/3.D0*(-32-t215-t21
     #6)*t9+128.D0/3.D0*(114+t220-t106)*t10+256.D0/3.D0*(2-27*s1l1+t1)*t
     #24+1280.D0/3.D0+t228+t229+t234+2560.D0/3.D0*t25-128.D0/3.D0*(208+t
     #236+t237)*t11)*t104+(-t246+512.D0/3.D0*(10+t247)*t11+1024.D0/3.D0+
     #t251+t252-256.D0/3.D0*(25*s1l1+t254+66)*t10-128.D0/3.D0*(-148-122*
     #s1l1-t259)*t24+t266)*l1+(128.D0/3.D0*(28+t119+t269)*t9-8192.D0/3.D
     #0*t25-128.D0/3.D0*(-t231-208-t274)*t11+128.D0/3.D0*(t278-24-24*s2l
     #1)*t17-256.D0/3.D0*(22-22*s1l1+t284)*t24-128.D0/3.D0*(104+t220-t24
     #3)*t10)*t60+(128.D0/3.D0*(90+t293-t120)*t10+128.D0/3.D0*(-140-t284
     #+t29)*t11-2560.D0/3.D0*t300*t24-128.D0/3.D0*(22+t293-t20)*t17+t307
     #)*t147
      s1 = s2+(-128.D0/3.D0*(-t62-t7-12)*t17-t26)*t79+(t318-t321+256.D0/
     #3.D0*(42+t7+35*s1l1)*t11-256.D0/3.D0*(t247+25)*t17-t331+128.D0/3.D
     #0*(38-t236+48*s2l1)*t10)*t38+(-128.D0/3.D0*(-t1-100+t29)*t11+128.D
     #0/3.D0*(10+t274-t341)*t17-t266+t246-1024.D0/3.D0-t251-t252+128.D0/
     #3.D0*(-148-102*s1l1-78*s2l1)*t24-2560*t25-128.D0/3.D0*(-42+t179-72
     #*s2l1)*t10)*t39+256.D0/3.D0*(t81+20+20*s2l1)*t24
      t405 = s1-1280.D0/3.D0+(128.D0/3.D0*t329*t11-128.D0/3.D0*(44+t29-t
     #1)*t17+5632.D0/3.D0*t25)*t156-t234+(128.D0/3.D0*(-128-t369-32*s1l1
     #)*t11+t331+t307+128.D0/3.D0*(-6+t374-36*s2l1)*t10+t321-t318-128.D0
     #/3.D0*(-t284+10+t379)*t17)*t132-t228-t229+(128.D0/3.D0*(60+t379-t3
     #41)*t17-128.D0/3.D0*(32+t215+t106)*t10-t307-128.D0/3.D0*(-44+38*s1
     #l1-t237)*t11)*t75+1280.D0/3.D0*t185*t10+128.D0/3.D0*(-60*s2l1-t230
     #-60)*t9
      t410 = -t190
      t415 = 2*l6-1
      t416 = t415**2
      t419 = l5-1
      t421 = z*l1
      t423 = (1-l1+t421)**2
      t424 = 1/t423
      t429 = dsqrt(z)
      t431 = dsqrt(l5)
      t433 = dsqrt(-t419)
      t436 = dsqrt(l1)
      t440 = l5*z*l1-l1+1-l5+l5*l1
      t442 = dsqrt(t190)
      t447 = t440**2
      t450 = dsqrt(4*l1*t416*z*l5*t419*t410*t424-4*t415*t429*t431*t433*t
     #436*t440*t424*t442+t447*t424)
      t455 = 11*s1l1
      t461 = 25*s2l1
      t475 = 15*s2l1
      t476 = 107*s1l1
      t481 = 28*s1l1
      t485 = 15616.D0/3.D0*t25
      t502 = 131*s1l1
      t508 = 17*s2l1
      t516 = 19*s2l1
      t520 = 40*s1l1
      t524 = 252*s1l1
      t538 = 128.D0/3.D0*(1-t455+t516)*t3
      t546 = 128.D0/3.D0*(-55-66*s1l1+t508)*t9
      t552 = 128.D0/3.D0*(-126-t481-t216)*t24
      t558 = 128.D0/3.D0*(30+52*s1l1+t7)*t3
      t568 = 128*s1l1
      t569 = 128*s2l1
      t579 = 128.D0/3.D0*(148+169*s1l1+53*s2l1)*t9
      t585 = 48*s1l1
      t598 = 128.D0/3.D0*(73+51*s1l1+73*s2l1)*t3
      t599 = 2816.D0/3.D0*s1l1
      t600 = 3200.D0/3.D0*s2l1
      t602 = 136*s1l1
      s2 = 128*t410*(1+l1+t38)*t450+(3328.D0/3.D0*t25-128.D0/3.D0*(26+11
     #*s2l1+t455)*t17)*t79+(128.D0/3.D0*(7-t81+t461)*t11-128.D0/3.D0*(-t
     #247-t278-77)*t17-3584*t25)*t156+(-128.D0/3.D0*(-52-19*s1l1-t461)*t
     #10+128.D0/3.D0*(-80+t475-t476)*t17-128.D0/3.D0*(70*s2l1+94-t481)*t
     #11+t485)*t75
      s3 = s2+(128.D0/3.D0*(-107-90*s1l1+37*s2l1)*t10+128.D0/3.D0*(209+s
     #1l1+41*s2l1)*t11-t485-128.D0/3.D0*(11-13*s1l1+47*s2l1)*t24-128.D0/
     #3.D0*(t1-t502-31)*t17)*t147+(128.D0/3.D0*(-26-t81-t508)*t9-128.D0/
     #3.D0*(311+t502+69*s2l1)*t11+128.D0/3.D0*(t516-t374+27)*t17-256.D0/
     #3.D0*(-19-t369+t520)*t24-128.D0/3.D0*(-156-t524+t369)*t10+14848.D0
     #/3.D0*t25)*t60
      s1 = s3+(-t485+128.D0/3.D0*(86*s2l1+238+84*s1l1)*t11+t538+128.D0/3
     #.D0*(-46-155*s1l1+t475)*t10-t546-128.D0/3.D0*(t254-t476)*t17+t552)
     #*t132+(t558+128.D0/3.D0*(t508-5-120*s1l1)*t17+128.D0/3.D0*(281+187
     #*s1l1+167*s2l1)*t24-t568+t569-128.D0/3.D0*(109-t67+141*s2l1)*t10+4
     #096*t25-t579-128.D0/3.D0*(145+s1l1+t211)*t11)*t39+(256.D0/3.D0*(13
     #-t7+t585)*t24-128.D0/3.D0*(73+43*s1l1+82*s2l1)*t9+t598-3200.D0/3.D
     #0-t599-t600-128.D0/3.D0*(-304-44*s2l1-t602)*t11-4096.D0/3.D0*t25+1
     #28.D0/3.D0*(-172-t524+t269)*t10-256.D0/3.D0*(-t179+t106+52)*t17)*t
     #104
      t651 = s1+3200.D0/3.D0+(256.D0/3.D0*(-72-t47-t379)*t11-t538-t552+1
     #28.D0/3.D0*(-t231-6+t602)*t10-256.D0/3.D0*(-t47-40)*t17+t546)*t38+
     #(512.D0/3.D0*(-16-t47)*t11-256.D0/3.D0*(-t520-52*s2l1-108)*t10-128
     #.D0/3.D0*(270+120*s2l1+200*s1l1)*t24-t558+t579+t568-t569)*l1+128.D
     #0/3.D0*(99+99*s2l1+t585)*t9+256.D0/3.D0*(-t55-32-32*s2l1)*t24-2048
     #.D0/3.D0*t185*t10-t598+t599+t600
      t655 = (l1+z-t421)**2
      taux20 = (t188*(-t189*t192/2-(t195-t196)*t191+t199-t200-t201*t203/
     #2+t206*t202)+t405*(t189*t191+t195-t196+t201*t202)-t651*t189)/t655*
     #t424

CCC Coefficients

      tcoff37 = taux19
      tcoff38 = taux19
      tcoff39 = taux20
      tcoff40 = taux20

CCC Jet functions

      tFjet37 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff38)
      tFjet39 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff40)

CCC Sum

      qgRVnum10e2 = tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tFje
     #t40*tcoff40
      RETURN
      END

      double precision function qgRVnum10e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = s1l1*s2l1
      t2 = 256.D0/3.D0*t1
      t3 = -1+l1
      t4 = z**2
      t5 = t4**2
      t6 = t3*t5
      t7 = z*l1
      t9 = (1-l1+t7)**2
      t10 = 1/t9
      t11 = l6**2
      t12 = t10*t11
      t14 = 1024.D0/3.D0*t6*t12
      t16 = 2*t7
      t17 = l1**2
      t18 = 2*t17
      t19 = t4*t17
      t21 = t4*z
      t22 = t21*t17
      t24 = t17*l1
      t25 = z*t24
      t27 = t4*t24
      t29 = t21*t24
      t31 = t5*t24
      t32 = z*t17
      t34 = t4*l1
      t36 = -4*z-t16+l1-t18-14*t19+6*t22-4*t25+6*t27-4*t29+t31+10*t32+9*
     #t34+t24
      t40 = t10*l6
      t42 = 1024.D0/3.D0*t6*t40
      t44 = l5**2
      t45 = (-t14-64.D0/3.D0*t36*t21*t10+t42)*t44
      t46 = dsqrt(z)
      t51 = dsqrt(l1)
      t52 = t46*t21*(t34+3*z-t16+l1-1)*t51
      t53 = 1-l5
      t54 = dsqrt(t53)
      t55 = -t3
      t56 = dsqrt(t55)
      t57 = t54*t56
      t59 = t52*t57*t40
      t61 = t57*t10
      t62 = t52*t61
      t65 = dsqrt(l5)
      t66 = t65*l5
      t67 = (512.D0/3.D0*t59-256.D0/3.D0*t62)*t66
      t68 = s1l1**2
      t69 = 512.D0/3.D0*t68
      t70 = 3*s2l1
      t73 = 1-z
      t74 = t73**2
      t75 = t74**2
      t76 = t75*t74
      t79 = 576*s1l1
      t80 = s2l1**2
      t81 = 512.D0/3.D0*t80
      t82 = 448.D0/3.D0*t68
      t83 = 832.D0/3.D0*s2l1
      t84 = 128.D0/3.D0*t1
      t87 = 256.D0/3.D0*s1l1
      t88 = 320.D0/3.D0*t80
      t89 = 320.D0/3.D0*t68
      t90 = 256.D0/3.D0*s2l1
      t91 = 128*t1
      t94 = 1280.D0/3.D0*s1l1
      t95 = 320*t68
      t96 = 256*t1
      t98 = t74*t73
      t100 = 128*s2l1
      t101 = 64*t68
      t102 = 448.D0/3.D0*s1l1
      t104 = t75*t73
      t106 = t75*t98
      t107 = 128*t106
      t109 = t17**2
      t110 = t109*t17
      t112 = 64.D0/3.D0*t80
      t113 = 64.D0/3.D0*t68
      t114 = 2*s2l1
      t115 = 7*s1l1
      t119 = 1408.D0/3.D0*s1l1
      t120 = 896.D0/3.D0*t80
      t121 = 1024.D0/3.D0*t68
      t122 = 1472.D0/3.D0*s2l1
      t123 = t119-t120+t121-t122+t84-64.D0/3.D0
      t125 = 192*t80
      t126 = 128*t68
      t127 = 512.D0/3.D0*t1
      t128 = 512.D0/3.D0*s1l1
      t129 = 704.D0/3.D0*s2l1
      t132 = 448.D0/3.D0*t80
      t133 = 256.D0/3.D0*t68
      t135 = 640.D0/3.D0*t1
      t138 = 64.D0/3.D0*s1l1
      t139 = 256*t80
      t140 = 512*s1l1
      t141 = 192*s2l1
      t144 = 64.D0/3.D0*s2l1
      t145 = 64*t80
      t146 = 192*t68
      t147 = t129+t145-t146-t87+256.D0/3.D0
      t149 = t112+t113-64.D0/3.D0+64.D0/3.D0*(1-t114+t115)*t76+t123*t74-
     #t107+(320.D0/3.D0-t125-t126-t127-t128-t129)*t104+(192+t128+t132+t1
     #33+1088.D0/3.D0*s2l1+t135)*t75-t138+(t139-t126-704.D0/3.D0-t140+t1
     #41-t2)*t98-t144+t147*t73
      t151 = 640.D0/3.D0*t80
      t152 = 320*t80
      t153 = 320.D0/3.D0*s2l1
      t156 = 192*s1l1
      t157 = 256.D0/3.D0*t80
      t163 = 256*s1l1
      t164 = 832.D0/3.D0*s1l1
      t165 = 128*t80
      t166 = 320*s2l1
      t169 = 640.D0/3.D0*s2l1
      t173 = 704.D0/3.D0*t80
      t176 = 1280.D0/3.D0*s2l1
      t177 = t96+t151+640.D0/3.D0+t101+(-832.D0/3.D0-t102+t152+t69+t153+
     #t127)*t98+(-t100+128+t156-t157)*t104-64.D0/3.D0*(s2l1+5)*t76+t163+
     #(-640.D0/3.D0-t164-t165-t89-t166-t127)*t73+(448.D0/3.D0+t169-t69-t
     #163-t2)*t75+(192+1088.D0/3.D0*s1l1-t173+t133-t153-t84)*t74+t176
      t179 = 128.D0/3.D0*t24
      t180 = 64*t17
      t183 = 320.D0/3.D0*t17
      t188 = 2*l6-1
      t189 = t188**2
      t204 = l5*z*l1-l1+1-l5+l5*l1
      t206 = t10*t56
      t210 = t204**2
      t212 = -4*l1*t189*z*l5*t53*t3*t10-4*t188*t46*t65*t54*t51*t204*t206
     #+t210*t10
      t213 = dsqrt(t212)
      t215 = 320.D0/3.D0*s1l1
      t218 = 5*s1l1
      t222 = 1024.D0/3.D0*t80
      t223 = 448.D0/3.D0*s2l1
      t228 = 448.D0/3.D0*t106
      t235 = t109*l1
      t239 = 704.D0/3.D0*s1l1
      t246 = 1792.D0/3.D0*s1l1
      t251 = 64*s1l1
      t259 = 128.D0/3.D0*t80
      t260 = 128.D0/3.D0*t68
      t264 = t109*t24
      t272 = 256.D0/3.D0*t106
      t274 = t109**2
      t275 = t274*l1
      s1 = -64-t2+t45+t67-t69+(64.D0/3.D0*(-3-t70+4*s1l1)*t76+(-192-t79+
     #t81-t82+t83-t84)*t75+(-t87-t88-t89-t90-t91)*t74+(t94+t88+t95+t96)*
     #t98+(384-t100-t81-t101+t102-t2)*t104-t107)*t110
      s2 = s1+t149*t109+t177*t17+((128.D0/3.D0-t179-t180)*s1l1-t183+(-t1
     #83-128.D0/3.D0+t179)*s2l1)*t213
      t277 = s2+((-t153-t157+t133+t215)*t73-64.D0/3.D0*(1+t218-t70)*t76+
     #(256+t102-t222-t69-t223-t127)*t98+(t132+t101-512.D0/3.D0-t138+t223
     #+t2)*t104+t228+(t222-t133-256.D0/3.D0-320*s1l1+t166+t2)*t74+(-128-
     #t145+t89+t156-t83)*t75)*t235+(-t112-t113+t144+(320.D0/3.D0+t112+t9
     #0)*t104+64.D0/3.D0+(-t239-896.D0/3.D0+t145-t84-t90)*t75+t138-t147*
     #t73-t123*t74+(704.D0/3.D0+t246-t120+t69-t83+t2)*t98)*l1+(t228+(320
     #.D0/3.D0-t151-t133-t127+t251-t83)*t75-64.D0/3.D0*(t115-t114+2)*t76
     #+(-640.D0/3.D0+t128+t81+t126+t223+t127)*t104+(-t87+t259-t260+t90)*
     #t98)*t264+((64.D0/3.D0*s1l1-64.D0/3.D0*s2l1)*t104-64.D0/3.D0*(4+s1
     #l1-s2l1)*t76+t272)*t275
      t278 = 384*s2l1
      t281 = 1408.D0/3.D0*s2l1
      t282 = 1472.D0/3.D0*s1l1
      t285 = t17*t212
      t294 = (t14+128.D0/3.D0*(t19+3*t7-2*t32+t17+l1-2)*t5*t10-t42)*l5
      t296 = t46*t5*t51
      t299 = t296*t54*t206*l6
      t301 = t296*t61
      t304 = (-512.D0/3.D0*t299+256.D0/3.D0*t301)*t65
      t305 = -t106+t76
      t307 = t274*t17
      t311 = t5*z*l1*t10
      t312 = 64.D0/3.D0*t311
      t318 = 64*s2l1
      t323 = 1+s2l1
      t328 = 640.D0/3.D0*t68
      t329 = 512.D0/3.D0*s2l1
      t330 = 256+t96-t157+t328+t329+t79
      t332 = 128.D0/3.D0*t106
      t337 = 1664.D0/3.D0*s2l1
      t345 = 1024.D0/3.D0*t1
      t350 = 64+t2-t165+t69+t332-t318+64.D0/3.D0*(t114+7-3*s1l1)*t76+(38
     #4+t79+t125+t146+t337+384*t1)*t74-t330*t73+t164+(256+t79-t165+t82-t
     #329+t84)*t75+(-256-t79-t173-t95-t278-t345)*t98+(-384-t128+t81+t101
     #+t223+t2)*t104
      t352 = (t278+t165+t2+256+t102)*t98+(-384-t281-t157-t133-t282-t96)*
     #t74+64.D0/3.D0*t285+t294+t304+64.D0/3.D0*t305*t307-t312+t165+(64.D
     #0/3.D0*(6+t218-t114)*t76+(128.D0/3.D0-t144-t145-t101-t128-t2)*t104
     #-t228+(-64.D0/3.D0+t251+t145+t101+t318+t2)*t75)*t274+t318-t164-64.
     #D0/3.D0*t323*(t114+3)*t75+t330*t73+t350*t24
      t354 = sqe1(l1)
      t355 = zqe1(l1,z,l3)
      t356 = t354-t355
      t357 = dlog(t55)
      t358 = t357**2
      t361 = sqe2(l1)
      t362 = zqe2(l1,z,l3)
      t365 = sbqe1(l1)
      t366 = dlog(z)
      t367 = t366**2
      t370 = sbqe2(l1)
      t374 = 896.D0/3.D0*t1
      t377 = 896.D0/3.D0*s2l1
      t378 = 384*s1l1
      t381 = 1408.D0/3.D0*t80
      t382 = 3776.D0/3.D0*s1l1
      t385 = t2+384-t381+t382-1984.D0/3.D0*s2l1+2048.D0/3.D0*t68
      t388 = 1216.D0/3.D0*s2l1
      t390 = -128.D0/3.D0-1600.D0/3.D0*s1l1+t81+t388-1664.D0/3.D0*t68
      t394 = 256*t68
      t402 = 2560.D0/3.D0*s2l1
      t403 = 1024.D0/3.D0*s1l1
      t407 = 1280.D0/3.D0*t80
      t408 = 640*t68
      t413 = 2176.D0/3.D0*s2l1
      t416 = 896.D0/3.D0*t68
      t420 = -64.D0/3.D0*(-14-t114)*t76-t402+(t165-1472.D0/3.D0-t403+t14
     #1)*t104-512*t1-t126-t407+(-t135+1600.D0/3.D0-t120+t318-t164-t408)*
     #t98+(t2-1408.D0/3.D0+t69+t151-t140-t90)*t74-1280.D0/3.D0-t140+(t34
     #5+1408.D0/3.D0+t120+t246+t413+t69)*t73+(t84-256.D0/3.D0+t416-t157-
     #t176+2368.D0/3.D0*s1l1)*t75
      t430 = 1408.D0/3.D0*t106
      t432 = 10*s1l1
      t436 = 1024.D0/3.D0*s2l1
      t439 = 1664.D0/3.D0*s1l1
      t447 = 4*s2l1
      t456 = 1216*s1l1
      t467 = 640+640*t1+t157+1088*s1l1+t413+t121
      t474 = 1280.D0/3.D0*t1
      t477 = 2176.D0/3.D0*t68
      t480 = t157-640.D0/3.D0+(t127-2560.D0/3.D0+t151-t456+t83-t260)*t75
     #-t100-t79-t107+64.D0/3.D0*(-18-t447+6*s1l1)*t76+(-t2+3328.D0/3.D0-
     #t151-t129+t102-t126)*t104+t467*t73-t374+(-2432.D0/3.D0*t1-896-t381
     #-3968.D0/3.D0*s2l1-t382-512*t68)*t74-t416+(t474+704+t222+2368.D0/3
     #.D0*s2l1+4736.D0/3.D0*s1l1+t477)*t98
      t483 = 2*t55*t5
      t485 = 1024.D0/3.D0*t483*t12
      t487 = 1024.D0/3.D0*t483*t40
      t490 = 2*l1
      t501 = 8*z+4*t7-t490+4*t17+28*t19-12*t22+8*t25-12*t27+8*t29-2*t31-
     #20*t32-18*t34-2*t24
      t507 = 896.D0/3.D0*t106
      t508 = 128.D0/3.D0*s1l1
      t518 = 640.D0/3.D0+t374+(t251+(-t329-t259-256)*t104+(-t259+t377+t3
     #78+2432.D0/3.D0+t2)*t75+t126+t385*t74-t318+t390*t73+(-t127-896+384
     #*t80+t141-3520.D0/3.D0*s1l1-t394)*t98)*l1+t420*t17+t416+(1024.D0/3
     #.D0*t299-512.D0/3.D0*t301)*t65+((-128.D0/3.D0-t163+t165+t100-t126)
     #*t73+(1984.D0/3.D0+t2-t260-t165-t215-t100)*t104-t430-64.D0/3.D0*(-
     #4+10*s2l1-t432)*t76+(256-t127-t69-t407+t94-t436)*t74+(t91-t439+t81
     #+1600.D0/3.D0*s2l1-t133)*t75+(t135-1472.D0/3.D0+t222+t144+t164+t40
     #8)*t98)*t235-64.D0/3.D0*t323*(-6-t447)*t75+(-1024.D0/3.D0*t59+512.
     #D0/3.D0*t62)*t66+t480*t24+(-t485+t487-64.D0/3.D0*t501*t21*t10)*t44
     #+(-t507+(t260+t259+t90+t2-t508)*t104-64.D0/3.D0*(-14+t447-2*s1l1)*
     #t76)*t275-128.D0/3.D0*t285
      t522 = 1664.D0/3.D0*t1
      t524 = 3584.D0/3.D0*s2l1
      t533 = 384*t106
      t534 = 6*s2l1
      t535 = 14*s1l1
      t539 = 2944.D0/3.D0*s1l1
      t542 = 1408.D0/3.D0*t68
      t553 = 1408.D0/3.D0*t1
      t571 = -t305
      t578 = 640.D0/3.D0*t17
      t581 = 64.D0/3.D0*t24
      t589 = 8*s2l1
      t596 = 128.D0/3.D0*s2l1
      t599 = 384*t68
      s1 = (-832*s2l1-t163-t120-1600.D0/3.D0-t127)*t98+(2560.D0/3.D0+t52
     #2+t222+2624.D0/3.D0*s1l1+t126+t524)*t74-t467*t73+((-1088.D0/3.D0+t
     #2+t151+t246+t388+t328)*t104-t390*t73+t533+64.D0/3.D0*(t534-4-t535)
     #*t76-t251+(2624.D0/3.D0-t139+t374+t539+t394)*t98-t126+(-t522-1408.
     #D0/3.D0-t542-t120-t282-t402)*t75-t385*t74+t318)*t109+128.D0/3.D0*t
     #311+(-64.D0/3.D0*(t534-8-t535)*t76+(-1024.D0/3.D0+t553+t222+t215+t
     #337+t542)*t75+(-t2+1856.D0/3.D0-t81-t246-t129-t328)*t104+(-t91+64.
     #D0/3.D0+t156-t165-t141)*t98-t430)*t264+(t485-t487+128.D0/3.D0*(-2*
     #t19-6*t7+4*t32-t18-t490+4)*t5*t10)*l5
      t613 = s1+256.D0/3.D0*t571*t307+((256.D0/3.D0*t24-256.D0/3.D0+128*
     #t17)*s1l1+t578+(-t179+t578+128.D0/3.D0)*s2l1-64.D0/3.D0+t581)*t213
     #-t157+t100+((-512.D0/3.D0-t2+t260+448*s1l1-t318)*t104+t430+(-t127+
     #256.D0/3.D0-t157-t153-t239-t328)*t75+64.D0/3.D0*(-18+t589-t432)*t7
     #6)*t274+t79+((-t96-512.D0/3.D0-3968.D0/3.D0*s1l1-t259-t477+t596)*t
     #98+(t96+128.D0/3.D0+t165+t100+t378+t599)*t74+1280.D0/3.D0*t106+(-t
     #127+2176.D0/3.D0-t120+t260+t456-t122)*t75+(-3328.D0/3.D0-t215+t81+
     #t223+t133)*t104+64.D0/3.D0*(t589+4-8*s1l1)*t76)*t110
      t615 = zqe0(l1,z,l3)
      t619 = sbqe0(l1)
      t631 = 512.D0/3.D0*t76
      t651 = 1088.D0/3.D0*t106
      t652 = 256*s2l1
      t663 = 832.D0/3.D0*t68
      t664 = 1088.D0/3.D0*t80
      t678 = (-t87-t260-t2)*t73+t294-128.D0/3.D0*t323**2*t75+(-t631+(166
     #4.D0/3.D0*t80+t96+960+t140+t133+1152*s2l1)*t98-t82+(-t84-t508+1024
     #.D0/3.D0-t157-t329)*t104+(-t599-t139-t553-896*s2l1-768-t539)*t74+(
     #-t132-t169-1664.D0/3.D0+t84+t508)*t75+t132+(-t81+t599+192+t135+t24
     #6)*t73+t329-t128)*t17+320.D0/3.D0*t305*t307+t87+t67+((t112-t84+t59
     #6+t260+1856.D0/3.D0-t508)*t75+t651+(-192-t112-t135-t652-t403-t146)
     #*t73+(-576*t80-t96-t524-t101-t119-960)*t98+(t260+t139+t127+t329-19
     #84.D0/3.D0+t508)*t104-64.D0/3.D0*(-3-t589)*t76+(t663+t664+t553+320
     #0.D0/3.D0*s2l1+768+2432.D0/3.D0*s1l1)*t74)*t235+(-64.D0/3.D0*l1-64
     #.D0/3.D0)*t212+(128+t165+t508+t84+t652)*t98-t312+(-t88-320.D0/3.D0
     #-t84-t169-t508)*t74
      t690 = 704.D0/3.D0*t68
      t698 = 640.D0/3.D0*s1l1
      t704 = (t169+t151+t96-1600.D0/3.D0+t128)*t98-t260+(-t439-t553+t112
     #+1280.D0/3.D0-t169-t690)*t74+t329+64.D0/3.D0*(-1-t447)*t76+448.D0/
     #3.D0+t165+(-t139+t135-t652-1216.D0/3.D0+t163+t328)*t73+(-t113+t151
     #+t84+t413+576+t698)*t75-256*t106+(-t127+64-t87-t337-t222-t260)*t10
     #4
      t707 = 128.D0/3.D0*t17
      t708 = 320.D0/3.D0*l1
      t711 = 64*l1
      t714 = 64*t24
      t719 = 128*s1l1
      t723 = 896.D0/3.D0*s1l1
      t726 = 512*s2l1
      t731 = -640.D0/3.D0+t272+(t112+640.D0/3.D0+t96+t329+t140+t394)*t73
     #-t663+(t173+t91+t337-256+t719+t113)*t104-t145-t140+(-1792.D0/3.D0*
     #t80+t84-t337+t690+t723-896.D0/3.D0)*t98+(t664-t128+t726-t146+192)*
     #t74+t631+(t112-t101+128.D0/3.D0-t127-t94-t726)*t75-t345-t652
      t740 = 832.D0/3.D0*t80
      t777 = -t176-t96+(t145-832.D0/3.D0+t100+t84+t508)*t75+(t112+320.D0
     #/3.D0+t596)*t104-t94+(-t112+t126-1600.D0/3.D0+t474+t403+t329)*t74+
     #(-t125-t374-t281+1088.D0/3.D0-t723-t260)*t98-1088.D0/3.D0+(t176+t7
     #40+1856.D0/3.D0-t84+t128)*t73-t69-t125
      t779 = ((-t81+320-t91-t260)*t104-t651+(t165+t329-t260-64)*t75+64.D
     #0/3.D0*(5-t589)*t76)*t274+t45+t704*t109+((320.D0/3.D0*t24+t707+t70
     #8)*s1l1+128.D0/3.D0+(128.D0/3.D0+t711+t581-t707)*s2l1+t714+t708)*t
     #213+t731*t24+(t507-64.D0/3.D0*(11-t447)*t76+(t84-64+t112-t90+t113)
     #*t104)*t275+((-t2-t163-t402-t740-896.D0/3.D0+t113)*t75-64.D0/3.D0*
     #(4-4*s2l1)*t76+(-t112+192+t84+t652+t101+t128)*t98+t651+(t152+t260-
     #512.D0/3.D0+t127+t87+t726)*t104)*t264+((t381-t2+t377+192-t723-t328
     #)*t98+(-t139-t436+t133-192)*t74-t533+(t112+t94+t127+t101+1792.D0/3
     #.D0*s2l1)*t75+(-t139+320-t127-t260-t719-t281)*t104+64.D0/3.D0*(3-t
     #447)*t76)*t110+t777*l1+t304+t2+t101
      t800 = -t84-t145+t508+128-t90
      t802 = -t91-t112+t508-t89-512+t90
      t806 = 192*t76+(t84+t125-t508-t89-192+t278)*t74+(t84-t508+256-t100
     #)*t98+t800*t104+t802*t75+(-t2-t125-t508+t89+192-t278)*t73+t84+t145
     #+t508-t113-64+t100
      t816 = t2+t88-t128+1280.D0/3.D0+t90
      t830 = 64.D0/3.D0*t76
      t849 = 256*t76
      t856 = t2+t145-t508+t113+832.D0/3.D0+t596
      s2 = ((t135-t259+t246+t95+448.D0/3.D0+t329)*t73+(t145+t163-1664.D0
     #/3.D0+t652)*t98-t96-t145-t94-t328-640.D0/3.D0-t652+(-t2+t145-t246-
     #t328+832.D0/3.D0-t329)*t74+(-t145+256-t100)*t75)*t17+t806*t109+((t
     #96+t151+t119+t89+64.D0/3.D0+t337)*t73+(-t2-t125-t128+256.D0/3.D0-t
     #176)*t74+(t145-128.D0/3.D0+t100)*t98-t374-t165-t119-t146-448.D0/3.
     #D0-t278)*l1+(t816*t75+(-t125+t163+t328+256.D0/3.D0-t652)*t98+(-t2+
     #t157-t163-t328-704.D0/3.D0+t100)*t74+(t91-t145+t378+t328+128)*t73+
     #(t112-1024.D0/3.D0+t596)*t104-t91+t145-t723-t146-256.D0/3.D0)*t24+
     #((320.D0/3.D0+t84+t113+t112)*t104-t830-t272)*t275
      s1 = s2+(-t112-64.D0/3.D0-t596)*t74+(t84+t157+t508+256.D0/3.D0+t32
     #9)*t73+(t830+(-t84-t259-t113+704.D0/3.D0-t596)*t104+t272-t816*t75+
     #(t84-t145+t508+t89+448.D0/3.D0-t100)*t74+(t125-t698-t146-64+t652)*
     #t98)*t110+((-t84+t145-t128-t89+64)*t73-t849+(t84-t145+t378+t89-384
     #+t100)*t74+(-t145-t698+t113+576-t652)*t98+(t84+t259+t508+t101-704.
     #D0/3.D0+t90)*t75-t332+t856*t104)*t235+(-t332-t800*t104+(-t84-t145+
     #t87+t113-192)*t98-t802*t75-448.D0/3.D0*t76)*t264
      t886 = s1-128+(t849+t332+(-t84+t112-t508-t101-64.D0/3.D0+t596)*t75
     #-t856*t104)*t274+128.D0/3.D0*t571*t307+(-64.D0/3.D0-64.D0/3.D0*l1-
     #64.D0/3.D0*t17)*t212-t91-t165+((128.D0/3.D0+t714+t183+t708)*s1l1+(
     #256.D0/3.D0-t581+t180+t711)*s2l1+256.D0/3.D0+t581+t183+t708)*t213-
     #t719-t113-t652
      t890 = (l1+z-t7)**2
      taux19 = ((t277+t352)*(t356*t358/2-(t361-t362)*t357-t365*t367/2+t3
     #70*t366)+(t518+t613)*(-t615*t358/2-t356*t357+t361-t362-t619*t367/2
     #+t365*t366)+(t678+t779)*(t615*t357+t354-t355+t619*t366)-t886*t615)
     #/t890*t10
      t1 = 6*s2l1
      t3 = 1-z
      t5 = 128.D0/3.D0*(-4-t1)*t3
      t6 = 4*s1l1
      t7 = 2*s2l1
      t9 = t3**2
      t10 = t9**2
      t11 = t10*t3
      t14 = 14*s1l1
      t15 = 4*s2l1
      t17 = t10*t9
      t20 = 14*s2l1
      t24 = t9*t3
      t25 = t10*t24
      t26 = 512*t25
      t27 = 256.D0/3.D0*s1l1
      t28 = 256.D0/3.D0*s2l1
      t29 = 16*s1l1
      t32 = 128.D0/3.D0*(8+t29-t7)*t9
      t33 = 20*s1l1
      t38 = l1**2
      t39 = t38**2
      t47 = 8*s2l1
      t51 = 34*s1l1
      t55 = 8*s1l1
      t60 = t39*t38
      t62 = 2*s1l1
      t66 = 1792.D0/3.D0*t25
      t67 = 10*s1l1
      t75 = t39**2
      t79 = t75*t38
      t81 = 5*s1l1
      t99 = 6*s1l1
      t102 = 128.D0/3.D0*(8+t99+t47)*t3
      t104 = t38*l1
      t106 = 12*s2l1
      t115 = 128.D0/3.D0*(-12-t14-t7)*t9
      t119 = 18*s1l1
      t120 = 10*s2l1
      t126 = 128.D0/3.D0*(-28-t67-t20)*t24
      t127 = -s1l1+s2l1
      t128 = 2*t127
      t130 = 128.D0/3.D0*t128*t3
      t132 = t39*l1
      t147 = t39*t104
      t156 = t75*l1
      t179 = 24*s1l1
      t185 = -1-s2l1
      s1 = 256.D0/3.D0+(t5-128.D0/3.D0*(20-t6+t7)*t11+128.D0/3.D0*(-t14-
     #2+t15)*t17-128.D0/3.D0*(8-t6+t20)*t10+t26+256.D0/3.D0+t27+t28-t32+
     #128.D0/3.D0*(28+t33+t20)*t24)*t39+(128.D0/3.D0*(-6-t6-t15)*t9-256.
     #D0/3.D0*(-5+t6-t1)*t24-128.D0/3.D0*(42+t14+t47)*t11-128.D0/3.D0*(-
     #20-t51+t1)*t10+t26+128.D0/3.D0*(6-t55+t1)*t17)*t60+(-128.D0/3.D0*(
     #-6-t62-t7)*t10+t66+128.D0/3.D0*(-12-t67+t15)*t17-128.D0/3.D0*(-t55
     #+t1+8)*t11)*t75+t27+t28+(256.D0/3.D0*t25-256.D0/3.D0*t17)*t79+(256
     #.D0/3.D0*(-1+t81-t7)*t24-256.D0/3.D0*(t7+7-3*s1l1)*t17-128.D0/3.D0
     #*(6+t62+t47)*t9-128.D0/3.D0*(-t1-t29-42)*t11-512.D0/3.D0*t25-256.D
     #0/3.D0-t27-t28+128.D0/3.D0*(-22-t51+t15)*t10+t102)*t104
      t188 = s1+128.D0/3.D0*(t99+t106+12)*t9+(-t66+128.D0/3.D0*(t99+t1+2
     #4)*t11-t115-128.D0/3.D0*(-2-t67+t1)*t17+128.D0/3.D0*(4-t119+t120)*
     #t10+t126+t130)*t132+(128.D0/3.D0*(-18-t14+t7)*t10+128.D0/3.D0*(28-
     #t6+t1)*t11-512.D0/3.D0*t127*t24-128.D0/3.D0*(t15-4-t14)*t17-t66)*t
     #147+(128.D0/3.D0*t128*t11-128.D0/3.D0*(-t62-8+t7)*t17-1024.D0/3.D0
     #*t25)*t156+(256.D0/3.D0*(-7*s1l1-8)*t11-256.D0/3.D0*(-s2l1-5)*t17-
     #t126+128.D0/3.D0*(t29-t106-10)*t10-t130+t115)*t38+256.D0/3.D0*(-s1
     #l1-4-t15)*t24+(512.D0/3.D0*(-2-s2l1)*t11-t5+t32-256.D0/3.D0*(-13-t
     #1-t81)*t10-256.D0/3.D0-t27-t28-128.D0/3.D0*(t179+t120+28)*t24)*l1-
     #t102-256.D0/3.D0*t185*t10
      t189 = sqe1(l1)
      t190 = zqe1(l1,z,l3)
      t191 = t189-t190
      t192 = 1-l1
      t193 = dlog(t192)
      t194 = t193**2
      t197 = sqe2(l1)
      t198 = zqe2(l1,z,l3)
      t201 = sbqe1(l1)
      t202 = dlog(z)
      t203 = t202**2
      t206 = sbqe2(l1)
      t211 = 9*s2l1
      t215 = 12*s1l1
      t216 = 42*s2l1
      t220 = 174*s1l1
      t228 = 1280.D0/3.D0*s1l1
      t229 = 1280.D0/3.D0*s2l1
      t230 = 30*s1l1
      t231 = 40*s2l1
      t234 = 128.D0/3.D0*(-40-t230-t231)*t3
      t236 = 82*s1l1
      t237 = 30*s2l1
      t243 = 22*s2l1
      t246 = 128.D0/3.D0*(12-t99+t243)*t3
      t247 = 5*s2l1
      t251 = 1024.D0/3.D0*s1l1
      t252 = 1024.D0/3.D0*s2l1
      t253 = 31*s2l1
      t259 = 58*s2l1
      t266 = 128.D0/3.D0*(-52-86*s1l1-t7)*t9
      t269 = 18*s2l1
      t273 = 8192.D0/3.D0*t25
      t274 = 72*s1l1
      t278 = 46*s1l1
      t283 = 22*s1l1
      t284 = 26*s2l1
      t293 = 74*s1l1
      t300 = -t127
      t307 = 3072*t25
      t318 = 128.D0/3.D0*(56+68*s1l1+t1)*t9
      t321 = 128.D0/3.D0*(128+t278+t259)*t24
      t322 = 35*s1l1
      t329 = 10*t300
      t331 = 128.D0/3.D0*t329*t3
      t341 = 16*s2l1
      t351 = 72*s2l1
      t357 = 20*s2l1
      t369 = 34*s2l1
      t374 = 94*s1l1
      t379 = 56*s1l1
      t397 = -t185
      t401 = 60*s2l1
      s2 = (-256.D0/3.D0*(-34+15*s1l1-t211)*t17-128.D0/3.D0*(-32-t215-t2
     #16)*t9+128.D0/3.D0*(114+t220-t106)*t10+256.D0/3.D0*(2-27*s1l1+t1)*
     #t24+1280.D0/3.D0+t228+t229+t234+2560.D0/3.D0*t25-128.D0/3.D0*(208+
     #t236+t237)*t11)*t104+(-t246+512.D0/3.D0*(10+t247)*t11+1024.D0/3.D0
     #+t251+t252-256.D0/3.D0*(t253+66+25*s1l1)*t10-128.D0/3.D0*(-148-122
     #*s1l1-t259)*t24+t266)*l1+(128.D0/3.D0*(28+t119+t269)*t9-t273-128.D
     #0/3.D0*(-t231-208-t274)*t11+128.D0/3.D0*(t278-24-24*s2l1)*t17-256.
     #D0/3.D0*(22-t283+t284)*t24-128.D0/3.D0*(104+t220-t243)*t10)*t60+(1
     #28.D0/3.D0*(90+t293-t120)*t10+128.D0/3.D0*(-140-t284+t29)*t11-2560
     #.D0/3.D0*t300*t24-128.D0/3.D0*(22+t293-t20)*t17+t307)*t147
      s1 = s2+(-128.D0/3.D0*(-t62-t7-12)*t17-t26)*t79+(t318-t321+256.D0/
     #3.D0*(t322+42+t7)*t11-256.D0/3.D0*(t247+25)*t17-t331+128.D0/3.D0*(
     #38-t236+48*s2l1)*t10)*t38+(-128.D0/3.D0*(-t1-100+t29)*t11+128.D0/3
     #.D0*(10+t274-t341)*t17-t266+t246-1024.D0/3.D0-t251-t252+128.D0/3.D
     #0*(-148-102*s1l1-78*s2l1)*t24-2560*t25-128.D0/3.D0*(-42+t179-t351)
     #*t10)*t39+256.D0/3.D0*(t81+20+t357)*t24
      t405 = s1-1280.D0/3.D0+(128.D0/3.D0*t329*t11-128.D0/3.D0*(44+t29-t
     #1)*t17+5632.D0/3.D0*t25)*t156-t234+(128.D0/3.D0*(-128-t369-32*s1l1
     #)*t11+t331+t307+128.D0/3.D0*(-6+t374-36*s2l1)*t10+t321-t318-128.D0
     #/3.D0*(-t284+10+t379)*t17)*t132-t228-t229+(128.D0/3.D0*(60+t379-t3
     #41)*t17-128.D0/3.D0*(32+t215+t106)*t10-t307-128.D0/3.D0*(-44+38*s1
     #l1-t237)*t11)*t75-1280.D0/3.D0*t397*t10+128.D0/3.D0*(-t401-t230-60
     #)*t9
      t406 = zqe0(l1,z,l3)
      t410 = sbqe0(l1)
      t419 = 2816.D0/3.D0*s1l1
      t420 = 3200.D0/3.D0*s2l1
      t421 = 252*s1l1
      t425 = 48*s1l1
      t436 = 136*s1l1
      t441 = 73*s2l1
      t444 = 128.D0/3.D0*(73+51*s1l1+t441)*t3
      t447 = 52*s2l1
      t448 = 40*s1l1
      t455 = 128.D0/3.D0*(30+52*s1l1+t7)*t3
      t456 = 128*s1l1
      t457 = 128*s2l1
      t470 = 128.D0/3.D0*(148+169*s1l1+53*s2l1)*t9
      t477 = 17*s2l1
      t496 = 28*s1l1
      t499 = 128.D0/3.D0*(-126-t496-t216)*t24
      t506 = 11*s1l1
      t507 = 19*s2l1
      t510 = 128.D0/3.D0*(1-t506+t507)*t3
      t511 = 66*s1l1
      t514 = 128.D0/3.D0*(-55-t511+t477)*t9
      t520 = -t192
      t522 = 1+l1+t38
      t525 = 2*l6-1
      t526 = t525**2
      t529 = l5-1
      t531 = z*l1
      t533 = (1-l1+t531)**2
      t534 = 1/t533
      t539 = dsqrt(z)
      t541 = dsqrt(l5)
      t543 = dsqrt(-t529)
      t546 = dsqrt(l1)
      t550 = l5*z*l1-l1+1-l5+l5*l1
      t552 = dsqrt(t192)
      t557 = t550**2
      t560 = dsqrt(4*l1*t526*z*l5*t529*t520*t534-4*t525*t539*t541*t543*t
     #546*t550*t534*t552+t557*t534)
      t577 = 131*s1l1
      t599 = 25*s2l1
      t603 = 15*s2l1
      t604 = 107*s1l1
      t608 = 70*s2l1
      t612 = 15616.D0/3.D0*t25
      t625 = 84*s1l1
      t648 = 41*s2l1
      s1 = (-256.D0/3.D0*(52-t179+t106)*t17-3200.D0/3.D0-t419-t420+128.D
     #0/3.D0*(-172-t421+t269)*t10+256.D0/3.D0*(13-t7+t425)*t24-4096.D0/3
     #.D0*t25-128.D0/3.D0*(73+43*s1l1+82*s2l1)*t9-128.D0/3.D0*(-304-44*s
     #2l1-t436)*t11+t444)*t104+(-256.D0/3.D0*(-t447-108-t448)*t10-t455+t
     #456-t457-128.D0/3.D0*(270+120*s2l1+200*s1l1)*t24+512.D0/3.D0*(-16-
     #t47)*t11+t470)*l1+(-128.D0/3.D0*(145+s1l1+t211)*t11-t456+t457-t470
     #+t455+4096*t25+128.D0/3.D0*(t477-5-120*s1l1)*t17-128.D0/3.D0*(109-
     #t67+141*s2l1)*t10+128.D0/3.D0*(281+187*s1l1+167*s2l1)*t24)*t39-204
     #8.D0/3.D0*t185*t10+(-t499+256.D0/3.D0*(-t379-72-t47)*t11-256.D0/3.
     #D0*(-t47-40)*t17-t510+t514+128.D0/3.D0*(-t231-6+t436)*t10)*t38+128
     #*t520*t522*t560-t444+(-128.D0/3.D0*(26+11*s2l1+t506)*t17+3328.D0/3
     #.D0*t25)*t79+128.D0/3.D0*(99+99*s2l1+t425)*t9
      t657 = s1+3200.D0/3.D0+(128.D0/3.D0*(-26-t81-t477)*t9-128.D0/3.D0*
     #(311+t577+69*s2l1)*t11+128.D0/3.D0*(t507-t374+27)*t17-256.D0/3.D0*
     #(-19-t369+t448)*t24-128.D0/3.D0*(-156-t421+t369)*t10+14848.D0/3.D0
     #*t25)*t60+256.D0/3.D0*(-t55-32-32*s2l1)*t24+(-128.D0/3.D0*(-52-19*
     #s1l1-t599)*t10+128.D0/3.D0*(-80+t603-t604)*t17-128.D0/3.D0*(t608+9
     #4-t496)*t11+t612)*t75+t419+t420+(128.D0/3.D0*(7-t81+t599)*t11-128.
     #D0/3.D0*(-t247-t278-77)*t17-3584*t25)*t156+(-t612+128.D0/3.D0*(86*
     #s2l1+238+t625)*t11+t510+128.D0/3.D0*(-46-155*s1l1+t603)*t10-t514-1
     #28.D0/3.D0*(t253-t604)*t17+t499)*t132+(128.D0/3.D0*(-107-90*s1l1+3
     #7*s2l1)*t10-128.D0/3.D0*(11-13*s1l1+47*s2l1)*t24+128.D0/3.D0*(209+
     #s1l1+t648)*t11-128.D0/3.D0*(t1-t577-31)*t17-t612)*t147
      t673 = 128.D0/3.D0*(-126*s1l1-125-t441)*t3
      t683 = 5120.D0/3.D0*s1l1
      t684 = 768*s2l1
      t695 = 128.D0/3.D0*(-150*s1l1-211-116*s2l1)*t9
      t696 = 65*s1l1
      t709 = 128.D0/3.D0*(-t322-65-65*s2l1)*t3
      t725 = 1280*s1l1
      t726 = 5248.D0/3.D0*s2l1
      t736 = 28*s2l1
      t749 = 128.D0/3.D0*(-t506-t507-4)*t3
      t757 = 128.D0/3.D0*(-92*s1l1-60-t608)*t9
      t758 = 76*s1l1
      t765 = 128.D0/3.D0*(-108*s1l1-119-t351)*t24
      s2 = 1408.D0/3.D0*t192*t522*t560+128.D0/3.D0*(-59-59*s2l1-t179)*t9
     #+(t673+128.D0/3.D0*(-34-t15+64*s1l1)*t17+128.D0/3.D0*(-159*s1l1-23
     #6-172*s2l1)*t24+3712.D0/3.D0+t683+t684-128.D0/3.D0*(-106*s1l1-175-
     #139*s2l1)*t10-2048*t25-t695-128.D0/3.D0*(-t696-28-27*s2l1)*t11)*t3
     #9+256.D0/3.D0*(t6+16+t341)*t24
      s1 = s2-t709+(256.D0/3.D0*(-26*s1l1-39+s2l1)*t24-256.D0/3.D0*(-24+
     #t215-t15)*t17-128.D0/3.D0*(-t625-93-75*s2l1)*t9+128.D0/3.D0*(t119+
     #27-t447)*t10+2048.D0/3.D0*t25+t725+5248.D0/3.D0+t726+t709-128.D0/3
     #.D0*(t211+71+t274)*t11)*t104+(128.D0/3.D0*(-26-t648-33*s1l1)*t11-1
     #28.D0/3.D0*(38+t736+t379)*t17+t273)*t156+(t307+128.D0/3.D0*(-183-1
     #23*s2l1-138*s1l1)*t11+t749+128.D0/3.D0*(t274+122+51*s2l1)*t10-t757
     #-128.D0/3.D0*(-52-t15+t758)*t17+t765)*t132+(-128.D0/3.D0*(19+t243)
     #*t10+128.D0/3.D0*(-12+t15+t758)*t17-128.D0/3.D0*(-103-107*s2l1-t23
     #6)*t11-t307)*t75
      t854 = s1-5248.D0/3.D0+(512.D0/3.D0*(8+t15)*t11-128.D0/3.D0*(-195-
     #101*s2l1-104*s1l1)*t24-256.D0/3.D0*(t736+56+t33)*t10-t673+t695-371
     #2.D0/3.D0-t683-t684)*l1+(-1024*t25-128.D0/3.D0*(-24-t357-t33)*t17)
     #*t79+(256.D0/3.D0*(t496+40+t47)*t11-t749-256.D0/3.D0*(t15+20)*t17-
     #t765+128.D0/3.D0*(-29*s2l1-103-t274)*t10+t757)*t38+(128.D0/3.D0*(-
     #60*s1l1-34-t341)*t9-128.D0/3.D0*(-97-50*s2l1-105*s1l1)*t11+128.D0/
     #3.D0*(-10+t357+80*s1l1)*t17-256.D0/3.D0*(-t283-23+t477)*t24-128.D0
     #/3.D0*(t119+19-t401)*t10-10240.D0/3.D0*t25)*t60-1024.D0/3.D0*t397*
     #t10+(128.D0/3.D0*(-t511-63-83*s2l1)*t10+128.D0/3.D0*(-t696-60-43*s
     #2l1)*t11+t307-128.D0/3.D0*(-55*s1l1-41-71*s2l1)*t24-128.D0/3.D0*(-
     #10+t341+t625)*t17)*t147-t725-t726
      t858 = (l1+z-t531)**2
      taux20 = (t188*(t191*t194/2-(t197-t198)*t193-t201*t203/2+t206*t202
     #)+t405*(-t406*t194/2-t191*t193+t197-t198-t410*t203/2+t201*t202)+t6
     #57*(t406*t193+t189-t190+t410*t202)-t854*t406)/t858*t534

CCC Coefficients

      tcoff37 = taux19
      tcoff38 = taux19
      tcoff39 = taux20
      tcoff40 = taux20

CCC Jet functions

      tFjet37 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff37)
      tFjet38 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff38)
      tFjet39 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff39)
      tFjet40 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff40)

CCC Sum

      qgRVnum10e3 = tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tFje
     #t40*tcoff40
      RETURN
      END

      double precision function qgRVnum11e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux21 = 0
      taux22 = 0

CCC Coefficients

      tcoff41 = taux21
      tcoff42 = taux21
      tcoff43 = taux22
      tcoff44 = taux22

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff41)
      tFjet42 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff42)
      tFjet43 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff43)
      tFjet44 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff44)

CCC Sum

      qgRVnum11e0 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFje
     #t44*tcoff44
      RETURN
      END

      double precision function qgRVnum11e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -1+l1
      t2 = t1**2
      t3 = 1-z
      t4 = t3**2
      t7 = l1**2
      t8 = 768*t7
      t13 = 1536*l1
      t21 = s1l1**2
      t26 = s2l1**2
      t28 = -t1
      t35 = t7*l1
      t36 = 384*t35
      t70 = dlog(1-l3+l3*l1)
      t71 = z*l3
      t75 = dlog(l1+z-t71-z*l1+t71*l1)
      t76 = dlog(z)
      taux21 = ((-384*t2*t4-t8+1152*t1**2*t3+t13-768)*s2l1-576*t7+1152*l
     #1-576+((t8-t13+768)*t3-768-t8+t13)*t21+(384*t2*t3-384*t2)*t26+((-7
     #68*t1*t28*t3+t13-768-t8)*s2l1+(1920-t36+2688*t7-4224*l1)*t3+192*t2
     #8*(-2*t7+6*l1-4)*t4-1152*t7-1152+2304*l1)*s1l1-192*t2*(-t7+2*l1-2)
     #*t4*t3+(1344+2112*t7-t36-3072*l1)*t3+192*t28*(-5*t7+10*l1-6+t35)*t
     #4)*(t70-t75+t76)
      t1 = -1+l1
      t2 = l1**2
      t3 = 4*t2
      t7 = 1-z
      t8 = t7**2
      t14 = -t1
      t15 = 4*t14
      t20 = t1**2
      t54 = dlog(1-l3+l3*l1)
      t55 = z*l3
      t59 = dlog(l1+z-t55-z*l1+t55*l1)
      t60 = dlog(z)
      taux22 = ((384*t1*(-t3-8+12*l1)*t8-384*t1*(-12-t3+16*l1)*t7-384*t1
     #*t15)*s1l1+768*t20*(-t2+2*l1-2)*t8*t7+384*t1*(-12-10*t2+20*l1+2*t2
     #*l1)*t8-768*t1*t14+(1536*t20*t8+384*t1*t15*t7)*s2l1-384*t1*(-t3-10
     #+14*l1)*t7)*(t54-t59+t60)

CCC Coefficients

      tcoff41 = taux21
      tcoff42 = taux21
      tcoff43 = taux22
      tcoff44 = taux22

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff41)
      tFjet42 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff42)
      tFjet43 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff43)
      tFjet44 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff44)

CCC Sum

      qgRVnum11e1 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFje
     #t44*tcoff44
      RETURN
      END

      double precision function qgRVnum11e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -1+l1
      t2 = t1**2
      t3 = 1-z
      t4 = t3**2
      t7 = l1**2
      t8 = 768*t7
      t13 = 1536*l1
      t21 = s1l1**2
      t26 = s2l1**2
      t28 = -t1
      t29 = 2*t28
      t35 = t7*l1
      t36 = 384*t35
      t37 = 2688*t7
      t43 = -2*t7+6*l1-4
      t47 = 1152*t7
      t51 = 2*l1
      t53 = t2*(-t7+t51-2)
      t54 = t4*t3
      t58 = 3072*l1
      t63 = -5*t7+10*l1-6+t35
      t70 = dlog(1-l3+l3*l1)
      t71 = t70**2
      t73 = z*l3
      t74 = z*l1
      t77 = dlog(l1+z-t73-t74+t73*l1)
      t78 = t77**2
      t80 = dlog(z)
      t82 = t80**2
      t87 = 2*l6-1
      t88 = t87**2
      t91 = l5-1
      t94 = (1-l1+t74)**2
      t95 = 1/t94
      t98 = l1*t88*z*l5*t91*t1*t95
      t100 = dsqrt(z)
      t102 = dsqrt(l5)
      t104 = dsqrt(-t91)
      t107 = dsqrt(l1)
      t111 = l5*z*l1-l1+1-l5+l5*l1
      t113 = dsqrt(t28)
      t116 = t87*t100*t102*t104*t107*t111*t95*t113
      t118 = t111**2
      t119 = t118*t95
      t133 = dsqrt(4*t98-4*t116+t119)
      t139 = 2688*l1
      t153 = -t29
      t207 = 1536*t98-1536*t116+384*t119+(192*(-5*l1-1)*s1l1-192*(3-l1)*
     #s2l1-384*l1-768)*t133+((-1728*t7+t58-1536-192*t35)*t3-t139+1536+23
     #04*t7)*t21+((-384*t1*(t7+4*l1-4)*t3+1920*t7-t139+1536)*s2l1+t37-28
     #80*l1+2496+192*(t28*t7+t153*t43)*t4+(960*t35-5952*t7+8064*l1-3840)
     #*t3)*s1l1+(-192*t2*(5+l1)*t3-384*t1*(-t51+3))*t26+(-192*(-t2*l1-4*
     #t2)*t4+192*t1*(-t7-9*l1+13)*t3+t47-t139+2304)*s2l1+1920-2112*l1-19
     #2*(t2*t7-3*t53)*t54+192*(t28*(7*t7-8*l1+4-2*t35)+t153*t63)*t4+1344
     #*t7+(-3456+1344*t35+7488*l1-5760*t7)*t3
      taux21 = ((-384*t2*t4-t8+1152*t1**2*t3+t13-768)*s2l1-576*t7+1152*l
     #1-576+((t8-t13+768)*t3-768-t8+t13)*t21+(384*t2*t3-384*t2)*t26+((-3
     #84*t1*t29*t3+t13-768-t8)*s2l1+(1920-t36+t37-4224*l1)*t3+192*t28*t4
     #3*t4-t47-1152+2304*l1)*s1l1-192*t53*t54+(1344+2112*t7-t36-t58)*t3+
     #192*t28*t63*t4)*(t71/2-t78/2+t80*t77-t82/2)+t207*(t70-t77+t80)
      t1 = -1+l1
      t2 = l1**2
      t3 = 4*t2
      t4 = 12*l1
      t5 = -t3-8+t4
      t7 = 1-z
      t8 = t7**2
      t11 = -12-t3+16*l1
      t14 = -t1
      t15 = 4*t14
      t20 = t1**2
      t23 = t20*(-t2+2*l1-2)
      t24 = t8*t7
      t29 = t2*l1
      t30 = 2*t29
      t31 = -12-10*t2+20*l1+t30
      t35 = 2*t14
      t40 = -t15
      t41 = t1*t40
      t47 = -t3-10+14*l1
      t54 = dlog(1-l3+l3*l1)
      t55 = t54**2
      t57 = z*l3
      t58 = z*l1
      t61 = dlog(l1+z-t57-t58+t57*l1)
      t62 = t61**2
      t64 = dlog(z)
      t66 = t64**2
      t72 = 2*l6-1
      t73 = t72**2
      t76 = l5-1
      t79 = (1-l1+t58)**2
      t80 = 1/t79
      t85 = dsqrt(z)
      t87 = dsqrt(l5)
      t89 = dsqrt(-t76)
      t92 = dsqrt(l1)
      t96 = l5*z*l1-l1+1-l5+l5*l1
      t98 = dsqrt(t14)
      t103 = t96**2
      t106 = dsqrt(4*l1*t73*z*l5*t76*t1*t80-4*t72*t85*t87*t89*t92*t96*t8
     #0*t98+t103*t80)
      t109 = 9*t2
      t115 = 18*l1
      t124 = t15**2
      taux22 = ((384*t1*t5*t8-384*t1*t11*t7-384*t1*t15)*s1l1+768*t23*t24
     #+384*t1*t31*t8-384*t1*t35+(1536*t20*t8-384*t41*t7)*s2l1-384*t1*t47
     #*t7)*(t55/2-t62/2+t64*t61-t66/2)+(1152*t14*t106+(384*(t1*(8+t109-t
     #4-t30)+t15*t5)*t8-384*(t1*(12-t115+t109)+t15*t11)*t7-384*t1*(-9+6*
     #l1)-384*t124)*s1l1+(384*(t20*(-4-2*t2+l1)-16*t20)*t8-384*(t1*(5-7*
     #l1+5*t2)+t40*t15)*t7-384*t41)*s2l1-768*(-t20*t2+5*t23)*t24+384*(t1
     #*(-21*l1-4*t29+12+16*t2)+t15*t31)*t8-384*(t1*(13+8*t2-t115)+t15*t4
     #7)*t7-384*t1*(-8+5*l1)-384*t15*t35)*(t54-t61+t64)

CCC Coefficients

      tcoff41 = taux21
      tcoff42 = taux21
      tcoff43 = taux22
      tcoff44 = taux22

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff41)
      tFjet42 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff42)
      tFjet43 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff43)
      tFjet44 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff44)

CCC Sum

      qgRVnum11e2 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFje
     #t44*tcoff44
      RETURN
      END

      double precision function qgRVnum11e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -1+l1
      t2 = t1**2
      t3 = 1-z
      t4 = t3**2
      t7 = l1**2
      t8 = 768*t7
      t13 = 1536*l1
      t16 = 576*t7
      t21 = s1l1**2
      t26 = s2l1**2
      t28 = -t1
      t29 = 2*t28
      t35 = t7*l1
      t36 = 384*t35
      t37 = 2688*t7
      t41 = 2*t7
      t43 = -t41+6*l1-4
      t47 = 1152*t7
      t51 = 2*l1
      t53 = t2*(-t7+t51-2)
      t54 = t4*t3
      t57 = 2112*t7
      t58 = 3072*l1
      t63 = -5*t7+10*l1-6+t35
      t68 = dlog(z)
      t69 = t68**2
      t70 = z*l3
      t71 = z*l1
      t74 = dlog(l1+z-t70-t71+t70*l1)
      t77 = t74**2
      t86 = dlog(1-l3+l3*l1)
      t87 = t86**2
      t92 = 4*l1
      t97 = t2*l1
      t100 = 3456*l1
      t104 = 2*l6-1
      t105 = t104**2
      t108 = l5-1
      t111 = (1-l1+t71)**2
      t112 = 1/t111
      t115 = l1*t105*z*l5*t108*t1*t112
      t116 = 1536*t115
      t117 = dsqrt(z)
      t119 = dsqrt(l5)
      t121 = dsqrt(-t108)
      t124 = dsqrt(l1)
      t128 = l5*z*l1-l1+1-l5+l5*l1
      t130 = dsqrt(t28)
      t133 = t104*t117*t119*t121*t124*t128*t112*t130
      t134 = 1536*t133
      t135 = t128**2
      t136 = t135*t112
      t137 = 384*t136
      t139 = 576*t35
      t142 = 1920*l1
      t145 = t2*t7
      t151 = 3*l1
      t152 = t151-5
      t161 = -t29
      t165 = 384*l1
      t182 = 7*t7-8*l1+4-2*t35
      t187 = 5+l1
      t197 = dsqrt(4*t115-4*t133+t136)
      t201 = (192*t1*(t92+t41-8)*t3-384*t97*t4-3072-t47+t100)*s2l1-t116+
     #t134-t137+((960*t7-t139)*t3+t8-t142)*t21-t16-2496-192*(-3*t145+2*t
     #53)*t54+(-192*t2*t152*t3-1536*t1**2)*t26+(t47+192*t161*t7*t4-t142+
     #(-1152+t165-384*t1*(-3*l1+3*t7)*t3)*s2l1+(-t8-t36+t142)*t3-1536)*s
     #1l1+192*(t28*(t35+t41-t151)+t161*t182)*t4+(192*t187*s1l1-192*t152*
     #s2l1+1344-192*l1)*t197+(-t100-t139+2304+t57)*t3+t142
      t213 = 2688*l1
      t276 = (-192*(-t97-4*t2)*t4+t47+2304+192*t1*(-t7+13-9*l1)*t3-t213)
     #*s2l1+(2304*t7+(t58-1536-1728*t7-192*t35)*t3+1536-t213)*t21+(192*(
     #-5*l1-1)*s1l1-192*(3-l1)*s2l1-t165-768)*t197+t116-t134+t137-192*(t
     #145-3*t53)*t54+(-192*t2*t187*t3-384*t1*(3-t51))*t26+((8064*l1+960*
     #t35-5952*t7-3840)*t3+192*(t28*t7+t161*t43)*t4-2880*l1+t37+(1920*t7
     #-384*t1*(t7-4+t92)*t3+1536-t213)*s2l1+2496)*s1l1+1920+192*(t28*t18
     #2+t161*t63)*t4-2112*l1+1344*t7+(1344*t35+7488*l1-5760*t7-3456)*t3
      taux21 = ((-384*t2*t4-t8+1152*t1**2*t3+t13-768)*s2l1-t16+1152*l1-5
     #76+((t8-t13+768)*t3-768-t8+t13)*t21+(384*t2*t3-384*t2)*t26+((-384*
     #t1*t29*t3+t13-768-t8)*s2l1+(1920-t36+t37-4224*l1)*t3+192*t28*t43*t
     #4-t47-1152+2304*l1)*s1l1-192*t53*t54+(1344+t57-t36-t58)*t3+192*t28
     #*t63*t4)*(-t69*t74/2-t77*t74/6+t69*t68/6+t68*t77/2+t87*t86/6)+t201
     #*(t86-t74+t68)+t276*(t87/2-t77/2+t68*t74-t69/2)
      t1 = -1+l1
      t2 = l1**2
      t3 = 4*t2
      t4 = 12*l1
      t5 = -t3-8+t4
      t7 = 1-z
      t8 = t7**2
      t11 = -12-t3+16*l1
      t14 = -t1
      t15 = 4*t14
      t20 = t1**2
      t21 = 2*l1
      t23 = t20*(-t2+t21-2)
      t24 = t8*t7
      t29 = t2*l1
      t30 = 2*t29
      t31 = -12-10*t2+20*l1+t30
      t35 = 2*t14
      t40 = -t15
      t41 = t1*t40
      t47 = -t3-10+14*l1
      t52 = dlog(z)
      t53 = t52**2
      t54 = z*l3
      t55 = z*l1
      t58 = dlog(l1+z-t54-t55+t54*l1)
      t61 = t58**2
      t70 = dlog(1-l3+l3*l1)
      t71 = t70**2
      t78 = 2*l6-1
      t79 = t78**2
      t82 = l5-1
      t85 = (1-l1+t55)**2
      t86 = 1/t85
      t91 = dsqrt(z)
      t93 = dsqrt(l5)
      t95 = dsqrt(-t82)
      t98 = dsqrt(l1)
      t102 = l5*z*l1-l1+1-l5+l5*l1
      t104 = dsqrt(t14)
      t109 = t102**2
      t112 = dsqrt(4*l1*t79*z*l5*t82*t1*t86-4*t78*t91*t93*t95*t98*t102*t
     #86*t104+t109*t86)
      t115 = 9*t2
      t116 = 8+t115-t4-t30
      t121 = 18*l1
      t122 = 12-t121+t115
      t127 = 6*l1
      t128 = -9+t127
      t130 = t15**2
      t136 = t20*(-4-2*t2+l1)
      t137 = 16*t20
      t141 = 5*t2
      t142 = 5-7*l1+t141
      t144 = t40*t15
      t150 = t20*t2
      t158 = -21*l1-4*t29+12+16*t2
      t165 = 13+8*t2-t121
      t172 = -8+5*l1
      t238 = t40**2
      s1 = ((384*t1*t5*t8-384*t1*t11*t7-384*t1*t15)*s1l1+768*t23*t24+384
     #*t1*t31*t8-384*t1*t35+(1536*t20*t8-384*t41*t7)*s2l1-384*t1*t47*t7)
     #*(-t53*t58/2-t61*t58/6+t53*t52/6+t52*t61/2+t71*t70/6)
      s3 = (1152*t14*t112+(384*(t1*t116+t15*t5)*t8-384*(t1*t122+t15*t11)
     #*t7-384*t1*t128-384*t130)*s1l1+(384*(t136-t137)*t8-384*(t1*t142+t1
     #44)*t7-384*t41)*s2l1-768*(-t150+5*t23)*t24+384*(t1*t158+t15*t31)*t
     #8-384*(t1*t165+t15*t47)*t7-384*t1*t172-384*t15*t35)*(t71/2-t61/2+t
     #52*t58-t53/2)
      s4 = (384*(t1*(t30-t127+3*t2)+t15*t158+t40*t31)*t8+4992*t1*t112+(-
     #384*t1*(l1+4)-384*t15*t128-384*t144-384*(t1*(-t4+7*t2)+t15*t122+t4
     #0*t11)*t7+384*(t1*(-t3+3*t29)+t15*t116+t40*t5)*t8)*s1l1-768*(5*t15
     #0-8*t23)*t24-384*(t1*(-2-8*l1+t141)+t15*t165+t40*t47)*t7-384*t1*(3
     #+t21)-384*t15*t172-384*t40*t35+(-384*(t1*(t141-t127)+t15*t142+t238
     #)*t7-384*t1*(1+3*l1)-384*t144+384*(3*t150-4*t136+t137)*t8)*s2l1)*(
     #t70-t58+t52)
      s2 = s3+s4
      taux22 = s1+s2

CCC Coefficients

      tcoff41 = taux21
      tcoff42 = taux21
      tcoff43 = taux22
      tcoff44 = taux22

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff41)
      tFjet42 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff42)
      tFjet43 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff43)
      tFjet44 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff44)

CCC Sum

      qgRVnum11e3 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFje
     #t44*tcoff44
      RETURN
      END

      double precision function qgRVnum12e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -1+l1
      t2 = l1**2
      t3 = t2**2
      t4 = t3*t2
      t5 = 3*t4
      t6 = t2*l1
      t8 = t3*l1
      t12 = t3*t6
      t13 = 9*t12
      t16 = 1-z
      t17 = t16**2
      t18 = t17**2
      t21 = 2944.D0/3.D0*t6
      t22 = 256*t4
      t24 = 3008.D0/3.D0*t3
      t26 = 13312*l1
      t31 = 320*t8
      t32 = 36992.D0/3.D0*l1
      t36 = 64*t3
      t37 = t1**2
      t38 = t37*t1
      t39 = t38*t2
      t40 = 3*t3
      t45 = t18*t17
      t49 = t37*l1
      t56 = t18*t16
      t61 = 36*t3
      t62 = 18*t8
      t67 = t17*t16
      t82 = 9344.D0/3.D0*l1
      t96 = 21248.D0/3.D0*l1
      t100 = 11264.D0/3.D0*l1
      t108 = 640.D0/3.D0*l1
      t110 = 256*t12
      t113 = 12416.D0/3.D0*l1
      t114 = 320.D0/3.D0*t2
      t115 = -3328+64.D0/3.D0*t1*(72+t5+213*t6+39*t8+432*l1-730*t2-t13)*
     #t18-t21+(t22-35072.D0/3.D0-14912.D0/3.D0*t6-t24+1088.D0/3.D0*t8+t2
     #6+11200.D0/3.D0*t2)*t17+(320*t6+t24+27392.D0/3.D0-t31-t32+6592.D0/
     #3.D0*t2)*t16+t36-64.D0/3.D0*t39*(-t40+13*t6-30*t2+72*l1-72)*t45-64
     #.D0/3.D0*t49*(t5-25*t8+49*t3-121*t6+234*t2-36*l1-144)*t56+(-2560+1
     #28.D0/3.D0*t1*(-36-82*t6+t61-t62-134*l1+214*t2)*t67-256*t6-128.D0/
     #3.D0*t39*(-6*t6+8*t2-18*l1+36)*t56+(384*t4-20864.D0/3.D0*t6+7552.D
     #0/3.D0*t3-1152*t8-14336.D0/3.D0+t82+20608.D0/3.D0*t2)*t17-128.D0/3
     #.D0*t49*(72+17*t6-14*t3+6*t8+4*l1-45*t2)*t18+(15872.D0/3.D0+768*t6
     #-t96+1024*t2)*t16+t100-2816.D0/3.D0*t2)*s2l1+(6912-5632.D0/3.D0*t4
     #+68800.D0/3.D0*t6-24256.D0/3.D0*t3+3648*t8+t108-24000*t2+t110)*t67
     #+t113+t114
      t118 = 54*l1
      t123 = t18*t67
      t124 = t2*t123
      t127 = 256*t3
      t128 = 26176.D0/3.D0*l1
      t132 = 1600.D0/3.D0*t6
      t134 = 256*t8
      t135 = 15872.D0/3.D0*l1
      t136 = 6400.D0/3.D0*t2
      t139 = 9472.D0/3.D0*l1
      t140 = 6*t4
      t142 = 68*t3
      t144 = 147*l1
      t176 = s2l1**2
      t180 = 69*t3
      t184 = 3*t12
      t192 = 234*l1
      t200 = t3**2
      t202 = 530*t2
      t215 = 3840*l1
      t217 = 128*t12
      t224 = 10496.D0/3.D0*l1
      t230 = 6208*l1
      t235 = 6208.D0/3.D0*l1
      t237 = t1*l1
      t242 = 144*l1
      t248 = t37*t2
      t252 = 108*l1
      t258 = s1l1**2
      t268 = 39680.D0/3.D0*l1
      t272 = 128*t6
      t274 = -t272+832.D0/3.D0*l1-256+t114
      t276 = 2240.D0/3.D0*l1
      t277 = 320*t2
      t280 = t272-1280.D0/3.D0-1856.D0/3.D0*t2+2752.D0/3.D0*l1
      t284 = 2*l6-1
      t285 = t284**2
      t288 = l5-1
      t290 = z*l1
      t292 = (1-l1+t290)**2
      t293 = 1/t292
      t298 = dsqrt(z)
      t300 = dsqrt(l5)
      t302 = dsqrt(-t288)
      t305 = dsqrt(l1)
      t309 = l5*z*l1-l1+1-l5+l5*l1
      t312 = dsqrt(-t1)
      t317 = t309**2
      t319 = 4*l1*t285*z*l5*t288*t1*t293-4*t284*t298*t300*t302*t305*t309
     #*t293*t312+t317*t293
      t320 = dsqrt(t319)
      t334 = 448.D0/3.D0*l1-256.D0/3.D0-64*t2
      t338 = 44*t2
      t339 = 3*t8
      s1 = t115*s1l1+64.D0/3.D0*t38*(36+85*t6+t62-t118-20*t2-63*t3)*t124
     #+(t21-t127+17536.D0/3.D0-t128+6464.D0/3.D0*t2)*t16+((t132-2368.D0/
     #3.D0*t3+9472.D0/3.D0+t134-t135+t136)*t16+t132-4864.D0/3.D0-t36+t13
     #9+64.D0/3.D0*t1*(-40-t140-178*t6+t142-7*t8-t144+290*t2)*t67+64.D0/
     #3.D0*t1*(128+142*t6-t142+15*t8+19*l1-176*t2)*t17-6016.D0/3.D0*t2-6
     #4.D0/3.D0*t39*(-9*t6+18*t2-29*l1+40)*t56-64.D0/3.D0*t49*(80+49*t6-
     #33*t3+9*t8-10*l1-55*t2)*t18)*t176+64.D0/3.D0*t1*(180-19*t4+362*t6-
     #t180+38*t8+246*l1-691*t2+t184)*t18-64.D0/3.D0*t1*(36-20*t4+189*t6+
     #t61-5*t8+t192-466*t2+6*t12)*t56+64.D0/3.D0*t1*(36-438*t8-51*t200+t
     #192-t202-309*t4-223*t6+266*t12+1017*t3)*t56
      t344 = s1+(-1792+(1536-3776.D0/3.D0*t4+11776*t6-15104.D0/3.D0*t3+8
     #768.D0/3.D0*t8+t215-41728.D0/3.D0*t2+t217)*t67+(-4608+320*t4-11008
     #.D0/3.D0*t6+1024.D0/3.D0*t3-1472.D0/3.D0*t8+t224+4608*t2)*t17+(-t1
     #34+3520.D0/3.D0*t3+14080.D0/3.D0+4928.D0/3.D0*t2-t230-3136.D0/3.D0
     #*t6)*t16-1088.D0/3.D0*t6-t36+t235+448.D0/3.D0*t2+64.D0/3.D0*t237*(
     #144-9*t4+92*t6-31*t3+24*t8-t242-76*t2)*t18-64.D0/3.D0*t248*(-9*t3+
     #21*t6-48*t2+t252-72)*t56)*t258+(21184.D0/3.D0*t2-1088*t8-18560*t6+
     #29120.D0/3.D0*t3-17536.D0/3.D0+t128)*t16+(-3456*t6-27776.D0/3.D0+1
     #280.D0/3.D0*t3+t134+t268-3584.D0/3.D0*t2)*t17+(t274*s1l1+t276-t277
     #-1280.D0/3.D0+t280*s2l1)*t320+(3840-157376.D0/3.D0*t8+1408*l1-768*
     #t200-80704.D0/3.D0*t2+24128*t4+64768.D0/3.D0*t6-1856*t12+93056.D0/
     #3.D0*t3)*t18+t334*t319+64.D0/3.D0*t38*(-10*t3+23*t6-t338+t118-36+t
     #339)*t124
      t351 = l1*t45
      t356 = 2176*t12
      t357 = 24448.D0/3.D0*l1
      t388 = 51*t4
      t392 = 290*l1
      t418 = 6912*l1
      t424 = 10368*l1
      t428 = 16640.D0/3.D0*l1
      t430 = 64.D0/3.D0*t1*(-188-12*t4-247*t6+21*t3+16*t8-336*l1+646*t2)
     #*t67-64.D0/3.D0*t1*(-36-t388+88*t6-226*t3+127*t8-t392+359*t2+t13)*
     #t18+128.D0/3.D0*t6+t36+64.D0/3.D0*t49*(72+t5+80*t6-23*t3+2*t8+88*l
     #1-182*t2)*t56-64.D0/3.D0*t39*(t40-4*t6+3*t2+14*l1-36)*t45-9728.D0/
     #3.D0+(t22-36800.D0/3.D0*t6+16768.D0/3.D0*t3-23296.D0/3.D0-4672.D0/
     #3.D0*t8+t418+8832*t2)*t17+(21248.D0/3.D0+7744.D0/3.D0*t6-4352.D0/3
     #.D0*t3+t31-t424+5504.D0/3.D0*t2)*t16+t428-7232.D0/3.D0*t2
      t435 = 54*t6
      t456 = t37**2
      t460 = -4+3*l1
      t461 = t460*t123
      t463 = 64.D0/3.D0*t456*t1*t3*t461
      t491 = 64.D0/3.D0*t39*(t339+81*t2-4*t3+6*l1-t435-36)*t45+(-21248.D
     #0/3.D0-61312.D0/3.D0*t6+9216*t3-384*t8+t424+8320*t2)*t16-t428+64.D
     #0/3.D0*t1*(-36+156*t8-1055*t3+303*t4+602*t6-t392-99*t12+415*t2)*t1
     #8-t463-64.D0/3.D0*t49*(72+t140-55*t6+372*t3-309*t2-170*t8+92*l1)*t
     #56-1472.D0/3.D0*t2+15104.D0/3.D0*t6+(384*t12-89408.D0/3.D0*t3+1118
     #08.D0/3.D0*t8-49600.D0/3.D0*t6+83840.D0/3.D0*t2-12096*t4-t139-1203
     #2.D0/3.D0)*t67-2240*t3+(4160*t4+29248*t6+576*t3+23296.D0/3.D0-t418
     #-65408.D0/3.D0*t2-39104.D0/3.D0*t8)*t17+9728.D0/3.D0
      t552 = 128*t8
      t556 = 2560-128.D0/3.D0*t456*t8*t461-t100-5888.D0/3.D0*t2-128.D0/3
     #.D0*t248*(131*t6-58*l1+t339+36-t36-47*t2)*t56+(-1536+28544*t8-7667
     #2.D0/3.D0*t3-27520.D0/3.D0*t4-18304.D0/3.D0*t6-12544.D0/3.D0*l1+t1
     #10+53248.D0/3.D0*t2)*t67+128.D0/3.D0*t237*(-72-t388-258*t6-39*t3+1
     #93*t2+158*t8+68*l1)*t18+16000.D0/3.D0*t6-2176*t3+(4096*t4+50176.D0
     #/3.D0*t6+19456.D0/3.D0*t3-40192.D0/3.D0*t8-t82-46720.D0/3.D0*t2+14
     #336.D0/3.D0)*t17-128.D0/3.D0*t38*t3*t460*(-t2+2)*t45+(-48512.D0/3.
     #D0*t6+20480.D0/3.D0*t3-15872.D0/3.D0-t552+t96+7680*t2)*t16
      t614 = 3328+t556*s2l1-64.D0/3.D0*t38*t8*t460*(1+l1)*t123-64.D0/3.D
     #0*t237*(-144+t252-423*t4-654*t6-340*t3+867*t8+54*t12+t202)*t56+844
     #8*t6-3264*t3-t113+(1536-3136*t200+141440.D0/3.D0*t3+24128.D0/3.D0*
     #t4+21376.D0/3.D0*t6-149504.D0/3.D0*t8+27136.D0/3.D0*t12+7680*l1-82
     #816.D0/3.D0*t2)*t18+(6208*t4+35072.D0/3.D0+115712.D0/3.D0*t6-128.D
     #0/3.D0*t3-t26-74944.D0/3.D0*t2-54400.D0/3.D0*t8)*t17+(-87808.D0/3.
     #D0*t6-27392.D0/3.D0+40640.D0/3.D0*t3-896*t8+t32+40256.D0/3.D0*t2)*
     #t16+64.D0/3.D0*t248*(72+270*t6+50*t8-22*t2-t242+t5-225*t3)*t45-131
     #20.D0/3.D0*t2+(-6912+58624*t8-118016.D0/3.D0*t3-70336.D0/3.D0*t4-8
     #2816.D0/3.D0*t6-t108+t356+110144.D0/3.D0*t2)*t67
      t667 = 1792-8896.D0/3.D0*t2+(-41728.D0/3.D0*t6+16384.D0/3.D0*t3+t2
     #30+20416.D0/3.D0*t2-14080.D0/3.D0+t552)*t16+64.D0/3.D0*t1*t2*(t338
     #+287*t6+72-338*t3-t140+122*t8-180*l1)*t56+64.D0/3.D0*l1*(144-96*t1
     #2-383*t3+392*t4+792*t6-373*t8-288*l1-188*t2)*t18-64.D0/3.D0*t37*t8
     #*t460*(-t2-l1+4)*t45+(-1536-70592.D0/3.D0*t3-t215-24704.D0/3.D0*t4
     #-5952*t6+t217+26240*t8+50176.D0/3.D0*t2)*t67+15872.D0/3.D0*t6-2048
     #*t3-64.D0/3.D0*t38*t4*t461+(4608+3968*t4+40000.D0/3.D0*t6+8000*t3-
     #13056*t2-40064.D0/3.D0*t8-t224)*t17-t235
      t716 = 4864.D0/3.D0+(1920*t4+33536.D0/3.D0*t6+2048.D0/3.D0*t3-5760
     #*t8-6976.D0/3.D0*l1-25280.D0/3.D0*t2+8192.D0/3.D0)*t17-t139+6016.D
     #0/3.D0*t6-64.D0/3.D0*t456*t6*t460*(-t2+2+l1)*t45-1024*t3+(-24448.D
     #0/3.D0*t6+12160.D0/3.D0*t3-9472.D0/3.D0-t134+t135+6656.D0/3.D0*t2)
     #*t16+1664.D0/3.D0*t2-64.D0/3.D0*t39*(-t435+16*l1-44+78*t2)*t56-t46
     #3+64.D0/3.D0*t49*(80+70*t6+106*t3+2*l1-204*t2-48*t8)*t18-64.D0/3.D
     #0*t1*(-40-t140+76*t6-438*t3-t144+346*t2+206*t8)*t67
      t718 = 64.D0/3.D0*t37*(-72+13*t4-139*t6+t180-30*t8+142*t2-t184)*t3
     #51+(-10880*t3-24448.D0/3.D0+34112*t8+t356+t357+72320.D0/3.D0*t2-98
     #944.D0/3.D0*t6-49600.D0/3.D0*t4)*t67-64.D0/3.D0*t37*(-72+327*t8+27
     #0*t2-132*t4-159*t6-t127+18*t12)*t351+(-t274*s1l1-t280*s2l1-t276+12
     #80.D0/3.D0+t277)*t320+t430*s2l1+t491*s2l1-t334*t319+(39808.D0/3.D0
     #*t6+24448.D0/3.D0-3840*t3-t357-28288.D0/3.D0*t2)*t67+(27776.D0/3.D
     #0-16960.D0/3.D0*t8+2752*t4+99584.D0/3.D0*t6-33920.D0/3.D0*t3-45056
     #.D0/3.D0*t2-t268)*t17+t614*s1l1+t667*t258-t136+11200.D0/3.D0*t6-16
     #00*t3+t716*t176
      t721 = (l1+z-t290)**2
      taux23 = (t344+t718)/t721*t293
      t1 = -1+l1
      t2 = t1**2
      t3 = t2*t1
      t4 = l1**2
      t5 = t3*t4
      t6 = t4*l1
      t7 = 26*t6
      t8 = t4**2
      t9 = 6*t8
      t13 = 1-z
      t14 = t13**2
      t15 = t14**2
      t16 = t15*t14
      t19 = t2*l1
      t20 = t8*t4
      t21 = 6*t20
      t24 = t8*l1
      t25 = 26*t24
      t29 = t15*t13
      t34 = 184*t8
      t36 = 576*l1
      t38 = t8*t6
      t39 = 6*t38
      t43 = 12*t20
      t47 = 116*l1
      t51 = t14*t13
      t55 = 12*t24
      t56 = 372*l1
      t61 = 8*t6
      t62 = 228*l1
      t67 = 6*t6
      t68 = 28*l1
      t69 = 14*t4
      t80 = 22*t8
      t89 = 60*t8
      t90 = 260*l1
      t96 = 90*t8
      t98 = 68*l1
      t104 = 60*l1
      t105 = 44*t4
      t109 = 4*l1
      t110 = 48*t4
      t114 = 8*l1
      t116 = t1*(-t67-t114+t69)
      t126 = t15*t51
      t127 = t4*t126
      t130 = 26*t20
      t133 = 60*t24
      t137 = l1*t16
      t144 = 468*l1
      t155 = 332*l1
      t165 = 310*l1
      t171 = 54*t6
      t173 = 474*l1
      t180 = 166*l1
      t186 = 10*l1
      t199 = t8**2
      t236 = -2*(-4+3*l1)*t1*t126
      t294 = 460*t4
      t320 = t2**2
      t337 = 6*t24
      t339 = 20*t8
      s3 = (128.D0/3.D0*t5*(-144+t7-t9+144*l1-60*t4)*t16-128.D0/3.D0*t19
     #*(288-t21+166*t6-54*t8+t25-72*l1-268*t4)*t29+128.D0/3.D0*t1*(-144+
     #50*t20-646*t6+t34-122*t24-t36+1220*t4-t39)*t15-128.D0/3.D0*t1*(-50
     #4+t43-392*t6+74*t8-40*t24-t47+810*t4)*t51-128.D0/3.D0*t1*(648-32*t
     #6+58*t8-t55-t56-62*t4)*t14+128.D0/3.D0*t1*(360-t61+t9-t62+18*t4)*t
     #13+128.D0/3.D0*t1*(-72-t67+t68+t69))*s1l1+(128.D0/3.D0*t5*(-72-t61
     #+t9+t68+6*t4)*t16-128.D0/3.D0*t19*(144+t21+16*t6+t80-20*t24+16*l1-
     #104*t4)*t29+128.D0/3.D0*t1*(-72+2*t20-60*t6-t89+t25-t90+390*t4-t39
     #)*t15-128.D0/3.D0*t1*(-216-t43+14*t6-t96+52*t24-t98+200*t4)*t51-12
     #8.D0/3.D0*t1*(216-2*t6+t80-t55-t104-t105)*t14+128.D0/3.D0*t1*(72+t
     #7-t9-t109-t110)*t13+128.D0/3.D0*t116)*s2l1
      s2 = s3-256.D0/3.D0*t3*(-10*t8+23*t6-t105+54*l1-36+3*t24)*t127+128
     #.D0/3.D0*t2*(144-t130+278*t6-138*t8+t133-284*t4+t39)*t137-128.D0/3
     #.D0*t1*(-72+80*t20-580*t6+108*t8-76*t24-t144+1012*t4-24*t38)*t29
      s1 = s2+128.D0/3.D0*t1*(-360+78*t20-854*t6+302*t8-150*t24-t155+123
     #4*t4-18*t38)*t15-128.D0/3.D0*t1*(-684-182*t6-52*t8+30*t24+t165+380
     #*t4)*t51-128.D0/3.D0*t1*(612-t171+t89-18*t24-t173+68*t4)*t14+128.D
     #0/3.D0*t1*(252+12*t6-t180-4*t4)*t13+128.D0/3.D0*t1*(-36-t67+t186+t
     #69)
      s2 = s1+128.D0/3.D0*t1*(-684+t165-180*t20+1644*t4+876*t24-1114*t6-
     #872*t8)*t51+128.D0/3.D0*t1*(-72+90*t199+138*t20+1148*t4-t144+254*t
     #6+1312*t24-354*t38-2050*t8)*t29+128.D0/3.D0*t1*(612-t133-1052*t4+1
     #154*t6-160*t8-t173)*t14-128.D0/3.D0*t1*(-36+78*t4+t186-t171)
      s3 = s2-128.D0/3.D0*t3*(72-108*l1-t110-154*t8+48*t24+192*t6)*t127+
     #(128.D0/3.D0*t3*t24*t236-128.D0/3.D0*t1*(-144+1484*t4-30*t20-t36+3
     #24*t6-90*t38-2348*t8+1376*t24)*t15-128.D0/3.D0*t1*(-72+142*t4+t68-
     #102*t6)+128.D0/3.D0*t1*(-504-1674*t8-168*t20+1142*t24-426*t6-t47+1
     #730*t4)*t51+128.D0/3.D0*t1*(648-198*t24-1238*t4-t56+962*t6+222*t8)
     #*t14-128.D0/3.D0*t2*t4*(144+104*t24-t21-t105+476*t6-386*t8-288*l1)
     #*t16+128.D0/3.D0*t1*l1*(-288+360*l1-538*t20+84*t38-1224*t6+716*t4+
     #970*t24-80*t8)*t29-128.D0/3.D0*t1*(360-84*t8-t62-630*t4+598*t6)*t1
     #3)*s1l1
      t363 = s3-128.D0/3.D0*t1*(252-t96-t180-t294+474*t6)*t13-128.D0/3.D
     #0*t1*(-360+2064*t24-578*t20-t155-430*t6+1906*t4+42*t38-2322*t8)*t1
     #5-128.D0/3.D0*t2*(144+354*t6+292*t20-652*t24-556*t4-48*t38+470*t8)
     #*t137+(128.D0/3.D0*t320*t8*t236-128.D0/3.D0*t1*(-72+550*t4-t130-t9
     #0+340*t6+t39-986*t8+444*t24)*t15-128.D0/3.D0*t116+128.D0/3.D0*t1*(
     #-216+632*t4+210*t24-62*t6-508*t8-t98)*t51+128.D0/3.D0*t1*(216-t337
     #-t294+302*t6+t339-t104)*t14-128.D0/3.D0*t5*(-72+44*l1-112*t6-t337+
     #t339+122*t4)*t16+128.D0/3.D0*t19*(144-84*t24-394*t4+t34+150*t6+t11
     #4)*t29-128.D0/3.D0*t1*(72-t109-184*t4+120*t6)*t13)*s2l1
      t364 = z*l1
      t366 = (l1+z-t364)**2
      t370 = (1-l1+t364)**2
      taux24 = t363/t366/t370

CCC Coefficients

      tcoff45 = taux23
      tcoff46 = taux23
      tcoff47 = taux24
      tcoff48 = taux24

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff45)
      tFjet46 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff46)
      tFjet47 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff47)
      tFjet48 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff48)

CCC Sum

      qgRVnum12e0 = tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tFje
     #t48*tcoff48
      RETURN
      END

      double precision function qgRVnum12e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = l1**2
      t2 = 64*t1
      t4 = -1+l1
      t5 = t4**2
      t6 = t5*l1
      t7 = 13*l1
      t8 = t1**2
      t9 = t8*l1
      t11 = t1*l1
      t12 = t8*t11
      t13 = 12*t12
      t17 = t8*t1
      t18 = 79*t17
      t20 = 1-z
      t21 = t20**2
      t22 = t21**2
      t27 = t5*t4
      t28 = t27*t1
      t30 = 12*t9
      t31 = 76*t8
      t33 = t22*t20
      t39 = 8192*l1
      t40 = 10432.D0/3.D0*t8
      t41 = 1024*t12
      t42 = 14848.D0/3.D0*t17
      t46 = t8**2
      t47 = 256*t46
      t57 = 2624.D0/3.D0*t9
      t59 = t5**2
      t60 = t59*t11
      t61 = 3*l1
      t62 = -4+t61
      t63 = 4*t1
      t67 = t22*t21
      t73 = t21*t20
      t74 = t22*t73
      t75 = t62*t74
      t76 = t59*t4*t8*t75
      t78 = 256*t17
      t82 = 24*t46
      t85 = 74*t12
      t91 = -t2+4992*l1-21248.D0/3.D0+64.D0/3.D0*t6*(-t7-90*t9-t13-217*t
     #8-192+200*t1+204*t11+t18)*t22+24832.D0/3.D0*t11-64.D0/3.D0*t28*(t7
     #+90+t11-147*t1-t30+t31)*t33+(17024.D0/3.D0*t9-4224*t1+28160.D0/3.D
     #0-t39+t40+t41-t42-6592.D0/3.D0*t11)*t20+(-20992.D0/3.D0+t47+2176*l
     #1+17408.D0/3.D0*t8-44672.D0/3.D0*t9-13376*t11-11200.D0/3.D0*t12+17
     #472*t1+39872.D0/3.D0*t17)*t21+t57-21760.D0/3.D0*t8-64.D0/3.D0*t60*
     #t62*(t63-4*l1-3)*t67+128.D0/3.D0*t76+t78-64.D0/3.D0*t4*(-103*t8+39
     #5*l1-685*t1+98+t82+461*t9+44*t11-t85-139*t17)*t73
      t92 = s2l1**2
      t96 = 36*l1
      t103 = l1*t67
      t106 = 1088*t9
      t108 = 36*t12
      t126 = t46*l1
      t132 = 6*l1
      t143 = 12*t126
      t149 = 163*t46
      t164 = 1536*t12
      t166 = 9344*t17
      t170 = 768*t17
      t172 = 3072*l1
      t173 = 46528.D0/3.D0*t8
      t174 = 6144*t12
      t181 = -44800.D0/3.D0+t106-64.D0/3.D0*t6*(240*t17-t108+342*t1-109*
     #t8-352*t9-210*l1-144+214*t11)*t33-2944.D0/3.D0*t1-17088*t8+(25216.
     #D0/3.D0+135616.D0/3.D0*t9-134144.D0/3.D0*t1-54016.D0/3.D0*t8+3520*
     #t46+33344.D0/3.D0*l1+73600.D0/3.D0*t11-768*t126-94976.D0/3.D0*t17+
     #2368*t12)*t73+64.D0/3.D0*t28*(72-t132-33*t9-46*t1+128*t8-91*t11)*t
     #67+64.D0/3.D0*t4*(747*t9+626*l1-t143-102*t8-430*t12-453*t11+72+42*
     #t17-623*t1+t149)*t22-64.D0/3.D0*t59*t8*t62*(2-t61)*t74+(14848*t9-4
     #544*t1+60928.D0/3.D0-25216.D0/3.D0*l1-78784.D0/3.D0*t11+t164+35584
     #.D0/3.D0*t8-t166)*t20+11392.D0/3.D0*l1+t170+(-53504.D0/3.D0+79808.
     #D0/3.D0*t1+t172+t173-t174-46272*t9-14464.D0/3.D0*t11+29888*t17)*t2
     #1+82048.D0/3.D0*t11
      t189 = 19328.D0/3.D0*t12
      t190 = 512*t46
      t198 = 21760.D0/3.D0*t17
      t203 = 3584.D0/3.D0*l1
      t205 = 2176.D0/3.D0*t9
      t206 = 512*t17
      t207 = 346*t1
      t214 = 12*t46
      t233 = (28672.D0/3.D0-3712*t1-107392.D0/3.D0*t8-7168.D0/3.D0*l1+40
     #96.D0/3.D0*t11+155392.D0/3.D0*t9-26752*t17+t189-t190)*t21+9728+(-3
     #6352.D0/3.D0-39808.D0/3.D0*t1-14720.D0/3.D0*t8+4352*l1+106496.D0/3
     #.D0*t11-47360.D0/3.D0*t9+t198-t41)*t20+6912*t1+44032.D0/3.D0*t8+t2
     #03-93824.D0/3.D0*t11-t205-t206+128.D0/3.D0*t4*(72-t207-655*t8+268*
     #l1+354*t11+441*t9-63*t17-43*t12+t214)*t73-128.D0/3.D0*t6*(-144+52*
     #t1+116*t8-8*l1-t11-162*t9+t18-t13)*t22-128.D0/3.D0*t28*(t30-t31+86
     #*t11-26*t1+t96-72)*t33
      t238 = 256*t126
      t243 = 5056.D0/3.D0*t46
      t249 = 192*t9
      t262 = 24*t9
      t265 = 78*t1
      t266 = 144*l1
      t287 = 4608*t12
      t294 = 20096.D0/3.D0*t17
      t295 = 512*t12
      t298 = 12416+t233*s2l1+(-13824+33792*t1+257728.D0/3.D0*t8-1792*l1+
     #t238-152128.D0/3.D0*t11-219776.D0/3.D0*t9+20096*t17+1216*t12-t243)
     #*t73+27200.D0/3.D0*t1+t173+5568*l1-41984*t11+t249-t170+64.D0/3.D0*
     #t4*(-144+1368*t1+689*t8-864*l1+t143-984*t11+353*t9-763*t17+444*t12
     #-151*t46)*t22-64.D0/3.D0*t28*(-t262+35*t8-49*t11+t265-t266+144)*t6
     #7-64.D0/3.D0*t6*(t108-186*t17+334*t9-458*t8+398*t11-404*t1+72*l1+2
     #88)*t33+(71296.D0/3.D0+59072.D0/3.D0*t1-154048.D0/3.D0*t8-52288.D0
     #/3.D0*l1-22080*t11+208256.D0/3.D0*t9-26624*t17+t287)*t21+(-62848.D
     #0/3.D0-82496.D0/3.D0*t1-t40+6976*l1+59584*t11-62464.D0/3.D0*t9+t29
     #4-t295)*t20
      t324 = 132*t17
      t326 = 36*t46
      t360 = 44800.D0/3.D0+(53504.D0/3.D0-29888.D0/3.D0*t1-84416.D0/3.D0
     #*t8-4096.D0/3.D0*l1-53504.D0/3.D0*t11+210176.D0/3.D0*t9-36736*t17+
     #t174)*t21+12544.D0/3.D0*t1+(-60928.D0/3.D0-8256*t1-68416.D0/3.D0*t
     #8+20096.D0/3.D0*l1+52096*t11-15232*t9+t166-t164)*t20+20608*t8-2944
     #*l1-104768.D0/3.D0*t11-t106-t170+64.D0/3.D0*t4*(394-691*t1-2418*t8
     #+955*l1+559*t11+1706*t9-t324-129*t12+t326)*t73+64.D0/3.D0*t6*(-144
     #+269*t1+485*t8-230*l1-192*t11-452*t9+180*t17-t108)*t33-64.D0/3.D0*
     #t28*(-t262+41*t8-36*t11+17*t1-10*l1+72)*t67-64.D0/3.D0*t4*(72-489*
     #t1-615*t8+634*l1-t143+30*t11+5*t9+891*t17-619*t12+t149)*t22
      t369 = 256*t46*t1
      t372 = -39296.D0/3.D0*t17-10752-8960*t46+20288*t8+35264.D0/3.D0*t1
     #2-t203+11264*t9+9664.D0/3.D0*t126-t369-68608*t11+169024.D0/3.D0*t1
      t377 = 108*l1
      t378 = 55*t11
      t383 = 44*t1
      t384 = 54*l1
      t385 = 3*t9
      t393 = t27*t17*t75
      t399 = t59*t9*t75
      t408 = t27*t8
      t414 = t5*t1
      t415 = 12*t17
      t442 = l1*t4
      t447 = 91*t12
      t453 = -14080.D0/3.D0*t1-512*l1+t205-36352.D0/3.D0*t8+256.D0/3.D0*
     #t399+(15488.D0/3.D0*t1-4608*t8+50048.D0/3.D0*t9-5376*l1-t198+36352
     #.D0/3.D0-53248.D0/3.D0*t11+t41)*t20+t206-128.D0/3.D0*t408*t62*(t63
     #-t61-3)*t67-9728-128.D0/3.D0*t414*(-t415+t377-72-94*t11+85*t9+42*t
     #1-62*t8)*t33+(45824.D0/3.D0*t9+t39-512*t126-14080*t8+4480.D0/3.D0*
     #t12+26240*t11+3072-37120.D0/3.D0*t17-88960.D0/3.D0*t1+7040.D0/3.D0
     #*t46)*t73+77440.D0/3.D0*t11+(-103552.D0/3.D0*t9+t172+28160*t8-t189
     #+t190-42752.D0/3.D0*t11-28672.D0/3.D0+20608*t17+12416*t1)*t21+128.
     #D0/3.D0*t442*(-157*t17+188*t11-170*t1+71*t8+t447-2*t9-t266+144-t21
     #4)*t22
      t480 = -3072+148736.D0/3.D0*t1-t369-27520.D0/3.D0*t46+10432.D0/3.D
     #0*t126+99328.D0/3.D0*t9-12608.D0/3.D0*t8-123712.D0/3.D0*t11-15360*
     #l1+14912*t12-27776*t17
      t487 = 288*l1
      t496 = 216*l1
      t515 = 64*t393+t453*s2l1+(14784*t1-35648.D0/3.D0*t8-8000*l1-t294+6
     #2848.D0/3.D0-97216.D0/3.D0*t11+t295+68224.D0/3.D0*t9)*t20+(46336.D
     #0/3.D0*t9+4864.D0/3.D0*l1-t238-134336.D0/3.D0*t8-2880*t12+166208.D
     #0/3.D0*t11+13824+10496.D0/3.D0*t17-130688.D0/3.D0*t1+t243)*t73+t48
     #0*t22-37120.D0/3.D0*t8+64.D0/3.D0*t414*(38*t9-100*t11-33*t17-144+9
     #2*t8+t487-154*t1)*t67-t249+105536.D0/3.D0*t11-64.D0/3.D0*t442*(-62
     #8*t1+1008*t11-t496+180*t12+288-413*t8-125*t17-t326-44*t9)*t33-1241
     #6-6080*t1+(-71296.D0/3.D0-7360.D0/3.D0*t1+18112*l1+53248*t8-153344
     #.D0/3.D0*t9-25664.D0/3.D0*t11-t287+19136*t17)*t21-14656.D0/3.D0*l1
     #+t170
      t534 = t1*t74
      t540 = 15808.D0/3.D0*t12
      t552 = 94*t17
      t578 = 64*t8
      t594 = 21248.D0/3.D0+64.D0/3.D0*t4*(98-524*t1-860*t8+415*l1+356*t1
     #1+719*t9-t552-t85+t82)*t73+64.D0/3.D0*t4*(-328+183*t1+668*t8-266*l
     #1+64*t11-652*t9+163*t17-t13)*t21+2304*t1+(-28160.D0/3.D0-3328*t1-1
     #0560*t8+22016.D0/3.D0*l1+52672.D0/3.D0*t11-16640.D0/3.D0*t9+t42-t4
     #1)*t20+28480.D0/3.D0*t8-13696.D0/3.D0*l1-13184*t11-t57-t78-64.D0/3
     #.D0*t28*(t30-t578+121*t11-76*t1+45*l1-98)*t33-64.D0/3.D0*t6*(-196+
     #t383+173*t8-42*l1+59*t11-225*t9+t18-t13)*t22
      t596 = 256*t8
      t599 = t596-1600.D0/3.D0*t1+2816.D0/3.D0-1984.D0/3.D0*t11
      t601 = 2*l6-1
      t602 = t601**2
      t605 = l5-1
      t607 = z*l1
      t609 = (1-l1+t607)**2
      t610 = 1/t609
      t615 = dsqrt(z)
      t617 = dsqrt(l5)
      t619 = dsqrt(-t605)
      t622 = dsqrt(l1)
      t626 = l5*z*l1-l1+1-l5+l5*l1
      t629 = dsqrt(-t4)
      t634 = t626**2
      t636 = 4*l1*t602*z*l5*t605*t4*t610-4*t601*t615*t617*t619*t622*t626
     #*t610*t629+t634*t610
      t640 = 9344.D0/3.D0*t17
      t646 = t91*t92+64.D0/3.D0*t5*(216-340*t1-325*t8-t96+460*t11+194*t9
     #-149*t17+60*t12)*t103+t181*s2l1+t298*s1l1+t360*s2l1+t372*t22-64.D0
     #/3.D0*(-t27*(106*t1-9*t9+24*t8-t377+72-t378)+t27*(-10*t8+23*t11-t3
     #83+t384-36+t385))*t1*t74+t515*s1l1-64.D0/3.D0*t5*(-135*t9+6*t12+52
     #0*t11+64*t17-74*t8-t96+216-532*t1)*t103+64.D0/3.D0*t27*(-54*t1-42*
     #t9-108+101*t8-78*t11+162*l1)*t534+(-73984.D0/3.D0+22528*t8+2816.D0
     #/3.D0*t1-9856*t11-t540+27776*l1-88064.D0/3.D0*t9+52928.D0/3.D0*t17
     #+t47)*t21+1280.D0/3.D0*l1+t594*t92+t599*t636+(53504.D0/3.D0+8512.D
     #0/3.D0*t1+4416*t8+t295-t640-48896.D0/3.D0*t11+17792.D0/3.D0*t9-363
     #52.D0/3.D0*l1)*t20
      t662 = 630*l1
      t663 = 60*t126
      t673 = 1792*t8
      t675 = 6208.D0/3.D0*l1
      t680 = 512*t9
      t681 = 7808.D0/3.D0*t8
      t688 = 7360.D0/3.D0*t8
      t689 = 256*t9
      t690 = 512*t1
      t692 = dsqrt(t636)
      t695 = dlog(z)
      t696 = 128*t11
      t763 = 2176*t12
      t814 = 128*t9
      t818 = 2560-11264.D0/3.D0*l1-5888.D0/3.D0*t1+128.D0/3.D0*t442*(-72
     #-39*t8-51*t17+193*t1+158*t9-258*t11+68*l1)*t22+(-1536-18304.D0/3.D
     #0*t11-27520.D0/3.D0*t17-76672.D0/3.D0*t8-12544.D0/3.D0*l1+53248.D0
     #/3.D0*t1+256*t12+28544*t9)*t73-128.D0/3.D0*t399-128.D0/3.D0*t414*(
     #131*t11-58*l1+t385+36-t578-47*t1)*t33+(50176.D0/3.D0*t11+4096*t17+
     #19456.D0/3.D0*t8-9344.D0/3.D0*l1-46720.D0/3.D0*t1+14336.D0/3.D0-40
     #192.D0/3.D0*t9)*t21-128.D0/3.D0*t408*t62*(-t1+2)*t67-2176*t8+16000
     #.D0/3.D0*t11+(-48512.D0/3.D0*t11+20480.D0/3.D0*t8+7680*t1-15872.D0
     #/3.D0-t814+21248.D0/3.D0*l1)*t20
      t820 = 530*t1
      t835 = 7680*l1
      t850 = 3328+8448*t11-64.D0/3.D0*t27*t9*t62*(1+l1)*t74-13120.D0/3.D
     #0*t1+64.D0/3.D0*t414*(72+270*t11+50*t9-22*t1-t266+3*t17-225*t8)*t6
     #7+(-87808.D0/3.D0*t11-27392.D0/3.D0+40640.D0/3.D0*t8-896*t9+40256.
     #D0/3.D0*t1+36992.D0/3.D0*l1)*t20+(-6912-82816.D0/3.D0*t11-70336.D0
     #/3.D0*t17-118016.D0/3.D0*t8-640.D0/3.D0*l1+110144.D0/3.D0*t1+t763+
     #58624*t9)*t73+t818*s2l1-64.D0/3.D0*t442*(-144+t820-423*t17-654*t11
     #-340*t8+t377+867*t9+54*t12)*t33-3264*t8+(1536+141440.D0/3.D0*t8+24
     #128.D0/3.D0*t17-149504.D0/3.D0*t9+27136.D0/3.D0*t12+t835-3136*t46+
     #21376.D0/3.D0*t11-82816.D0/3.D0*t1)*t22-12416.D0/3.D0*l1+(6208*t17
     #-128.D0/3.D0*t8-74944.D0/3.D0*t1+35072.D0/3.D0+115712.D0/3.D0*t11-
     #13312*l1-54400.D0/3.D0*t9)*t21
      t882 = t4*t1
      t885 = 6*t17
      t892 = t5*t9
      t925 = 1792+(4608+8000*t8+3968*t17-13056*t1-40064.D0/3.D0*t9+40000
     #.D0/3.D0*t11-10496.D0/3.D0*l1)*t21-64.D0/3.D0*t393-2048*t8+64.D0/3
     #.D0*t882*(t383+287*t11+72-338*t8-t885+122*t9-180*l1)*t33-64.D0/3.D
     #0*t892*t62*(-t1-l1+4)*t67+(16384.D0/3.D0*t8-14080.D0/3.D0+6208*l1+
     #t814-41728.D0/3.D0*t11+20416.D0/3.D0*t1)*t20+64.D0/3.D0*l1*(144-38
     #3*t8+392*t17-188*t1-96*t12-373*t9+792*t11-t487)*t22-t675-8896.D0/3
     #.D0*t1+(-1536-3840*l1-24704.D0/3.D0*t17+128*t12-70592.D0/3.D0*t8+2
     #6240*t9-5952*t11+50176.D0/3.D0*t1)*t73+15872.D0/3.D0*t11
      t926 = s1l1**2
      t934 = 2*l1
      t939 = 64.D0/3.D0*t76
      t954 = 9472.D0/3.D0*l1
      t956 = 54*t11
      t976 = 6016.D0/3.D0*t11+64.D0/3.D0*t6*(80+106*t8-204*t1-48*t9+70*t
     #11+t934)*t22-t939-64.D0/3.D0*t4*(-40+76*t11-t885+t207-147*l1+206*t
     #9-438*t8)*t73+(-24448.D0/3.D0*t11+12160.D0/3.D0*t8+15872.D0/3.D0*l
     #1-9472.D0/3.D0-t689+6656.D0/3.D0*t1)*t20-t954-1024*t8-64.D0/3.D0*t
     #28*(-t956+16*l1-44+t265)*t33+(-25280.D0/3.D0*t1+1920*t17+2048.D0/3
     #.D0*t8-6976.D0/3.D0*l1+8192.D0/3.D0+33536.D0/3.D0*t11-5760*t9)*t21
     #+1664.D0/3.D0*t1-64.D0/3.D0*t60*t62*(-t1+2+l1)*t67+4864.D0/3.D0
      t1050 = -t939+9728.D0/3.D0-1472.D0/3.D0*t1+(384*t12-12096*t17-8940
     #8.D0/3.D0*t8+83840.D0/3.D0*t1-12032.D0/3.D0-t954+111808.D0/3.D0*t9
     #-49600.D0/3.D0*t11)*t73+(-21248.D0/3.D0-61312.D0/3.D0*t11+9216*t8-
     #384*t9+8320*t1+10368*l1)*t20+64.D0/3.D0*t28*(t385+81*t1-4*t8+t132-
     #t956-36)*t67-64.D0/3.D0*t6*(72+372*t8+t885-170*t9-t378-309*t1+92*l
     #1)*t33-2240*t8-16640.D0/3.D0*l1+64.D0/3.D0*t4*(-36+602*t11+303*t17
     #-1055*t8-290*l1+415*t1-99*t12+156*t9)*t22+15104.D0/3.D0*t11+(4160*
     #t17+576*t8+23296.D0/3.D0-39104.D0/3.D0*t9+29248*t11-65408.D0/3.D0*
     #t1-6912*l1)*t21
      t1052 = ((-t696+832.D0/3.D0*l1-256+320.D0/3.D0*t1)*s1l1+2240.D0/3.
     #D0*l1-320*t1-1280.D0/3.D0+(t696-1280.D0/3.D0-1856.D0/3.D0*t1+2752.
     #D0/3.D0*l1)*s2l1)*t692+(3840-80704.D0/3.D0*t1+24128*t17+93056.D0/3
     #.D0*t8+1408*l1-157376.D0/3.D0*t9-1856*t12-768*t46+64768.D0/3.D0*t1
     #1)*t22+(-45056.D0/3.D0*t1+2752*t17-39680.D0/3.D0*l1-33920.D0/3.D0*
     #t8+99584.D0/3.D0*t11+27776.D0/3.D0-16960.D0/3.D0*t9)*t21+(-18560*t
     #11-t106+29120.D0/3.D0*t8+21184.D0/3.D0*t1-17536.D0/3.D0+26176.D0/3
     #.D0*l1)*t20-1664*t8+t850*s1l1+12224.D0/3.D0*t11+(34112*t9-49600.D0
     #/3.D0*t17+72320.D0/3.D0*t1+t763-98944.D0/3.D0*t11-10880*t8+24448.D
     #0/3.D0*l1-24448.D0/3.D0)*t73+64.D0/3.D0*t4*(36-t820-309*t17+266*t1
     #2-223*t11+1017*t8+234*l1-438*t9-51*t46)*t33+t925*t926-8896.D0/3.D0
     #*l1+t976*t92+6016.D0/3.D0+64.D0/3.D0*t27*(36+85*t11+18*t9-t384-20*
     #t1-63*t8)*t534+(t2-448.D0/3.D0*l1+256.D0/3.D0)*t636-4352.D0/3.D0*t
     #1-64.D0/3.D0*t5*(-72-159*t11-t324-t596+327*t9+270*t1+18*t12)*t103+
     #t1050*s2l1
      t1077 = 1856.D0/3.D0*t46
      t1102 = 8128.D0/3.D0*t12
      t1107 = 448.D0/3.D0*t9
      t1123 = 2304*t17
      t1137 = 128.D0/3.D0*t393+64.D0/3.D0*l1*(-288+850*t8-56*t1-146*t12+
     #576*l1-t143-124*t9-838*t11+103*t46-65*t17)*t22-64.D0/3.D0*t892*t62
     #*(-t934+t63-5)*t67+(11200.D0/3.D0*t17-45824.D0/3.D0*t11+16064.D0/3
     #.D0*t1+69056.D0/3.D0*t8-t1102-10880*t9+17152.D0/3.D0*l1-9216+t47)*
     #t21-t1107-5120+(-82688.D0/3.D0*t1+704.D0/3.D0*t12+83584.D0/3.D0*t1
     #1-36032.D0/3.D0*t9+3072+t835+7040*t17-19712.D0/3.D0*t8+t47)*t73+t7
     #8+14528*t11-4480*t1+(11136*t9-21440.D0/3.D0*l1-14144*t8+29312.D0/3
     #.D0+28672.D0/3.D0*t1-20608.D0/3.D0*t11-t1123)*t20-6272.D0/3.D0*l1-
     #2944*t8+64.D0/3.D0*t882*(41*t9+360*l1+207*t8+t13-t552-113*t11-144-
     #274*t1)*t33
      t1194 = 5120+19520.D0/3.D0*t1+4800*t8+7040.D0/3.D0*l1-56000.D0/3.D
     #0*t11+t1107-t78+64.D0/3.D0*t442*(-288+88*t1+272*t8+t487-274*t11-16
     #7*t9+160*t17-t447+t214)*t22-64.D0/3.D0*t414*(t415-100*t9+121*t8-75
     #*t11+114*t1-t496+144)*t33+(-3072+25856*t1+87680.D0/3.D0*t8-t835-94
     #592.D0/3.D0*t11-44224.D0/3.D0*t9+2624*t17-1280.D0/3.D0*t12-t47)*t7
     #3+(9216+2176.D0/3.D0*t1-87040.D0/3.D0*t8-5632*l1+20480.D0/3.D0*t11
     #+72064.D0/3.D0*t9-25792.D0/3.D0*t17+t1102-t47)*t21+(-10240*t9-4832
     #0.D0/3.D0*t1+20416.D0/3.D0*t8+60352.D0/3.D0*t11-29312.D0/3.D0+2067
     #2.D0/3.D0*l1+t1123)*t20
      s1 = 8000.D0/3.D0*t1+64.D0/3.D0*t4*(-504+1588*t1+1339*t8-564*l1+t1
     #43-1551*t11-415*t9-233*t17+287*t12-139*t46)*t22-64.D0/3.D0*t4*(-10
     #8+1338*t1+193*t8-t662+t663-908*t11+369*t9-577*t17+501*t12-278*t46)
     #*t33+((-t673-2432+20992.D0/3.D0*t11-t675-704*t1)*s1l1+(1984.D0/3.D
     #0*t1+t680-t681-10624.D0/3.D0+3392.D0/3.D0*l1+3840*t11)*s2l1-11776.
     #D0/3.D0+6080*t11-1408.D0/3.D0*l1-t688+t689+t690)*t692+2560*t8-t695
     #*t1052-16960.D0/3.D0*t11+(-53504.D0/3.D0-45760.D0/3.D0*t1-51712.D0
     #/3.D0*t8+11264*l1+124096.D0/3.D0*t11-14720.D0/3.D0*t9+t640-t295)*t
     #20
      t1221 = s1+(73984.D0/3.D0+21120*t1-37696.D0/3.D0*t8-80768.D0/3.D0*
     #l1-29440*t11+39872*t9-65216.D0/3.D0*t17+t540-t47)*t21+(-65024.D0/3
     #.D0+20672.D0/3.D0*t1+164032.D0/3.D0*t8+65024.D0/3.D0*l1+t238-66944
     #.D0/3.D0*t11-195392.D0/3.D0*t9+101696.D0/3.D0*t17-8896*t12+t1077)*
     #t73+t1137*t926-t599*t636+(82240.D0/3.D0*t9+65024.D0/3.D0-105344.D0
     #/3.D0*t8-t1077-66304.D0/3.D0*l1+8192*t12-t238-20032*t17+48320*t11-
     #82432.D0/3.D0*t1)*t73+t1194*t926+((2432+1984*l1+t673-21184.D0/3.D0
     #*t11+2560.D0/3.D0*t1)*s1l1-6144*t11-t689+11776.D0/3.D0+(t681-3904*
     #t11-t680+10624.D0/3.D0-t690-1216*l1)*s2l1+t688-1088.D0/3.D0*t1+384
     #*l1)*t692+64.D0/3.D0*t4*(-401*t11+1446*t1-876*t8-108-t662-254*t46+
     #860*t9+t663+433*t12-544*t17)*t33
      t1224 = (l1+z-t607)**2
      taux23 = (t646+t1221)/t1224*t610
      t1 = -1+l1
      t2 = t1**2
      t3 = t2*t1
      t4 = l1**2
      t5 = t3*t4
      t6 = t4*l1
      t8 = t4**2
      t11 = t8*l1
      t12 = 54*t11
      t15 = 1-z
      t16 = t15**2
      t17 = t16**2
      t18 = t17*t16
      t21 = t2*l1
      t22 = t8*t6
      t23 = 72*t22
      t28 = t8*t4
      t30 = 360*l1
      t32 = t17*t15
      t39 = t8**2
      t40 = 302*t39
      t43 = 2880*l1
      t44 = t39*l1
      t45 = 24*t44
      t49 = 58*t22
      t53 = 24*t39
      t56 = 724*l1
      t59 = t16*t15
      t61 = 48*t22
      t66 = 604*t28
      t67 = 1635*l1
      t74 = 100*t11
      t75 = 48*t28
      t76 = 1023*l1
      t81 = 302*t8
      t83 = 24*t11
      t84 = 77*l1
      t109 = 230*t39
      t112 = 1318*l1
      t121 = 439*l1
      t129 = 460*t28
      t130 = 192*l1
      t137 = 2*l1
      t142 = 230*t8
      t144 = 58*l1
      t153 = 108*l1
      t159 = 44*t4
      t167 = t17*t59
      t178 = l1*t18
      t188 = 2340*l1
      t189 = 120*t44
      t198 = 278*t39
      t201 = 1750*l1
      t206 = 154*t22
      t212 = 1307*l1
      t221 = 556*t28
      t222 = 2154*l1
      t229 = 142*t4
      t230 = 52*t11
      t231 = 740*l1
      t237 = 278*t8
      t239 = 4*l1
      t265 = t2**2
      t266 = t265*t8
      t268 = -4+3*l1
      t269 = 2*t1
      t270 = t268*t269
      t273 = -t268*t269
      t341 = t3*t11
      t351 = t1*l1
      t363 = t2*t4
      t367 = 30*t28
      t431 = t3*(72-t153-48*t4-154*t8+48*t11+192*t6)
      t454 = dlog(z)
      t495 = t273*t167
      t556 = 460*t4
      t593 = 8*l1
      t606 = 20*t8
      t607 = 6*t11
      s2 = -128.D0/3.D0*t1*(t201-697*t22-t45+1800-8485*t4+5963*t8-4909*t
     #11+t198+1611*t28+2803*t6)*t17+128.D0/3.D0*t1*(-5486*t4+304*t39+649
     #9*t8-t189+83*t6+t188+869*t28-144*t22-4685*t11+360)*t32
      s3 = s2+(128.D0/3.D0*(t266*t270-4*t266*t273)*t167-128.D0/3.D0*t1*(
     #-24*t4+180+t144-t83+t142-390*t6)+128.D0/3.D0*t21*(-720-243*t11-t23
     #+249*t28-721*t6-152*t8-30*l1+1605*t4)*t32-128.D0/3.D0*t5*(360-31*t
     #8-359*t4-246*l1+30*t11+300*t6)*t18-128.D0/3.D0*t1*(475*t4-576+t75-
     #t74-495*t8+586*t6+t137)*t15+128.D0/3.D0*t1*(1098+1085*t11+21*t6-25
     #38*t4+t121+686*t8+t49+t53-783*t28)*t59+128.D0/3.D0*t1*(-930*t11+2*
     #t6-73*t8-1188-t61+1485*t4+t129+t130)*t16-128.D0/3.D0*t1*(t112-594*
     #t22-t45+360-2563*t4+2713*t8-808*t11+t109+362*t28-964*t6)*t17)*s2l1
      s1 = s3+128.D0/3.D0*t1*(-t212-229*t11+3437*t6+3474-6351*t4+1735*t8
     #+t206+t53-767*t28)*t59+(128.D0/3.D0*(t341*t270-4*t341*t273)*t167-1
     #28.D0/3.D0*t1*(-537*t4+558-t84-t83+t81-162*t6)+128.D0/3.D0*t351*(4
     #938*t6+1440-2652*t4-1800*l1-789*t8+69*t22-2129*t11-72*t39+1005*t28
     #)*t32-128.D0/3.D0*t363*(-262*t4+1318*t8-720-403*t11+t367-1421*t6+1
     #440*l1)*t18-128.D0/3.D0*t1*(-2034-375*t8+t75+2167*t4-t74+t76-879*t
     #6)*t15+128.D0/3.D0*t1*(-1557*t11+1810*t6+2520-7202*t4+t56+4208*t8+
     #t49+t53-435*t28)*t59+128.D0/3.D0*t1*(-642*t11-3330-2123*t6-701*t8-
     #t61+t66+4405*t4+t67)*t16-128.D0/3.D0*t1*(t43-510*t22-t45+720-6920*
     #t4+6794*t8-3931*t11+t40+653*t28+96*t6)*t17)*s1l1
      s2 = s1+128.D0/3.D0*t1*(-896*t11-2346*t6+158*t8-3276-t61+3518*t4+t
     #221+t222)*t16-128.D0/3.D0*t1*(-285*t4+396+t239-t83+t237-319*t6)
      s3 = s2+128.D0/3.D0*(-t3*(-216-120*t11-222*t6+286*t8-84*t4+324*l1)
     #+2*t431)*t4*t167
      s4 = s3-128.D0/3.D0*t2*(-720+2334*t4-1741*t6+120*t22+1813*t11-856*
     #t28-1006*t8)*t178
      s5 = s4
      s7 = -128.D0/3.D0*t1*(-1584-472*t8+t75+1522*t4-t230+t231-322*t6)*t
     #15
      s9 = -1
      s11 = t454
      s13 = 128.D0/3.D0*t1*(-684-1114*t6-180*t28+876*t11+310*l1-872*t8+1
     #644*t4)*t59+128.D0/3.D0*t1*(-72+1148*t4+138*t28+254*t6-2050*t8-354
     #*t22+1312*t11+90*t39-468*l1)*t32+128.D0/3.D0*t1*(612-160*t8-60*t11
     #+1154*t6-1052*t4-474*l1)*t16-128.D0/3.D0*t1*(-36+10*l1+78*t4-54*t6
     #)-128.D0/3.D0*t431*t4*t167
      s14 = s13+(128.D0/3.D0*t341*t495-128.D0/3.D0*t1*(-144+1484*t4-t367
     #+324*t6-2348*t8-90*t22-576*l1+1376*t11)*t17-128.D0/3.D0*t1*(-72+28
     #*l1+t229-102*t6)+128.D0/3.D0*t1*(-504-1674*t8-168*t28+1142*t11-426
     #*t6+1730*t4-116*l1)*t59+128.D0/3.D0*t1*(648+222*t8-198*t11+962*t6-
     #1238*t4-372*l1)*t16-128.D0/3.D0*t363*(144-386*t8-6*t28+104*t11+476
     #*t6-t159-288*l1)*t18+128.D0/3.D0*t351*(-288-80*t8-538*t28+970*t11-
     #1224*t6+716*t4+84*t22+t30)*t32-128.D0/3.D0*t1*(360-84*t8-228*l1-63
     #0*t4+598*t6)*t15)*s1l1-128.D0/3.D0*t1*(252-90*t8-166*l1-t556+474*t
     #6)*t15
      s12 = s14-128.D0/3.D0*t1*(-360+1906*t4-578*t28-430*t6-2322*t8+42*t
     #22-332*l1+2064*t11)*t17-128.D0/3.D0*t2*(144+470*t8+292*t28-652*t11
     #+354*t6-556*t4-t61)*t178+(128.D0/3.D0*t266*t495-128.D0/3.D0*t1*(-7
     #2+550*t4-26*t28+340*t6-986*t8+6*t22-260*l1+444*t11)*t17-128.D0/3.D
     #0*t1*(-t593+14*t4-6*t6)+128.D0/3.D0*t1*(-216+210*t11-62*t6-508*t8-
     #68*l1+632*t4)*t59+128.D0/3.D0*t1*(216+t606-t607+302*t6-t556-60*l1)
     #*t16-128.D0/3.D0*t5*(-72+122*t4-t607-112*t6+t606+44*l1)*t18+128.D0
     #/3.D0*t21*(144-394*t4-84*t11+150*t6+184*t8+t593)*t32-128.D0/3.D0*t
     #1*(72-t239-184*t4+120*t6)*t15)*s2l1
      s10 = s11*s12
      s8 = s9*s10
      s6 = s7+s8
      t636 = s5+s6
      t638 = z*l1
      t640 = (l1+z-t638)**2
      t644 = (1-l1+t638)**2
      s4 = (128.D0/3.D0*t5*(720-193*t6+75*t8+390*t4-t12-720*l1)*t18-128.
     #D0/3.D0*t21*(-1440-t23-884*t6+463*t8+1196*t4-260*t11+201*t28+t30)*
     #t32+128.D0/3.D0*t1*(720-654*t22+3536*t6-1612*t8-5992*t4+t40+727*t1
     #1+335*t28+t43-t45)*t17-128.D0/3.D0*t1*(2520+t49+2082*t6-1690*t8-39
     #06*t4+t53+1771*t11-771*t28+t56)*t59-128.D0/3.D0*t1*(-3330-t61+869*
     #t6+301*t8+85*t4-1254*t11+t66+t67)*t16+128.D0/3.D0*t1*(-2034+1145*t
     #6-507*t8-297*t4-t74+t75+t76)*t15+128.D0/3.D0*t1*(558-546*t6+t81-25
     #*t4-t83-t84))*s1l1
      s5 = (128.D0/3.D0*t5*(360-52*t6+57*t8+87*t4-t12-158*l1)*t18-128.D0
     #/3.D0*t21*(-720-t23-121*t6+30*t8+349*t4-t11+117*t28-62*l1)*t32+128
     #.D0/3.D0*t1*(360-306*t22+748*t6-515*t8-1923*t4+t109+356*t11-4*t28+
     #t112-t45)*t17-128.D0/3.D0*t1*(1098+t49+53*t6-338*t8-874*t4+t53+701
     #*t11-459*t28+t121)*t59-128.D0/3.D0*t1*(-1188-t61+570*t6+25*t8+13*t
     #4-702*t11+t129+t130)*t16+128.D0/3.D0*t1*(-576+562*t6-279*t8+123*t4
     #-t74+t75+t137)*t15+128.D0/3.D0*t1*(180-342*t6+t142-88*t4-t83+t144)
     #)*s2l1
      s3 = s4+s5
      s2 = s3-256.D0/3.D0*(t3*(106*t4-9*t11+24*t8-t153+72-55*t6)-3*t3*(-
     #10*t8+23*t6-t159+54*l1-36+3*t11))*t4*t167+128.D0/3.D0*t2*(-720-132
     #*t22-1543*t6+974*t8+1366*t4-535*t11+350*t28)*t178-128.D0/3.D0*t1*(
     #360-786*t22+3179*t6-925*t8-5006*t4+556*t39+163*t11+359*t28+t188-t1
     #89)*t32
      s1 = s2+128.D0/3.D0*t1*(1800-295*t22+4699*t6-2923*t8-6053*t4+t198+
     #2103*t11-755*t28+t201-t45)*t17-128.D0/3.D0*t1*(3474+t206+669*t6-10
     #15*t8-1711*t4+t53+2001*t11-1193*t28-t212)*t59-128.D0/3.D0*t1*(-327
     #6-t61+1502*t6-382*t8-610*t4-914*t11+t221+t222)*t16+128.D0/3.D0*t1*
     #(-1584+1158*t6-646*t8-t229-t230+t75+t231)*t15+128.D0/3.D0*t1*(396-
     #487*t6+t237-61*t4-t83+t239)+t636
      s2 = 1/t640/t644
      taux24 = s1*s2

CCC Coefficients

      tcoff45 = taux23
      tcoff46 = taux23
      tcoff47 = taux24
      tcoff48 = taux24

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff45)
      tFjet46 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff46)
      tFjet47 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff47)
      tFjet48 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff48)

CCC Sum

      qgRVnum12e1 = tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tFje
     #t48*tcoff48
      RETURN
      END

      double precision function qgRVnum12e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -1+l1
      t2 = t1**2
      t3 = l1**2
      t4 = t3**2
      t6 = t3*l1
      t8 = t4*t3
      t10 = t4*t6
      t12 = t4*l1
      t14 = 72*l1
      t18 = 1-z
      t19 = t18**2
      t20 = t19**2
      t21 = t20*t19
      t22 = l1*t21
      t25 = t2*t1
      t30 = 108*l1
      t33 = t19*t18
      t34 = t20*t33
      t35 = t3*t34
      t39 = 159*t6
      t43 = 22*t3
      t48 = 12928.D0/3.D0*t12
      t49 = t2*t12
      t50 = 3*l1
      t51 = -4+t50
      t61 = t4**2
      t62 = t61*l1
      t63 = 12*t62
      t71 = 256*t8
      t79 = 512*t62
      t84 = t1*t3
      t88 = 42*t3
      t89 = 12*t10
      t92 = t20*t18
      t97 = 3968.D0/3.D0*t10
      t104 = 256*t61
      t108 = t51*t34
      t109 = t25*t8*t108
      t110 = 64.D0/3.D0*t109
      t114 = 6080*t3
      t115 = 1024*t10
      t120 = 4864+t48-64.D0/3.D0*t49*t51*(-1+t50)*t21+64.D0/3.D0*l1*(5*t
     #10+112*t3-22*t6+195*t4+t63-636*t12-58*t61+410*t8)*t20-t71+21760.D0
     #/3.D0*t6-52928.D0/3.D0*t4+(-16448*t10-8704*t3+17024.D0/3.D0*t61+34
     #240.D0/3.D0*t6-76736.D0/3.D0*t4-t79+80960.D0/3.D0*t12+4864*t8)*t33
     #+64.D0/3.D0*t84*(240*t6+7*t8-269*t4-t88-t89+99*t12)*t92+7232*l1+(-
     #t97+62528.D0/3.D0*t8+22784.D0/3.D0*l1+6976.D0/3.D0*t6-39424*t12+31
     #36*t3+12864*t4-t104)*t19-t110+(-18944*t6-41216.D0/3.D0*l1+28096.D0
     #/3.D0*t12+t114+t115-25600.D0/3.D0*t8+19072*t4-6016.D0/3.D0)*t18
      t121 = s1l1**2
      t123 = 1280.D0/3.D0*l1
      t124 = 6*t12
      t125 = 38*t6
      t127 = 36*t3
      t133 = 55*t6
      t147 = 324*l1
      t177 = dlog(z)
      t198 = 256*t12
      t201 = 512*t12
      t211 = 2*l6-1
      t212 = t211**2
      t215 = l5-1
      t217 = z*l1
      t219 = (1-l1+t217)**2
      t220 = 1/t219
      t225 = dsqrt(z)
      t227 = dsqrt(l5)
      t229 = dsqrt(-t215)
      t232 = dsqrt(l1)
      t236 = l5*z*l1-l1+1-l5+l5*l1
      t239 = dsqrt(-t1)
      t244 = t236**2
      t246 = 4*l1*t212*z*l5*t215*t1*t220-4*t211*t225*t227*t229*t232*t236
     #*t220*t239+t244*t220
      t247 = dsqrt(t246)
      t249 = 256*t62
      t255 = 8192*t10
      t260 = 532*t3
      t261 = 6*t10
      t264 = 74*t4
      t265 = 36*l1
      t274 = 256*t4
      t278 = t25*t3
      t279 = 6*l1
      t288 = 1088*t12
      t289 = 768*t8
      t295 = 3072*l1
      t299 = t2*l1
      t301 = 36*t10
      t317 = 768*t62
      t337 = 1536*t10
      t345 = t2**2
      t346 = t345*t4
      t353 = -17088*t4-44800.D0/3.D0+64.D0/3.D0*t278*(72-t279-33*t12-46*
     #t3+128*t4-91*t6)*t21+t288+t289+(-53504.D0/3.D0+46528.D0/3.D0*t4+79
     #808.D0/3.D0*t3-6144*t10-14464.D0/3.D0*t6+29888*t8+t295-46272*t12)*
     #t19-64.D0/3.D0*t299*(-144+342*t3-t301-210*l1-352*t12+214*t6-109*t4
     #+240*t8)*t92+(-134144.D0/3.D0*t3+25216.D0/3.D0+33344.D0/3.D0*l1+23
     #68*t10+73600.D0/3.D0*t6+3520*t61+135616.D0/3.D0*t12-t317-54016.D0/
     #3.D0*t4-94976.D0/3.D0*t8)*t33+64.D0/3.D0*t1*(-t63-430*t10-623*t3+7
     #47*t12-453*t6+163*t61+72-102*t4+626*l1+42*t8)*t20+11392.D0/3.D0*l1
     #-2944.D0/3.D0*t3+(-4544*t3+60928.D0/3.D0+t337-78784.D0/3.D0*t6+355
     #84.D0/3.D0*t4-25216.D0/3.D0*l1+14848*t12-9344*t8)*t18-64.D0/3.D0*t
     #346*t51*(2-t50)*t34+82048.D0/3.D0*t6
      t359 = t2*t3
      t360 = 12*t8
      t368 = 512*t61
      t379 = t25*t4
      t380 = 4*t3
      t395 = 512*t8
      t396 = 512*l1
      t397 = t1*l1
      t399 = 144*l1
      t403 = 12*t61
      t404 = 2*t12
      t416 = 8192*l1
      t422 = t345*t12*t108
      t424 = -128.D0/3.D0*t359*(-t360+t30-72-94*t6+85*t12+t88-62*t4)*t92
     #-9728+(t368-19328.D0/3.D0*t10-103552.D0/3.D0*t12+20608*t8-42752.D0
     #/3.D0*t6-28672.D0/3.D0+12416*t3+t295+28160*t4)*t19+77440.D0/3.D0*t
     #6-36352.D0/3.D0*t4-128.D0/3.D0*t379*t51*(t380-t50-3)*t21+(36352.D0
     #/3.D0+15488.D0/3.D0*t3+t115-53248.D0/3.D0*t6-4608*t4-5376*l1+50048
     #.D0/3.D0*t12-21760.D0/3.D0*t8)*t18+2176.D0/3.D0*t12+t395-t396+128.
     #D0/3.D0*t397*(144-157*t8-t399+91*t10-170*t3+188*t6-t403-t404+71*t4
     #)*t20+(-t79+4480.D0/3.D0*t10-88960.D0/3.D0*t3+45824.D0/3.D0*t12+26
     #240*t6+7040.D0/3.D0*t61+3072-14080*t4+t416-37120.D0/3.D0*t8)*t33-1
     #4080.D0/3.D0*t3+256.D0/3.D0*t422
      t426 = 7360.D0/3.D0*t3
      t439 = 216*l1
      t451 = 256*t61*t3
      t459 = -3072+99328.D0/3.D0*t12+14912*t10-t451+10432.D0/3.D0*t62-12
     #3712.D0/3.D0*t6-12608.D0/3.D0*t4+148736.D0/3.D0*t3-15360*l1-27520.
     #D0/3.D0*t61-27776*t8
      t462 = 100*t6
      t465 = 288*l1
      t478 = 10496.D0/3.D0*t8
      t483 = 512*t10
      t491 = -192*t12-12416+64*t109+105536.D0/3.D0*t6+t424*s2l1+(-t426-7
     #1296.D0/3.D0+53248*t4-4608*t10-25664.D0/3.D0*t6+19136*t8+18112*l1-
     #153344.D0/3.D0*t12)*t19-14656.D0/3.D0*l1-64.D0/3.D0*t397*(288-125*
     #t8+1008*t6-44*t12-t439+180*t10-628*t3-413*t4-36*t61)*t92+t459*t20+
     #64.D0/3.D0*t359*(38*t12-t462-33*t8-144+92*t4+t465-154*t3)*t21+(-t2
     #49-2880*t10-130688.D0/3.D0*t3+46336.D0/3.D0*t12+166208.D0/3.D0*t6+
     #5056.D0/3.D0*t61+13824-134336.D0/3.D0*t4+4864.D0/3.D0*l1+t478)*t33
     #+t289-37120.D0/3.D0*t4-t114+(62848.D0/3.D0+14784*t3+t483-97216.D0/
     #3.D0*t6-35648.D0/3.D0*t4-8000*l1+68224.D0/3.D0*t12-20096.D0/3.D0*t
     #8)*t18
      t494 = 2176.D0/3.D0*l1
      t508 = t345*t6
      t515 = 64*t3
      t519 = 13*l1
      t528 = 24*t61
      t538 = 4224*t3
      t546 = 12*t12
      t554 = t345*t1*t4*t108
      t557 = (39872.D0/3.D0*t8-11200.D0/3.D0*t10-20992.D0/3.D0+17472*t3-
     #13376*t6+t104+17408.D0/3.D0*t4-44672.D0/3.D0*t12+2176*l1)*t19-2176
     #0.D0/3.D0*t4-21248.D0/3.D0+24832.D0/3.D0*t6+4992*l1-64.D0/3.D0*t50
     #8*t51*(t380-4*l1-3)*t21-t515+64.D0/3.D0*t299*(200*t3-192-t89+204*t
     #6-217*t4-t519-90*t12+79*t8)*t20+t71-64.D0/3.D0*t1*(98-685*t3-74*t1
     #0+t528+395*l1-139*t8-103*t4+461*t12+44*t6)*t33+(28160.D0/3.D0-t538
     #+t115-6592.D0/3.D0*t6+10432.D0/3.D0*t4-t416+17024.D0/3.D0*t12-1484
     #8.D0/3.D0*t8)*t18-64.D0/3.D0*t278*(t519+90+t6-147*t3-t546+76*t4)*t
     #92+128.D0/3.D0*t554+2624.D0/3.D0*t12
      t558 = s2l1**2
      t603 = 45824.D0/3.D0*t6
      t604 = 10880*t12
      t624 = 7680*l1
      t629 = 2*l1
      t646 = 14528*t6+64.D0/3.D0*l1*(-838*t6-288-146*t10-56*t3-124*t12+8
     #50*t4+103*t61-t63-65*t8+576*l1)*t20+(-9216+11200.D0/3.D0*t8+17152.
     #D0/3.D0*l1-8128.D0/3.D0*t10+16064.D0/3.D0*t3-t603+t104-t604+69056.
     #D0/3.D0*t4)*t19-448.D0/3.D0*t12-2944*t4+(-14144*t4+28672.D0/3.D0*t
     #3+11136*t12+29312.D0/3.D0-21440.D0/3.D0*l1-2304*t8-20608.D0/3.D0*t
     #6)*t18-5120+(3072+83584.D0/3.D0*t6-36032.D0/3.D0*t12+704.D0/3.D0*t
     #10+t104-19712.D0/3.D0*t4-82688.D0/3.D0*t3+7040*t8+t624)*t33-6272.D
     #0/3.D0*l1+t71-4480*t3-64.D0/3.D0*t49*t51*(-t629+t380-5)*t21+64.D0/
     #3.D0*t84*(41*t12+360*l1+207*t4+t89-94*t8-113*t6-144-274*t3)*t92+12
     #8.D0/3.D0*t109
      t657 = 169024.D0/3.D0*t3+20288*t4-t451+35264.D0/3.D0*t10-39296.D0/
     #3.D0*t8-10752-3584.D0/3.D0*l1+9664.D0/3.D0*t62-8960*t61-68608*t6+1
     #1264*t12
      s1 = -28160.D0/3.D0+t71+64.D0/3.D0*t1*(433*t10+1446*t3+60*t62-876*
     #t4-254*t61-544*t8-401*t6-630*l1-108+860*t12)*t92+((2432+1984*l1+17
     #92*t4-21184.D0/3.D0*t6+2560.D0/3.D0*t3)*s1l1-6144*t6-t198+11776.D0
     #/3.D0+(7808.D0/3.D0*t4-3904*t6-t201+10624.D0/3.D0-512*t3-1216*l1)*
     #s2l1+7360.D0/3.D0*t4-1088.D0/3.D0*t3+384*l1)*t247+(-t249+65024.D0/
     #3.D0-20032*t8-66304.D0/3.D0*l1-82432.D0/3.D0*t3-1856.D0/3.D0*t61+4
     #8320*t6+t255-105344.D0/3.D0*t4+82240.D0/3.D0*t12)*t33-64.D0/3.D0*t
     #2*(-t260+t261+520*t6+64*t8-t264-t265+216-135*t12)*t22-9856*t4+(198
     #4.D0/3.D0*t6-2816.D0/3.D0+1600.D0/3.D0*t3-t274)*t246+t353*s2l1+416
     #0.D0/3.D0*t12
      t659 = s1+t491*s1l1-2432.D0/3.D0*t3+t494+53056.D0/3.D0*t6+t557*t55
     #8+64.D0/3.D0*t25*(-54*t3-42*t12-108+101*t4-78*t6+162*l1)*t35+(-880
     #64.D0/3.D0*t12-73984.D0/3.D0+2816.D0/3.D0*t3-9856*t6+52928.D0/3.D0
     #*t8-15808.D0/3.D0*t10+t104+27776*l1+22528*t4)*t19+(4416*t4+t483+53
     #504.D0/3.D0+17792.D0/3.D0*t12-36352.D0/3.D0*l1-9344.D0/3.D0*t8+851
     #2.D0/3.D0*t3-48896.D0/3.D0*t6)*t18+t646*t121+t657*t20
      t678 = 58*t8
      t684 = 18*t10
      t709 = 64.D0/3.D0*t554
      t718 = -17536.D0/3.D0*l1-8320*t6-4544.D0/3.D0*t4+8320.D0/3.D0*t12-
     #t71+25600.D0/3.D0+20800.D0/3.D0*t3-64.D0/3.D0*t278*(-90+37*t3-99*t
     #6+83*l1+t546+29*t4)*t92+64.D0/3.D0*t299*(t89+234*t6-7*t12+12*t4-73
     #*l1+184-t678-238*t3)*t20-64.D0/3.D0*t1*(t684-364*t6+493*t12-366*l1
     #-94-171*t8-384*t4+809*t3)*t33-64.D0/3.D0*t508*t51*(t50+1)*t21+(-83
     #2*t8-42560.D0/3.D0*t3-54464.D0/3.D0*t4+18944.D0/3.D0*t12+2176.D0/3
     #.D0*t10-t104-4288*l1+22016.D0/3.D0+71168.D0/3.D0*t6)*t19-t709+(-58
     #24*t6-38464.D0/3.D0*t12+17536*t4+11264*l1-32512.D0/3.D0+2560*t8-10
     #240.D0/3.D0*t3)*t18
      t723 = t177**2
      t725 = 18*t12
      t727 = 20*t3
      t734 = 128*t6
      t761 = 530*t3
      t783 = 9472.D0/3.D0*l1
      t785 = 54*t6
      t805 = 6*t8
      t815 = 64.D0/3.D0*t299*(80-204*t3+70*t6+106*t4-48*t12+t629)*t20+(8
     #192.D0/3.D0+1920*t8+33536.D0/3.D0*t6+2048.D0/3.D0*t4-25280.D0/3.D0
     #*t3-5760*t12-6976.D0/3.D0*l1)*t19+4864.D0/3.D0-t783-t709+6016.D0/3
     #.D0*t6-64.D0/3.D0*t278*(-t785+16*l1-44+78*t3)*t92-1024*t4-64.D0/3.
     #D0*t508*t51*(-t3+2+l1)*t21+(6656.D0/3.D0*t3-24448.D0/3.D0*t6+12160
     #.D0/3.D0*t4-9472.D0/3.D0-t198+15872.D0/3.D0*l1)*t18-64.D0/3.D0*t1*
     #(-40-147*l1-t805+76*t6-438*t4+206*t12+346*t3)*t33+1664.D0/3.D0*t3
      t866 = 3*t12
      t873 = (-12032.D0/3.D0-12096*t8-49600.D0/3.D0*t6-89408.D0/3.D0*t4+
     #111808.D0/3.D0*t12-t783+83840.D0/3.D0*t3+384*t10)*t33+15104.D0/3.D
     #0*t6-2240*t4+9728.D0/3.D0+(-6912*l1+4160*t8+29248*t6+576*t4+23296.
     #D0/3.D0-65408.D0/3.D0*t3-39104.D0/3.D0*t12)*t19-1472.D0/3.D0*t3+64
     #.D0/3.D0*t1*(-36+415*t3+303*t8+602*t6-1055*t4-99*t10+156*t12-290*l
     #1)*t20+(-21248.D0/3.D0+8320*t3-61312.D0/3.D0*t6+9216*t4-384*t12+10
     #368*l1)*t18-16640.D0/3.D0*l1-64.D0/3.D0*t299*(72+t805-t133+372*t4-
     #309*t3-170*t12+92*l1)*t92+64.D0/3.D0*t278*(t866+81*t3-4*t4+t279-t7
     #85-36)*t21-t709
      t882 = 4352.D0/3.D0*t3
      t883 = t25*t12
      t893 = 2176*t4
      t902 = 128*t12
      t903 = 21248.D0/3.D0*l1
      t906 = 128.D0/3.D0*t422
      t908 = 58*l1
      t916 = 51*t8
      t943 = 2560-5888.D0/3.D0*t3+16000.D0/3.D0*t6-t893-128.D0/3.D0*t379
     #*t51*(-t3+2)*t21+(-48512.D0/3.D0*t6+20480.D0/3.D0*t4-15872.D0/3.D0
     #+7680*t3-t902+t903)*t18-t906-128.D0/3.D0*t359*(131*t6-t908+t866+36
     #-64*t4-47*t3)*t92-11264.D0/3.D0*l1+128.D0/3.D0*t397*(-72-t916-258*
     #t6-39*t4+158*t12+193*t3+68*l1)*t20+(-46720.D0/3.D0*t3+4096*t8+5017
     #6.D0/3.D0*t6+19456.D0/3.D0*t4-9344.D0/3.D0*l1-40192.D0/3.D0*t12+14
     #336.D0/3.D0)*t19+(-1536+53248.D0/3.D0*t3-27520.D0/3.D0*t8-18304.D0
     #/3.D0*t6-76672.D0/3.D0*t4+256*t10+28544*t12-12544.D0/3.D0*l1)*t33
      t957 = 2176*t10
      t963 = 50*t12
      t964 = 3*t8
      t997 = 3328-64.D0/3.D0*t883*t51*(1+l1)*t34+8448*t6-3264*t4+t943*s2
     #l1+(-128.D0/3.D0*t4+6208*t8+35072.D0/3.D0+115712.D0/3.D0*t6-54400.
     #D0/3.D0*t12-13312*l1-74944.D0/3.D0*t3)*t19+(-6912+110144.D0/3.D0*t
     #3-70336.D0/3.D0*t8-82816.D0/3.D0*t6-118016.D0/3.D0*t4+t957+58624*t
     #12-640.D0/3.D0*l1)*t33+64.D0/3.D0*t359*(72+270*t6+t963-t43-t399+t9
     #64-225*t4)*t21+(1536-82816.D0/3.D0*t3+24128.D0/3.D0*t8+21376.D0/3.
     #D0*t6+141440.D0/3.D0*t4+t624-149504.D0/3.D0*t12+27136.D0/3.D0*t10-
     #3136*t61)*t20-13120.D0/3.D0*t3-64.D0/3.D0*t397*(-144-423*t8-654*t6
     #-340*t4+t30+867*t12+t761+54*t10)*t92+(-87808.D0/3.D0*t6-27392.D0/3
     #.D0+40640.D0/3.D0*t4+40256.D0/3.D0*t3-896*t12+36992.D0/3.D0*l1)*t1
     #8-12416.D0/3.D0*l1
      t1002 = 24448.D0/3.D0*l1
      t1013 = 1408*l1
      t1033 = 41728.D0/3.D0*t6
      t1038 = 3968*t8
      t1062 = 24704.D0/3.D0*t8
      t1070 = 1792+64.D0/3.D0*t84*(44*t3+287*t6+72-338*t4-t805+122*t12-1
     #80*l1)*t92+15872.D0/3.D0*t6-2048*t4-6208.D0/3.D0*l1-8896.D0/3.D0*t
     #3+(6208*l1-t1033+16384.D0/3.D0*t4+20416.D0/3.D0*t3+t902-14080.D0/3
     #.D0)*t18+(4608+t1038+40000.D0/3.D0*t6+8000*t4-40064.D0/3.D0*t12-13
     #056*t3-10496.D0/3.D0*l1)*t19-64.D0/3.D0*t49*t51*(-t3-l1+4)*t21+64.
     #D0/3.D0*l1*(144-188*t3+392*t8+792*t6-383*t4-96*t10-t465-373*t12)*t
     #20-t110+(-1536+50176.D0/3.D0*t3-t1062-5952*t6-70592.D0/3.D0*t4-384
     #0*l1+26240*t12+128*t10)*t33
      t1078 = 64.D0/3.D0*t25*(36+85*t6+t725-54*l1-t727-63*t4)*t35+12224.
     #D0/3.D0*t6+((-t734+832.D0/3.D0*l1-256+320.D0/3.D0*t3)*s1l1+2240.D0
     #/3.D0*l1-320*t3-1280.D0/3.D0+(t734-1280.D0/3.D0-1856.D0/3.D0*t3+27
     #52.D0/3.D0*l1)*s2l1)*t247+(-33920.D0/3.D0*t4+2752*t8+27776.D0/3.D0
     #+99584.D0/3.D0*t6-45056.D0/3.D0*t3-16960.D0/3.D0*t12-39680.D0/3.D0
     #*l1)*t19+64.D0/3.D0*t1*(36+266*t10-309*t8-223*t6+1017*t4-438*t12+2
     #34*l1-t761-51*t61)*t92+t815*t558+(t515-448.D0/3.D0*l1+256.D0/3.D0)
     #*t246+t873*s2l1-64.D0/3.D0*t2*(-72+270*t3-132*t8-t39-t274+327*t12+
     #t684)*t22-t882+t997*s1l1+(-10880*t4-49600.D0/3.D0*t8-24448.D0/3.D0
     #-98944.D0/3.D0*t6+t1002+72320.D0/3.D0*t3+34112*t12+t957)*t33-1664*
     #t4+(3840-768*t61+24128*t8+64768.D0/3.D0*t6+93056.D0/3.D0*t4-1856*t
     #10+t1013-80704.D0/3.D0*t3-157376.D0/3.D0*t12)*t20-8896.D0/3.D0*l1+
     #6016.D0/3.D0+t1070*t121+(-17536.D0/3.D0+26176.D0/3.D0*l1-18560*t6+
     #29120.D0/3.D0*t4-t288+21184.D0/3.D0*t3)*t18
      t1081 = -64.D0/3.D0*t2*(27*t4-308*t6-96*t8-9*t10+198*t12+t14-144+2
     #10*t3)*t22+64.D0/3.D0*t25*(72+50*t3+51*t12-122*t4+84*t6-t30)*t35+6
     #4.D0/3.D0*t2*(-87*t10-144-t39+125*t4-68*t12+149*t8+t43+t14)*t22+t1
     #20*t121-t123-64.D0/3.D0*(-t25*(t124+t125-23*t4-t127)+t25*(106*t3-9
     #*t12+24*t4-t30+72-t133))*t3*t34+64.D0/3.D0*t1*(-89*t8-710*t3+586*t
     #4+171*t61+72+205*t6-456*t12+14*t10+t147-72*t62)*t92-64.D0/3.D0*t1*
     #(59*t10+72+642*t6-522*t4+205*t61-84*t62+299*t12-340*t8-610*t3+t147
     #)*t92+64.D0/3.D0*t1*(294*t10+288+322*t6-270*t4+34*t61-t63-13*t12-3
     #94*t8-369*t3+350*l1)*t20-t177*t659+t718*t558+12800.D0/3.D0*t6-1465
     #6.D0/3.D0*t4+t198+8128.D0/3.D0*t3+t723*t1078/2
      t1085 = 5248.D0/3.D0*t6-2560.D0/3.D0-2368.D0/3.D0*t3-t274-704.D0/3
     #.D0*l1
      t1093 = 7552.D0/3.D0*t3
      t1104 = 58*t10
      t1121 = 18*t3
      t1136 = 24*t12
      t1161 = -55808.D0/3.D0+14784*t6+7232*t4-8320*t12+t289+64.D0/3.D0*t
     #1*(-t1104-152+166*t6-385*t4+t403+560*t12-220*t8+834*t3-970*l1)*t33
     #-64.D0/3.D0*t1*(-97*t10+807*t6-817*t4-22*t61+t63+311*t12+98*t8+t11
     #21-340*l1)*t20-10880.D0/3.D0*t3+t783+64.D0/3.D0*t299*(t301+t125+17
     #8*t4-107*t12-t916-258*t3+224*l1)*t92-64.D0/3.D0*t278*(t1136-5*t4-5
     #1*t6+38*t3-t265)*t21+(-t483+65024.D0/3.D0+2688*t6-26176*t4+64064.D
     #0/3.D0*t12-4096.D0/3.D0*t8-36608.D0/3.D0*t3-7168.D0/3.D0*l1)*t18+(
     #4096*t10-39424.D0/3.D0-68992.D0/3.D0*t6+20544*t4-24128.D0/3.D0*t12
     #-26368.D0/3.D0*t8+115072.D0/3.D0*t3-32384.D0/3.D0*l1)*t19
      t1165 = 2*t51*t1*t34
      t1217 = -128.D0/3.D0*t359*(t360+73*t6-50*l1-75*t4+14*t3+11*t12)*t9
     #2+5376*t3-51328.D0/3.D0*t4+(t115-55808.D0/3.D0*t6-3712*t12+90496.D
     #0/3.D0*t4-21376.D0/3.D0*l1-27392.D0/3.D0-17920.D0/3.D0*t8+17024.D0
     #/3.D0*t3)*t18+21248.D0/3.D0*t12-7168.D0/3.D0*t6-128*t883*t51*t21-t
     #395-t906+(14720.D0/3.D0*t3+15872.D0/3.D0*t61-11136*t10-7936.D0/3.D
     #0*t6-2304*l1-54784.D0/3.D0*t4-t79+83584.D0/3.D0*t12-t1038)*t33+268
     #8*l1+(46208.D0/3.D0*t8-t893+t903-t368-68608.D0/3.D0*t12-1792.D0/3.
     #D0*t10+41344.D0/3.D0*t6-26240.D0/3.D0*t3+7424.D0/3.D0)*t19+12544+1
     #28.D0/3.D0*t397*(-t1104+41*t8+t908-74*t6+83*t12-t727-51*t4+t403)*t
     #20
      t1236 = 9984*l1
      t1240 = 20*t6
      t1247 = 7552*t61
      t1248 = 7552*t10
      t1270 = 12544.D0/3.D0*t6
      t1271 = -64.D0/3.D0*t883*t1165+t604+t1217*s2l1+(-t255+94208.D0/3.D
     #0*t6-168704.D0/3.D0*t12+144896.D0/3.D0*t8-50816.D0/3.D0*t4+82304.D
     #0/3.D0*l1+9728.D0/3.D0-52480.D0/3.D0*t3)*t19+(t451+25856.D0/3.D0*t
     #3+98816.D0/3.D0*t10-55552.D0/3.D0*t6-119296.D0/3.D0*t12-14336*t61+
     #116992.D0/3.D0*t4+2432.D0/3.D0*t62-t1062)*t20-t289+t1236-95488.D0/
     #3.D0*t4-64.D0/3.D0*t397*(28*t10-568*t4-t1240+1002*t12+t528-t380-t3
     #95)*t92+(-16000.D0/3.D0*t3+t1247-t1248+t603-t1236-156416.D0/3.D0*t
     #4-t317+282880.D0/3.D0*t12-45952*t8)*t33+t538+(t337-114560.D0/3.D0*
     #t6-5888*t12+180736.D0/3.D0*t4-77824.D0/3.D0*l1-34304.D0/3.D0-10496
     #*t8+44416.D0/3.D0*t3)*t18+64.D0/3.D0*t359*(238*t6-160*t4+t124-92*t
     #3+18*t8)*t21+t1270+14848
      t1344 = 35584.D0/3.D0*t3+(-5632*t10-17152.D0/3.D0*t12+138496.D0/3.
     #D0*t6-113408.D0/3.D0*t4+9088*l1+39424.D0/3.D0+62464.D0/3.D0*t8-362
     #24*t3)*t19-64.D0/3.D0*t299*(-104*t6+494*t4+24*t10-274*t12+t439-20*
     #t8-276*t3)*t92+64.D0/3.D0*t1*(t515+14*t61-340*t10+892*t6-332*l1-13
     #38*t4+t63+440*t12+502*t8)*t20+9088*t12-13184*t4-12032.D0/3.D0*l1+6
     #4.D0/3.D0*t278*(-6*t6-42*t4-44*l1+t725+t515)*t21-t289-t1033+55808.
     #D0/3.D0+(7296*t10-23424*t6+4480.D0/3.D0*t61-75904.D0/3.D0*t4+18342
     #4.D0/3.D0*t12-49792.D0/3.D0*l1-9728.D0/3.D0-41600*t8+119168.D0/3.D
     #0*t3-t249)*t33-64.D0/3.D0*t346*t1165+(t483-18176*t6-75776.D0/3.D0*
     #t12+144896.D0/3.D0*t4+4096*l1-65024.D0/3.D0+1792.D0/3.D0*t8+11776.
     #D0/3.D0*t3)*t18
      t1366 = 68416.D0/3.D0*t4+42496.D0/3.D0*t10-73472.D0/3.D0*t3+6144-6
     #080*t61+9920*t8+3712.D0/3.D0*l1+t451+20032*t6-42880*t12-2176.D0/3.
     #D0*t62
      t1395 = 1024.D0/3.D0*t4
      t1436 = 28*t3
      t1455 = 8*t4
      t1461 = -12544+4992*t6+32768.D0/3.D0*t4+(-1280.D0/3.D0*t10-7424.D0
     #/3.D0+5120.D0/3.D0*t6-640.D0/3.D0*t4+t368+7040*t12-20608.D0/3.D0*t
     #8+23168.D0/3.D0*t3-23296.D0/3.D0*l1)*t19-19712.D0/3.D0*t12+t395+12
     #8.D0/3.D0*t1*(-124*t10-26*t6-t264+t403+47*t12+174*t8+t1436-t908)*t
     #33+(-t115+27392.D0/3.D0+7552.D0/3.D0*t6-41600.D0/3.D0*t4+8192.D0/3
     #.D0*t12+16384.D0/3.D0*t8-30080.D0/3.D0*t3+t1002)*t18+t882-10112.D0
     #/3.D0*l1-128.D0/3.D0*t299*(t89+47*t6+21*t4+t404-t678+t43-t908)*t20
     #-128.D0/3.D0*t278*(-t1136+t1455+41*t6+t1121-t908)*t92
      t1504 = -14848+2496*t6+23360*t4-10624*t12+t289+64.D0/3.D0*t1*(441*
     #t10+t462-567*t4-38*t61-t63-15*t12-366*t8+t260)*t20+t1461*s2l1+(t12
     #48-99008.D0/3.D0*t6+13696.D0/3.D0*t4-t1247+t317-44224.D0/3.D0*t12+
     #17152*t8+15872*t3+30464.D0/3.D0*l1)*t33+5504.D0/3.D0*t3-32000.D0/3
     #.D0*l1+(-t337+34304.D0/3.D0+448.D0/3.D0*t6-30592*t4+5312*t12+9728*
     #t8-41344.D0/3.D0*t3+80896.D0/3.D0*l1)*t18-64.D0/3.D0*t278*(t1136+7
     #*t4+2*t6-t127)*t21-64.D0/3.D0*t299*(-t301+t964+114*t12+105*t4+16*t
     #6-244*t3)*t92+(7680*t10-9728.D0/3.D0+61888.D0/3.D0*t6+10240.D0/3.D
     #0*t4+75200.D0/3.D0*t12-34560*t8+15872.D0/3.D0*t3-84352.D0/3.D0*l1)
     #*t19
      t1552 = -25600.D0/3.D0+8512*t6-1088*t4-6784.D0/3.D0*t12+t71-3456*t
     #3+16256.D0/3.D0*l1+64.D0/3.D0*t1*(t261-94-368*t6-32*t4+t963-14*t8+
     #808*t3-386*l1)*t33+64.D0/3.D0*t1*(t89+344+274*t6-298*t4+219*t12-70
     #*t8-520*t3+183*l1)*t19-64.D0/3.D0*t278*(-t1136+t1455+t125+t1436-11
     #4*l1+94)*t92-64.D0/3.D0*t299*(t89+188+171*t6-25*t4+73*t12-82*t8-20
     #5*t3-t14)*t20+(32512.D0/3.D0-256*t6-23296.D0/3.D0*t4+30464.D0/3.D0
     #*t12-2816*t8+640*t3-31232.D0/3.D0*l1)*t18
      t1599 = -4864-3776.D0/3.D0*t6+37568.D0/3.D0*t4+(39424.D0/3.D0*t10-
     #6656*t6-26240.D0/3.D0*t4-17792.D0/3.D0*t61+t79+5888.D0/3.D0*t12-52
     #48*t8+34048.D0/3.D0*t3)*t33+(t97+40448.D0/3.D0*t6+t1395+t104+43712
     #.D0/3.D0*t12-42176.D0/3.D0*t8-30464.D0/3.D0*t3-t624)*t19-t48+t71+3
     #904.D0/3.D0*t3-7488*l1+64.D0/3.D0*t397*(46*t10-208*t6-155*t4-t403-
     #54*t12+59*t8+360*t3)*t20-64.D0/3.D0*t359*(-24*t8+32*t12+48*t4+t124
     #0-94*t3)*t92+(-22144.D0/3.D0*t4+6016.D0/3.D0-1664*t3-t115+24832.D0
     #/3.D0*t8-2944*t6-7424*t12+41984.D0/3.D0*l1)*t18
      s1 = ((-5888.D0/3.D0*t6+t201+t1013+3328.D0/3.D0*t3+3712-8576.D0/3.
     #D0*t4)*s1l1+(t893-20992.D0/3.D0*t6-t123+t1093+11648.D0/3.D0)*s2l1-
     #18688.D0/3.D0*t6+t494+7808.D0/3.D0*t3+t198+13952.D0/3.D0-256.D0/3.
     #D0*t4)*t247+t1161*s2l1+t1271*s1l1+(-10752*t4+14848.D0/3.D0*t61+140
     #480.D0/3.D0*t12+25280*t3-7360*t10-64640.D0/3.D0*t8-37120.D0/3.D0-8
     #8576.D0/3.D0*t6+7360.D0/3.D0*l1-t249)*t33+t1344*s2l1+(-20416*t3-80
     #512.D0/3.D0*t4+38400*t6+53120.D0/3.D0-3328*l1-t104-5504.D0/3.D0*t1
     #0-15168*t12+52480.D0/3.D0*t8)*t19+t1366*t20
      t1608 = s1+(-54400.D0/3.D0-6720*t12-11264.D0/3.D0*t8+15616.D0/3.D0
     #*t3+t483+93760.D0/3.D0*t4+2624*l1-56000.D0/3.D0*t6)*t18+(-t483+544
     #00.D0/3.D0-19328.D0/3.D0*t6-34688.D0/3.D0*t4+15040.D0/3.D0*t12+t47
     #8-t1093-5312.D0/3.D0*l1)*t18+(26752.D0/3.D0*t10+37120.D0/3.D0-3904
     #.D0/3.D0*t6-5248.D0/3.D0*t4-15616.D0/3.D0*t61+t249-5056*t12+10432.
     #D0/3.D0*t8-27904.D0/3.D0*t3-6080.D0/3.D0*l1)*t33+(-640*l1-13952.D0
     #/3.D0+6848*t6-11840.D0/3.D0*t3+t1395+(7744.D0/3.D0*t6-3968.D0/3.D0
     #*l1-t426-t201+9344.D0/3.D0*t4-3712)*s1l1-t198+(-1920*t4-11584.D0/3
     #.D0*t3+7616*t6-11648.D0/3.D0+t396)*s2l1)*t247+t1504*s1l1+t1552*t55
     #8+t1599*t121+(t97-53120.D0/3.D0+t1270+8640*t4+t104+1408*t12-29696.
     #D0/3.D0*t8+22208.D0/3.D0*t3+7424.D0/3.D0*l1)*t19
      t1611 = (l1+z-t217)**2
      taux23 = (t1081+t1608)/t1611*t220
      t1 = -1+l1
      t2 = t1**2
      t3 = t2*t1
      t4 = l1**2
      t5 = t3*t4
      t6 = t4*l1
      t8 = t4**2
      t10 = t8*l1
      t11 = 228*t10
      t15 = 1-z
      t16 = t15**2
      t17 = t16**2
      t18 = t17*t16
      t21 = t2*l1
      t24 = t8*t6
      t25 = 312*t24
      t27 = t8*t4
      t28 = 575*t27
      t30 = 576*l1
      t32 = t17*t15
      t38 = t8**2
      t40 = t38*l1
      t41 = 72*t40
      t45 = 4608*l1
      t51 = 410*t24
      t52 = 48*t38
      t56 = 1864*l1
      t59 = t16*t15
      t63 = 144*t24
      t67 = 1410*l1
      t73 = 164*t10
      t74 = 96*t27
      t76 = 975*l1
      t82 = 72*t10
      t84 = 124*l1
      t97 = 150*t6
      t112 = 2152*l1
      t121 = 1084*l1
      t130 = 315*l1
      t137 = 355*l1
      t144 = 253*l1
      t150 = 6*t10
      t158 = 24*t8
      t159 = 108*l1
      t166 = 44*t4
      t174 = t17*t59
      t185 = l1*t18
      t196 = 3744*l1
      t208 = 3196*l1
      t219 = 950*l1
      t229 = 2277*l1
      t238 = 599*l1
      t249 = dlog(z)
      t255 = 24*t40
      t275 = t2**2
      t276 = t275*t8
      t278 = -4+3*l1
      t279 = 2*t1
      t280 = t278*t279
      t281 = t276*t280
      t283 = -t278*t279
      t285 = 4*t276*t283
      t288 = 24*t10
      t299 = 152*t8
      t316 = 100*t10
      t318 = 48*t27
      t325 = 24*t38
      t327 = 58*t24
      t333 = 48*t24
      t363 = 2*l6-1
      t364 = t363**2
      t367 = l5-1
      t369 = z*l1
      t371 = (1-l1+t369)**2
      t372 = 1/t371
      t377 = dsqrt(z)
      t379 = dsqrt(l5)
      t381 = dsqrt(-t367)
      t384 = dsqrt(l1)
      t388 = l5*z*l1-l1+1-l5+l5*l1
      t391 = dsqrt(-t1)
      t396 = t388**2
      t399 = dsqrt(4*l1*t364*z*l5*t367*t1*t372-4*t363*t377*t379*t381*t38
     #4*t388*t372*t391+t396*t372)
      t413 = t3*t10
      t414 = t413*t280
      t416 = 4*t413*t283
      t425 = t1*l1
      t437 = t2*t4
      t442 = 30*t27
      t475 = 2880*l1
      t497 = 4*l1
      t509 = t3*(-216-120*t10-222*t6+286*t8-84*t4+324*l1)
      t515 = t3*(72-t159-48*t4-154*t8+48*t10+192*t6)
      s1 = -128.D0/3.D0*t1*(1800+2803*t6+5963*t8-697*t24-4909*t10+1750*l
     #1-t255+1611*t27+278*t38-8485*t4)*t17+128.D0/3.D0*t1*(-t63-5486*t4+
     #6499*t8+2340*l1+83*t6+360-120*t40+304*t38-4685*t10+869*t27)*t32+(1
     #28.D0/3.D0*(t281-t285)*t174-128.D0/3.D0*t1*(-t288+230*t8+58*l1-390
     #*t6-24*t4+180)+128.D0/3.D0*t21*(-721*t6-72*t24-720+1605*t4-30*l1-t
     #299-243*t10+249*t27)*t32-128.D0/3.D0*t5*(360+300*t6-31*t8-359*t4+3
     #0*t10-246*l1)*t18-128.D0/3.D0*t1*(2*l1+586*t6-576+475*t4-t316-495*
     #t8+t318)*t15+128.D0/3.D0*t1*(21*t6-2538*t4+439*l1+t325-783*t27+109
     #8+t327+1085*t10+686*t8)*t59+128.D0/3.D0*t1*(-t333+192*l1+2*t6-1188
     #+1485*t4-73*t8-930*t10+460*t27)*t16-128.D0/3.D0*t1*(360-594*t24-96
     #4*t6-2563*t4+1318*l1+230*t38-t255+2713*t8-808*t10+362*t27)*t17)*s2
     #l1+128.D0/3.D0*t1*(t158+27*l1-9*t4-86*t6+54)*t399+128.D0/3.D0*t1*(
     #3474+1735*t8-6351*t4+3437*t6+t325-1307*l1-767*t27-229*t10+154*t24)
     #*t59
      s2 = s1+(128.D0/3.D0*(t414-t416)*t174-128.D0/3.D0*t1*(-t288+302*t8
     #-77*l1-162*t6-537*t4+558)+128.D0/3.D0*t425*(1440+4938*t6+1005*t27-
     #2652*t4-1800*l1-72*t38-2129*t10+69*t24-789*t8)*t32-128.D0/3.D0*t43
     #7*(-720+1440*l1-403*t10+1318*t8-262*t4+t442-1421*t6)*t18-128.D0/3.
     #D0*t1*(-2034-375*t8+2167*t4-879*t6-t316+t318+1023*l1)*t15+128.D0/3
     #.D0*t1*(1810*t6-7202*t4+724*l1+t325-435*t27+2520+t327-1557*t10+420
     #8*t8)*t59+128.D0/3.D0*t1*(-t333-3330-2123*t6+4405*t4+604*t27-701*t
     #8+1635*l1-642*t10)*t16-128.D0/3.D0*t1*(720-510*t24+96*t6-6920*t4+t
     #475+302*t38-t255+6794*t8-3931*t10+653*t27)*t17)*s1l1+128.D0/3.D0*t
     #1*(-t333+2154*l1-2346*t6-3276+3518*t4+158*t8-896*t10+556*t27)*t16
      t540 = s2-128.D0/3.D0*t1*(-t288+278*t8+t497-319*t6-285*t4+396)+128
     #.D0/3.D0*(-t509+2*t515)*t4*t174-128.D0/3.D0*t2*(120*t24-856*t27-72
     #0-1006*t8+2334*t4-1741*t6+1813*t10)*t185-128.D0/3.D0*t1*(-1584-472
     #*t8+1522*t4-322*t6-52*t10+t318+740*l1)*t15
      s3 = (128.D0/3.D0*t5*(454*t6-199*t8+t11-1152-912*t4+1152*l1)*t18-1
     #28.D0/3.D0*t21*(1218*t6-913*t8+t25+296*t10+2304-t28-1208*t4-t30)*t
     #32+128.D0/3.D0*t1*(-5776*t6+3774*t8+1763*t24-876*t38+t41-931*t10-1
     #152-1447*t27+8752*t4-t45)*t17-128.D0/3.D0*t1*(-2994*t6+4905*t8-t51
     #-t52-4873*t10-4032+2282*t27+5468*t4-t56)*t59-128.D0/3.D0*t1*(-2711
     #*t6-1699*t8+t63+3774*t10+5616-1704*t27+167*t4-t67)*t16+128.D0/3.D0
     #*t1*(-3815*t6+1443*t8+t73+3978-t74+726*t4-t76)*t15+128.D0/3.D0*t1*
     #(1831*t6-828*t8+t82-1476+87*t4-t84))*s1l1
      s4 = (128.D0/3.D0*t5*(193*t6-271*t8+t11-576-370*t4+296*l1)*t18-128
     #.D0/3.D0*t21*(t97+188*t8+t25-299*t10+1152-419*t27-140*t4+56*l1)*t3
     #2+128.D0/3.D0*t1*(-1594*t6+1496*t8+521*t24-660*t38+t41-1010*t10-57
     #6+t28+2828*t4-t112)*t17-128.D0/3.D0*t1*(191*t6+675*t8-t51-t52-1631
     #*t10-1800+1538*t27+1104*t4-t121)*t59-128.D0/3.D0*t1*(-1728*t6+173*
     #t8+t63+1722*t10+2178-1272*t27-103*t4+t130)*t16+128.D0/3.D0*t1*(-17
     #56*t6+855*t8+t73+1494-t74-500*t4+t137)*t15+128.D0/3.D0*t1*(1021*t6
     #-612*t8+t82-774+441*t4-t144))*s2l1
      s2 = s3+s4
      s1 = s2-256.D0/3.D0*(t3*(t150+38*t6-23*t8-36*t4)-3*t3*(106*t4-9*t1
     #0+t158-t159+72-55*t6)+2*t3*(-10*t8+23*t6-t166+54*l1-36+3*t10))*t4*
     #t174+128.D0/3.D0*t2*(2508*t6-1955*t8+462*t24+1257*t10+1152-1048*t2
     #7-1876*t4)*t185-128.D0/3.D0*t1*(-5436*t6+2620*t8+2388*t24-1868*t38
     #+456*t40-558*t10-576-1124*t27+7592*t4-t196)*t32
      t542 = s1+128.D0/3.D0*t1*(-7364*t6+6680*t8+565*t24-792*t38+t41-552
     #2*t10-2880+2326*t27+8914*t4-t208)*t17-128.D0/3.D0*t1*(89*t6+3148*t
     #8-698*t24-t52-6075*t10-5724+3788*t27+2340*t4+t219)*t59-128.D0/3.D0
     #*t1*(-4740*t6+988*t8+t63+2536*t10+5886-1536*t27+1060*t4-t229)*t16+
     #128.D0/3.D0*t1*(-3878*t6+1968*t8+20*t10+3474-t74+134*t4-t238)*t15+
     #128.D0/3.D0*t1*(1557*t6-744*t8+t82-1242+312*t4-t144)-t249*t540
      t626 = 362*t24
      t643 = 188*t10
      t675 = 8*l1
      t725 = 6*t6
      t751 = t249**2
      t792 = t283*t174
      t855 = 460*t4
      t902 = 20*t8
      s1 = 128.D0/3.D0*t1*(1064*t8-t229+5886-1488*t27-3588*t4+t63+2685*t
     #10-1734*t6)*t16-128.D0/3.D0*t1*(-176*t27-768*t38-333*t10-4722*t6-2
     #880+t41+11514*t4-t208+1379*t24-1240*t8)*t17+128.D0/3.D0*t1*(t219+3
     #688*t27-5334*t10-1541*t6-5724+580*t8+7436*t4-650*t24-t52)*t59+128.
     #D0/3.D0*t1*(-2072*t6-1710*t38-3498*t27+8088*t4-4434*t8+4463*t10+29
     #70*t24-t196+432*t40-576)*t32-128.D0/3.D0*t1*(-2552*t6+44*t10-t74-1
     #794*t4+2026*t8-t238+3474)*t15
      s3 = s1
      s5 = (-128.D0/3.D0*t437*(774*t6-2304*l1+1152-1305*t8-63*t27+1248*t
     #4+570*t10)*t18-128.D0/3.D0*t1*(-804*t8-1476+1259*t6+807*t4-t84+t82
     #)+128.D0/3.D0*t425*(-1035*t24+t475+288*t38+1153*t27+2456*t4+3057*t
     #8-5610*t6-2304-938*t10)*t32+128.D0/3.D0*(3*t3*t24*t278-4*t414+t416
     #)*t174+128.D0/3.D0*t1*(2677*t27-671*t8-4032-3463*t10-t626-1740*t6-
     #t52-t56+8916*t4)*t59-128.D0/3.D0*t1*(1700*t10-852*t38+t41-1152-327
     #6*t8+1875*t24-2316*t6-t45+9696*t4-1405*t27)*t17-128.D0/3.D0*t1*(-1
     #885*t6+t643-t74-2274*t4+1651*t8-t76+3978)*t15+128.D0/3.D0*t1*(667*
     #t8+2975*t10-971*t6-4617*t4+t63-1656*t27-t67+5616)*t16)*s1l1
      s6 = (-128.D0/3.D0*t5*(-141*t6-576+45*t8+440*l1+138*t4-81*t10)*t18
     #-128.D0/3.D0*t1*(-588*t8-774+1121*t6+265*t4-t144+t82)+128.D0/3.D0*
     #t21*(1299*t10+994*t6+t675-825*t27+288*t24-861*t8+1152-1768*t4)*t32
     #+128.D0/3.D0*(3*t275*t27*t278-4*t281+t285)*t174+128.D0/3.D0*t1*(24
     #67*t27+1219*t8-1800-4470*t10-t626+689*t6-t52-t121+3064*t4)*t59-128
     #.D0/3.D0*t1*(-793*t10-636*t38+t41-576-962*t8+1455*t24+464*t6-t112+
     #3636*t4-617*t27)*t17-128.D0/3.D0*t1*(-2502*t6+1459*t8-t74-636*t4+t
     #643+1494+t137)*t15+128.D0/3.D0*t1*(427*t8+t130+2178-1224*t27-1615*
     #t4+t63+2399*t10-2250*t6)*t16)*s2l1
      s4 = s5+s6
      s2 = s3+s4
      s3 = s2+128.D0/3.D0*(-t3*(204*t4-t725+144+87*t10-216*l1-t299)+2*t5
     #09)*t4*t174
      s4 = s3-128.D0/3.D0*t1*(-720*t8-1242+1321*t6+584*t4-t144+t82)
      s5 = s4
      s7 = -128.D0/3.D0*t2*(-505*t8-2884*t4+2696*t6-707*t10+1152+421*t27
     #+33*t24)*t185
      s9 = 1.D0/2.D0
      s11 = t751
      s13 = 128.D0/3.D0*t1*(-684-872*t8-180*t27+1644*t4+310*l1-1114*t6+8
     #76*t10)*t59+128.D0/3.D0*t1*(-72-354*t24+138*t27+1312*t10+90*t38+25
     #4*t6-468*l1+1148*t4-2050*t8)*t32+128.D0/3.D0*t1*(612-60*t10-1052*t
     #4-160*t8-474*l1+1154*t6)*t16-128.D0/3.D0*t1*(-36+10*l1-54*t6+78*t4
     #)-128.D0/3.D0*t515*t4*t174
      s14 = s13+(128.D0/3.D0*t413*t792-128.D0/3.D0*t1*(-144-90*t24-t442+
     #1484*t4+1376*t10-t30-2348*t8+324*t6)*t17-128.D0/3.D0*t1*(-72+28*l1
     #-102*t6+142*t4)+128.D0/3.D0*t1*(-504+1142*t10-168*t27-116*l1-1674*
     #t8+1730*t4-426*t6)*t59+128.D0/3.D0*t1*(648-1238*t4+962*t6+222*t8-1
     #98*t10-372*l1)*t16-128.D0/3.D0*t437*(144+104*t10-6*t27-288*l1-386*
     #t8+476*t6-t166)*t18+128.D0/3.D0*t425*(-288-538*t27+360*l1+84*t24+7
     #16*t4+970*t10-80*t8-1224*t6)*t32-128.D0/3.D0*t1*(360+598*t6-228*l1
     #-84*t8-630*t4)*t15)*s1l1-128.D0/3.D0*t1*(252+474*t6-166*l1-90*t8-t
     #855)*t15
      s12 = s14-128.D0/3.D0*t1*(-360+42*t24-578*t27-332*l1+2064*t10-430*
     #t6+1906*t4-2322*t8)*t17-128.D0/3.D0*t2*(144-t333+292*t27-652*t10+4
     #70*t8-556*t4+354*t6)*t185+(128.D0/3.D0*t276*t792-128.D0/3.D0*t1*(-
     #72+6*t24-26*t27+550*t4+444*t10-260*l1-986*t8+340*t6)*t17-128.D0/3.
     #D0*t1*(-t675-t725+14*t4)+128.D0/3.D0*t1*(-216+632*t4-68*l1-508*t8-
     #62*t6+210*t10)*t59+128.D0/3.D0*t1*(216-t150-t855+t902-60*l1+302*t6
     #)*t16-128.D0/3.D0*t5*(-72-112*t6+122*t4-t150+44*l1+t902)*t18+128.D
     #0/3.D0*t21*(144+t97-394*t4+184*t8+t675-84*t10)*t32-128.D0/3.D0*t1*
     #(72-184*t4-t497+120*t6)*t15)*s2l1
      s10 = s11*s12
      s8 = s9*s10
      s6 = s7+s8
      t931 = s5+s6
      t934 = (l1+z-t369)**2
      taux24 = (t542+t931)/t934*t372

CCC Coefficients

      tcoff45 = taux23
      tcoff46 = taux23
      tcoff47 = taux24
      tcoff48 = taux24

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff45)
      tFjet46 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff46)
      tFjet47 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff47)
      tFjet48 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff48)

CCC Sum

      qgRVnum12e2 = tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tFje
     #t48*tcoff48
      RETURN
      END

      double precision function qgRVnum12e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -1+l1
      t2 = l1**2
      t3 = t2**2
      t4 = t3*t2
      t6 = t2*l1
      t7 = t3*t6
      t8 = 24*t7
      t9 = 80*t2
      t11 = t3*l1
      t12 = 50*t11
      t16 = 1-z
      t17 = t16**2
      t18 = t17*t16
      t22 = 4544.D0/3.D0*t3
      t23 = 1024.D0/3.D0*l1
      t31 = 512*l1
      t32 = t1**2
      t33 = t32*l1
      t34 = 27*t4
      t35 = 66*t6
      t38 = 21*t11
      t40 = t17**2
      t45 = 256.D0/3.D0*l1
      t49 = 2560.D0/3.D0*t4
      t50 = 6976.D0/3.D0*t6
      t53 = t32**2
      t54 = t53*t11
      t55 = 3*l1
      t56 = -4+t55
      t57 = t40*t17
      t58 = t56*t57
      t62 = 448*t11
      t63 = t32*t1
      t64 = t63*t2
      t65 = 27*t3
      t68 = t40*t16
      t72 = -64.D0/3.D0*t1*(-182*t4+t8-t9-336*t6+t12+497*t3)*t18+(1024*t
     #4+t22-t23-19712.D0/3.D0*t2-14272.D0/3.D0*t11+27392.D0/3.D0*t6)*t16
     #+5312.D0/3.D0*t3+8576.D0/3.D0*t2+t31+64.D0/3.D0*t33*(t34+t35-211*t
     #3+152*t2+t38)*t40+(16960.D0/3.D0*t2+t45-1088*t7+12288*t11-13760*t3
     #-t49-t50)*t17-64.D0/3.D0*t54*t58-16768.D0/3.D0*t6+t62-64.D0/3.D0*t
     #64*(-t9+t65+24*t6)*t68
      t73 = s2l1**2
      t75 = dlog(z)
      t76 = t75**2
      t79 = t3**2
      t83 = t79*l1
      t89 = 256*t79*t2
      t90 = -10752-39296.D0/3.D0*t4+35264.D0/3.D0*t7-8960*t79-3584.D0/3.
     #D0*l1+169024.D0/3.D0*t2+9664.D0/3.D0*t83+20288*t3-68608*t6+11264*t
     #11-t89
      t92 = 256*t4
      t100 = t40*t18
      t101 = t2*t100
      t106 = 74*t3
      t113 = l1*t57
      t119 = t32*t11
      t120 = 2*l1
      t121 = 4*t2
      t127 = 12*t83
      t140 = 69056.D0/3.D0*t3
      t141 = 45824.D0/3.D0*t6
      t145 = 256*t79
      t147 = 10880*t11
      t159 = 7680*l1
      t161 = 19712.D0/3.D0*t3
      t167 = t1*t2
      t168 = 41*t11
      t171 = 12*t7
      t172 = 94*t4
      t179 = 14528*t6
      t181 = t63*t4
      t182 = t56*t100
      t183 = t181*t182
      t185 = -448.D0/3.D0*t11+t92-2944*t3-4480*t2-64.D0/3.D0*t119*t56*(-
     #t120+t121-5)*t57+64.D0/3.D0*l1*(-288-t127-65*t4-56*t2+576*l1+850*t
     #3-838*t6+103*t79-146*t7-124*t11)*t40-5120+(-9216+t140-t141-8128.D0
     #/3.D0*t7+17152.D0/3.D0*l1+16064.D0/3.D0*t2+t145+11200.D0/3.D0*t4-t
     #147)*t17+(-14144*t3+29312.D0/3.D0-2304*t4+28672.D0/3.D0*t2-20608.D
     #0/3.D0*t6-21440.D0/3.D0*l1+11136*t11)*t16+(704.D0/3.D0*t7+t159+704
     #0*t4+t145-t161+3072-82688.D0/3.D0*t2+83584.D0/3.D0*t6-36032.D0/3.D
     #0*t11)*t18+64.D0/3.D0*t167*(t168+360*l1+207*t3+t171-t172-113*t6-14
     #4-274*t2)*t68+t179-6272.D0/3.D0*l1+128.D0/3.D0*t183
      t186 = s1l1**2
      t206 = -3072-27776*t4-t89-27520.D0/3.D0*t79-12608.D0/3.D0*t3+14873
     #6.D0/3.D0*t2+10432.D0/3.D0*t83-15360*l1+14912*t7-123712.D0/3.D0*t6
     #+99328.D0/3.D0*t11
      t208 = 192*t11
      t210 = l1*t1
      t213 = 216*l1
      t228 = 7360.D0/3.D0*t2
      t238 = 512*t7
      t244 = t63*t3
      t252 = 8192*l1
      t256 = 512*t83
      t261 = 512*t4
      t262 = t32*t2
      t263 = 12*t4
      t264 = 108*l1
      t265 = 94*t6
      t267 = 42*t2
      t277 = 3072*l1
      t279 = 512*t79
      t290 = 1024*t7
      t295 = 71*t3
      t298 = 144*l1
      t300 = 12*t79
      t307 = t54*t182
      t309 = 2176.D0/3.D0*t11-128.D0/3.D0*t244*t56*(t121-t55-3)*t57-9728
     #+(45824.D0/3.D0*t11+26240*t6+3072+t252-37120.D0/3.D0*t4+4480.D0/3.
     #D0*t7-88960.D0/3.D0*t2-t256-14080*t3+7040.D0/3.D0*t79)*t18+t261-12
     #8.D0/3.D0*t262*(-t263+t264-72-t265+85*t11+t267-62*t3)*t68-36352.D0
     #/3.D0*t3-14080.D0/3.D0*t2+77440.D0/3.D0*t6+(20608*t4-28672.D0/3.D0
     #+t277-42752.D0/3.D0*t6+t279+12416*t2-19328.D0/3.D0*t7-103552.D0/3.
     #D0*t11+28160*t3)*t17+(-21760.D0/3.D0*t4+15488.D0/3.D0*t2-53248.D0/
     #3.D0*t6-5376*l1+t290+50048.D0/3.D0*t11+36352.D0/3.D0-4608*t3)*t16+
     #128.D0/3.D0*t210*(144+t295+188*t6+91*t7-t298-170*t2-t300-157*t4-2*
     #t11)*t40-t31+256.D0/3.D0*t307
      t311 = 6080*t2
      t318 = 256*t83
      t328 = 288*l1
      t334 = 768*t4
      t335 = t206*t40-t208+105536.D0/3.D0*t6-64.D0/3.D0*t210*(180*t7-125
     #*t4-t213-413*t3+288-36*t79+1008*t6-628*t2-44*t11)*t68-12416+(19136
     #*t4-71296.D0/3.D0-4608*t7+18112*l1-153344.D0/3.D0*t11+53248*t3-t22
     #8-25664.D0/3.D0*t6)*t17-14656.D0/3.D0*l1+64*t183+(-20096.D0/3.D0*t
     #4+14784*t2-97216.D0/3.D0*t6-8000*l1+t238+68224.D0/3.D0*t11+62848.D
     #0/3.D0-35648.D0/3.D0*t3)*t16+t309*s2l1-t311+(46336.D0/3.D0*t11+166
     #208.D0/3.D0*t6+13824+4864.D0/3.D0*l1+10496.D0/3.D0*t4-2880*t7-1306
     #88.D0/3.D0*t2-t318-134336.D0/3.D0*t3+5056.D0/3.D0*t79)*t18-37120.D
     #0/3.D0*t3+64.D0/3.D0*t262*(38*t11-100*t6-33*t4-144+92*t3+t328-154*
     #t2)*t57+t334
      t344 = 256*t11
      t347 = 512*t11
      t348 = 512*t2
      t353 = 1088.D0/3.D0*t2
      t357 = 2*l6-1
      t358 = t357**2
      t361 = l5-1
      t363 = z*l1
      t365 = (1-l1+t363)**2
      t366 = 1/t365
      t371 = dsqrt(z)
      t373 = dsqrt(l5)
      t375 = dsqrt(-t361)
      t378 = dsqrt(l1)
      t382 = l5*z*l1-l1+1-l5+l5*l1
      t385 = dsqrt(-t1)
      t390 = t382**2
      t392 = 4*l1*t358*z*l5*t361*t1*t366-4*t357*t371*t373*t375*t378*t382
     #*t366*t385+t390*t366
      t393 = dsqrt(t392)
      t397 = 256*t3
      t404 = 36*t7
      t434 = 6*l1
      t437 = 128*t3
      t449 = 1536*t7
      t452 = t53*t3
      t463 = 768*t83
      t470 = 1088*t11
      t473 = -64.D0/3.D0*t33*(-144-109*t3+214*t6-t404-210*l1+342*t2-352*
     #t11+240*t4)*t68-44800.D0/3.D0+64.D0/3.D0*t1*(747*t11-453*t6+72+626
     #*l1+42*t4-430*t7-623*t2-t127-102*t3+163*t79)*t40+(29888*t4-46272*t
     #11-53504.D0/3.D0+79808.D0/3.D0*t2-14464.D0/3.D0*t6+t277-6144*t7+46
     #528.D0/3.D0*t3)*t17+t334-17088*t3+64.D0/3.D0*t64*(72-t434-33*t11-4
     #6*t2+t437-91*t6)*t57+(60928.D0/3.D0-9344*t4-78784.D0/3.D0*t6-25216
     #.D0/3.D0*l1-4544*t2+14848*t11+35584.D0/3.D0*t3+t449)*t16-64.D0/3.D
     #0*t452*t56*(2-t55)*t100+(135616.D0/3.D0*t11-94976.D0/3.D0*t4+25216
     #.D0/3.D0-134144.D0/3.D0*t2+2368*t7+73600.D0/3.D0*t6-t463-54016.D0/
     #3.D0*t3+3520*t79+33344.D0/3.D0*l1)*t18+11392.D0/3.D0*l1+t470-2944.
     #D0/3.D0*t2+82048.D0/3.D0*t6
      t483 = 2432.D0/3.D0*t2
      t485 = t53*t6
      t486 = 4*l1
      t493 = 64*t2
      t496 = 4224*t2
      t505 = 24*t79
      t524 = 13*l1
      t535 = 12*t11
      t543 = t53*t1*t3*t182
      t545 = 4992*l1-64.D0/3.D0*t485*t56*(t121-t486-3)*t57+2624.D0/3.D0*
     #t11-t493-21248.D0/3.D0+(-6592.D0/3.D0*t6-14848.D0/3.D0*t4-t252-t49
     #6+28160.D0/3.D0+10432.D0/3.D0*t3+t290+17024.D0/3.D0*t11)*t16-64.D0
     #/3.D0*t1*(98-139*t4-74*t7+395*l1+461*t11+t505-103*t3-685*t2+44*t6)
     #*t18+(-44672.D0/3.D0*t11+39872.D0/3.D0*t4-13376*t6+17472*t2+t145+2
     #176*l1-11200.D0/3.D0*t7-20992.D0/3.D0+17408.D0/3.D0*t3)*t17+24832.
     #D0/3.D0*t6+64.D0/3.D0*t33*(-t171-192-217*t3-t524+79*t4+200*t2+204*
     #t6-90*t11)*t40-21760.D0/3.D0*t3-64.D0/3.D0*t64*(t524+90+t6-147*t2-
     #t535+76*t3)*t68+t92+128.D0/3.D0*t543
      t568 = 8192*t7
      s1 = -28160.D0/3.D0+t90*t40+t92+64.D0/3.D0*t63*(-54*t2-42*t11-108+
     #101*t3-78*t6+162*l1)*t101-64.D0/3.D0*t32*(64*t4-135*t11+216-t106-3
     #6*l1+520*t6+6*t7-532*t2)*t113+t185*t186+(-88064.D0/3.D0*t11+52928.
     #D0/3.D0*t4-73984.D0/3.D0-9856*t6+t145-15808.D0/3.D0*t7+22528*t3+28
     #16.D0/3.D0*t2+27776*l1)*t17+t335*s1l1+((2432+1984*l1+1792*t3-21184
     #.D0/3.D0*t6+2560.D0/3.D0*t2)*s1l1-6144*t6-t344+11776.D0/3.D0+(7808
     #.D0/3.D0*t3-3904*t6-t347+10624.D0/3.D0-t348-1216*l1)*s2l1+7360.D0/
     #3.D0*t3-t353+384*l1)*t393+(1984.D0/3.D0*t6-2816.D0/3.D0+1600.D0/3.
     #D0*t2-t397)*t392
      t572 = s1+2176.D0/3.D0*l1+4160.D0/3.D0*t11+t473*s2l1+(-48896.D0/3.
     #D0*t6+53504.D0/3.D0+t238+4416*t3+8512.D0/3.D0*t2-36352.D0/3.D0*l1+
     #17792.D0/3.D0*t11-9344.D0/3.D0*t4)*t16-t483+t545*t73+53056.D0/3.D0
     #*t6+64.D0/3.D0*t1*(1446*t2-108+60*t83+860*t11-544*t4-401*t6-254*t7
     #9-876*t3+433*t7-630*l1)*t68+(-1856.D0/3.D0*t79+48320*t6-105344.D0/
     #3.D0*t3-82432.D0/3.D0*t2+65024.D0/3.D0-66304.D0/3.D0*l1-t318-20032
     #*t4+82240.D0/3.D0*t11+t568)*t18-9856*t3
      t579 = 8832*l1
      t586 = 448.D0/3.D0*t2-t45-64*t6
      t600 = 20032*t6+6144+3712.D0/3.D0*l1+42496.D0/3.D0*t7+68416.D0/3.D
     #0*t3-73472.D0/3.D0*t2+t89+9920*t4-2176.D0/3.D0*t83-6080*t79-42880*
     #t11
      t660 = 12928.D0/3.D0*t11
      t670 = 64.D0/3.D0*t183
      t671 = 64.D0/3.D0*l1*(t127+112*t2+195*t3-22*t6-58*t79+410*t4-636*t
     #11+5*t7)*t40-64.D0/3.D0*t119*t56*(-1+t55)*t57+7232*l1+(62528.D0/3.
     #D0*t4+3136*t2-t145+t50-3968.D0/3.D0*t7+22784.D0/3.D0*l1-39424*t11+
     #12864*t3)*t17+(t311-25600.D0/3.D0*t4-18944*t6-6016.D0/3.D0+28096.D
     #0/3.D0*t11+t290-41216.D0/3.D0*l1+19072*t3)*t16+21760.D0/3.D0*t6+64
     #.D0/3.D0*t167*(240*t6+7*t4-269*t3-t267-t171+99*t11)*t68-52928.D0/3
     #.D0*t3+t660+(-t256+80960.D0/3.D0*t11+34240.D0/3.D0*t6+17024.D0/3.D
     #0*t79-76736.D0/3.D0*t3+4864*t4-16448*t7-8704*t2)*t18-t92-t670+4864
      t687 = 9*t7
      t693 = 832.D0/3.D0*l1
      t711 = 1920*t3
      t726 = 41728.D0/3.D0*t6
      t739 = 2*t56*t1*t100
      t745 = 104*t6
      t786 = 18*t11
      t791 = -t334-t726+35584.D0/3.D0*t2+(138496.D0/3.D0*t6+39424.D0/3.D
     #0+62464.D0/3.D0*t4-17152.D0/3.D0*t11+9088*l1-5632*t7-36224*t2-1134
     #08.D0/3.D0*t3)*t17-64.D0/3.D0*t452*t739-64.D0/3.D0*t33*(-20*t4-276
     #*t2+494*t3+t213-t745+t8-274*t11)*t68-13184*t3+(-18176*t6+1792.D0/3
     #.D0*t4+4096*l1+11776.D0/3.D0*t2+t238-65024.D0/3.D0-75776.D0/3.D0*t
     #11+144896.D0/3.D0*t3)*t16+(-9728.D0/3.D0-41600*t4-23424*t6-75904.D
     #0/3.D0*t3+4480.D0/3.D0*t79+119168.D0/3.D0*t2-t318-49792.D0/3.D0*l1
     #+7296*t7+183424.D0/3.D0*t11)*t18+9088*t11+64.D0/3.D0*t1*(-332*l1-1
     #338*t3+14*t79+t127+t493-340*t7+502*t4+440*t11+892*t6)*t40-12032.D0
     #/3.D0*l1+64.D0/3.D0*t64*(-6*t6-42*t3-44*l1+t786+t493)*t57+55808.D0
     #/3.D0
      t806 = 18*t7
      t812 = 64.D0/3.D0*t543
      t841 = 12*t3
      t849 = (-32512.D0/3.D0+11264*l1+17536*t3-38464.D0/3.D0*t11-5824*t6
     #-10240.D0/3.D0*t2+2560*t4)*t16-t92-t22-64.D0/3.D0*t1*(-171*t4-364*
     #t6+809*t2-366*l1-94+493*t11+t806-384*t3)*t18-t812-64.D0/3.D0*t485*
     #t56*(t55+1)*t57-64.D0/3.D0*t64*(-90+37*t2-99*t6+83*l1+t535+29*t3)*
     #t68+(22016.D0/3.D0-4288*l1-832*t4-54464.D0/3.D0*t3+2176.D0/3.D0*t7
     #+71168.D0/3.D0*t6-t145+18944.D0/3.D0*t11-42560.D0/3.D0*t2)*t17-175
     #36.D0/3.D0*l1+64.D0/3.D0*t33*(-58*t4+234*t6+184-7*t11-238*t2-73*l1
     #+t841+t171)*t40+20800.D0/3.D0*t2+25600.D0/3.D0-8320*t6+8320.D0/3.D
     #0*t11
      t872 = 9984*l1
      t891 = t63*t11
      t894 = 20*t6
      t913 = 20*t2
      t916 = 58*l1
      t928 = 3968*t4
      t930 = 7936.D0/3.D0*t6
      t933 = 7168.D0/3.D0*t6
      t935 = 21248.D0/3.D0*l1
      t938 = 2176*t3
      t956 = 128.D0/3.D0*t307
      t957 = -51328.D0/3.D0*t3-128.D0/3.D0*t262*(t263+73*t6-50*l1-75*t3+
     #14*t2+11*t11)*t68+128.D0/3.D0*t210*(41*t4-t913+t300-74*t6-58*t7+t9
     #16+83*t11-51*t3)*t40+(-2304*l1-54784.D0/3.D0*t3+15872.D0/3.D0*t79-
     #t256+14720.D0/3.D0*t2-11136*t7-t928+83584.D0/3.D0*t11-t930)*t18-t9
     #33+(-t279+46208.D0/3.D0*t4+t935+7424.D0/3.D0-1792.D0/3.D0*t7-26240
     #.D0/3.D0*t2-t938+41344.D0/3.D0*t6-68608.D0/3.D0*t11)*t17+(-17920.D
     #0/3.D0*t4-27392.D0/3.D0-21376.D0/3.D0*l1-55808.D0/3.D0*t6+t290+904
     #96.D0/3.D0*t3+17024.D0/3.D0*t2-3712*t11)*t16+2688*l1-t261+21248.D0
     #/3.D0*t11-128*t891*t58+5376*t2-t956+12544
      t961 = 6*t11
      t974 = 24704.D0/3.D0*t4
      t978 = (-10496*t4-34304.D0/3.D0-77824.D0/3.D0*l1-114560.D0/3.D0*t6
     #+t449+180736.D0/3.D0*t3+44416.D0/3.D0*t2-5888*t11)*t16+(-t872-1564
     #16.D0/3.D0*t3+7552*t79-t463-16000.D0/3.D0*t2-7552*t7-45952*t4+2828
     #80.D0/3.D0*t11+t141)*t18+12544.D0/3.D0*t6-95488.D0/3.D0*t3+(144896
     #.D0/3.D0*t4+9728.D0/3.D0+94208.D0/3.D0*t6-52480.D0/3.D0*t2+82304.D
     #0/3.D0*l1-t568-168704.D0/3.D0*t11-50816.D0/3.D0*t3)*t17+14848-64.D
     #0/3.D0*t891*t739-64.D0/3.D0*t210*(-t894+t505+1002*t11-t261-568*t3+
     #28*t7-t121)*t68+t957*s2l1+t872+t496+64.D0/3.D0*t262*(238*t6-160*t3
     #+t961-92*t2+18*t4)*t57-t334+(2432.D0/3.D0*t83-55552.D0/3.D0*t6-143
     #36*t79-119296.D0/3.D0*t11+116992.D0/3.D0*t3+t89+25856.D0/3.D0*t2-t
     #974+98816.D0/3.D0*t7)*t40+t147
      s1 = 35968.D0/3.D0+4736*t2+t600*t40-t92+(2368.D0/3.D0*t2-5248.D0/3
     #.D0*t6+2560.D0/3.D0+t397+704.D0/3.D0*l1)*t392-13184.D0/3.D0*t6-310
     #40.D0/3.D0*t3+(-3328*l1+38400*t6-80512.D0/3.D0*t3+53120.D0/3.D0-20
     #416*t2+52480.D0/3.D0*t4-5504.D0/3.D0*t7-15168*t11-t145)*t17+t671*t
     #186+(-10752*t3+25280*t2+14848.D0/3.D0*t79+140480.D0/3.D0*t11-t318-
     #7360*t7+7360.D0/3.D0*l1-64640.D0/3.D0*t4-37120.D0/3.D0-88576.D0/3.
     #D0*t6)*t18
      t980 = s1-64.D0/3.D0*t32*(210*t2+72*l1-96*t4-308*t6+t65-144-t687+1
     #98*t11)*t113-t693+64.D0/3.D0*t63*(72+50*t2+51*t11-122*t3+84*t6-t26
     #4)*t101+(-640*l1-13952.D0/3.D0+6848*t6-11840.D0/3.D0*t2+1024.D0/3.
     #D0*t3+(7744.D0/3.D0*t6-3968.D0/3.D0*l1-t228-t347+9344.D0/3.D0*t3-3
     #712)*s1l1-t344+(-t711-11584.D0/3.D0*t2+7616*t6-11648.D0/3.D0+t31)*
     #s2l1)*t393+(15616.D0/3.D0*t2+93760.D0/3.D0*t3-6720*t11-56000.D0/3.
     #D0*t6+2624*l1-54400.D0/3.D0+t238-11264.D0/3.D0*t4)*t16+t791*s2l1+t
     #849*t73+64.D0/3.D0*t1*(205*t6+171*t79-710*t2-89*t4+14*t7+324*l1+58
     #6*t3-456*t11-72*t83+72)*t68+t660+t978*s1l1
      t986 = 133*t11
      t988 = 338*t3
      t1006 = 3584*t11
      t1010 = 6016.D0/3.D0*t6
      t1018 = 52*t2
      t1021 = 9*t11
      t1034 = 51*t4
      t1037 = 108*t6
      t1048 = 14720.D0/3.D0*l1
      t1055 = 128*t4
      t1056 = 4352.D0/3.D0*t2
      t1057 = 4480.D0/3.D0*l1
      t1062 = 1408.D0/3.D0*t2
      t1064 = 512.D0/3.D0*l1
      t1068 = 256*t2
      t1072 = 64*t3
      t1075 = 1024*l1
      t1084 = 76*t6
      t1091 = 8448*t6
      t1093 = 1408*t4
      t1095 = 2048.D0/3.D0*l1
      t1108 = 5888.D0/3.D0*t2
      t1120 = 51*t79
      t1127 = t1075+(t1064-2048*t7+5120.D0/3.D0*t4+11776.D0/3.D0*t6+4940
     #8.D0/3.D0*t11+20096.D0/3.D0*t2-80768.D0/3.D0*t3)*t17+64.D0/3.D0*t6
     #4*(-t1084+41*t3+t535)*t57+(22528.D0/3.D0*t3+t1091-24064.D0/3.D0*t2
     #+t1093-25984.D0/3.D0*t11-t1095)*t16-22144.D0/3.D0*t6-64.D0/3.D0*t3
     #3*(228*t6+t171-t295-203*t11+80*t4)*t68+6080.D0/3.D0*t3+10496.D0/3.
     #D0*t2+(-t1108-320*t79+25280.D0/3.D0*t7+t140-59072.D0/3.D0*t4-11584
     #*t6+2112*t11)*t18+832*t11+64.D0/3.D0*t1*(-345*t4-122*t7+868*t11+t1
     #120-239*t3-236*t6)*t40
      t1129 = 64*t79
      s1 = t72*t73+t76*t572/2+(59584.D0/3.D0*t3+19264.D0/3.D0*t4-62528.D
     #0/3.D0*t11+5760*t6+t579-20480*t2+192*t79-t238+768)*t18-t75*t980+64
     #.D0/3.D0*t1*(452*t6-141*t4-2*t79-224*t2+t986+114*t7-t988-t127)*t68
     #-64.D0/3.D0*t1*(-t300+36*t4-288*t2+487*t3+200*t6+88*t7-521*t11)*t6
     #8+(-t483-t31+(-t334+7424.D0/3.D0*t2+t23-8704.D0/3.D0*t3-8192.D0/3.
     #D0*t6+t1006)*t16-704.D0/3.D0*t3+t1010-64.D0/3.D0*t64*(t535-40*t3+1
     #4*t6+24*t2)*t68-64.D0/3.D0*t33*(t34+t1018-111*t3+61*t6-t1021)*t40+
     #64.D0/3.D0*t1*(-192*t4-32*t2+t437-136*t6+t404+186*t11)*t18+64.D0/3
     #.D0*t1*(t1034+109*t2+t486-224*t3+t1037-8*t11)*t17-t62)*t73
      t1138 = s1+4352.D0/3.D0*t3+t344+(13568.D0/3.D0*t4-11712*t2+t1048+7
     #744.D0/3.D0*t3+11648.D0/3.D0*t6+2304-1472*t7-15040.D0/3.D0*t11)*t1
     #7+(-t1055+t1056+t1057+8768.D0/3.D0*t3-t50-2304-3328.D0/3.D0*t11)*t
     #16+((-t1062+1280.D0/3.D0*t6+t1064-t437)*s1l1+(-1024.D0/3.D0*t6+t23
     #-t1068+t397)*s2l1+t23-t348+320.D0/3.D0*t6+t1072)*t393+t1127*s2l1+(
     #t1129-4864.D0/3.D0*t4+31744.D0/3.D0*t2-t579-41024.D0/3.D0*t3+25984
     #.D0/3.D0*t6-768+2240.D0/3.D0*t7+14528.D0/3.D0*t11)*t18
      t1139 = 960*t7
      t1144 = 63872.D0/3.D0*t6
      t1163 = 180*l1
      t1196 = 1024.D0/3.D0*t2
      t1197 = 13312.D0/3.D0*t3
      t1198 = 3328.D0/3.D0*t6
      t1199 = 128*t7
      t1203 = 1664.D0/3.D0*t2
      t1214 = 27*t7
      t1220 = 3*t4
      t1223 = 17*t11
      t1228 = 28*t3
      t1234 = 128*t11
      t1238 = 832*t3
      t1252 = 9*t4
      t1262 = 128*t2
      t1268 = ((2944.D0/3.D0*t4-t1196+t1197-t1198+t1199-4096*t11)*t17+(-
     #t334+t1203+t1064-4736.D0/3.D0*t3-2048.D0/3.D0*t6+2304*t11)*t16+640
     #.D0/3.D0*t2-t31-2432.D0/3.D0*t3+2944.D0/3.D0*t6+128.D0/3.D0*t1*(-9
     #3*t4-5*t3-t894+t1214+91*t11)*t18-128.D0/3.D0*t33*(t1220-52*t3+32*t
     #6+t1223)*t40-128.D0/3.D0*t64*(t535-t1228+16*t6)*t68+t1234)*s2l1+(-
     #t1129-10304.D0/3.D0*t4-t1238-256*t6+3520.D0/3.D0*t7+10240.D0/3.D0*
     #t11)*t18+t483-t31-3328.D0/3.D0*t3+1088.D0/3.D0*t6+64.D0/3.D0*t1*(9
     #*t79+5*t4-t1228-t1214+t168)*t40-64.D0/3.D0*t33*(t1252-t38+t841)*t6
     #8+(4736.D0/3.D0*t4-t1068-t711+1728*t6-t238-1856.D0/3.D0*t11)*t17+t
     #62+(t1055+t1262+t1064+3328*t3-2240*t6-4544.D0/3.D0*t11)*t16
      t1270 = 512.D0/3.D0*t6
      t1279 = 320*t3
      t1286 = 640*t3
      t1287 = 2560.D0/3.D0*t6
      t1294 = 16*t3
      t1308 = 2048.D0/3.D0*t3
      t1315 = 15*t7
      t1328 = 1280*t4
      t1344 = 384*t11
      t1351 = 128.D0/3.D0*t210*(4*t6+t1315-t1223-23*t4+25*t3)*t40+(-1152
     #*t7+t1196+8192.D0/3.D0*t11-t1197+t1093+t1198)*t17+(12416.D0/3.D0*t
     #3+512*t6+t1328-3712.D0/3.D0*t2-13568.D0/3.D0*t11-t1064)*t16+(-1868
     #8.D0/3.D0*t4+11392.D0/3.D0*t7+2944*t11-640*t79+896.D0/3.D0*t3-t127
     #0)*t18-2688*t6-128.D0/3.D0*t181*t58+1152*t2+t1286+t31+t1344-128.D0
     #/3.D0*t262*(-8*t6-6*t3+t1021)*t68
      t1355 = 640*t4
      t1361 = 4*t3
      t1367 = t31+(-16000.D0/3.D0*t4+5696.D0/3.D0*t7+t1006+t1129+t1286-t
     #1287)*t18+1792.D0/3.D0*t3-64.D0/3.D0*t210*(t300-4*t7+72*t11-61*t4-
     #t1294)*t68+(5632.D0/3.D0*t4-t290+1792.D0/3.D0*t2+3776.D0/3.D0*t6+4
     #672.D0/3.D0*t11-12800.D0/3.D0*t3)*t17+896*t2+(576*t83+t1308+7616.D
     #0/3.D0*t7-2432*t79+t49-6656.D0/3.D0*t11)*t40+t1351*s2l1+(8384.D0/3
     #.D0*t3+3392.D0/3.D0*t6+t1355-4480.D0/3.D0*t2-8704.D0/3.D0*t11-t106
     #4)*t16-t50+320*t11+64.D0/3.D0*t262*(t263+t1361-19*t11)*t57
      t1406 = -64.D0/3.D0*t32*t7*t58+(1792.D0/3.D0*t6-64*t7+t347+256.D0/
     #3.D0*t4-3392.D0/3.D0*t3)*t17+64.D0/3.D0*t167*(t1294+t1252-24*t11)*
     #t68+1600.D0/3.D0*t3-2752.D0/3.D0*t6+(1664.D0/3.D0*t7-128*t79+1536*
     #t11-t1328-t1308)*t18+64.D0/3.D0*l1*(3*t79-19*t7+53*t4+t1228-65*t11
     #)*t40+t353+t45-64*t11+(t92+t1270+t1238-t1196-2752.D0/3.D0*t11)*t16
      t1413 = 3840*l1
      t1420 = 256*t7
      t1434 = 320*t2
      t1448 = 576*t11
      t1465 = 44*t2
      t1495 = (-2816*t4-2944*t2-t1064+12480*t3-t933+t449-5696*t11)*t17-1
     #280.D0/3.D0*t2-t1075+(-t1355+8704.D0/3.D0*t2+t1095-4800*t3-t930+13
     #504.D0/3.D0*t11)*t16-t1279+7040.D0/3.D0*t6+64.D0/3.D0*t1*(-254*t4-
     #t1465+119*t3-291*t6+t1315+405*t11)*t18-64.D0/3.D0*t1*(39*t79-169*t
     #4-162*t3-t745-94*t7+480*t11)*t40+64.D0/3.D0*t33*(69*t4-17*t3+t1037
     #-140*t11)*t68-64.D0/3.D0*t64*(30*t3-40*t6)*t57-t1448
      t1535 = 6*t4
      t1542 = 1792+(-1536+t1199-t974-5952*t6-70592.D0/3.D0*t3-t1413+2624
     #0*t11+50176.D0/3.D0*t2)*t18-8896.D0/3.D0*t2+64.D0/3.D0*l1*(144-383
     #*t3+392*t4-t328+792*t6-188*t2-373*t11-96*t7)*t40+15872.D0/3.D0*t6-
     #2048*t3-t670-64.D0/3.D0*t119*t56*(-t2-l1+4)*t57+(-t726+16384.D0/3.
     #D0*t3+6208*l1+20416.D0/3.D0*t2-14080.D0/3.D0+t1234)*t16+(4608+t928
     #+8000*t3+40000.D0/3.D0*t6-40064.D0/3.D0*t11-13056*t2-10496.D0/3.D0
     #*l1)*t17+64.D0/3.D0*t167*(t1465+287*t6+72-t988-t1535+122*t11-t1163
     #)*t68-6208.D0/3.D0*l1
      t1553 = 530*t2
      t1602 = 3*t11
      t1630 = 2560+(-1536-76672.D0/3.D0*t3-27520.D0/3.D0*t4+t1420+28544*
     #t11-12544.D0/3.D0*l1-18304.D0/3.D0*t6+53248.D0/3.D0*t2)*t18-128.D0
     #/3.D0*t244*t56*(-t2+2)*t57-t1108+16000.D0/3.D0*t6-t938-11264.D0/3.
     #D0*l1-t956-128.D0/3.D0*t262*(131*t6-t916+t1602+36-t1072-47*t2)*t68
     #+128.D0/3.D0*t210*(-72-t1034-39*t3-258*t6+158*t11+193*t2+68*l1)*t4
     #0+(-9344.D0/3.D0*l1+4096*t4+19456.D0/3.D0*t3-40192.D0/3.D0*t11-467
     #20.D0/3.D0*t2+50176.D0/3.D0*t6+14336.D0/3.D0)*t17+(-48512.D0/3.D0*
     #t6+20480.D0/3.D0*t3-15872.D0/3.D0+7680*t2+t935-t1234)*t16
      t1635 = 2176*t7
      t1649 = 3328-64.D0/3.D0*t891*t56*(1+l1)*t100-64.D0/3.D0*t210*(-144
     #+t264-423*t4-654*t6-340*t3+t1553+867*t11+54*t7)*t68-3264*t3+t1091+
     #(-13312*l1+6208*t4+35072.D0/3.D0-128.D0/3.D0*t3-54400.D0/3.D0*t11-
     #74944.D0/3.D0*t2+115712.D0/3.D0*t6)*t17+(1536+141440.D0/3.D0*t3+24
     #128.D0/3.D0*t4+t159+21376.D0/3.D0*t6-3136*t79-82816.D0/3.D0*t2-149
     #504.D0/3.D0*t11+27136.D0/3.D0*t7)*t40-13120.D0/3.D0*t2+(-87808.D0/
     #3.D0*t6-27392.D0/3.D0+40640.D0/3.D0*t3+40256.D0/3.D0*t2+36992.D0/3
     #.D0*l1-896*t11)*t16+t1630*s2l1-12416.D0/3.D0*l1+(-6912-118016.D0/3
     #.D0*t3-70336.D0/3.D0*t4+t1635+58624*t11-640.D0/3.D0*l1-82816.D0/3.
     #D0*t6+110144.D0/3.D0*t2)*t18+64.D0/3.D0*t262*(72+270*t6+t12-22*t2-
     #t298+t1220-225*t3)*t57
      t1676 = 128*t6
      t1727 = 54*t6
      t1738 = 9472.D0/3.D0*l1
      t1758 = 64.D0/3.D0*t1*(-36-1055*t3+303*t4-99*t7+156*t11-290*l1+602
     #*t6+415*t2)*t40-t812-64.D0/3.D0*t33*(72+t1535+372*t3-55*t6-170*t11
     #-309*t2+92*l1)*t68-1472.D0/3.D0*t2+64.D0/3.D0*t64*(t1602+81*t2-t13
     #61+t434-t1727-36)*t57+15104.D0/3.D0*t6+(-12032.D0/3.D0-12096*t4+83
     #840.D0/3.D0*t2-49600.D0/3.D0*t6-89408.D0/3.D0*t3+384*t7-t1738+1118
     #08.D0/3.D0*t11)*t18-2240*t3+9728.D0/3.D0-16640.D0/3.D0*l1+(-6912*l
     #1+4160*t4+576*t3+23296.D0/3.D0-39104.D0/3.D0*t11+29248*t6-65408.D0
     #/3.D0*t2)*t17+(-21248.D0/3.D0-61312.D0/3.D0*t6+9216*t3+8320*t2+103
     #68*l1-t1344)*t16
      t1804 = -64.D0/3.D0*t485*t56*(-t2+2+l1)*t57-t812-t1738+(8192.D0/3.
     #D0+1920*t4+t1308-6976.D0/3.D0*l1+33536.D0/3.D0*t6-5760*t11-25280.D
     #0/3.D0*t2)*t17+(-24448.D0/3.D0*t6+12160.D0/3.D0*t3-9472.D0/3.D0+15
     #872.D0/3.D0*l1+6656.D0/3.D0*t2-t344)*t16-64.D0/3.D0*t1*(-40+346*t2
     #-t1535-438*t3-147*l1+t1084+206*t11)*t18+t1010+t1203+64.D0/3.D0*t33
     #*(80+70*t6+106*t3+t120-204*t2-48*t11)*t40-64.D0/3.D0*t64*(-t1727+1
     #6*l1-44+78*t2)*t68+4864.D0/3.D0-1024*t3
      t1822 = t1542*t186-8896.D0/3.D0*l1+t1649*s1l1+(t1635-49600.D0/3.D0
     #*t4-10880*t3-24448.D0/3.D0+34112*t11-98944.D0/3.D0*t6+72320.D0/3.D
     #0*t2+24448.D0/3.D0*l1)*t18-64.D0/3.D0*t32*(-72+t806-132*t4-159*t6-
     #t397+327*t11+270*t2)*t113+(-45056.D0/3.D0*t2+2752*t4-16960.D0/3.D0
     #*t11+27776.D0/3.D0-33920.D0/3.D0*t3+99584.D0/3.D0*t6-39680.D0/3.D0
     #*l1)*t17-t1056-1664*t3+((-t1676+t693-256+320.D0/3.D0*t2)*s1l1+2240
     #.D0/3.D0*l1-t1434-1280.D0/3.D0+(t1676-1280.D0/3.D0-1856.D0/3.D0*t2
     #+2752.D0/3.D0*l1)*s2l1)*t393+64.D0/3.D0*t63*(36+85*t6+t786-54*l1-t
     #913-63*t3)*t101+12224.D0/3.D0*t6+(3840+64768.D0/3.D0*t6+24128*t4-1
     #856*t7-768*t79+1408*l1-80704.D0/3.D0*t2-157376.D0/3.D0*t11+93056.D
     #0/3.D0*t3)*t40+t1758*s2l1+(t493-448.D0/3.D0*l1+256.D0/3.D0)*t392+6
     #016.D0/3.D0+t1804*t73+(-17536.D0/3.D0-t470-18560*t6+29120.D0/3.D0*
     #t3+21184.D0/3.D0*t2+26176.D0/3.D0*l1)*t16+64.D0/3.D0*t1*(36-223*t6
     #-309*t4+234*l1-t1553+266*t7+1017*t3-t1120-438*t11)*t68
      s1 = (t1139-5824*t4+67904.D0/3.D0*t2-2304-704*t3+34304.D0/3.D0*t11
     #-t1048-t1144)*t17-64.D0/3.D0*(-t63*(-20*t3+15*t11)+t63*(t961+38*t6
     #-23*t3-36*t2))*t2*t100+64.D0/3.D0*t1*(33*t79+190*t4-414*t2-t1163-9
     #*t3+776*t6-83*t7-382*t11)*t40+64.D0/3.D0*t32*(-39*t4+180*t2-t106-2
     #02*t6-t687+164*t11)*t113-64.D0/3.D0*t32*(-57*t7-t35+118*t4-74*t11+
     #t1018+48*t3)*t113+64.D0/3.D0*t63*(-45*t11+108*t3+28*t2-t265)*t101+
     #t1268*s1l1+(-t23+(t1262+t1270-t1064-t437)*s1l1+1472.D0/3.D0*t6+(28
     #16.D0/3.D0*t6-512*t3-256.D0/3.D0*t2-t23)*s2l1-t1279+512.D0/3.D0*t2
     #)*t393
      t1826 = s1-13312.D0/3.D0*t6+t1367*s1l1+(384*t4-t161+t179-8960*t2-t
     #208+2304-t1057)*t16+t1406*t186+(-448*t83-5248.D0/3.D0*t4-8384.D0/3
     #.D0*t7+t1144+t1062-106816.D0/3.D0*t3+3904.D0/3.D0*t79-t1413+21376*
     #t11)*t40+((384*t79-6784.D0/3.D0*t4-5120.D0/3.D0*t3-t1420+3840*t11)
     #*t18+(3520*t4+640.D0/3.D0*t3+t1287-t1139-10880.D0/3.D0*t11)*t17+(-
     #1280*t11-5120.D0/3.D0*t6+8960.D0/3.D0*t3)*t16+t1434-t45-1216*t3+12
     #16.D0/3.D0*t6+64.D0/3.D0*t210*(-t172-60*t3+21*t7+t986)*t40-64.D0/3
     #.D0*t262*(t263-28*t11+t1294)*t68+t1448)*t186+t1495*s2l1-t76*t75*t1
     #822/6+8192.D0/3.D0*t2
      t1829 = (l1+z-t363)**2
      taux23 = (t1138+t1826)/t1829*t366
      t1 = -1+l1
      t2 = t1*l1
      t3 = l1**2
      t4 = t3**2
      t5 = t4**2
      t7 = t3*l1
      t8 = t4*t7
      t10 = t4*l1
      t12 = t4*t3
      t17 = 1440*l1
      t19 = 1-z
      t20 = t19**2
      t21 = t20**2
      t22 = t21*t19
      t25 = t1**2
      t26 = t25*t3
      t34 = t21*t20
      t37 = 24*t5
      t40 = 590*t8
      t43 = 2480*l1
      t47 = t20*t19
      t49 = t5*l1
      t50 = 24*t49
      t56 = 152*t5
      t57 = 2304*l1
      t62 = 516*l1
      t63 = 24*t10
      t64 = 56*t4
      t69 = 48*t8
      t71 = 208*t12
      t72 = 2076*l1
      t75 = 158*t4
      t79 = t25*t1
      t82 = -4+3*l1
      t83 = t79*t8*t82
      t85 = t79*t10
      t86 = 2*t1
      t87 = t82*t86
      t88 = t85*t87
      t89 = 4*t88
      t91 = t21*t47
      t93 = 48*t12
      t95 = 1076*l1
      t96 = 40*t10
      t105 = 632*l1
      t106 = 20*t4
      t117 = 204*t3
      t118 = 1964*l1
      t123 = t79*t3
      t124 = 297*t10
      t128 = 576*l1
      t132 = t25*l1
      t135 = 420*t8
      t139 = 288*l1
      t155 = 734*t8
      t170 = 112*t10
      t203 = 1637*t7
      t207 = 1244*l1
      t218 = 1112*l1
      t223 = 560*l1
      t226 = 475*t3
      t231 = 1640*l1
      t238 = t25**2
      t240 = t238*t12*t82
      t242 = t238*t4
      t243 = t242*t87
      t244 = 4*t243
      t250 = 1304*l1
      t269 = 8*l1
      t297 = 184*t4
      t302 = 680*t7
      t315 = 6*t7
      t318 = 152*t4
      t320 = t79*(t117-t315+144+87*t10-216*l1-t318)
      t327 = 1160*l1
      t341 = 1872*l1
      t351 = 6*t10
      t360 = 24*t4
      t361 = 108*l1
      t377 = l1*t34
      s2 = (128.D0/3.D0*t2*(-348*t5+1382*t8+2229*t10-2214*t12-1476*t4+11
     #52+396*t7+416*t3-t17)*t22-128.D0/3.D0*t26*(206*t7+60*t12-1104*t3-5
     #76+1152*l1-392*t10+554*t4)*t34+128.D0/3.D0*t1*(-t37-1436*t7-1608*t
     #12+t40-2264*t3-655*t4+t43+2100*t10+2016)*t47-128.D0/3.D0*t1*(-t50-
     #2568*t3+1214*t10+230*t8+576-732*t4-904*t12+t56+t57+270*t7)*t21-128
     #.D0/3.D0*t1*(t62-t63+t64-956*t7+1512-586*t3)+128.D0/3.D0*t1*(-t69-
     #474*t10+t71-t72+2176*t7-3240+2056*t3-t75)*t20+128.D0/3.D0*(-11*t83
     #+t89)*t91-128.D0/3.D0*t1*(-t93+1137*t3-t95-3024+t96+1992*t7-220*t4
     #)*t19)*s1l1-128.D0/3.D0*t1*(t105-t63+t106-917*t7+1368-607*t3)
      s3 = s2-128.D0/3.D0*t1*(3226*t10-2007*t12+806*t8-t37-2588*t7-2184*
     #t4+3096-t117+t118)*t47
      s4 = s3
      s6 = (128.D0/3.D0*t123*(-t124-440*t7+213*t4+576+888*t3-t128)*t34-1
     #28.D0/3.D0*t132*(309*t10+452*t12-t135+452*t7+287*t4-1152-944*t3+t1
     #39)*t22+128.D0/3.D0*t1*(123*t10-t50+580*t12+70*t8+200*t5+1612*t7-2
     #933*t4+576-2144*t3+t57)*t21-128.D0/3.D0*t1*(2224*t10-1029*t12+t155
     #-t37-385*t7-3841*t4+2016-732*t3+t43)*t47-128.D0/3.D0*t1*(-1571*t10
     #+304*t12-t69+2384*t7+2301*t4-3240-296*t3-t72)*t20+128.D0/3.D0*t1*(
     #t170-t93+3209*t7-214*t4-3024-667*t3-t95)*t19+128.D0/3.D0*t1*(-t63-
     #1440*t7+104*t4+1512-26*t3+t62))*s1l1
      s7 = (128.D0/3.D0*t132*(-1225*t10-576+48*l1+154*t7-82*t3+756*t4+89
     #0*t12-348*t8)*t22-128.D0/3.D0*t123*(288+214*t3-147*t4-58*t7+126*t1
     #0-240*l1)*t34+128.D0/3.D0*t1*(-t37-t203-1674*t12+t40-476*t3-1067*t
     #4+t207+2713*t10+936)*t47-128.D0/3.D0*t1*(-t50-948*t3+1698*t10-12*t
     #8+288-1256*t4-571*t12+104*t5+t218-202*t7)*t21-128.D0/3.D0*t1*(t223
     #-t63+8*t4-676*t7+936-t226)+128.D0/3.D0*t1*(-t69-490*t10+112*t12-t2
     #31+2092*t7-1368+989*t3-425*t4)*t20+128.D0/3.D0*(-11*t240+t244)*t91
     #-128.D0/3.D0*t1*(-502*t4+871*t3+1717*t7-t93-1368-t250+t96)*t19)*s2
     #l1
      s5 = s6+s7
      s1 = s4+s5
      s2 = s1+(128.D0/3.D0*t123*(-t124-265*t7+426*t4+288+432*t3-184*l1)*
     #t34-128.D0/3.D0*t132*(722*t10+518*t12-t135+662*t7-1038*t4-576-676*
     #t3+t269)*t22+128.D0/3.D0*t1*(1235*t10-t50-2075*t12+926*t8+t56+167*
     #t7-853*t4+288-528*t3+t218)*t21-128.D0/3.D0*t1*(-218*t10-1065*t12+t
     #155-t37-1822*t7+1027*t4+936+464*t3+t207)*t47-128.D0/3.D0*t1*(45*t1
     #0+t71-t69+t203-574*t4-1368+289*t3-t231)*t20+128.D0/3.D0*t1*(t170-t
     #93+1210*t7-t297-1368+811*t3-t250)*t19+128.D0/3.D0*t1*(-t63-t302+t6
     #4+936-555*t3+t223))*s2l1+128.D0/3.D0*(-t79*(-33*t10+92*t4-8*t3-76*
     #t7)+2*t320)*t3*t91
      t380 = s2-128.D0/3.D0*t1*(-t93+1275*t3-t327-2952+76*t10+2277*t7-52
     #9*t4)*t19+128.D0/3.D0*t1*(1218*t7-384*t49-2912*t3-694*t10+1265*t5-
     #522*t4+t341+1895*t12-1883*t8+288)*t22-256.D0/3.D0*(t79*(-t106+15*t
     #10)-3*t79*(t351+38*t7-23*t4-36*t3)+2*t79*(106*t3-9*t10+t360-t361+7
     #2-55*t7))*t3*t91+128.D0/3.D0*t25*(-1180*t10+1049*t12-474*t8-t302+1
     #261*t4-576+200*t3)*t377
      t386 = 1568*l1
      t410 = 2408*l1
      t469 = dlog(z)
      t474 = 144*t8
      t483 = 72*t49
      t498 = 48*t5
      t524 = 96*t12
      t537 = 72*t10
      t548 = 2880*l1
      t557 = -t82*t86
      t559 = 4*t85*t557
      t565 = 362*t8
      t587 = 188*t10
      t611 = 253*l1
      t628 = 4*t242*t557
      t677 = 2*l6-1
      t678 = t677**2
      t681 = l5-1
      t683 = z*l1
      t685 = (1-l1+t683)**2
      t686 = 1/t685
      t691 = dsqrt(z)
      t693 = dsqrt(l5)
      t695 = dsqrt(-t681)
      t698 = dsqrt(l1)
      t702 = l5*z*l1-l1+1-l5+l5*l1
      t705 = dsqrt(-t1)
      t710 = t702**2
      t713 = dsqrt(4*l1*t678*z*l5*t681*t1*t686-4*t677*t691*t693*t695*t69
     #8*t702*t686*t705+t710*t686)
      t722 = t79*(-216-120*t10-222*t7+286*t4-84*t3+324*l1)
      s1 = 128.D0/3.D0*t1*(1064*t4+5886-1734*t7+2685*t10-2277*l1+t474-14
     #88*t12-3588*t3)*t20-128.D0/3.D0*t1*(-3196*l1+11514*t3-2880+t483-47
     #22*t7-176*t12+1379*t8-768*t5-333*t10-1240*t4)*t21+128.D0/3.D0*t1*(
     #3688*t12-5724+580*t4-1541*t7+950*l1-t498+7436*t3-5334*t10-650*t8)*
     #t47+128.D0/3.D0*t1*(-1710*t5+2970*t8-3744*l1-4434*t4-3498*t12-2072
     #*t7-576+432*t49+4463*t10+8088*t3)*t22-128.D0/3.D0*t1*(2026*t4-599*
     #l1-1794*t3-2552*t7+44*t10+3474-t524)*t19
      s3 = s1
      s5 = (-128.D0/3.D0*t26*(-63*t12+774*t7-1305*t4+1152+1248*t3-t57+57
     #0*t10)*t34-128.D0/3.D0*t1*(t537-124*l1-1476+807*t3-804*t4+1259*t7)
     #+128.D0/3.D0*t2*(-938*t10+288*t5+1153*t12+3057*t4+t548-5610*t7-230
     #4+2456*t3-1035*t8)*t22+128.D0/3.D0*(3*t83-t89+t559)*t91+128.D0/3.D
     #0*t1*(-4032-671*t4-3463*t10-t498+2677*t12-t565-1864*l1+8916*t3-174
     #0*t7)*t47-128.D0/3.D0*t1*(-4608*l1-1152+t483-2316*t7-1405*t12+9696
     #*t3+1875*t8+1700*t10-3276*t4-852*t5)*t21-128.D0/3.D0*t1*(1651*t4-9
     #75*l1-2274*t3-1885*t7+t587+3978-t524)*t19+128.D0/3.D0*t1*(667*t4+5
     #616-971*t7+2975*t10-1410*l1+t474-1656*t12-4617*t3)*t20)*s1l1
      s6 = (-128.D0/3.D0*t123*(-81*t10+138*t3-141*t7+440*l1+45*t4-576)*t
     #34-128.D0/3.D0*t1*(t537-t611-774+265*t3-588*t4+1121*t7)+128.D0/3.D
     #0*t132*(1299*t10+288*t8-861*t4-1768*t3+t269+994*t7-825*t12+1152)*t
     #22+128.D0/3.D0*(3*t240-t244+t628)*t91+128.D0/3.D0*t1*(-1800+1219*t
     #4-4470*t10-t498+2467*t12-t565-1084*l1+3064*t3+689*t7)*t47-128.D0/3
     #.D0*t1*(-2152*l1-576+t483+464*t7-617*t12+3636*t3+1455*t8-793*t10-9
     #62*t4-636*t5)*t21-128.D0/3.D0*t1*(1459*t4+1494+t587-2502*t7-t524+3
     #55*l1-636*t3)*t19+128.D0/3.D0*t1*(427*t4+2178-2250*t7+2399*t10+315
     #*l1+t474-1224*t12-1615*t3)*t20)*s2l1
      s4 = s5+s6
      s2 = s3+s4
      t744 = s2+128.D0/3.D0*t1*(-48*t4+262*t7+96*t3-129*l1-234)*t713+128
     #.D0/3.D0*(-t320+2*t722)*t3*t91-128.D0/3.D0*t1*(t537-t611-1242+584*
     #t3-720*t4+1321*t7)-128.D0/3.D0*t25*(-707*t10+1152+2696*t7-505*t4+3
     #3*t8+421*t12-2884*t3)*t377
      t746 = t469**2
      t797 = 100*t10
      t809 = 58*t8
      t876 = 30*t12
      t929 = 4*l1
      t940 = t79*(72-t361-48*t3-154*t4+48*t10+192*t7)
      s1 = -128.D0/3.D0*t1*(1800-697*t8+1750*l1+2803*t7+5963*t4-8485*t3+
     #278*t5-4909*t10+1611*t12-t50)*t21+128.D0/3.D0*t1*(360+83*t7-5486*t
     #3+869*t12-t474+2340*l1+304*t5-120*t49+6499*t4-4685*t10)*t22+(128.D
     #0/3.D0*(t243-t628)*t91-128.D0/3.D0*t1*(230*t4+180-t63+58*l1-24*t3-
     #390*t7)+128.D0/3.D0*t132*(-72*t8-243*t10+1605*t3-721*t7+249*t12-72
     #0-30*l1-t318)*t22-128.D0/3.D0*t123*(-359*t3+30*t10+360-31*t4-246*l
     #1+300*t7)*t34-128.D0/3.D0*t1*(2*l1-576-t797+t226+586*t7-495*t4+t93
     #)*t19+128.D0/3.D0*t1*(1098-783*t12+t37+686*t4-2538*t3+1085*t10+21*
     #t7+439*l1+t809)*t47+128.D0/3.D0*t1*(460*t12-t69-1188+1485*t3-930*t
     #10+192*l1-73*t4+2*t7)*t20-128.D0/3.D0*t1*(360-594*t8+1318*l1+2713*
     #t4-808*t10-2563*t3+230*t5-964*t7+362*t12-t50)*t21)*s2l1+128.D0/3.D
     #0*t1*(t360+27*l1-9*t3-86*t7+54)*t713+128.D0/3.D0*t1*(-6351*t3-767*
     #t12+3474+t37+1735*t4-1307*l1+3437*t7+154*t8-229*t10)*t47
      t965 = s1+(128.D0/3.D0*(t88-t559)*t91-128.D0/3.D0*t1*(302*t4+558-t
     #63-77*l1-537*t3-162*t7)+128.D0/3.D0*t2*(1440-72*t5-789*t4+1005*t12
     #+4938*t7-2129*t10-2652*t3-1800*l1+69*t8)*t22-128.D0/3.D0*t26*(-403
     #*t10-262*t3-720+1318*t4+t876-1421*t7+t17)*t34-128.D0/3.D0*t1*(2167
     #*t3+1023*l1-375*t4-2034-t797+t93-879*t7)*t19+128.D0/3.D0*t1*(2520-
     #435*t12+t37+4208*t4-7202*t3-1557*t10+1810*t7+724*l1+t809)*t47+128.
     #D0/3.D0*t1*(604*t12-t69-3330+4405*t3-642*t10+1635*l1-701*t4-2123*t
     #7)*t20-128.D0/3.D0*t1*(720-510*t8+t548+6794*t4-3931*t10-6920*t3+30
     #2*t5+96*t7+653*t12-t50)*t21)*s1l1+128.D0/3.D0*t1*(556*t12-t69-3276
     #+3518*t3-896*t10+2154*l1+t75-2346*t7)*t20-128.D0/3.D0*t1*(278*t4+3
     #96-t63+t929-285*t3-319*t7)+128.D0/3.D0*(-t722+2*t940)*t3*t91-128.D
     #0/3.D0*t25*(120*t8-720+1813*t10-1741*t7-856*t12+2334*t3-1006*t4)*t
     #377-128.D0/3.D0*t1*(1522*t3+740*l1-472*t4-1584-52*t10+t93-322*t7)*
     #t19
      t1009 = t557*t91
      t1070 = 460*t3
      s1 = -128.D0/3.D0*t1*(-692*t10+232*t12-t69+4354*t7-243*t4-3816-451
     #*t3-t386)*t20-128.D0/3.D0*t1*(2034*t10-456*t49-45*t12-1262*t8+1433
     #*t5+2572*t7-3556*t4+288-2680*t3+t341)*t22+128.D0/3.D0*t1*(4225*t10
     #-t50-2483*t12+847*t8+164*t5+1295*t7-4532*t4+1440-2320*t3+t410)*t21
     #-128.D0/3.D0*t25*(-415*t10-1018*t7-576+1054*t4+656*t3+162*t12-153*
     #t8)*t377+128.D0/3.D0*t1*(148*t10-t93+3194*t7-631*t4-2952+63*t3-t32
     #7)*t19-128.D0/3.D0*t1*(-3508*t3+1440+2300*t10-t50+371*t8+116*t5+t4
     #10-1422*t4-1235*t12+210*t7)*t21
      s2 = s1+128.D0/3.D0*t1*(-t37-2616*t3-1737*t12+662*t8+t118-2091*t7+
     #2777*t10-740*t4+3096)*t47+128.D0/3.D0*t1*(-t63-1161*t7+68*t4+1368-
     #367*t3+t105)
      s3 = s2+128.D0/3.D0*t1*(-t69-603*t10+136*t12-t386+2997*t7-3816+197
     #7*t3-539*t4)*t20
      s4 = s3-t469*t744
      s5 = s4
      s7 = t746*t965/2
      s9 = -t746/6
      s11 = t469
      s13 = 128.D0/3.D0*t1*(-684+876*t10-180*t12+310*l1+1644*t3-872*t4-1
     #114*t7)*t47+128.D0/3.D0*t1*(-72-354*t8+138*t12-2050*t4-468*l1+254*
     #t7+1312*t10+1148*t3+90*t5)*t22+128.D0/3.D0*t1*(612-474*l1-1052*t3-
     #160*t4-60*t10+1154*t7)*t20-128.D0/3.D0*t1*(-36-54*t7+10*l1+78*t3)-
     #128.D0/3.D0*t940*t3*t91
      s14 = s13+(128.D0/3.D0*t85*t1009-128.D0/3.D0*t1*(-144+324*t7-t876-
     #2348*t4-90*t8+1376*t10+1484*t3-t128)*t21-128.D0/3.D0*t1*(-72-102*t
     #7+28*l1+142*t3)+128.D0/3.D0*t1*(-504-426*t7-168*t12+1730*t3+1142*t
     #10-116*l1-1674*t4)*t47+128.D0/3.D0*t1*(648-372*l1-1238*t3+222*t4-1
     #98*t10+962*t7)*t20-128.D0/3.D0*t26*(144+104*t10-6*t12-44*t3+476*t7
     #-386*t4-t139)*t34+128.D0/3.D0*t2*(-288+84*t8-538*t12+970*t10-1224*
     #t7+360*l1+716*t3-80*t4)*t22-128.D0/3.D0*t1*(360+598*t7-630*t3-84*t
     #4-228*l1)*t19)*s1l1-128.D0/3.D0*t1*(252+474*t7-t1070-90*t4-166*l1)
     #*t19
      s12 = s14-128.D0/3.D0*t1*(-360+42*t8-578*t12+2064*t10-2322*t4-430*
     #t7-332*l1+1906*t3)*t21-128.D0/3.D0*t25*(144-652*t10+292*t12+470*t4
     #-t69+354*t7-556*t3)*t377+(128.D0/3.D0*t242*t1009-128.D0/3.D0*t1*(-
     #72+340*t7-26*t12-986*t4+6*t8+444*t10+550*t3-260*l1)*t21-128.D0/3.D
     #0*t1*(-t315-t269+14*t3)+128.D0/3.D0*t1*(-216-68*l1+632*t3+210*t10-
     #62*t7-508*t4)*t47+128.D0/3.D0*t1*(216-60*l1-t1070+t106-t351+302*t7
     #)*t20-128.D0/3.D0*t123*(-72-t351+122*t3-112*t7+t106+44*l1)*t34+128
     #.D0/3.D0*t132*(144+t269-394*t3+t297-84*t10+150*t7)*t22-128.D0/3.D0
     #*t1*(72+120*t7-184*t3-t929)*t19)*s2l1
      s10 = s11*s12
      s8 = s9*s10
      s6 = s7+s8
      t1147 = s5+s6
      t1150 = (l1+z-t683)**2
      taux24 = (t380+t1147)/t1150*t686

CCC Coefficients

      tcoff45 = taux23
      tcoff46 = taux23
      tcoff47 = taux24
      tcoff48 = taux24

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff45)
      tFjet46 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff46)
      tFjet47 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff47)
      tFjet48 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff48)

CCC Sum

      qgRVnum12e3 = tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tFje
     #t48*tcoff48
      RETURN
      END

      double precision function qgRVnum13e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -2+l1
      t2 = -1+l1
      t3 = t2**2
      t4 = t1*t3
      t5 = l1**2
      t6 = 8*l1
      t7 = 4*t5
      t10 = 1-z
      t11 = t10**2
      t12 = t11**2
      t13 = t12*t10
      t17 = t2*t1
      t18 = t5**2
      t19 = t5*l1
      t32 = t11*t10
      t35 = 4*t18
      t37 = 6*t5
      t40 = t1*(6-t35+12*t19-t6-t37)*t11
      t43 = 12*t5
      t48 = 2*l1
      t49 = -1+t48-t5
      t53 = s1l1**2
      t55 = t3*t2
      t56 = t1*t55
      t57 = -t2
      t62 = 2*t19
      t67 = 8*t19
      t68 = 16*t5
      t69 = 6*l1
      t74 = 12*l1-6-t37
      t75 = t1*t74
      t77 = 2*t5
      t80 = t1*(-t77-2+4*l1)
      t84 = t5*t1
      t87 = t12*t11
      t100 = t18*l1
      t109 = 10*l1
      t110 = 4*t100
      t134 = t3**2
      t144 = 8*t5
      t161 = s2l1**2
      t184 = 28*l1
      t210 = 7*t18
      s2 = (-192*t4*t5*(t6-t7-4)*t13+768*t17*l1*(-t5+2-t18+3*t19-3*l1)*t
     #12+192*t1*(t6+48*t19-16*t18+4-44*t5)*t32-384*t40-192*t1*(-12+24*l1
     #-t43)*t10+768*t1*t49)*t53+((-768*t56*t5*t57*t13-384*t4*l1*(t62-t7-
     #t48+4)*t12+384*t17*(-t67+t68-t69-2)*t32-384*t40-384*t75*t10+384*t8
     #0)*s2l1-192*t84*t55*(-t77+t69-4)*t87-192*t4*l1*(18*t5+2*t18-8-14*t
     #19+t48)*t13+192*t17*(4+22*l1-6*t100-56*t5+26*t19+10*t18)*t12+192*t
     #1*(18-t109+t110-70*t5-44*t18+102*t19)*t32-192*t1*(30-54*l1-12*t18+
     #30*t19+t37)*t11-192*t1*(46*l1+t62-26*t5-22)*t10+192*t75)*s1l1+(384
     #*t1*t134*t5*t13-384*t56*l1*(t5-l1-2)*t12+192*t4*(-t144+t6+2)*t32+3
     #84*t17*(t62-t7-l1+3)*t11-192*t1*t74*t10+384*t1*t49)*t161
      s1 = s2+(-384*t84*t134*t87+192*t56*l1*(t37-t69-4)*t13-192*t4*(t35-
     #t67-t43+16*l1+2)*t12-192*t17*(20*t19-40*t5+t109+10)*t32+192*t1*(8*
     #t18-18-24*t19+t37+t184)*t11+192*t1*(-t184+14*t5+14)*t10+384*t80)*s
     #2l1-192*t134*(-t5+t48-2)*t84*t12*t32-192*l1*t1*t55*(-5*t19+t144+t1
     #8-t48-4)*t87
      t247 = s1-192*t17*(-10-t6+t210-40*t19+t100+50*t5)*t12-192*t1*(-13+
     #t184+t62-17*t5)*t10+192*t1*(-3+t69-3*t5)+192*t1*(21-34*l1+54*t19-t
     #68+t110-29*t18)*t32-192*t4*(-2-14*l1-t210+20*t5+2*t100)*t13-192*t1
     #*(23-48*l1+22*t5+t67-5*t18)*t11
      t248 = zqe0(t57,z,l3)
      t250 = z*l1
      t252 = (l1+z-t250)**2
      t255 = (1-l1+t250)**2
      taux25 = -t247*t248/t252/t255
      t1 = l1**2
      t2 = -2+l1
      t3 = t2*t1
      t4 = -1+l1
      t5 = t4**2
      t6 = t5*t4
      t7 = 4*t1
      t8 = 12*l1
      t11 = 1-z
      t12 = t11**2
      t13 = t12**2
      t14 = t13*t12
      t17 = t2*t5
      t18 = t1*l1
      t21 = t1**2
      t22 = 4*t21
      t25 = t13*t11
      t28 = t4*t2
      t30 = 28*l1
      t31 = t21*l1
      t32 = 4*t31
      t33 = 88*t1
      t38 = 8*l1
      t40 = 8*t21
      t42 = t12*t11
      t46 = 8*t18
      t54 = -t4
      t55 = 4*t54
      t60 = t5**2
      t64 = t2*t6
      t65 = 4*l1
      t71 = 16*t1
      t141 = zqe0(t54,z,l3)
      t143 = z*l1
      t145 = (l1+z-t143)**2
      t148 = (1-l1+t143)**2
      taux26 = -((384*t3*t6*(-t7-8+t8)*t14-384*t17*l1*(-t8+16+20*t18-20*
     #t1-t22)*t25-384*t28*(-t22+60*t18+8+t30-t32-t33)*t13+384*t28*(28-60
     #*t1-t38+48*t18-t40)*t42+384*t28*(t7-36+40*l1-t46)*t12-384*t28*(-20
     #+24*l1-t7)*t11-384*t28*t55)*s1l1+(1536*t3*t60*t14-384*t64*l1*(t7-t
     #65-8)*t25-384*t17*(t71-16*l1-4)*t13+384*t28*(-t71-t65+12+t46)*t42+
     #4608*t28*t4*t12+384*t28*t55*t11)*s2l1+768*t60*(-t1+2*l1-2)*t3*t13*
     #t42+384*t64*l1*(-8-t65+t71+2*t21-10*t18)*t14-384*t17*(4+t46+t30-t3
     #2+10*t21-44*t1)*t25-384*t28*(20+t38+2*t31-26*t21-t33+84*t18)*t13+3
     #84*t28*(38+34*t18-38*l1-t40-26*t1)*t42+384*t28*(-34+46*l1-2*t18-10
     #*t1)*t12-384*t28*(-14+18*l1-t7)*t11-768*t28*t54)*t141/t145/t148

CCC Coefficients

      tcoff49 = taux25
      tcoff50 = taux25
      tcoff51 = taux26
      tcoff52 = taux26

CCC Jet functions

      tFjet49 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff50)
      tFjet51 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff52)

CCC Sum

      qgRVnum13e0 = tFjet49*tcoff49+tFjet50*tcoff50+tFjet51*tcoff51+tFje
     #t52*tcoff52
      RETURN
      END

      double precision function qgRVnum13e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -2+l1
      t2 = -1+l1
      t3 = t2**2
      t4 = t1*t3
      t5 = l1**2
      t6 = 8*l1
      t7 = 4*t5
      t10 = 1-z
      t11 = t10**2
      t12 = t11**2
      t13 = t12*t10
      t17 = t2*t1
      t18 = t5**2
      t19 = t5*l1
      t20 = 3*t19
      t21 = 3*l1
      t32 = t11*t10
      t35 = 4*t18
      t37 = 6*t5
      t40 = t1*(6-t35+12*t19-t6-t37)*t11
      t43 = 12*t5
      t48 = 2*l1
      t49 = -1+t48-t5
      t53 = s1l1**2
      t55 = t3*t2
      t56 = t1*t55
      t57 = -t2
      t62 = 2*t19
      t67 = 8*t19
      t68 = 16*t5
      t69 = 6*l1
      t73 = 12*l1
      t74 = t73-6-t37
      t75 = t1*t74
      t77 = 2*t5
      t78 = 4*l1
      t80 = t1*(-t77-2+t78)
      t84 = t5*t1
      t86 = t55*(-t77+t69-4)
      t87 = t12*t11
      t92 = 2*t18
      t100 = t18*l1
      t109 = 10*l1
      t110 = 4*t100
      t134 = t3**2
      t135 = t1*t134
      t144 = 8*t5
      t161 = s2l1**2
      t171 = 16*l1
      t176 = 20*t19
      t182 = 8*t18
      t184 = 28*l1
      t198 = t134*(-t5+t48-2)
      t199 = t12*t32
      t204 = 5*t19
      t205 = -t204+t144+t18-t48-4
      t210 = 7*t18
      t222 = 3*t5
      t233 = 14*l1
      t242 = 5*t18
      s2 = (-192*t4*t5*(t6-t7-4)*t13+768*t17*l1*(-t5+2-t18+t20-t21)*t12+
     #192*t1*(t6+48*t19-16*t18+4-44*t5)*t32-384*t40-192*t1*(-12+24*l1-t4
     #3)*t10+768*t1*t49)*t53+((-768*t56*t5*t57*t13-384*t4*l1*(t62-t7-t48
     #+4)*t12+384*t17*(-t67+t68-t69-2)*t32-384*t40-384*t75*t10+384*t80)*
     #s2l1-192*t84*t86*t87-192*t4*l1*(18*t5+t92-8-14*t19+t48)*t13+192*t1
     #7*(4+22*l1-6*t100-56*t5+26*t19+10*t18)*t12+192*t1*(18-t109+t110-70
     #*t5-44*t18+102*t19)*t32-192*t1*(30-54*l1-12*t18+30*t19+t37)*t11-19
     #2*t1*(46*l1+t62-26*t5-22)*t10+192*t75)*s1l1+(384*t135*t5*t13-384*t
     #56*l1*(t5-l1-2)*t12+192*t4*(-t144+t6+2)*t32+384*t17*(t62-t7-l1+3)*
     #t11-192*t1*t74*t10+384*t1*t49)*t161
      s1 = s2+(-384*t84*t134*t87+192*t56*l1*(t37-t69-4)*t13-192*t4*(t35-
     #t67-t43+t171+2)*t12-192*t17*(t176-40*t5+t109+10)*t32+192*t1*(t182-
     #18-24*t19+t37+t184)*t11+192*t1*(-t184+14*t5+14)*t10+384*t80)*s2l1-
     #192*t198*t84*t199-192*l1*t1*t55*t205*t87
      t247 = s1-192*t17*(-10-t6+t210-40*t19+t100+50*t5)*t12-192*t1*(-13+
     #t184+t62-17*t5)*t10+192*t1*(-3+t69-t222)+192*t1*(21-34*l1+54*t19-t
     #68+t110-29*t18)*t32-192*t4*(-2-t233-t210+20*t5+2*t100)*t13-192*t1*
     #(23-48*l1+22*t5+t67-t242)*t11
      t248 = zqe0(t57,z,l3)
      t249 = dlog(l1)
      t251 = sqe1(t57)
      t252 = zqe1(t57,z,l3)
      t253 = sbqe0(t57)
      t254 = dlog(z)
      t259 = 2*l6-1
      t260 = t259**2
      t263 = l5-1
      t265 = z*l1
      t267 = (1-l1+t265)**2
      t268 = 1/t267
      t273 = dsqrt(z)
      t275 = dsqrt(l5)
      t277 = dsqrt(-t263)
      t280 = dsqrt(l1)
      t284 = l5*z*l1-l1+1-l5+l5*l1
      t286 = dsqrt(t57)
      t291 = t284**2
      t293 = -4*l1*t260*z*l5*t263*t2*t268+4*t259*t273*t275*t277*t280*t28
     #4*t268*t286-t291*t268
      t308 = dsqrt(-t293)
      t310 = 9*t5
      t316 = 7*t19
      t362 = 18*t19
      t415 = 47*l1
      t471 = 35*l1
      t477 = 68*t5
      t482 = 27*l1
      t484 = 36*t5
      s2 = -384*t293*t1+(192*t1*(-5*l1-5)*s1l1-192*t1*(-l1+7)*s2l1-192*t
     #1*(t48+8))*t308+(-192*t4*t5*(t19+t310-t171+8)*t13+768*t17*l1*(-t31
     #6+3*t18+t5-4+t69)*t12+192*t1*(49*t18-t171-8+82*t5-106*t19-3*t100)*
     #t32-384*t1*(-t176+t77+t233+11*t18-12)*t11-192*t1*(-9*t19+37*t5-33*
     #l1+25)*t10+768*t1*(t77+3))*t53
      s3 = s2+((-384*t56*t5*(t5+t78-4)*t13-384*t4*l1*(-t204+t310+t78-8)*
     #t12+384*t17*(16*t19-29*t5+t73+4)*t32-384*t1*(t310+t182-12+11*l1-t3
     #62)*t11-384*t1*(t310-t20-t73+14)*t10+384*t1*(10+t5+l1))*s2l1+192*(
     #-t18*t1*t55+2*t84*t86)*t87-192*t4*l1*(-t78+31*t19-32*t5+16-t242)*t
     #13+192*t17*(-8+106*t5-68*t19-44*l1-9*t18+14*t100)*t12+192*t1*(-36-
     #186*t19+t171+109*t18-17*t100+110*t5)*t32-192*t1*(-61-21*t19+83*l1+
     #24*t18-45*t5)*t11-192*t1*(51-19*t19+55*t5-t415)*t10+192*t1*(27+7*l
     #1+t37))*s1l1
      s1 = s3+(-192*t135*t5*(l1+5)*t13-384*t56*l1*(-t77+t78+5)*t12+192*t
     #4*(t20+13*t5-23*l1-5)*t32+384*t17*(-t20+11*t5-t21-8)*t11+192*t1*(2
     #0*l1-20-t20+t5)*t10-384*t1*(-7+t78))*t161+(-192*(-t19*t1*t134-4*t8
     #4*t134)*t87+192*t56*l1*(-t19-t310+15*l1+8)*t13-192*t4*(-6*t18+t362
     #+t43-t471-4)*t12-192*t17*(-t18-22*t19+t477-t471-21)*t32+192*t1*(-t
     #482+33*t19-t182-t484+42)*t11+192*t1*(-44+t19+t482)*t10+384*t1*(15-
     #t48-t5))*s2l1-192*(t134*t5-3*t198)*t1*t5*t199
      t552 = s1+192*(-t56*l1*(4+t316-t92-t7)+2*t56*l1*t205)*t87-192*t17*
     #(28+t171+117*t19-38*t18+t100-119*t5)*t12-192*t1*(42-11*t19-t415+t4
     #84)*t10+192*t1*(20+t222-t21)+192*t1*(-56+10*t5+89*l1-105*t19+76*t1
     #8-16*t100)*t32-192*t4*(6+38*l1+15*t19+14*t18-7*t100-58*t5)*t13-192
     #*t1*(-62+105*l1+t316+t182-t477)*t11
      t556 = (l1+z-t265)**2
      taux25 = (t247*(t248*t249+t251-t252+t253*t254)-t552*t248)/t556*t26
     #8
      t1 = -1+l1
      t2 = -2+l1
      t3 = t2*t1
      t5 = l1**2
      t6 = t5*l1
      t7 = 2*t6
      t9 = -34+46*l1-t7-10*t5
      t10 = 1-z
      t11 = t10**2
      t15 = t1**2
      t16 = t15**2
      t18 = -t5+2*l1-2
      t20 = t2*t5
      t21 = t11*t10
      t22 = t11**2
      t23 = t22*t21
      t27 = t15*t1
      t28 = t2*t27
      t29 = 4*l1
      t30 = 16*t5
      t31 = t5**2
      t33 = 10*t6
      t35 = l1*(-8-t29+t30+2*t31-t33)
      t36 = t22*t11
      t40 = 8*l1
      t41 = t31*l1
      t42 = 2*t41
      t44 = 88*t5
      t46 = 20+t40+t42-26*t31-t44+84*t6
      t50 = 4*t5
      t51 = 12*l1
      t53 = t27*(-t50-8+t51)
      t56 = t2*t15
      t59 = 4*t31
      t61 = l1*(-t51+16+20*t6-20*t5-t59)
      t62 = t22*t10
      t65 = 60*t5
      t67 = 8*t31
      t68 = 28-t65-t40+48*t6-t67
      t71 = 40*l1
      t72 = 8*t6
      t73 = t50-36+t71-t72
      t76 = -t1
      t77 = 4*t76
      t78 = t3*t77
      t80 = 28*l1
      t81 = 4*t41
      t82 = -t59+60*t6+8+t80-t81-t44
      t85 = 24*l1
      t86 = -20+t85-t50
      t93 = -14+18*l1-t50
      t98 = 2*t3*t76
      t100 = 10*t31
      t102 = 4+t72+t80-t81+t100-44*t5
      t108 = 26*t5
      t109 = 38+34*t6-38*l1-t67-t108
      t117 = l1*(t50-t29-8)
      t121 = -t30-t29+12+t72
      t125 = 12*t1
      t129 = -t77
      t134 = t30-16*l1-4
      t141 = zqe0(t76,z,l3)
      t142 = dlog(l1)
      t144 = sqe1(t76)
      t145 = zqe1(t76,z,l3)
      t146 = sbqe0(t76)
      t147 = dlog(z)
      t151 = 3*l1
      t154 = 2*l6-1
      t155 = t154**2
      t158 = l5-1
      t160 = z*l1
      t162 = (1-l1+t160)**2
      t163 = 1/t162
      t168 = dsqrt(z)
      t170 = dsqrt(l5)
      t172 = dsqrt(-t158)
      t175 = dsqrt(l1)
      t179 = l5*z*l1-l1+1-l5+l5*l1
      t181 = dsqrt(t76)
      t186 = t179**2
      t189 = dsqrt(4*l1*t155*z*l5*t158*t1*t163-4*t154*t168*t170*t172*t17
     #5*t179*t163*t181+t186*t163)
      t193 = 9*t5
      t240 = 3*t5
      t255 = 2*t5
      t259 = t20*t16
      t287 = 6*t6
      t288 = 9*l1
      s2 = 384*(6-t151)*t189*t1+(384*(t20*t27*(8+t193-t51-t7)-4*t20*t53)
     #*t36-384*(t56*l1*(-23*t6-16+12*t5+t51+9*t31)-4*t56*t61)*t62-384*(t
     #3*(t42-8-t80-70*t6+19*t31+82*t5)-4*t3*t82)*t22+384*(t3*(15*t31-56*
     #t6-28+t65)-4*t3*t68)*t21+384*(t3*(-t193-30*l1+41+t33)-4*t3*t73)*t1
     #1-384*(t3*(33-t85+t240)-4*t3*t86)*t10-384*t3*(6*l1-15)+1536*t78)*s
     #1l1+(384*(t20*t16*(-4-t255+l1)-16*t259)*t36-384*(t28*l1*(t151-t240
     #+8+5*t6)-4*t28*t117)*t62-384*(t56*(4-18*t5+t72+17*l1)-4*t56*t134)*
     #t22+384*(t3*(-22*t6+7*t31-13-l1+t108)-4*t3*t121)*t21+384*(t3*(t287
     #-t288+18-t193)-4*t3*t125)*t11-384*(t3*(-t5+16-t288)-4*t3*t129)*t10
     #-384*t3*(t29-10))*s2l1
      s1 = s2-768*(-t259+5*t2*t16*t18)*t5*t23+384*(t28*l1*(8-t59+16*t6-1
     #3*t5+t29)-4*t28*t35)*t36-384*(t56*(-4+41*t5-t80-t100-17*t6+8*t41)-
     #4*t56*t102)*t62
      t375 = s1-384*(t3*(-20-13*l1-93*t6+t44-3*t41+38*t31)-4*t3*t46)*t22
     #+384*(t3*(-44*t6+37*t5+26*l1+13*t31-41)-4*t3*t109)*t21+384*(t3*(46
     #-t71+t287)-4*t3*t9)*t11-384*(t3*(32+t255-22*l1)-4*t3*t93)*t10-384*
     #t3*(5*l1-14)+1536*t98
      t379 = (l1+z-t160)**2
      taux26 = ((384*t3*t9*t11+768*t16*t18*t20*t23+384*t28*t35*t36-384*t
     #3*t46*t22+(384*t20*t53*t36-384*t56*t61*t62+384*t3*t68*t21+384*t3*t
     #73*t11-384*t78-384*t3*t82*t22-384*t3*t86*t10)*s1l1-384*t3*t93*t10-
     #384*t98-384*t56*t102*t62+384*t3*t109*t21+(1536*t20*t16*t36-384*t28
     #*t117*t62+384*t3*t121*t21+384*t3*t125*t11-384*t3*t129*t10-384*t56*
     #t134*t22)*s2l1)*(t141*t142+t144-t145+t146*t147)-t375*t141)/t379*t1
     #63

CCC Coefficients

      tcoff49 = taux25
      tcoff50 = taux25
      tcoff51 = taux26
      tcoff52 = taux26

CCC Jet functions

      tFjet49 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff50)
      tFjet51 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff52)

CCC Sum

      qgRVnum13e1 = tFjet49*tcoff49+tFjet50*tcoff50+tFjet51*tcoff51+tFje
     #t52*tcoff52
      RETURN
      END

      double precision function qgRVnum13e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -2+l1
      t2 = -1+l1
      t3 = t2**2
      t4 = t1*t3
      t5 = l1**2
      t6 = t5**2
      t7 = 7*t6
      t8 = 14*l1
      t10 = t6*l1
      t13 = 1-z
      t14 = t13**2
      t15 = t14**2
      t16 = t15*t13
      t20 = t5*l1
      t21 = 8*t20
      t22 = 48*l1
      t24 = 5*t6
      t29 = 2*t20
      t31 = 28*l1
      t36 = t3**2
      t37 = 2*l1
      t39 = t36*(-t5+t37-2)
      t40 = t1*t5
      t41 = t14*t13
      t42 = t15*t41
      t46 = 8*l1
      t47 = 4*t5
      t53 = 4*t6
      t55 = 6*t5
      t58 = t1*(6-t53+12*t20-t46-t55)*t14
      t60 = -1-t5+t37
      t63 = t2*t1
      t64 = 3*l1
      t65 = 3*t20
      t71 = 12*t5
      t77 = 16*t6
      t85 = s1l1**2
      t87 = t15*t14
      t88 = t36*t87
      t91 = 2*t5
      t92 = 4*l1
      t94 = t1*(-t91-2+t92)
      t96 = 14*t5
      t101 = l1*t1
      t102 = t3*t2
      t103 = 6*l1
      t109 = 8*t6
      t110 = 24*t20
      t115 = 10*l1
      t116 = 20*t20
      t122 = 16*l1
      t129 = 3*t5
      t133 = t1*t102
      t134 = -t2
      t143 = 12*l1
      t144 = -6+t143-t55
      t145 = t1*t144
      t147 = 16*t5
      t155 = 26*t5
      t162 = 30*t20
      t168 = t5*(-t91+t103-4)
      t173 = 6*t10
      t174 = 22*l1
      t175 = 26*t20
      t177 = 10*t6
      t182 = 2*t6
      t190 = 4*t10
      t200 = 5*t20
      t201 = 8*t5
      t203 = t102*(-t200+t201+t6-t37-4)
      t231 = s2l1**2
      t239 = 54*t20
      t240 = 34*l1
      s1 = -192*t4*(-2-t7-t8+20*t5+2*t10)*t16-192*t1*(23+t21-t22+22*t5-t
     #24)*t14-192*t1*(-13+t29-17*t5+t31)*t13-192*t39*t40*t42+(-192*t4*t5
     #*(-4+t46-t47)*t16-384*t58+768*t1*t60+768*t63*l1*(2-t6-t5-t64+t65)*
     #t15-192*t1*(-12-t71+24*l1)*t13+192*t1*(-t77-44*t5+48*t20+t46+4)*t4
     #1)*t85+(-384*t40*t88+384*t94+192*t1*(-t31+t96+14)*t13+192*t101*t10
     #2*(-4+t55-t103)*t16+192*t1*(t31+t109-18+t55-t110)*t14-192*t63*(t11
     #5+t116+10-40*t5)*t41-192*t4*(t53-t21+t122-t71+2)*t15)*s2l1
      t246 = s1+192*t1*(-3-t129+t103)+((-768*t133*t5*t134*t16-384*t58+38
     #4*t94-384*t4*l1*(t29-t47+4-t37)*t15-384*t145*t13+384*t63*(-t103-t2
     #1+t147-2)*t41)*s2l1-192*t1*(t29+46*l1-t155-22)*t13-192*t1*(30+t55-
     #54*l1-12*t6+t162)*t14-192*t133*t168*t87+192*t145+192*t63*(4-t173+t
     #174+t175-56*t5+t177)*t15-192*t4*l1*(t37-8+t182+18*t5-14*t20)*t16+1
     #92*t1*(18-t115+t190-70*t5-44*t6+102*t20)*t41)*s1l1-192*t101*t203*t
     #87+(384*t40*t36*t16+384*t63*(t29+3-t47-l1)*t14+384*t1*t60-384*t101
     #*t102*(-l1+t5-2)*t15-192*t1*t144*t13+192*t4*(t46-t201+2)*t41)*t231
     #-192*t63*(-10+t7+t10-t46+50*t5-40*t20)*t15+192*t1*(21+t239-t147-t2
     #40-29*t6+t190)*t41
      t247 = zqe0(t134,z,l3)
      t248 = dlog(l1)
      t249 = t248**2
      t252 = sqe1(t134)
      t253 = zqe1(t134,z,l3)
      t256 = sqe2(t134)
      t257 = zqe2(t134,z,l3)
      t258 = sbqe0(t134)
      t259 = dlog(z)
      t260 = t259**2
      t263 = sbqe1(t134)
      t271 = 7*t10
      t292 = 9*t5
      t300 = 18*t20
      t301 = 11*l1
      t310 = 9*t6
      t321 = 47*l1
      t347 = 24*t6
      t354 = t36*t5
      t358 = t5*t42
      t361 = 6*t6
      t362 = 35*l1
      t373 = t20*t1
      t381 = 68*t5
      t389 = 27*l1
      t394 = 36*t5
      t409 = 7*t20
      t412 = t101*t102*(4+t409-t182-t47)
      t424 = 3*t10
      t466 = 2*l6-1
      t467 = t466**2
      t470 = l5-1
      t472 = z*l1
      t474 = (1-l1+t472)**2
      t475 = 1/t474
      t480 = dsqrt(z)
      t482 = dsqrt(l5)
      t484 = dsqrt(-t470)
      t487 = dsqrt(l1)
      t491 = l5*z*l1-l1+1-l5+l5*l1
      t493 = dsqrt(t134)
      t498 = t491**2
      t500 = 4*l1*t467*z*l5*t470*t2*t475-4*t466*t480*t482*t484*t487*t491
     #*t475*t493+t498*t475
      t501 = dsqrt(t500)
      t521 = 23*l1
      t532 = 20*l1
      t549 = 10*t5
      s2 = -192*t4*(6+15*t20+14*t6-58*t5+38*l1-t271)*t16+(-192*t4*l1*(-t
     #92+31*t20-t24+16-32*t5)*t16+(-384*t133*t5*(t5-4+t92)*t16+384*t63*(
     #16*t20-29*t5+t143+4)*t41-384*t4*l1*(-t200+t92+t292-8)*t15-384*t1*(
     #14-t65-t143+t292)*t13-384*t1*(t109+t292-t300-12+t301)*t14+384*t1*(
     #t5+10+l1))*s2l1+192*t63*(-t310-8-68*t20+106*t5+14*t10-44*l1)*t15-1
     #92*t1*(-19*t20+55*t5+51-t321)*t13+192*(-t133*t6+2*t133*t168)*t87+1
     #92*t1*(-36+t122-186*t20+109*t6-17*t10+110*t5)*t41+192*t1*(7*l1+t55
     #+27)-192*t1*(-21*t20-61-45*t5+83*l1+t347)*t14)*s1l1-192*(t354-3*t3
     #9)*t1*t358
      s1 = s2+(-192*t4*(-t361-t362+t300-4+t71)*t15+192*t101*t102*(8-t20+
     #15*l1-t292)*t16-192*(-t373*t36-4*t40*t36)*t87-192*t63*(-22*t20-t36
     #2-t6-21+t381)*t41+384*t1*(15-t37-t5)+192*t1*(-44+t389+t20)*t13+192
     #*t1*(-t389-t109+42-t394+33*t20)*t14)*s2l1-192*t63*(-38*t6-119*t5+t
     #10+t122+117*t20+28)*t15+192*(-t412+2*t101*t203)*t87+(-192*t4*t5*(8
     #+t20+t292-t122)*t16+192*t1*(-8+82*t5-t122-t424+49*t6-106*t20)*t41+
     #768*t63*l1*(3*t6-4-t409+t5+t103)*t15-192*t1*(-9*t20+37*t5-33*l1+25
     #)*t13-384*t1*(t91+t8+11*t6-t116-12)*t14+768*t1*(t91+3))*t85
      t559 = s1+(192*t1*(-5*l1-5)*s1l1-192*t1*(-l1+7)*s2l1-192*t1*(t37+8
     #))*t501-192*t1*(-11*t20+t394+42-t321)*t13+192*t1*(20+t129-t64)-192
     #*t1*(-62-t381+t109+t409+105*l1)*t14+(-192*t40*t36*(l1+5)*t16+192*t
     #4*(t65-t521+13*t5-5)*t41-384*t101*t102*(-t91+t92+5)*t15+192*t1*(t5
     #32-t65+t5-20)*t13+384*t63*(-t65-t64+11*t5-8)*t14-384*t1*(-7+t92))*
     #t231+192*t1*(76*t6-56-16*t10+t549-105*t20+89*l1)*t41+384*t500*t1
      t577 = 5*t5
      t588 = 30*t5
      t597 = 15*t5
      t616 = 6*t20
      t620 = 9*l1
      t666 = 18*l1
      s2 = -384*t500*t1+(192*t1*(l1+9)*s1l1-192*t1*(-9+t64)*s2l1-192*t1*
     #(l1-11))*t501+(-192*t4*t5*(t65-t577)*t16+768*t63*l1*(-t200+t6+t577
     #)*t15+192*t1*(-t24+t271+t588-t162)*t41-384*t1*(t549-t116+t24+t115)
     #*t14-192*t1*(t597+t200-t362-5)*t13+768*t1*(t91-t92-3))*t85
      s3 = s2+((-384*t133*t5*(3*t5-3*l1)*t16-384*t4*l1*(t577-t64)*t15+38
     #4*t63*(t24-t616-t64)*t41-384*t1*(-t616+t620-t55+3+t182)*t14-384*t1
     #*(-t46+t20+t549-11)*t13+384*t1*(-l1-15+t47))*s2l1+384*t6*t1*t102*t
     #87-192*t4*l1*(t182-t96+4*t20)*t16+192*t63*(-t29+t588+t173-t347)*t1
     #5+192*t1*(26*l1+t394-90*t20+28*t6+t190)*t41-192*t1*(-46*t20+70*l1-
     #28*t5+8+t77)*t14-192*t1*(48*t5-62*l1+t29-28)*t13+192*t1*(t96-36-t6
     #66))*s1l1
      s1 = s3+(-192*t1*t36*t5*(t64-5)*t16-384*t133*l1*(t47-l1-5)*t15+192
     #*t4*(t65-27*t5+t521+5)*t41+384*t63*(t409-t597+l1+9)*t14+192*t1*(26
     #+7*t5+t65-t240)*t13-384*t1*(10+t91-t620))*t231+(-384*t373*t88+192*
     #t133*l1*(t29+t47-t143)*t16-192*t4*(t361-t175+t71+t666)*t15-192*t63
     #*(t109-54*t5+t22+8)*t41+192*t1*(90*t5+t109-t122-32-t239)*t14+192*t
     #1*(-t532+10*t20-t155+52)*t13+384*t1*(t5-22+t620))*s2l1-192*(-3*t35
     #4+2*t39)*t1*t358
      t776 = s1+192*(-t101*t102*(t6+t29-t577)+2*t412)*t87-192*t63*(-16-3
     #4*t20+35*t5-t301-t271+27*t6)*t15-192*t1*(t20+t129+21*l1-45)*t13+19
     #2*t1*(t5+t46-29)+192*t1*(32+t300-t174-21*t5+5*t10-t177)*t41-192*t4
     #*(-t110+t424-t532+25*t5-4+t310)*t16-192*t1*(t20-t389-t549+t6+45)*t
     #14
      t780 = (l1+z-t472)**2
      taux25 = (t246*(-t247*t249/2-(t252-t253)*t248+t256-t257-t258*t260/
     #2+t263*t259)+t559*(t247*t248+t252-t253+t258*t259)-t776*t247)/t780*
     #t475
      t1 = -1+l1
      t2 = -2+l1
      t3 = t1*t2
      t5 = l1**2
      t6 = t5*l1
      t7 = 2*t6
      t8 = 10*t5
      t9 = -34+46*l1-t7-t8
      t10 = 1-z
      t11 = t10**2
      t15 = t1**2
      t16 = t15**2
      t17 = 2*l1
      t18 = -t5+t17-2
      t20 = t2*t5
      t21 = t11*t10
      t22 = t11**2
      t23 = t22*t21
      t27 = t15*t1
      t28 = t2*t27
      t29 = 4*l1
      t30 = 16*t5
      t31 = t5**2
      t32 = 2*t31
      t33 = 10*t6
      t35 = l1*(-8-t29+t30+t32-t33)
      t36 = t22*t11
      t40 = 8*l1
      t41 = t31*l1
      t42 = 2*t41
      t44 = 88*t5
      t46 = 20+t40+t42-26*t31-t44+84*t6
      t50 = 4*t5
      t51 = 12*l1
      t53 = t27*(-t50-8+t51)
      t56 = t2*t15
      t58 = 20*t5
      t59 = 4*t31
      t61 = l1*(-t51+16+20*t6-t58-t59)
      t62 = t22*t10
      t65 = 60*t5
      t66 = 48*t6
      t67 = 8*t31
      t68 = 28-t65-t40+t66-t67
      t71 = 40*l1
      t72 = 8*t6
      t73 = t50-36+t71-t72
      t76 = -t1
      t77 = 4*t76
      t78 = t3*t77
      t80 = 28*l1
      t81 = 4*t41
      t82 = -t59+60*t6+8+t80-t81-t44
      t85 = 24*l1
      t86 = -20+t85-t50
      t93 = -14+18*l1-t50
      t98 = 2*t3*t76
      t100 = 10*t31
      t102 = 4+t72+t80-t81+t100-44*t5
      t108 = 26*t5
      t109 = 38+34*t6-38*l1-t67-t108
      t117 = l1*(t50-t29-8)
      t121 = -t30-t29+12+t72
      t125 = 12*t1
      t129 = -t77
      t133 = 16*l1
      t134 = t30-t133-4
      t141 = zqe0(t76,z,l3)
      t142 = dlog(l1)
      t143 = t142**2
      t146 = sqe1(t76)
      t147 = zqe1(t76,z,l3)
      t150 = sqe2(t76)
      t151 = zqe2(t76,z,l3)
      t152 = sbqe0(t76)
      t153 = dlog(z)
      t154 = t153**2
      t157 = sbqe1(t76)
      t161 = t20*t16
      t163 = t2*t16*t18
      t173 = t28*l1*(8-t59+16*t6-13*t5+t29)
      t175 = 4*t28*t35
      t183 = t56*(-4+41*t5-t80-t100-17*t6+8*t41)
      t185 = 4*t56*t102
      t189 = 3*l1
      t192 = 2*l6-1
      t193 = t192**2
      t196 = l5-1
      t198 = z*l1
      t200 = (1-l1+t198)**2
      t201 = 1/t200
      t206 = dsqrt(z)
      t208 = dsqrt(l5)
      t210 = dsqrt(-t196)
      t213 = dsqrt(l1)
      t217 = l5*z*l1-l1+1-l5+l5*l1
      t219 = dsqrt(t76)
      t224 = t217**2
      t227 = dsqrt(4*l1*t193*z*l5*t196*t1*t201-4*t192*t206*t208*t210*t21
     #3*t217*t201*t219+t224*t201)
      t231 = 6*t6
      t233 = t3*(46-t71+t231)
      t235 = 4*t3*t9
      t239 = 2*t5
      t242 = t3*(32+t239-22*l1)
      t244 = 4*t3*t93
      t250 = t3*(5*l1-14)
      t252 = 1536*t98
      t253 = 13*l1
      t258 = t3*(-20-t253-93*t6+t44-3*t41+38*t31)
      t260 = 4*t3*t46
      t264 = 9*t5
      t267 = t20*t27*(8+t264-t51-t7)
      t269 = 4*t20*t53
      t273 = 15*t31
      t276 = t3*(t273-56*t6-28+t65)
      t278 = 4*t3*t68
      t282 = 3*t5
      t284 = t3*(33-t85+t282)
      t286 = 4*t3*t86
      t291 = 12*t5
      t292 = 9*t31
      t295 = t56*l1*(-23*t6-16+t291+t51+t292)
      t297 = 4*t56*t61
      t305 = t3*(t42-8-t80-70*t6+19*t31+82*t5)
      t307 = 4*t3*t82
      t311 = 30*l1
      t313 = t3*(-t264-t311+41+t33)
      t315 = 4*t3*t73
      t321 = t3*(6*l1-15)
      t323 = 1536*t78
      t328 = t20*t16*(-4-t239+l1)
      t329 = 16*t161
      t333 = 7*t31
      t335 = t3*(-22*t6+t333-13-l1+t108)
      t337 = 4*t3*t121
      t341 = t3*(t29-10)
      t342 = 9*l1
      t344 = t3*(-t5+16-t342)
      t346 = 4*t3*t129
      t349 = 5*t6
      t352 = t28*l1*(t189-t282+8+t349)
      t354 = 4*t28*t117
      t357 = 18*t5
      t360 = t56*(4-t357+t72+17*l1)
      t362 = 4*t56*t134
      t366 = t3*(t231-t342+18-t264)
      t368 = 4*t3*t125
      t379 = t3*(-44*t6+37*t5+26*l1+13*t31-41)
      t381 = 4*t3*t109
      t385 = -768*(-t161+5*t163)*t5*t23+384*(t173-t175)*t36-384*(t183-t1
     #85)*t62+384*(6-t189)*t227*t1+384*(t233-t235)*t11-384*(t242-t244)*t
     #10-384*t250+t252-384*(t258-t260)*t22+(384*(t267-t269)*t36+384*(t27
     #6-t278)*t21-384*(t284-t286)*t10-384*(t295-t297)*t62-384*(t305-t307
     #)*t22+384*(t313-t315)*t11-384*t321+t323)*s1l1+(384*(t328-t329)*t36
     #+384*(t335-t337)*t21-384*t341-384*(t344-t346)*t10-384*(t352-t354)*
     #t62-384*(t360-t362)*t22+384*(t366-t368)*t11)*s2l1+384*(t379-t381)*
     #t21
      t394 = 3*t6
      t465 = 11*t31
      t482 = 7*t5
      t483 = 14*l1
      s2 = 384*(-26+t253)*t227*t1+(384*(t20*t27*(-t50+t394)-4*t267+t269)
     #*t36-384*(t56*l1*(-25*t6+t333+t58)-4*t295+t297)*t62-384*(t3*(-40*t
     #6+7*t41+32*t5+t32)-4*t305+t307)*t22+384*(t3*(-t66+t273+20*l1+t58)-
     #4*t276+t278)*t21+384*(t3*(-35*l1+t357+t349-4)-4*t313+t315)*t11-384
     #*(t3*(-9-t133+t264)-4*t284+t286)*t10-384*t3*(l1+6)+1536*t321-t323)
     #*s1l1+(384*(3*t31*t2*t16-4*t328+t329)*t36-384*(t28*l1*(t349-t291)-
     #4*t352+t354)*t62-384*(t56*(t292-t30+t7)-4*t360+t362)*t22+384*(t3*(
     #t40+t465+t357-36*t6)-4*t335+t337)*t21+384*(t3*(9*t6+5*t5-1-19*l1)-
     #4*t366+t368)*t11-384*(t3*(t482-t483-3)-4*t344+t346)*t10-384*t3*(3+
     #3*l1)+1536*t341)*s2l1
      s1 = s2-768*(5*t161-8*t163)*t5*t23+384*(t28*l1*(t394+t32-t8)-4*t17
     #3+t175)*t36-384*(t56*(t30+4*t6-21*t31+5*t41)-4*t183+t185)*t62
      t554 = s1-384*(t3*(28*t5-54*t6+t465+10*l1+6*t41)-4*t258+t260)*t22+
     #384*(t3*(t465+t239-38*t6+2+t311)-4*t379+t381)*t21+384*(t3*(-32*l1+
     #t231-7+17*t5)-4*t233+t235)*t11-384*(t3*(-9-t483+t482)-4*t242+t244)
     #*t10-384*t3*(t17+5)+1536*t250-t252
      t558 = (l1+z-t198)**2
      taux26 = ((384*t3*t9*t11+768*t16*t18*t20*t23+384*t28*t35*t36-384*t
     #3*t46*t22+(384*t20*t53*t36-384*t56*t61*t62+384*t3*t68*t21+384*t3*t
     #73*t11-384*t78-384*t3*t82*t22-384*t3*t86*t10)*s1l1-384*t3*t93*t10-
     #384*t98-384*t56*t102*t62+384*t3*t109*t21+(1536*t20*t16*t36-384*t28
     #*t117*t62+384*t3*t121*t21+384*t3*t125*t11-384*t3*t129*t10-384*t56*
     #t134*t22)*s2l1)*(-t141*t143/2-(t146-t147)*t142+t150-t151-t152*t154
     #/2+t157*t153)+t385*(t141*t142+t146-t147+t152*t153)-t554*t141)/t558
     #*t201

CCC Coefficients

      tcoff49 = taux25
      tcoff50 = taux25
      tcoff51 = taux26
      tcoff52 = taux26

CCC Jet functions

      tFjet49 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff50)
      tFjet51 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff52)

CCC Sum

      qgRVnum13e2 = tFjet49*tcoff49+tFjet50*tcoff50+tFjet51*tcoff51+tFje
     #t52*tcoff52
      RETURN
      END

      double precision function qgRVnum13e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -2+l1
      t2 = -1+l1
      t3 = t2**2
      t4 = t1*t3
      t5 = l1**2
      t6 = t5**2
      t7 = 7*t6
      t8 = 14*l1
      t10 = t6*l1
      t11 = 2*t10
      t13 = 1-z
      t14 = t13**2
      t15 = t14**2
      t16 = t15*t13
      t20 = t5*l1
      t21 = 8*t20
      t22 = 48*l1
      t24 = 5*t6
      t29 = 2*t20
      t31 = 28*l1
      t36 = t3**2
      t37 = 2*l1
      t39 = t36*(-t5+t37-2)
      t41 = t14*t13
      t42 = t15*t41
      t46 = 8*l1
      t47 = 4*t5
      t53 = 4*t6
      t55 = 6*t5
      t58 = t1*(6-t53+12*t20-t46-t55)*t14
      t60 = -1-t5+t37
      t63 = t2*t1
      t64 = 3*l1
      t65 = 3*t20
      t71 = 12*t5
      t77 = 16*t6
      t85 = s1l1**2
      t87 = t1*t36
      t88 = t15*t14
      t92 = 2*t5
      t93 = 4*l1
      t94 = -t92-2+t93
      t95 = t1*t94
      t97 = 14*t5
      t102 = t3*t2
      t103 = t1*t102
      t104 = 6*l1
      t110 = 8*t6
      t111 = 24*t20
      t116 = 10*l1
      t117 = 20*t20
      t118 = 40*t5
      t123 = 16*l1
      t130 = 3*t5
      t134 = -t2
      t143 = 12*l1
      t144 = -6+t143-t55
      t145 = t1*t144
      t147 = 16*t5
      t155 = 26*t5
      t162 = 30*t20
      t168 = t5*(-t92+t104-4)
      t173 = 6*t10
      t174 = 22*l1
      t175 = 26*t20
      t177 = 10*t6
      t182 = 2*t6
      t183 = 18*t5
      t190 = 4*t10
      t200 = 5*t20
      t201 = 8*t5
      t203 = l1*(-t200+t201+t6-t37-4)
      t231 = s2l1**2
      t239 = 54*t20
      t240 = 34*l1
      s1 = -192*t4*(-2-t7-t8+20*t5+t11)*t16-192*t1*(23+t21-t22+22*t5-t24
     #)*t14-192*t1*(-13+t29-17*t5+t31)*t13-192*t39*t1*t5*t42+(-192*t4*t5
     #*(-4+t46-t47)*t16-384*t58+768*t1*t60+768*t63*l1*(2-t6-t5-t64+t65)*
     #t15-192*t1*(-12-t71+24*l1)*t13+192*t1*(-t77-44*t5+48*t20+t46+4)*t4
     #1)*t85+(-384*t87*t5*t88+384*t95+192*t1*(-t31+t97+14)*t13+192*t103*
     #l1*(-4+t55-t104)*t16+192*t1*(t31+t110-18+t55-t111)*t14-192*t63*(t1
     #16+t117+10-t118)*t41-192*t4*(t53-t21+t123-t71+2)*t15)*s2l1
      t246 = s1+192*t1*(-3-t130+t104)+((-768*t103*t5*t134*t16-384*t58+38
     #4*t95-384*t4*l1*(t29-t47+4-t37)*t15-384*t145*t13+384*t63*(-t104-t2
     #1+t147-2)*t41)*s2l1-192*t1*(t29+46*l1-t155-22)*t13-192*t1*(30+t55-
     #54*l1-12*t6+t162)*t14-192*t103*t168*t88+192*t145+192*t63*(4-t173+t
     #174+t175-56*t5+t177)*t15-192*t4*l1*(t37-8+t182+t183-14*t20)*t16+19
     #2*t1*(18-t116+t190-70*t5-44*t6+102*t20)*t41)*s1l1-192*t103*t203*t8
     #8+(384*t87*t5*t16+384*t63*(t29+3-t47-l1)*t14+384*t1*t60-384*t103*l
     #1*(-l1+t5-2)*t15-192*t1*t144*t13+192*t4*(t46-t201+2)*t41)*t231-192
     #*t63*(-10+t7+t10-t46+50*t5-40*t20)*t15+192*t1*(21+t239-t147-t240-2
     #9*t6+t190)*t41
      t247 = sqe1(t134)
      t248 = zqe1(t134,z,l3)
      t249 = t247-t248
      t250 = dlog(l1)
      t251 = t250**2
      t254 = sqe2(t134)
      t255 = zqe2(t134,z,l3)
      t258 = sbqe1(t134)
      t259 = dlog(z)
      t260 = t259**2
      t263 = sbqe2(t134)
      t272 = 5*t5
      t273 = t6+t29-t272
      t293 = 6*t20
      t309 = zqe0(t134,z,l3)
      t315 = 7*t10
      t336 = 9*t5
      t344 = 18*t20
      t345 = 11*l1
      t354 = 9*t6
      t365 = 47*l1
      t391 = 24*t6
      t398 = t36*t5
      t402 = t5*t42
      t405 = 6*t6
      t406 = 35*l1
      t424 = 68*t5
      t432 = 27*l1
      t437 = 36*t5
      t452 = 7*t20
      t455 = t103*l1*(4+t452-t182-t47)
      t467 = 3*t10
      t509 = 2*l6-1
      t510 = t509**2
      t513 = l5-1
      t515 = z*l1
      t517 = (1-l1+t515)**2
      t518 = 1/t517
      t523 = dsqrt(z)
      t525 = dsqrt(l5)
      t527 = dsqrt(-t513)
      t530 = dsqrt(l1)
      t534 = l5*z*l1-l1+1-l5+l5*l1
      t536 = dsqrt(t134)
      t541 = t534**2
      t543 = 4*l1*t510*z*l5*t513*t2*t518-4*t509*t523*t525*t527*t530*t534
     #*t518*t536+t541*t518
      t544 = dsqrt(t543)
      t564 = 23*l1
      t575 = 20*l1
      t592 = 10*t5
      s2 = -192*t4*(6+15*t20+14*t6-58*t5+38*l1-t315)*t16+(-192*t4*l1*(-t
     #93+31*t20-t24+16-32*t5)*t16+(-384*t103*t5*(t5-4+t93)*t16+384*t63*(
     #16*t20-29*t5+t143+4)*t41-384*t4*l1*(-t200+t93+t336-8)*t15-384*t1*(
     #14-t65-t143+t336)*t13-384*t1*(t110+t336-t344-12+t345)*t14+384*t1*(
     #t5+10+l1))*s2l1+192*t63*(-t354-8-68*t20+106*t5+14*t10-44*l1)*t15-1
     #92*t1*(-19*t20+55*t5+51-t365)*t13+192*(-t103*t6+2*t103*t168)*t88+1
     #92*t1*(-36+t123-186*t20+109*t6-17*t10+110*t5)*t41+192*t1*(7*l1+t55
     #+27)-192*t1*(-21*t20-61-45*t5+83*l1+t391)*t14)*s1l1-192*(t398-3*t3
     #9)*t1*t402
      s1 = s2+(-192*t4*(-t405-t406+t344-4+t71)*t15+192*t103*l1*(8-t20+15
     #*l1-t336)*t16-192*(-t87*t20-4*t87*t5)*t88-192*t63*(-22*t20-t406-t6
     #-21+t424)*t41+384*t1*(15-t37-t5)+192*t1*(-44+t432+t20)*t13+192*t1*
     #(-t432-t110+42-t437+33*t20)*t14)*s2l1-192*t63*(-38*t6-119*t5+t10+t
     #123+117*t20+28)*t15+192*(-t455+2*t103*t203)*t88+(-192*t4*t5*(8+t20
     #+t336-t123)*t16+192*t1*(-8+82*t5-t123-t467+49*t6-106*t20)*t41+768*
     #t63*l1*(3*t6-4-t452+t5+t104)*t15-192*t1*(-9*t20+37*t5-33*l1+25)*t1
     #3-384*t1*(t92+t8+11*t6-t117-12)*t14+768*t1*(t92+3))*t85
      t602 = s1+(192*t1*(-5*l1-5)*s1l1-192*t1*(-l1+7)*s2l1-192*t1*(t37+8
     #))*t544-192*t1*(-11*t20+t437+42-t365)*t13+192*t1*(20+t130-t64)-192
     #*t1*(-62-t424+t110+t452+105*l1)*t14+(-192*t87*t5*(l1+5)*t16+192*t4
     #*(t65-t564+13*t5-5)*t41-384*t103*l1*(-t92+t93+5)*t15+192*t1*(t575-
     #t65+t5-20)*t13+384*t63*(-t65-t64+11*t5-8)*t14-384*t1*(-7+t93))*t23
     #1+192*t1*(76*t6-56-16*t10+t592-105*t20+89*l1)*t41+384*t543*t1
      t606 = sbqe0(t134)
      t635 = 18*l1
      t653 = 9*l1
      t659 = 30*t5
      t741 = 15*t5
      s2 = 192*(-t103*l1*t273+2*t455)*t88-192*t63*(-t345-t315+35*t5-16+2
     #7*t6-34*t20)*t15+(192*t1*(90*t5+t110-t123-32-t239)*t14-192*t63*(t1
     #10-54*t5+t22+8)*t41-192*t4*(t405+t71-t175+t635)*t15+192*t1*(-t155-
     #t575+10*t20+52)*t13+192*t103*l1*(t29+t47-t143)*t16-384*t87*t20*t88
     #+384*t1*(-22+t5+t653))*s2l1
      s1 = s2+(192*t63*(-t29+t659+t173-t391)*t15-192*t1*(70*l1-28*t5-46*
     #t20+t77+8)*t14+192*t1*(t437-90*t20+28*t6+t190+26*l1)*t41+(-384*t1*
     #(t653-t55-t293+3+t182)*t14-384*t4*l1*(t272-t64)*t15-384*t1*(-t46+t
     #20+t592-11)*t13-384*t103*t5*(-3*l1+3*t5)*t16+384*t63*(-t293-t64+t2
     #4)*t41+384*t1*(-l1+t47-15))*s2l1-192*t1*(48*t5-62*l1+t29-28)*t13-1
     #92*t4*l1*(-t97+4*t20+t182)*t16+384*t103*t6*t88+192*t1*(-t635-36+t9
     #7))*s1l1-192*t1*(t20+21*l1+t130-45)*t13-192*(-3*t398+2*t39)*t1*t40
     #2+(-384*t1*(-t117+t592+t24+t116)*t14+768*t63*l1*(-t200+t6+t272)*t1
     #5-192*t1*(t741+t200-5-t406)*t13-192*t4*t5*(t65-t272)*t16+192*t1*(-
     #t162+t659-t24+t315)*t41+768*t1*(t92-3-t93))*t85
      t820 = s1+(384*t63*(t452-t741+9+l1)*t14-384*t103*l1*(t47-l1-5)*t15
     #+192*t1*(26+7*t5+t65-t240)*t13-192*t87*t5*(t64-5)*t16+192*t4*(-27*
     #t5+t564+t65+5)*t41-384*t1*(t92-t653+10))*t231+192*t1*(t46+t5-29)+(
     #192*t1*(l1+9)*s1l1-192*t1*(-9+t64)*s2l1-192*t1*(l1-11))*t544+192*t
     #1*(-21*t5-t174-t177+5*t10+32+t344)*t41-384*t543*t1-192*t4*(t467-t5
     #75+25*t5-4-t111+t354)*t16-192*t1*(-t592+45-t432+t6+t20)*t14
      t827 = (l1+z-t515)**2
      taux25 = (t246*(t249*t251/2-(t254-t255)*t250-t258*t260/2+t263*t259
     #)-(-384*t36*t6*t1*t42+384*l1*t1*t102*t273*t88-192*t63*(t183-52*t20
     #+t116+22*t6+t11)*t15-192*t1*(-6*t5-6*l1+6*t20+6)*t13-192*t1*t94+19
     #2*t1*(-34*t6+t293-t174+t118+12*t10-2)*t41-192*t4*(-t29+t147-18*t6+
     #t173)*t16-192*t1*(-6-t175-t116+t118+t182)*t14)*t309+t602*(-t309*t2
     #51/2-t249*t250+t254-t255-t606*t260/2+t258*t259)+t820*(t309*t250+t2
     #47-t248+t606*t259))/t827*t518
      t1 = -1+l1
      t2 = -2+l1
      t3 = t1*t2
      t5 = l1**2
      t6 = t5*l1
      t7 = 2*t6
      t8 = 10*t5
      t9 = -34+46*l1-t7-t8
      t10 = 1-z
      t11 = t10**2
      t15 = t1**2
      t16 = t15**2
      t17 = 2*l1
      t18 = -t5+t17-2
      t20 = t2*t5
      t21 = t11*t10
      t22 = t11**2
      t23 = t22*t21
      t27 = t15*t1
      t28 = t2*t27
      t29 = 4*l1
      t30 = 16*t5
      t31 = t5**2
      t32 = 2*t31
      t33 = 10*t6
      t35 = l1*(-8-t29+t30+t32-t33)
      t36 = t22*t11
      t40 = 8*l1
      t41 = t31*l1
      t42 = 2*t41
      t44 = 88*t5
      t46 = 20+t40+t42-26*t31-t44+84*t6
      t50 = 4*t5
      t51 = 12*l1
      t53 = t27*(-t50-8+t51)
      t56 = t2*t15
      t58 = 20*t5
      t59 = 4*t31
      t61 = l1*(-t51+16+20*t6-t58-t59)
      t62 = t22*t10
      t65 = 60*t5
      t66 = 48*t6
      t67 = 8*t31
      t68 = 28-t65-t40+t66-t67
      t71 = 40*l1
      t72 = 8*t6
      t73 = t50-36+t71-t72
      t76 = -t1
      t77 = 4*t76
      t78 = t3*t77
      t80 = 28*l1
      t81 = 4*t41
      t82 = -t59+60*t6+8+t80-t81-t44
      t85 = 24*l1
      t86 = -20+t85-t50
      t93 = -14+18*l1-t50
      t98 = 2*t3*t76
      t100 = 10*t31
      t102 = 4+t72+t80-t81+t100-44*t5
      t108 = 26*t5
      t109 = 38+34*t6-38*l1-t67-t108
      t113 = t2*t16
      t118 = l1*(t50-t29-8)
      t122 = -t30-t29+12+t72
      t126 = 12*t1
      t130 = -t77
      t134 = 16*l1
      t135 = t30-t134-4
      t142 = sqe1(t76)
      t143 = zqe1(t76,z,l3)
      t144 = t142-t143
      t145 = dlog(l1)
      t146 = t145**2
      t149 = sqe2(t76)
      t150 = zqe2(t76,z,l3)
      t153 = sbqe1(t76)
      t154 = dlog(z)
      t155 = t154**2
      t158 = sbqe2(t76)
      t162 = t113*t5
      t163 = t113*t18
      t173 = t28*l1*(8-t59+16*t6-13*t5+t29)
      t175 = 4*t28*t35
      t183 = t56*(-4+41*t5-t80-t100-17*t6+8*t41)
      t185 = 4*t56*t102
      t189 = 3*l1
      t192 = 2*l6-1
      t193 = t192**2
      t196 = l5-1
      t198 = z*l1
      t200 = (1-l1+t198)**2
      t201 = 1/t200
      t206 = dsqrt(z)
      t208 = dsqrt(l5)
      t210 = dsqrt(-t196)
      t213 = dsqrt(l1)
      t217 = l5*z*l1-l1+1-l5+l5*l1
      t219 = dsqrt(t76)
      t224 = t217**2
      t227 = dsqrt(4*l1*t193*z*l5*t196*t1*t201-4*t192*t206*t208*t210*t21
     #3*t217*t201*t219+t224*t201)
      t231 = 6*t6
      t233 = t3*(46-t71+t231)
      t235 = 4*t3*t9
      t239 = 2*t5
      t242 = t3*(32+t239-22*l1)
      t244 = 4*t3*t93
      t250 = t3*(5*l1-14)
      t252 = 1536*t98
      t253 = 13*l1
      t258 = t3*(-20-t253-93*t6+t44-3*t41+38*t31)
      t260 = 4*t3*t46
      t264 = 9*t5
      t267 = t20*t27*(8+t264-t51-t7)
      t269 = 4*t20*t53
      t273 = 15*t31
      t276 = t3*(t273-56*t6-28+t65)
      t278 = 4*t3*t68
      t282 = 3*t5
      t284 = t3*(33-t85+t282)
      t286 = 4*t3*t86
      t291 = 12*t5
      t292 = 9*t31
      t295 = t56*l1*(-23*t6-16+t291+t51+t292)
      t297 = 4*t56*t61
      t305 = t3*(t42-8-t80-70*t6+19*t31+82*t5)
      t307 = 4*t3*t82
      t311 = 30*l1
      t313 = t3*(-t264-t311+41+t33)
      t315 = 4*t3*t73
      t321 = t3*(6*l1-15)
      t323 = 1536*t78
      t328 = t113*t5*(-4-t239+l1)
      t329 = 16*t162
      t333 = 7*t31
      t335 = t3*(-22*t6+t333-13-l1+t108)
      t337 = 4*t3*t122
      t341 = t3*(t29-10)
      t342 = 9*l1
      t344 = t3*(-t5+16-t342)
      t346 = 4*t3*t130
      t349 = 5*t6
      t352 = t28*l1*(t189-t282+8+t349)
      t354 = 4*t28*t118
      t357 = 18*t5
      t360 = t56*(4-t357+t72+17*l1)
      t362 = 4*t56*t135
      t366 = t3*(t231-t342+18-t264)
      t368 = 4*t3*t126
      t379 = t3*(-44*t6+37*t5+26*l1+13*t31-41)
      t381 = 4*t3*t109
      t385 = -768*(-t162+5*t163)*t5*t23+384*(t173-t175)*t36-384*(t183-t1
     #85)*t62+384*(6-t189)*t227*t1+384*(t233-t235)*t11-384*(t242-t244)*t
     #10-384*t250+t252-384*(t258-t260)*t22+(384*(t267-t269)*t36+384*(t27
     #6-t278)*t21-384*(t284-t286)*t10-384*(t295-t297)*t62-384*(t305-t307
     #)*t22+384*(t313-t315)*t11-384*t321+t323)*s1l1+(384*(t328-t329)*t36
     #+384*(t335-t337)*t21-384*t341-384*(t344-t346)*t10-384*(t352-t354)*
     #t62-384*(t360-t362)*t22+384*(t366-t368)*t11)*s2l1+384*(t379-t381)*
     #t21
      t386 = zqe0(t76,z,l3)
      t390 = sbqe0(t76)
      t402 = 3*t6
      t405 = t28*l1*(t402+t32-t8)
      t413 = t3*(-32*l1+t231-7+17*t5)
      t418 = 14*l1
      t419 = 7*t5
      t421 = t3*(-9-t418+t419)
      t427 = t3*(t17+5)
      t429 = 1536*t250
      t434 = t56*(t30+4*t6-21*t31+5*t41)
      t445 = t3*(-35*l1+t357+t349-4)
      t452 = t20*t27*(-t50+t402)
      t460 = t56*l1*(-25*t6+t333+t58)
      t467 = t3*(-t66+t273+20*l1+t58)
      t473 = t3*(-9-t134+t264)
      t479 = t3*(l1+6)
      t481 = 1536*t321
      t486 = t3*(-40*t6+7*t41+32*t5+t32)
      t497 = t3*(9*t6+5*t5-1-19*l1)
      t503 = t3*(t419-t418-3)
      t508 = t113*t31
      t510 = 4*t328
      t516 = t28*l1*(t349-t291)
      t521 = 11*t31
      t524 = t3*(t40+t521+t357-36*t6)
      t531 = t3*(3+3*l1)
      t533 = 1536*t341
      t535 = t56*(t292-t30+t7)
      t547 = t3*(28*t5-54*t6+t521+10*l1+6*t41)
      t554 = t3*(t521+t239-38*t6+2+t311)
      t559 = -768*(5*t162-8*t163)*t5*t23+384*(t405-4*t173+t175)*t36+384*
     #(t413-4*t233+t235)*t11-384*(t421-4*t242+t244)*t10-384*t427+t429-t2
     #52-384*(t434-4*t183+t185)*t62+384*(-26+t253)*t227*t1+(384*(t445-4*
     #t313+t315)*t11+384*(t452-4*t267+t269)*t36-384*(t460-4*t295+t297)*t
     #62+384*(t467-4*t276+t278)*t21-384*(t473-4*t284+t286)*t10-384*t479+
     #t481-t323-384*(t486-4*t305+t307)*t22)*s1l1+(384*(t497-4*t366+t368)
     #*t11-384*(t503-4*t344+t346)*t10+384*(3*t508-t510+t329)*t36-384*(t5
     #16-4*t352+t354)*t62+384*(t524-4*t335+t337)*t21-384*t531+t533-384*(
     #t535-4*t360+t362)*t22)*s2l1-384*(t547-4*t258+t260)*t22+384*(t554-4
     #*t379+t381)*t21
      t653 = 384*(32-t134)*t227*t1+(384*(-4*t452+4*t267)*t36-384*(-4*t46
     #0+4*t295)*t62-384*(-4*t486+4*t305)*t22+384*(-4*t467+4*t276)*t21+38
     #4*(-4*t445+4*t313)*t11-384*(-4*t473+4*t284)*t10+1536*t479-t481)*s1
     #l1+(384*(-12*t508+t510)*t36-384*(-4*t516+4*t352)*t62-384*(-4*t535+
     #4*t360)*t22+384*(-4*t524+4*t335)*t21+384*(-4*t497+4*t366)*t11-384*
     #(-4*t503+4*t344)*t10+1536*t531-t533)*s2l1-768*(-8*t162+4*t163)*t5*
     #t23+384*(-4*t405+4*t173)*t36-384*(-4*t434+4*t183)*t62-384*(-4*t547
     #+4*t258)*t22+384*(-4*t554+4*t379)*t21+384*(-4*t413+4*t233)*t11-384
     #*(-4*t421+4*t242)*t10+1536*t427-t429
      t657 = (l1+z-t198)**2
      taux26 = ((384*t3*t9*t11+768*t16*t18*t20*t23+384*t28*t35*t36-384*t
     #3*t46*t22+(384*t20*t53*t36-384*t56*t61*t62+384*t3*t68*t21+384*t3*t
     #73*t11-384*t78-384*t3*t82*t22-384*t3*t86*t10)*s1l1-384*t3*t93*t10-
     #384*t98-384*t56*t102*t62+384*t3*t109*t21+(1536*t113*t5*t36-384*t28
     #*t118*t62+384*t3*t122*t21+384*t3*t126*t11-384*t3*t130*t10-384*t56*
     #t135*t22)*s2l1)*(t144*t146/2-(t149-t150)*t145-t153*t155/2+t158*t15
     #4)+t385*(-t386*t146/2-t144*t145+t149-t150-t390*t155/2+t153*t154)+t
     #559*(t386*t145+t142-t143+t390*t154)-t653*t386)/t657*t201

CCC Coefficients

      tcoff49 = taux25
      tcoff50 = taux25
      tcoff51 = taux26
      tcoff52 = taux26

CCC Jet functions

      tFjet49 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff49)
      tFjet50 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff50)
      tFjet51 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff51)
      tFjet52 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff52)

CCC Sum

      qgRVnum13e3 = tFjet49*tcoff49+tFjet50*tcoff50+tFjet51*tcoff51+tFje
     #t52*tcoff52
      RETURN
      END

      double precision function qgRVnum14e0(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux27 = 0
      taux28 = 0

CCC Coefficients

      tcoff53 = taux27
      tcoff54 = taux27
      tcoff55 = taux28
      tcoff56 = taux28

CCC Jet functions

      tFjet53 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff53)
      tFjet54 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff54)
      tFjet55 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff55)
      tFjet56 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff56)

CCC Sum

      qgRVnum14e0 = tFjet53*tcoff53+tFjet54*tcoff54+tFjet55*tcoff55+tFje
     #t56*tcoff56
      RETURN
      END

      double precision function qgRVnum14e1(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -1+l1
      t4 = 2*l6-1
      t5 = t4**2
      t8 = l5-1
      t10 = z*l1
      t12 = (1-l1+t10)**2
      t13 = 1/t12
      t18 = dsqrt(z)
      t20 = dsqrt(l5)
      t22 = dsqrt(-t8)
      t25 = dsqrt(l1)
      t29 = l5*z*l1-l1+1-l5+l5*l1
      t31 = -t1
      t32 = dsqrt(t31)
      t37 = t29**2
      t39 = 4*l1*t5*z*l5*t8*t1*t13-4*t4*t18*t20*t22*t25*t29*t13*t32+t37*
     #t13
      t41 = l1**2
      t42 = 128*t41
      t43 = 64*l1
      t49 = 192*l1
      t51 = dsqrt(t39)
      t55 = 1-z
      t57 = 3*t41
      t62 = s1l1**2
      t64 = t41*t1
      t68 = 128*l1
      t69 = t41*l1
      t74 = t55**2
      t78 = 3*l1
      t84 = 192*t69
      t87 = 2*l1
      t96 = s2l1**2
      t98 = l1*t1
      t125 = 64*t1*t39+((64-t42+t43)*s1l1+(t42-320*l1+192)*s2l1+192-t49)
     #*t51+(192*t41*t31*t55-64*l1*(-t57+2+l1))*t62+((-256*t64*t55-256*t4
     #1+t68-128+256*t69)*s2l1-64*t64*(-3+l1)*t74+64*l1*(3+t69+t78-7*t41)
     #*t55-128-t49+t42+t84)*s1l1+(64*t1*(-t57+t87-1)*t55+64*t1*(t57-4*l1
     #+3))*t96+(-64*t98*(-t41-1)*t74-64*t1*(t69+t57+2)*t55+t84+576*l1-38
     #4*t41-384)*s2l1-192+t68+64*t69+(64+t43-128*t69)*t55+64*t64*(-t41+t
     #87-2)*t74*t55+64*t98*(-2*t41+t78+t69+1)*t74
      t128 = dlog(1-l3*l1)
      t132 = dlog(1-l1+t10-z*l3*l1)
      t133 = dlog(z)
      taux27 = t125*(t128-t132+t133)
      t1 = l1**2
      t2 = 2*l1
      t5 = -1+l1
      t6 = 1-z
      t7 = t6**2
      t8 = t5*t7
      t10 = 2*t1
      t14 = t5*t6
      t19 = -2*l1*t5**2
      t54 = dlog(1-l3*l1)
      t59 = dlog(1-l1+z*l1-z*l3*l1)
      t60 = dlog(z)
      taux28 = ((-128*t1*(6-t2)*t8+128*l1*(-t10+2+4*l1)*t14-128*t19)*s1l
     #1+(-128*l1*(2*t1+2)*t8+128*l1*(t10-t2+4)*t14-128*t19)*s2l1-256*t1*
     #t5*(-t1+t2-2)*t7*t6-128*l1*(2+6*l1-4*t1+2*t1*l1)*t8+512*l1*t5*t6-1
     #28*t19)*(t54-t59+t60)

CCC Coefficients

      tcoff53 = taux27
      tcoff54 = taux27
      tcoff55 = taux28
      tcoff56 = taux28

CCC Jet functions

      tFjet53 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff53)
      tFjet54 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff54)
      tFjet55 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff55)
      tFjet56 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff56)

CCC Sum

      qgRVnum14e1 = tFjet53*tcoff53+tFjet54*tcoff54+tFjet55*tcoff55+tFje
     #t56*tcoff56
      RETURN
      END

      double precision function qgRVnum14e2(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -1+l1
      t4 = 2*l6-1
      t5 = t4**2
      t8 = l5-1
      t10 = z*l1
      t12 = (1-l1+t10)**2
      t13 = 1/t12
      t18 = dsqrt(z)
      t20 = dsqrt(l5)
      t22 = dsqrt(-t8)
      t25 = dsqrt(l1)
      t29 = l5*z*l1-l1+1-l5+l5*l1
      t31 = -t1
      t32 = dsqrt(t31)
      t37 = t29**2
      t39 = 4*l1*t5*z*l5*t8*t1*t13-4*t4*t18*t20*t22*t25*t29*t13*t32+t37*
     #t13
      t41 = 128*l1
      t42 = 2*l1
      t43 = l1**2
      t44 = 3*t43
      t47 = 1-z
      t49 = 4*l1
      t54 = s2l1**2
      t56 = 64*l1
      t57 = 128*t43
      t60 = 192*l1
      t61 = 320*l1
      t65 = dsqrt(t39)
      t67 = t43*t1
      t68 = -t43+t42-2
      t69 = t47**2
      t70 = t69*t47
      t74 = t1*l1
      t75 = -t43-1
      t79 = t43*l1
      t84 = 192*t79
      t86 = 384*t43
      t96 = s1l1**2
      t98 = 2*t43
      t99 = 3*l1
      t100 = -t98+t99+t79+1
      t106 = 256*t43
      t107 = 256*t79
      t117 = t67*(-3+l1)*t69
      t125 = 64*t1*t39+t41+(64*t1*(t42-t44-1)*t47+64*t1*(3-t49+t44))*t54
     #+((64+t56-t57)*s1l1-t60+(-t61+t57+192)*s2l1+192)*t65+64*t67*t68*t7
     #0+(-64*t74*t75*t69-384-64*t1*(2+t79+t44)*t47+t84+576*l1-t86)*s2l1+
     #(192*t43*t31*t47-64*l1*(2+l1-t44))*t96+64*t74*t100*t69-192+((-256*
     #t67*t47-t106+t41-128+t107)*s2l1+t84+64*l1*(-7*t43+3+t99+t79)*t47+t
     #57-64*t117-128-t60)*s1l1+64*t79+(64-128*t79+t56)*t47
      t128 = dlog(1-l3*l1)
      t129 = t128**2
      t134 = dlog(1-l1+t10-z*l3*l1)
      t135 = t134**2
      t137 = dlog(z)
      t139 = t137**2
      t143 = t43**2
      t160 = 2*t79
      t173 = 512*t79
      t196 = 4*t43
      t225 = 64*(t143*t1-3*t67*t68)*t70-t41+128*t31*t39+((-128-t41+t106)
     #*s1l1+64*t43+(-t57-384+512*l1)*s2l1+t61-384)*t65-64*(-t74*(-t49+t9
     #8-t160)+2*t74*t100)*t69+((t86+128*t74*(2*t43+2*l1)*t47-t41+256-t17
     #3)*s2l1+64*l1*(19*t43-4-13*l1-t160)*t47+128*t117-704*t79+256+t57+t
     #61)*s1l1+(64*t43*(t49+t98-6)*t47-64*l1*(10*t43-8*l1-2))*t96-t107+(
     #64*t1*(t160+2-t42+t196)*t47+64*t1*(-6-t196+t49))*t54+(768+64*(-2*t
     #79*t1+2*t74*t75)*t69-320*t79+320*t43-64*t1*(-t99-4*t79-4-5*t43)*t4
     #7-768*l1)*s2l1+384+(-t60-128-192*t43+t173)*t47
      taux27 = t125*(t129/2-t135/2+t137*t134-t139/2)+t225*(t128-t134+t13
     #7)
      t1 = l1**2
      t2 = 2*l1
      t4 = t1*(6-t2)
      t5 = -1+l1
      t6 = 1-z
      t7 = t6**2
      t8 = t5*t7
      t10 = 2*t1
      t11 = 4*l1
      t13 = l1*(-t10+2+t11)
      t14 = t5*t6
      t18 = -2*l1*t5
      t19 = t18*t5
      t25 = l1*(2*t1+2)
      t28 = l1*(t10-t2+4)
      t34 = -t1+t2-2
      t35 = t7*t6
      t39 = 6*l1
      t40 = 4*t1
      t41 = t1*l1
      t44 = l1*(2+t39-t40+2*t41)
      t54 = dlog(1-l3*l1)
      t55 = t54**2
      t61 = dlog(1-l1+z*l1-z*l3*l1)
      t62 = t61**2
      t64 = dlog(z)
      t66 = t64**2
      t78 = (-t2-4*t18)*t5
      t121 = t1**2
      taux28 = ((-128*t4*t8+128*t13*t14-128*t19)*s1l1+(-128*t25*t8+128*t
     #28*t14-128*t19)*s2l1-256*t1*t5*t34*t35-128*t44*t8+512*l1*t5*t6-128
     #*t19)*(t55/2-t62/2+t64*t61-t66/2)+((128*(l1*(-2+t10-t39)-4*t13)*t5
     #*t6-128*t78-128*(t1*(t11-t10-6)-4*t4)*t5*t7)*s1l1+128*(l1*(-4*l1-4
     #)-16*l1)*t5*t6-128*(l1*(-2-8*l1-4*t41+t40)-4*t44)*t5*t7+(128*(l1*(
     #-4-t2-t10)-4*t28)*t5*t6-128*t78-128*(l1*(-2*t1-2*t41-2-2*l1)-4*t25
     #)*t5*t7)*s2l1+256*(-t121+5*t1*t34)*t5*t35-128*t78)*(t54-t61+t64)

CCC Coefficients

      tcoff53 = taux27
      tcoff54 = taux27
      tcoff55 = taux28
      tcoff56 = taux28

CCC Jet functions

      tFjet53 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff53)
      tFjet54 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff54)
      tFjet55 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff55)
      tFjet56 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff56)

CCC Sum

      qgRVnum14e2 = tFjet53*tcoff53+tFjet54*tcoff54+tFjet55*tcoff55+tFje
     #t56*tcoff56
      RETURN
      END

      double precision function qgRVnum14e3(var,ndim,z,ux,l1,
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
      double precision s13a,s23a,s1v,s2v
      double precision s13e0,s13e1,s13e2,s13e3
      external fillspin,s13e0,s13e1,s13e2,s13e3
      double precision sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      external sqe1,sqe2,zqe0,zqe1,zqe2,sbqe0,sbqe1,sbqe2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      s13a = -(1.0d0-l1)*(1.0d0-z)
      s23a = -l1*(1.0d0-z)
      s1v = s23a
      s2v = s13a
      call fillspin(z,s13a,0.0d0,s23a,0.0d0,0.0d0,s1v,s2v,1.0d0,z,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = -1+l1
      t4 = 2*l6-1
      t5 = t4**2
      t8 = l5-1
      t10 = z*l1
      t12 = (1-l1+t10)**2
      t13 = 1/t12
      t16 = l1*t5*z*l5*t8*t1*t13
      t18 = dsqrt(z)
      t20 = dsqrt(l5)
      t22 = dsqrt(-t8)
      t25 = dsqrt(l1)
      t29 = l5*z*l1-l1+1-l5+l5*l1
      t31 = -t1
      t32 = dsqrt(t31)
      t35 = t4*t18*t20*t22*t25*t29*t13*t32
      t37 = t29**2
      t38 = t37*t13
      t39 = 4*t16-4*t35+t38
      t41 = 128*l1
      t42 = 2*l1
      t43 = l1**2
      t44 = 3*t43
      t47 = 1-z
      t49 = 4*l1
      t54 = s2l1**2
      t56 = 64*l1
      t57 = 128*t43
      t60 = 192*l1
      t61 = 320*l1
      t65 = dsqrt(t39)
      t67 = t43*t1
      t68 = -t43+t42-2
      t69 = t47**2
      t70 = t69*t47
      t74 = l1*t1
      t75 = -t43-1
      t79 = t43*l1
      t84 = 192*t79
      t86 = 384*t43
      t96 = s1l1**2
      t98 = 2*t43
      t99 = 3*l1
      t100 = -t98+t99+t79+1
      t106 = 256*t43
      t107 = 256*t79
      t117 = t67*(-3+l1)*t69
      t125 = 64*t1*t39+t41+(64*t1*(t42-t44-1)*t47+64*t1*(3-t49+t44))*t54
     #+((64+t56-t57)*s1l1-t60+(-t61+t57+192)*s2l1+192)*t65+64*t67*t68*t7
     #0+(-64*t74*t75*t69-384-64*t1*(2+t79+t44)*t47+t84+576*l1-t86)*s2l1+
     #(192*t43*t31*t47-64*l1*(2+l1-t44))*t96+64*t74*t100*t69-192+((-256*
     #t67*t47-t106+t41-128+t107)*s2l1+t84+64*l1*(-7*t43+3+t99+t79)*t47+t
     #57-64*t117-128-t60)*s1l1+64*t79+(64-128*t79+t56)*t47
      t126 = dlog(z)
      t127 = t126**2
      t131 = dlog(1-l1+t10-z*l3*l1)
      t134 = t131**2
      t143 = dlog(1-l3*l1)
      t144 = t143**2
      t149 = t43**2
      t150 = t149*t1
      t151 = t67*t68
      t160 = 64*t43
      t161 = 512*l1
      t166 = 2*t79
      t168 = t74*(-t49+t98-t166)
      t179 = 512*t79
      t202 = 4*t43
      t211 = t79*t1
      t218 = 320*t43
      t219 = 4*t79
      t225 = 768*l1
      t228 = 192*t43
      t231 = 64*(t150-3*t151)*t70-t41+128*t31*t39+((-128-t41+t106)*s1l1+
     #t160+(-t57-384+t161)*s2l1+t61-384)*t65-64*(-t168+2*t74*t100)*t69+(
     #(t86+128*t74*(2*t43+2*l1)*t47-t41+256-t179)*s2l1+64*l1*(19*t43-4-1
     #3*l1-t166)*t47+128*t117-704*t79+256+t57+t61)*s1l1+(64*t43*(t49+t98
     #-6)*t47-64*l1*(10*t43-8*l1-2))*t96-t107+(64*t1*(t166+2-t42+t202)*t
     #47+64*t1*(-6-t202+t49))*t54+(768+64*(-2*t211+2*t74*t75)*t69-320*t7
     #9+t218-64*t1*(-t99-t219-4-5*t43)*t47-t225)*s2l1+384+(-t60-128-t228
     #+t179)*t47
      t247 = t43+1-l1
      t268 = 6*t43
      t305 = -256*t16+256*t35-64*t38+((-t60+64+t218)*s1l1+(t160-t60+192)
     #*s2l1-t60+t228+192)*t65+(64*t43*t247*t47-64*l1*(t99+t98-2))*t96+((
     #128*t74*t247*t47+t161-512*t43-128)*s2l1+64*l1*(t49-2)*t47+t161-768
     #*t43-128)*s1l1+(64*t1*(t79-t268+t42-1)*t47+64*t1*(t268-t99+3))*t54
     #+(256*t211*t69-64*t1*(t219+8*t43-t49+2)*t47+t225+t179-1024*t43-384
     #)*s2l1+448*t79+384*l1-832*t43+(-t60-t179+896*t43-192*t149+64)*t47-
     #192+64*(-3*t150+2*t151)*t70-64*(-t74*(t44+t79)+2*t168)*t69
      taux27 = t125*(-t127*t131/2-t134*t131/6+t127*t126/6+t126*t134/2+t1
     #44*t143/6)+t231*(t144/2-t134/2+t126*t131-t127/2)+t305*(t143-t131+t
     #126)
      t1 = l1**2
      t2 = 2*l1
      t4 = t1*(6-t2)
      t5 = -1+l1
      t6 = 1-z
      t7 = t6**2
      t8 = t5*t7
      t10 = 2*t1
      t11 = 4*l1
      t13 = l1*(-t10+2+t11)
      t14 = t5*t6
      t16 = -t5
      t18 = 2*l1*t16
      t19 = t18*t5
      t28 = l1*(2*t1+2)
      t31 = l1*(t10-t2+4)
      t37 = -t1+t2-2
      t38 = t7*t6
      t42 = 6*l1
      t43 = 4*t1
      t44 = l1*t1
      t45 = 2*t44
      t47 = l1*(2+t42-t43+t45)
      t52 = dlog(z)
      t53 = t52**2
      t54 = z*l1
      t58 = dlog(1-l1+t54-z*l3*l1)
      t61 = t58**2
      t70 = dlog(1-l3*l1)
      t71 = t70**2
      t77 = l1*(-2+t10-t42)
      t78 = 4*t13
      t82 = 4*t18
      t84 = (-t2-t82)*t5
      t86 = t1*(t11-t10-6)
      t87 = 4*t4
      t96 = l1*(-4*l1-4)
      t97 = 16*l1
      t102 = 8*l1
      t105 = l1*(-2-t102-4*t44+t43)
      t106 = 4*t47
      t112 = l1*(-4-t2-t10)
      t113 = 4*t31
      t119 = l1*(-2*t1-2*t44-2-2*l1)
      t120 = 4*t28
      t127 = t1**2
      t128 = t1*t37
      t143 = 2*l6-1
      t144 = t143**2
      t147 = l5-1
      t150 = (1-l1+t54)**2
      t151 = 1/t150
      t156 = dsqrt(z)
      t158 = dsqrt(l5)
      t160 = dsqrt(-t147)
      t163 = dsqrt(l1)
      t167 = l5*z*l1-l1+1-l5+l5*l1
      t169 = dsqrt(t16)
      t174 = t167**2
      t177 = dsqrt(4*l1*t144*z*l5*t147*t5*t151-4*t143*t156*t158*t160*t16
     #3*t167*t151*t169+t174*t151)
      t225 = 7*t1
      s1 = ((-128*t4*t8+128*t13*t14-128*t19)*s1l1+512*l1*t5*t6+(-128*t28
     #*t8+128*t31*t14-128*t19)*s2l1-256*t1*t5*t37*t38-128*t47*t8-128*t19
     #)*(-t53*t58/2-t61*t58/6+t53*t52/6+t52*t61/2+t71*t70/6)
      s2 = ((128*(t77-t78)*t5*t6-128*t84-128*(t86-t87)*t5*t7)*s1l1+128*(
     #t96-t97)*t5*t6-128*(t105-t106)*t5*t7+(128*(t112-t113)*t5*t6-128*t8
     #4-128*(t119-t120)*t5*t7)*s2l1+256*(-t127+5*t128)*t5*t38-128*t84)*(
     #t71/2-t61/2+t52*t58-t53/2)+(384*l1*t177*t5+(-128*(3*t127-4*t86+t87
     #)*t5*t7+128*(11*t44-4*t77+t78)*t5*t6-128*(11*t1+t102+t82)*t5)*s1l1
     #+(-128*(l1*(3*t44+t1)-4*t119+t120)*t5*t7+128*(l1*(t2+9*t1)-4*t112+
     #t113)*t5*t6-128*(l1*(5*l1+3)+t102+t82)*t5)*s2l1+256*(5*t127-8*t128
     #)*t5*t38-128*(l1*(t45+t225)-4*t105+t106)*t5*t7+128*(l1*(t102+t225)
     #-4*t96+t97)*t5*t6-128*(l1*(3+t102)+t102+t82)*t5)*(t70-t58+t52)
      taux28 = s1+s2

CCC Coefficients

      tcoff53 = taux27
      tcoff54 = taux27
      tcoff55 = taux28
      tcoff56 = taux28

CCC Jet functions

      tFjet53 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoeff
     #*tcoff53)
      tFjet54 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,5.D0,1.D0,1.D0,Tcoe
     #ff*tcoff54)
      tFjet55 = Fjet(var,ndim,z,ux,s13a,s23a,0.D0,0.D0,0.D0,s23a,s13a,1.
     #D0,z,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoeff
     #*tcoff55)
      tFjet56 = Fjet(var,ndim,z,1/ux,s23a,s13a,0.D0,0.D0,0.D0,s13a,s23a,
     #1.D0,z,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,7.D0,1.D0,1.D0,Tcoe
     #ff*tcoff56)

CCC Sum

      qgRVnum14e3 = tFjet53*tcoff53+tFjet54*tcoff54+tFjet55*tcoff55+tFje
     #t56*tcoff56
      RETURN
      END

