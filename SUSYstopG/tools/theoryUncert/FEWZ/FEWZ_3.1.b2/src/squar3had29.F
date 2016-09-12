CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function sq39e0int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision sq3num5e0
      external sq3num5e0
      double precision sq3num5e1
      external sq3num5e1
      double precision sq3num5e2
      external sq3num5e2
      double precision sq3num5e3
      external sq3num5e3
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = 0.3141592653589793D1**2
      t2 = 2*t1
      t3 = dlog(2.D0)
      t4 = t3**2
      t5 = 12*t4
      t6 = f2(0.D0,1.D0)
      t7 = dlog(t6)
      t8 = t7**2
      t9 = 3*t8
      t11 = dlog(1-l4)
      t12 = dlog(l1)
      t13 = t11*t12
      t14 = 6*t13
      t15 = f3(0.D0,1.D0)
      t17 = dlog(1-t15)
      t18 = t11*t17
      t19 = 6*t18
      t20 = t12*t17
      t21 = 6*t20
      t23 = 6*t11*t7
      t24 = dlog(l4)
      t25 = t24*t12
      t26 = 6*t25
      t28 = 6*t24*t7
      t29 = 1-z
      t30 = dlog(t29)
      t34 = 6*t12*t7
      t36 = 12*t17*t3
      t37 = t17*t7
      t38 = 6*t37
      t40 = 12*t3*t7
      t41 = -t2+t5+t9+t14+t19+t21-t23+t26-t28+24*t30*t17-t34+t36-t38-t40
      t42 = t24*t17
      t43 = 6*t42
      t46 = t12*t3
      t47 = 12*t46
      t50 = t24*t3
      t51 = 12*t50
      t54 = t11*t3
      t55 = 12*t54
      t56 = t11*t24
      t57 = 6*t56
      t60 = t17**2
      t61 = 3*t60
      t64 = t30**2
      t66 = t24**2
      t67 = 3*t66
      t68 = t11**2
      t69 = 3*t68
      t70 = t12**2
      t71 = 3*t70
      t72 = t43-24*t30*t7+t47+24*t30*t12+t51+24*t11*t30+t55+t57+48*t30*t
     #3+t61+24*t24*t30+48*t64+t67+t69+t71
      t75 = -1/t29
      t77 = Dsquar311(1.D0,l1,0.D0,0.D0,l4)
      t78 = 1/t77
      t80 = s2vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t81 = 1/t80
      t82 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t83 = 1/t82
      t85 = s1vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t86 = 1/t85
      t87 = f1(0.D0,1.D0,1.D0)
      t89 = t81*t83*t86*t87
      t92 = 1/l2
      t94 = 1/l3
      t95 = t94*t78
      t99 = -t2+t5+t9+t14+t19+t21-t23+t26-t28-t34+t36-t38-t40+t43
      t100 = dlog(l3)
      t101 = t100**2
      t115 = t47+t51+t55+t57+t61+t67+t69+t71+12*t101+12*t100*t17-12*t100
     #*t7+12*t11*t100+12*t24*t100+12*t100*t12+24*t100*t3
      t150 = t25*t17/8192+t60*t17/49152+t70*t12/49152-t54*t7/4096-t13*t7
     #/8192+t13*t3/4096+t56*t3/4096-t8*t7/49152+z3/3072-t1*t3/12288+t20*
     #t3/4096+t25*t3/4096+t56*t12/8192-t18*t7/8192+t68*t11/49152
      t183 = -t11*t1/24576+t11*t60/16384+t11*t8/16384+t24*t4/4096+t24*t7
     #0/16384+t66*t24/49152-t46*t7/4096+t42*t3/4096+t18*t3/4096-t25*t7/8
     #192-t50*t7/4096-t20*t7/8192+t4*t3/6144+t11*t4/4096+t66*t17/16384-t
     #7*t4/4096
      t217 = t11*t70/16384+t68*t17/16384+t24*t8/16384+t68*t24/16384+t68*
     #t3/8192-t68*t7/16384+t68*t12/16384+t12*t4/4096+t12*t8/16384-t66*t7
     #/16384+t70*t17/16384+t70*t3/8192-t70*t7/16384+t12*t60/16384+t66*t3
     #/8192+t66*t12/16384
      t250 = -t17*t1/24576-t12*t1/24576+t60*t3/8192+t17*t8/16384+t17*t4/
     #4096+t7*t1/24576+t8*t3/8192-t24*t1/24576+t56*t17/8192-t56*t7/8192+
     #t13*t17/8192+t24*t60/16384+t11*t66/16384-t60*t7/16384-t37*t3/4096-
     #t42*t7/8192
      t258 = dlog(l2)
      t259 = 2*t258
      t260 = 2*t3
      t261 = 2*t100
      t267 = 4*t30
      t270 = t75*t78
      t276 = t258**2
      t282 = -t2+t5+t9+t14+t19+t21-t23+24*t258*t3+12*t276+t26-t28-t34+12
     #*t11*t258+12*t24*t258
      t289 = 12*t258*t12+12*t258*t17-12*t258*t7+t36-t38-t40+t43+t47+t51+
     #t55+t57+t61+t67+t69+t71
      tncoff1 = (t41+t72)*t75*t78*t89/12288+t75*t92*t95*t89/512-(t99+t11
     #5)*t94*t78*t89/24576+(t150+t183+t217+t250)*t78*t81*t83*t86*t87+(t1
     #1+t259+t260+t17-t7+t12+t261+t24)*t92*t95*t89/2048-(t12+t260+t17-t7
     #+t259+t267+t11+t24)*t92*t270*t89/1024-(t282+t289)*t92*t78*t89/2457
     #6-(t261+t12+t260+t17-t7+t24+t267+t11)*t94*t270*t89/1024
      t1 = dlog(l3)
      t2 = 2*t1
      t3 = 2-l3
      t4 = dlog(t3)
      t5 = 2*t4
      t6 = l3/2
      t7 = f3(t6,1.D0)
      t9 = dlog(1-t7)
      t10 = dlog(l1)
      t11 = f2(t6,1.D0)
      t12 = dlog(t11)
      t13 = dlog(l4)
      t14 = 1-z
      t15 = dlog(t14)
      t18 = dlog(1-l4)
      t21 = -1/t3
      t23 = 1/l3
      t25 = -1/t14
      t28 = Dsquar311(1.D0,l1,0.D0,t6,l4)
      t29 = 1/t28
      t30 = s2vr2(1.D0,1.D0,1.D0,t6,l4)
      t31 = 1/t30
      t33 = s23r2(1.D0,1.D0,1.D0,t6,l4)
      t34 = 1/t33
      t35 = s1vr2(1.D0,1.D0,1.D0,t6,l4)
      t36 = 1/t35
      t38 = f1(t6,1.D0,1.D0)
      t40 = t29*t31*t34*t36*t38
      t43 = dlog(l2)
      t47 = 1/l2
      t52 = t18**2
      t54 = t13**2
      t56 = t10**2
      t60 = t9**2
      t62 = t4**2
      t64 = t12**2
      t66 = t1**2
      t80 = 3*t52+3*t54+3*t56+24*t1*t4+3*t60+12*t62+3*t64+12*t66+12*t13*
     #t4-6*t13*t12+12*t18*t4+12*t10*t4-6*t10*t12+12*t9*t4
      t109 = 0.3141592653589793D1**2
      t111 = 6*t10*t9+6*t18*t9-6*t18*t12-12*t12*t1-12*t12*t4+12*t1*t9+6*
     #t13*t9+6*t18*t13+12*t1*t10-6*t9*t12+6*t18*t10+6*t13*t10+12*t13*t1+
     #12*t18*t1-2*t109
      t116 = t31*t34
      tncoff2 = -(t2+t5+t9+t10-t12+t13+4*t15+t18)*t21*t23*t25*t40/512+(t
     #18+2*t43+t5-t12+t10+t9+t2+t13)*t21*t47*t23*t40/1024-(t80+t111)*t21
     #*t23*t29*t116*t36*t38/12288+t21*t47*t23*t36*t116*t29*t25*t38/256
      t1 = -1+z
      t2 = 1/t1
      t3 = 1/l2
      t6 = Dsquar311(1.D0,l1,l2,0.D0,l4)
      t7 = 1/t6
      t8 = 1/l3*t7
      t11 = 1-l1*l2
      t12 = 1-l2
      t13 = s2vr2(1.D0,t11,t12,0.D0,l4)
      t15 = s23r2(1.D0,t11,t12,0.D0,l4)
      t18 = s1vr2(1.D0,t11,t12,0.D0,l4)
      t20 = f1(0.D0,1.D0,t11)
      t22 = 1/t13/t15/t18*t20
      t25 = dlog(t12)
      t26 = t25**2
      t28 = f2(0.D0,1.D0)
      t29 = dlog(t28)
      t30 = t29**2
      t32 = dlog(l4)
      t33 = dlog(t11)
      t36 = dlog(l2)
      t37 = t36**2
      t39 = 0.3141592653589793D1**2
      t41 = dlog(2.D0)
      t42 = t41**2
      t45 = dlog(1-l4)
      t46 = t45**2
      t50 = f3(0.D0,1.D0)
      t54 = dlog(1-t50+t50*l1*l2)
      t55 = t54**2
      t57 = t33**2
      t61 = -3*t26-3*t30-6*t32*t33-12*t37+2*t39-12*t42-3*t46-6*t45*t32-3
     #*t55-3*t57-12*t32*t41
      t66 = t32**2
      t68 = dlog(l1)
      t69 = t68**2
      t87 = -12*t25*t41+6*t54*t29-3*t66-3*t69-12*t45*t41-6*t25*t54-12*t3
     #3*t41-6*t45*t54-6*t45*t25-6*t45*t33-6*t68*t33-12*t36*t33
      t111 = -6*t25*t33-6*t68*t54+6*t33*t29-12*t36*t54-12*t36*t25+6*t25*
     #t29-6*t32*t54-6*t32*t25-6*t54*t33-12*t54*t41-6*t68*t25
      t136 = -6*t45*t68+6*t45*t29-6*t32*t68+6*t32*t29+6*t68*t29+12*t41*t
     #29-12*t68*t41-24*t36*t41-12*t45*t36-12*t32*t36-12*t36*t68+12*t36*t
     #29
      t143 = 2*t41
      t144 = 2*t36
      t145 = dlog(l3)
      t153 = dlog(-t1)
      tncoff3 = -t2*t3*t8*t22/512-(t61+t87+t111+t136)*t3*t7*t22/24576+(-
     #t45-t143-t25-t144-t54+t29-t68-t33-2*t145-t32)*t3*t8*t22/2048-(-t54
     #-t25-t143-t68-t33+t29-t144-4*t153-t45-t32)*t3*t2*t7*t22/1024
      t1 = -2+l3
      t2 = 1/t1
      t3 = 1/l2
      t5 = 1/l3
      t7 = 1-l1*l2
      t8 = 1-l2
      t9 = l3/2
      t10 = s1vr2(1.D0,t7,t8,t9,l4)
      t11 = 1/t10
      t14 = s2vr2(1.D0,t7,t8,t9,l4)
      t15 = 1/t14
      t16 = s23r2(1.D0,t7,t8,t9,l4)
      t17 = 1/t16
      t19 = Dsquar311(1.D0,l1,l2,t9,l4)
      t20 = 1/t19
      t24 = f1(t9,1.D0,t7)
      t30 = dlog(1-l4)
      t31 = dlog(l1)
      t33 = dlog(-t1)
      t35 = f2(t9,1.D0)
      t36 = dlog(t35)
      t37 = dlog(t8)
      t38 = dlog(t7)
      t39 = dlog(l2)
      t41 = f3(t9,1.D0)
      t45 = dlog(1-t41+t41*l1*l2)
      t46 = dlog(l3)
      t48 = dlog(l4)
      tncoff4 = -t2*t3*t5*t11*t15*t17*t20/(-1+z)*t24/256-(t30+t31+2*t33-
     #t36+t37+t38+2*t39+t45+2*t46+t48)*t2*t3*t5*t20*t15*t17*t11*t24/1024
      t1 = f3(0.D0,z)
      t3 = dlog(1-t1)
      t4 = dlog(z)
      t5 = 2*t4
      t6 = dlog(2.D0)
      t7 = 2*t6
      t8 = dlog(l3)
      t10 = f2(0.D0,z)
      t11 = dlog(t10)
      t12 = dlog(l1)
      t13 = dlog(l4)
      t14 = 1-z
      t15 = dlog(t14)
      t16 = 4*t15
      t18 = dlog(1-l4)
      t20 = 1/l3
      t23 = -1/t14
      t24 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t25 = 1/t24
      t26 = t23*t25
      t28 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t30 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t33 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t35 = f1(0.D0,z,1.D0)
      t37 = 1/t28/t30/t33*t35
      t76 = -12*t13*t4-24*t15*t11+6*t18*t12+12*t11*t4+6*t13*t12-6*t12*t1
     #1+6*t12*t3-12*t12*t4+12*t3*t6+24*t13*t15+24*t15*t12-48*t15*t4-24*t
     #4*t6+48*t15*t6-12*t18*t4+6*t18*t13-6*t3*t11+24*t18*t15
      t77 = 0.3141592653589793D1**2
      t83 = t6**2
      t85 = t18**2
      t87 = t3**2
      t89 = t15**2
      t91 = t11**2
      t93 = t4**2
      t95 = t12**2
      t97 = t13**2
      t115 = -2*t77-12*t3*t4-6*t18*t11+12*t83+3*t85+3*t87+48*t89+3*t91+1
     #2*t93+3*t95+3*t97+12*t12*t6-12*t11*t6-6*t13*t11+6*t18*t3+6*t13*t3+
     #24*t15*t3+12*t13*t6+12*t18*t6
      t121 = 1/l2
      t127 = dlog(l2)
      tncoff5 = -(-t3+t5-t7-2*t8+t11-t12-t13-t16-t18)*t20*t26*t37/1024-(
     #t76+t115)*t23*t25*t37/12288-t23*t121*t20*t25*t37/512-(-2*t127-t12-
     #t3+t11-t7+t5-t16-t18-t13)*t121*t26*t37/1024
      t1 = 2-l3
      t2 = dlog(t1)
      t4 = l3/2
      t5 = f3(t4,z)
      t7 = dlog(1-t5)
      t8 = dlog(l1)
      t9 = dlog(l3)
      t11 = dlog(z)
      t13 = f2(t4,z)
      t14 = dlog(t13)
      t15 = dlog(l4)
      t16 = 1-z
      t17 = dlog(t16)
      t20 = dlog(1-l4)
      t23 = -1/t1
      t28 = -1/l3/t16
      t30 = Dsquar311(z,l1,0.D0,t4,l4)
      t31 = 1/t30
      t32 = s2vr2(z,1.D0,1.D0,t4,l4)
      t33 = 1/t32
      t35 = s23r2(z,1.D0,1.D0,t4,l4)
      t36 = 1/t35
      t37 = s1vr2(z,1.D0,1.D0,t4,l4)
      t38 = 1/t37
      t40 = f1(t4,z,1.D0)
      tncoff6 = (2*t2+t7+t8+2*t9-2*t11-t14+t15+4*t17+t20)*t23*t28*t31*t3
     #3*t36*t38*t40/512-t23*t38*t33*t36*t31/l2*t28*t40/256
      t1 = 1-l2
      t2 = dlog(t1)
      t3 = dlog(l1)
      t4 = f3(0.D0,z)
      t8 = dlog(1-t4+t4*l1*l2)
      t10 = 1-l1*l2
      t11 = dlog(t10)
      t12 = dlog(z)
      t14 = dlog(l2)
      t16 = dlog(2.D0)
      t18 = f2(0.D0,z)
      t19 = dlog(t18)
      t20 = 1-z
      t21 = dlog(t20)
      t24 = dlog(1-l4)
      t25 = dlog(l4)
      t26 = -t2-t3-t8-t11+2*t12-2*t14-2*t16+t19-4*t21-t24-t25
      t27 = 1/l2
      t30 = -1/t20
      t31 = Dsquar311(z,l1,l2,0.D0,l4)
      t32 = 1/t31
      t35 = s2vr2(z,t10,t1,0.D0,l4)
      t37 = s23r2(z,t10,t1,0.D0,l4)
      t40 = s1vr2(z,t10,t1,0.D0,l4)
      t42 = f1(0.D0,z,t10)
      t44 = 1/t35/t37/t40*t42
      tncoff7 = t26*t27*t30*t32*t44/1024+t30*t27/l3*t32*t44/512
      t1 = l3/2
      t3 = 1-l1*l2
      t4 = f1(t1,z,t3)
      t5 = 1-l2
      t6 = s1vr2(z,t3,t5,t1,l4)
      t9 = s2vr2(z,t3,t5,t1,l4)
      t15 = s23r2(z,t3,t5,t1,l4)
      t17 = Dsquar311(z,l1,l2,t1,l4)
      tncoff8 = t4/t6/t9/(-2+l3)/t15/t17/l2/l3/(-1+z)/256
      t1 = 1/l2
      t2 = 1/l3
      t4 = Dsquar311(1.D0,l1,0.D0,0.D0,l4)
      t5 = 1/t4
      t7 = s2vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t8 = 1/t7
      t9 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t10 = 1/t9
      t12 = s1vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t13 = 1/t12
      t14 = f1(0.D0,1.D0,1.D0)
      t16 = t8*t10*t13*t14
      t19 = -1+z
      t20 = 1/t19
      t25 = dlog(l2)
      t27 = f3(0.D0,1.D0)
      t29 = dlog(1-t27)
      t30 = 6*t29
      t31 = f2(0.D0,1.D0)
      t32 = dlog(t31)
      t33 = 6*t32
      t34 = dlog(2.D0)
      t35 = 12*t34
      t36 = dlog(l4)
      t37 = 6*t36
      t39 = dlog(1-l4)
      t40 = 6*t39
      t41 = dlog(l1)
      t42 = 6*t41
      t53 = dlog(-t19)
      t60 = dlog(l3)
      t67 = t36**2
      t69 = t39**2
      t71 = t41**2
      t89 = -t67/16384-t69/16384-t71/16384+t32*t41/8192-t41*t34/4096-t36
     #*t34/4096-t41*t29/8192+t32*t36/8192-t39*t36/8192-t29*t34/4096+t32*
     #t29/8192
      t92 = t29**2
      t94 = t32**2
      t96 = 0.3141592653589793D1**2
      t108 = t34**2
      t112 = t34*t32/4096-t92/16384-t94/16384+t96/24576-t36*t29/8192-t39
     #*t34/4096-t36*t41/8192+t39*t32/8192-t39*t29/8192-t108/4096-t39*t41
     #/8192
      tncoff9 = -t1*t2*t5*t16/2048+t1*t20*t5*t16/1024-(-12*t25-t30+t33-t
     #35-t37-t40-t42)*t1*t5*t16/24576+t2*t20*t5*t16/1024+(-t30+t33-t35-t
     #37-t40-24*t53-t42)*t20*t5*t16/12288-(-t30+t33-t35-t37-t40-t42-12*t
     #60)*t2*t5*t16/24576+(t89+t112)*t5*t8*t10*t13*t14
      t1 = dlog(l4)
      t4 = dlog(1-l4)
      t6 = dlog(l3)
      t8 = 2-l3
      t9 = dlog(t8)
      t11 = l3/2
      t12 = f2(t11,1.D0)
      t13 = dlog(t12)
      t15 = dlog(l1)
      t17 = f3(t11,1.D0)
      t19 = dlog(1-t17)
      t23 = -1/t8
      t25 = 1/l3
      t26 = Dsquar311(1.D0,l1,0.D0,t11,l4)
      t27 = 1/t26
      t28 = t25*t27
      t30 = s2vr2(1.D0,1.D0,1.D0,t11,l4)
      t32 = s23r2(1.D0,1.D0,1.D0,t11,l4)
      t35 = s1vr2(1.D0,1.D0,1.D0,t11,l4)
      t37 = f1(t11,1.D0,1.D0)
      t39 = 1/t30/t32/t35*t37
      tncoff10 = -(-6*t1-6*t4-12*t6-12*t9+6*t13-6*t15-6*t19)*t23*t28*t39
     #/12288+t23*t25/(-1+z)*t27*t39/512-t23/l2*t28*t39/1024
      t1 = dlog(2.D0)
      t4 = dlog(1-l4)
      t6 = dlog(l4)
      t8 = dlog(l2)
      t10 = f3(0.D0,1.D0)
      t14 = dlog(1-t10+t10*l1*l2)
      t16 = f2(0.D0,1.D0)
      t17 = dlog(t16)
      t19 = dlog(l1)
      t21 = 1-l2
      t22 = dlog(t21)
      t25 = 1-l1*l2
      t26 = dlog(t25)
      t29 = 1/l2
      t31 = Dsquar311(1.D0,l1,l2,0.D0,l4)
      t32 = 1/t31
      t34 = s2vr2(1.D0,t25,t21,0.D0,l4)
      t36 = s23r2(1.D0,t25,t21,0.D0,l4)
      t39 = s1vr2(1.D0,t25,t21,0.D0,l4)
      t41 = f1(0.D0,1.D0,t25)
      t43 = 1/t34/t36/t39*t41
      tncoff11 = -(12*t1+6*t4+6*t6+12*t8+6*t14-6*t17+6*t19+6*t22+6*t26)*
     #t29*t32*t43/24576+t29/l3*t32*t43/2048-t29/(-1+z)*t32*t43/1024
      t6 = l3/2
      t7 = Dsquar311(1.D0,l1,l2,t6,l4)
      t12 = 1-l1*l2
      t13 = 1-l2
      t14 = s2vr2(1.D0,t12,t13,t6,l4)
      t16 = s23r2(1.D0,t12,t13,t6,l4)
      t19 = s1vr2(1.D0,t12,t13,t6,l4)
      t21 = f1(t6,1.D0,t12)
      tncoff12 = 1/(-2+l3)/l2/l3/t7/t14/t16/t19*t21/1024
      t2 = -1+z
      t3 = 1/t2
      t5 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t6 = 1/t5
      t8 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t10 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t13 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t15 = f1(0.D0,z,1.D0)
      t17 = 1/t8/t10/t13*t15
      t21 = dlog(1-l4)
      t24 = dlog(-t2)
      t26 = dlog(l4)
      t28 = dlog(2.D0)
      t30 = f2(0.D0,z)
      t31 = dlog(t30)
      t33 = dlog(z)
      t35 = dlog(l1)
      t37 = f3(0.D0,z)
      t39 = dlog(1-t37)
      tncoff13 = -1/l3*t3*t6*t17/1024-(-6*t21-24*t24-6*t26-12*t28+6*t31+
     #12*t33-6*t35-6*t39)*t3*t6*t17/12288-1/l2*t3*t6*t17/1024
      t7 = l3/2
      t8 = Dsquar311(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,1.D0,1.D0,t7,l4)
      t14 = s23r2(z,1.D0,1.D0,t7,l4)
      t17 = s1vr2(z,1.D0,1.D0,t7,l4)
      t19 = f1(t7,z,1.D0)
      tncoff14 = -1/(-2+l3)/l3/(-1+z)/t8/t12/t14/t17*t19/512
      t5 = Dsquar311(z,l1,l2,0.D0,l4)
      t9 = 1-l1*l2
      t10 = 1-l2
      t11 = s2vr2(z,t9,t10,0.D0,l4)
      t13 = s23r2(z,t9,t10,0.D0,l4)
      t16 = s1vr2(z,t9,t10,0.D0,l4)
      t18 = f1(0.D0,z,t9)
      tncoff15 = 1/l2/(-1+z)/t5/t11/t13/t16*t18/1024
      t3 = Dsquar311(1.D0,l1,0.D0,0.D0,l4)
      t4 = 1/t3
      t6 = s2vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t7 = 1/t6
      t9 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t11 = s1vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t14 = f1(0.D0,1.D0,1.D0)
      t15 = 1/t9/t11*t14
      t18 = f2(0.D0,1.D0)
      t19 = dlog(t18)
      t21 = dlog(2.D0)
      t23 = dlog(l1)
      t25 = f3(0.D0,1.D0)
      t27 = dlog(1-t25)
      t29 = dlog(l4)
      t32 = dlog(1-l4)
      tncoff16 = 1/(-1+z)*t4*t7*t15/2048+(-t19/8192+t21/4096+t23/8192+t2
     #7/8192+t29/8192+t32/8192)*t4*t7*t15-1/l3*t4*t7*t15/4096-1/l2*t4*t7
     #*t15/4096
      t5 = l3/2
      t6 = Dsquar311(1.D0,l1,0.D0,t5,l4)
      t9 = s2vr2(1.D0,1.D0,1.D0,t5,l4)
      t11 = s23r2(1.D0,1.D0,1.D0,t5,l4)
      t14 = s1vr2(1.D0,1.D0,1.D0,t5,l4)
      t16 = f1(t5,1.D0,1.D0)
      tncoff17 = -1/(-2+l3)/l3/t6/t9/t11/t14*t16/2048
      t2 = Dsquar311(1.D0,l1,l2,0.D0,l4)
      t6 = 1-l1*l2
      t7 = 1-l2
      t8 = s2vr2(1.D0,t6,t7,0.D0,l4)
      t11 = s23r2(1.D0,t6,t7,0.D0,l4)
      t13 = s1vr2(1.D0,t6,t7,0.D0,l4)
      t16 = f1(0.D0,1.D0,t6)
      tncoff18 = 1/l2/t2/t8/t11/t13*t16/4096
      t3 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t6 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t9 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t11 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t14 = f1(0.D0,z,1.D0)
      tncoff19 = -1/(-1+z)/t3/t6/t9/t11*t14/2048
      t1 = Dsquar311(1.D0,l1,0.D0,0.D0,l4)
      t3 = s2vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t8 = s1vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t11 = f1(0.D0,1.D0,1.D0)
      tncoff20 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num5e0(var,ndim,1.D0,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num5e0(var,ndim,1.D0,ux,1-l1*l2,1-l2,l3/2,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num5e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num5e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num5e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num5e0(var,ndim,z,ux,1-l1*l2,1-l2,l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num5e1(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = sq3num5e1(var,ndim,1.D0,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = sq3num5e1(var,ndim,1.D0,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = sq3num5e1(var,ndim,1.D0,ux,1-l1*l2,1-l2,l3/2,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = sq3num5e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = sq3num5e1(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = sq3num5e1(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = sq3num5e2(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = sq3num5e2(var,ndim,1.D0,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = sq3num5e2(var,ndim,1.D0,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = sq3num5e2(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = sq3num5e3(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      sq39e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+
     #tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff1
     #6+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncof
     #f20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq39e1int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision sq3num5e0
      external sq3num5e0
      double precision sq3num5e1
      external sq3num5e1
      double precision sq3num5e2
      external sq3num5e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = dlog(l4)
      t2 = t1**2
      t5 = dlog(1-l4)
      t6 = t5**2
      t8 = dlog(l1)
      t9 = t8**2
      t11 = f2(0.D0,1.D0)
      t12 = dlog(t11)
      t15 = dlog(2.D0)
      t20 = f3(0.D0,1.D0)
      t22 = dlog(1-t20)
      t33 = -t2/16384-t6/16384-t9/16384+t8*t12/8192-t8*t15/4096-t1*t15/4
     #096-t8*t22/8192+t1*t12/8192-t5*t1/8192-t22*t15/4096+t22*t12/8192
      t36 = t22**2
      t38 = t12**2
      t40 = 0.3141592653589793D1**2
      t52 = t15**2
      t56 = t15*t12/4096-t36/16384-t38/16384+t40/24576-t1*t22/8192-t5*t1
     #5/4096-t1*t8/8192+t5*t12/8192-t5*t22/8192-t52/4096-t5*t8/8192
      t58 = Dsquar311(1.D0,l1,0.D0,0.D0,l4)
      t59 = 1/t58
      t61 = s2vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t62 = 1/t61
      t64 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t65 = 1/t64
      t66 = s1vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t67 = 1/t66
      t69 = f1(0.D0,1.D0,1.D0)
      t72 = 1-z
      t73 = dlog(t72)
      t75 = 2*t15
      t78 = -1/t72
      t83 = t62*t65*t67*t69
      t86 = dlog(l3)
      t89 = 1/l3
      t98 = dlog(l2)
      t101 = 1/l2
      tncoff1 = (t33+t56)*t59*t62*t65*t67*t69-(t1+4*t73+t5+t8+t75+t22-t1
     #2)*t78*t59*t83/2048+(t5+t1+2*t86+t8+t75+t22-t12)*t89*t59*t83/4096+
     #t89*t78*t59*t83/1024+(t5+t8+t22-t12+t75+2*t98+t1)*t101*t59*t83/409
     #6+t101*t78*t59*t83/1024-t101*t89*t59*t83/2048
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l3)
      t6 = dlog(l1)
      t7 = l3/2
      t8 = f2(t7,1.D0)
      t9 = dlog(t8)
      t10 = 2-l3
      t11 = dlog(t10)
      t13 = f3(t7,1.D0)
      t15 = dlog(1-t13)
      t18 = -1/t10
      t20 = 1/l3
      t21 = Dsquar311(1.D0,l1,0.D0,t7,l4)
      t22 = 1/t21
      t23 = t20*t22
      t25 = s2vr2(1.D0,1.D0,1.D0,t7,l4)
      t27 = s23r2(1.D0,1.D0,1.D0,t7,l4)
      t30 = s1vr2(1.D0,1.D0,1.D0,t7,l4)
      t32 = f1(t7,1.D0,1.D0)
      t34 = 1/t25/t27/t30*t32
      tncoff2 = (t2+t3+2*t4+t6-t9+2*t11+t15)*t18*t23*t34/2048+t18*t20/(-
     #1+z)*t22*t34/512-t18/l2*t23*t34/1024
      t2 = dlog(1-l4)
      t3 = dlog(l1)
      t4 = dlog(2.D0)
      t6 = 1-l2
      t7 = dlog(t6)
      t9 = 1-l1*l2
      t10 = dlog(t9)
      t11 = f2(0.D0,1.D0)
      t12 = dlog(t11)
      t13 = f3(0.D0,1.D0)
      t17 = dlog(1-t13+t13*l1*l2)
      t18 = dlog(l2)
      t20 = dlog(l4)
      t22 = 1/l2
      t24 = Dsquar311(1.D0,l1,l2,0.D0,l4)
      t25 = 1/t24
      t27 = s2vr2(1.D0,t9,t6,0.D0,l4)
      t29 = s23r2(1.D0,t9,t6,0.D0,l4)
      t32 = s1vr2(1.D0,t9,t6,0.D0,l4)
      t34 = f1(0.D0,1.D0,t9)
      t36 = 1/t27/t29/t32*t34
      tncoff3 = (-t2-t3-2*t4-t7-t10+t12-t17-2*t18-t20)*t22*t25*t36/4096-
     #t22/(-1+z)*t25*t36/1024+t22/l3*t25*t36/2048
      t6 = l3/2
      t7 = Dsquar311(1.D0,l1,l2,t6,l4)
      t12 = 1-l1*l2
      t13 = 1-l2
      t14 = s2vr2(1.D0,t12,t13,t6,l4)
      t16 = s23r2(1.D0,t12,t13,t6,l4)
      t19 = s1vr2(1.D0,t12,t13,t6,l4)
      t21 = f1(t6,1.D0,t12)
      tncoff4 = 1/(-2+l3)/l2/l3/t7/t14/t16/t19*t21/1024
      t1 = dlog(l4)
      t2 = 1-z
      t3 = dlog(t2)
      t6 = dlog(1-l4)
      t7 = dlog(l1)
      t8 = dlog(2.D0)
      t10 = f2(0.D0,z)
      t11 = dlog(t10)
      t12 = dlog(z)
      t14 = f3(0.D0,z)
      t16 = dlog(1-t14)
      t19 = -1/t2
      t21 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t22 = 1/t21
      t24 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t26 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t29 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t31 = f1(0.D0,z,1.D0)
      t33 = 1/t24/t26/t29*t31
      tncoff5 = -(-t1-4*t3-t6-t7-2*t8+t11+2*t12-t16)*t19*t22*t33/2048-1/
     #l3*t19*t22*t33/1024-1/l2*t19*t22*t33/1024
      t7 = l3/2
      t8 = Dsquar311(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,1.D0,1.D0,t7,l4)
      t14 = s23r2(z,1.D0,1.D0,t7,l4)
      t17 = s1vr2(z,1.D0,1.D0,t7,l4)
      t19 = f1(t7,z,1.D0)
      tncoff6 = -1/(-2+l3)/l3/(-1+z)/t8/t12/t14/t17*t19/512
      t5 = Dsquar311(z,l1,l2,0.D0,l4)
      t9 = 1-l1*l2
      t10 = 1-l2
      t11 = s2vr2(z,t9,t10,0.D0,l4)
      t13 = s23r2(z,t9,t10,0.D0,l4)
      t16 = s1vr2(z,t9,t10,0.D0,l4)
      t18 = f1(0.D0,z,t9)
      tncoff7 = 1/l2/(-1+z)/t5/t11/t13/t16*t18/1024
      t3 = Dsquar311(1.D0,l1,0.D0,0.D0,l4)
      t4 = 1/t3
      t6 = s2vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t7 = 1/t6
      t9 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t11 = s1vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t14 = f1(0.D0,1.D0,1.D0)
      t15 = 1/t9/t11*t14
      t18 = f2(0.D0,1.D0)
      t19 = dlog(t18)
      t21 = dlog(2.D0)
      t23 = dlog(l1)
      t25 = f3(0.D0,1.D0)
      t27 = dlog(1-t25)
      t29 = dlog(l4)
      t32 = dlog(1-l4)
      tncoff8 = 1/(-1+z)*t4*t7*t15/2048+(-t19/8192+t21/4096+t23/8192+t27
     #/8192+t29/8192+t32/8192)*t4*t7*t15-1/l3*t4*t7*t15/4096-1/l2*t4*t7*
     #t15/4096
      t5 = l3/2
      t6 = Dsquar311(1.D0,l1,0.D0,t5,l4)
      t9 = s2vr2(1.D0,1.D0,1.D0,t5,l4)
      t11 = s23r2(1.D0,1.D0,1.D0,t5,l4)
      t14 = s1vr2(1.D0,1.D0,1.D0,t5,l4)
      t16 = f1(t5,1.D0,1.D0)
      tncoff9 = -1/(-2+l3)/l3/t6/t9/t11/t14*t16/2048
      t2 = Dsquar311(1.D0,l1,l2,0.D0,l4)
      t6 = 1-l1*l2
      t7 = 1-l2
      t8 = s2vr2(1.D0,t6,t7,0.D0,l4)
      t11 = s23r2(1.D0,t6,t7,0.D0,l4)
      t13 = s1vr2(1.D0,t6,t7,0.D0,l4)
      t16 = f1(0.D0,1.D0,t6)
      tncoff10 = 1/l2/t2/t8/t11/t13*t16/4096
      t3 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t6 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t9 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t11 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t14 = f1(0.D0,z,1.D0)
      tncoff11 = -1/(-1+z)/t3/t6/t9/t11*t14/2048
      t1 = Dsquar311(1.D0,l1,0.D0,0.D0,l4)
      t3 = s2vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t8 = s1vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t11 = f1(0.D0,1.D0,1.D0)
      tncoff12 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num5e0(var,ndim,1.D0,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num5e0(var,ndim,1.D0,ux,1-l1*l2,1-l2,l3/2,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num5e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num5e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num5e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num5e1(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num5e1(var,ndim,1.D0,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = sq3num5e1(var,ndim,1.D0,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = sq3num5e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = sq3num5e2(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      sq39e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq39e2int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision sq3num5e0
      external sq3num5e0
      double precision sq3num5e1
      external sq3num5e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t3 = Dsquar311(1.D0,l1,0.D0,0.D0,l4)
      t4 = 1/t3
      t6 = s2vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t7 = 1/t6
      t9 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t11 = s1vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t14 = f1(0.D0,1.D0,1.D0)
      t15 = 1/t9/t11*t14
      t18 = f2(0.D0,1.D0)
      t19 = dlog(t18)
      t21 = dlog(2.D0)
      t23 = dlog(l1)
      t25 = f3(0.D0,1.D0)
      t27 = dlog(1-t25)
      t29 = dlog(l4)
      t32 = dlog(1-l4)
      tncoff1 = 1/(-1+z)*t4*t7*t15/2048+(-t19/8192+t21/4096+t23/8192+t27
     #/8192+t29/8192+t32/8192)*t4*t7*t15-1/l3*t4*t7*t15/4096-1/l2*t4*t7*
     #t15/4096
      t5 = l3/2
      t6 = Dsquar311(1.D0,l1,0.D0,t5,l4)
      t9 = s2vr2(1.D0,1.D0,1.D0,t5,l4)
      t11 = s23r2(1.D0,1.D0,1.D0,t5,l4)
      t14 = s1vr2(1.D0,1.D0,1.D0,t5,l4)
      t16 = f1(t5,1.D0,1.D0)
      tncoff2 = -1/(-2+l3)/l3/t6/t9/t11/t14*t16/2048
      t2 = Dsquar311(1.D0,l1,l2,0.D0,l4)
      t6 = 1-l1*l2
      t7 = 1-l2
      t8 = s2vr2(1.D0,t6,t7,0.D0,l4)
      t11 = s23r2(1.D0,t6,t7,0.D0,l4)
      t13 = s1vr2(1.D0,t6,t7,0.D0,l4)
      t16 = f1(0.D0,1.D0,t6)
      tncoff3 = 1/l2/t2/t8/t11/t13*t16/4096
      t3 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t6 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t9 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t11 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t14 = f1(0.D0,z,1.D0)
      tncoff4 = -1/(-1+z)/t3/t6/t9/t11*t14/2048
      t1 = Dsquar311(1.D0,l1,0.D0,0.D0,l4)
      t3 = s2vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t8 = s1vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t11 = f1(0.D0,1.D0,1.D0)
      tncoff5 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num5e0(var,ndim,1.D0,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num5e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num5e1(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      sq39e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq39e3int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision sq3num5e0
      external sq3num5e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = Dsquar311(1.D0,l1,0.D0,0.D0,l4)
      t3 = s2vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t8 = s1vr2(1.D0,1.D0,1.D0,0.D0,l4)
      t11 = f1(0.D0,1.D0,1.D0)
      tncoff1 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      sq39e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq39e4int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients


CCC Numerator functions


CCC Sum

      sq39e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

