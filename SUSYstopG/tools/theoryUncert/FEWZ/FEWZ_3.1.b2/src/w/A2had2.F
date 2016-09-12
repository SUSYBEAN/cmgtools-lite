CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function A22e0int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
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
      double precision A2num1e0
      external A2num1e0
      double precision A2num1e1
      external A2num1e1
      double precision A2num1e2
      external A2num1e2
      double precision A2num1e3
      external A2num1e3
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l5 = var(7)
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = var(9)
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = dlog(l4)
      t2 = dlog(l1)
      t3 = f2(1.D0,1.D0)
      t4 = dlog(t3)
      t5 = dlog(2.D0)
      t6 = 2*t5
      t8 = dlog(1-l2)
      t9 = dlog(l2)
      t10 = dlog(l3)
      t11 = 2*t10
      t13 = dlog(1-l4)
      t16 = 1/l3
      t17 = 1/l1
      t18 = t16*t17
      t20 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t21 = 1/t20
      t22 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t23 = t22**2
      t24 = 1/t23
      t26 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t27 = t26**2
      t28 = 1/t27
      t29 = f1(1.D0,1.D0,0.D0)
      t31 = t21*t24*t28*t29
      t34 = t5**2
      t35 = 12*t34
      t36 = 0.3141592653589793D1**2
      t37 = 2*t36
      t38 = 1-z
      t39 = dlog(t38)
      t47 = 12*t4*t5
      t54 = t4**2
      t55 = 3*t54
      t57 = 6*t9*t4
      t59 = 6*t8*t4
      t60 = t8*t5
      t61 = 12*t60
      t62 = t9*t8
      t63 = 6*t62
      t64 = t35-t37+24*t13*t39+24*t39*t9-24*t39*t4-t47+48*t39*t5+24*t39*
     #t8+24*t1*t39+t55-t57-t59+t61+t63
      t65 = t1*t8
      t66 = 6*t65
      t67 = t1*t5
      t68 = 12*t67
      t69 = t1*t9
      t70 = 6*t69
      t72 = 6*t1*t4
      t74 = 6*t13*t4
      t75 = t9*t5
      t76 = 12*t75
      t77 = t13*t5
      t78 = 12*t77
      t79 = t13*t1
      t80 = 6*t79
      t81 = t13*t8
      t82 = 6*t81
      t83 = t13*t9
      t84 = 6*t83
      t85 = t8**2
      t86 = 3*t85
      t87 = t1**2
      t88 = 3*t87
      t89 = t13**2
      t90 = 3*t89
      t91 = t9**2
      t92 = 3*t91
      t93 = t39**2
      t95 = t66+t68+t70-t72-t74+t76+t78+t80+t82+t84+t86+t88+t90+t92+48*t
     #93
      t99 = -1/t38
      t103 = 4*t39
      t139 = 16*z3+8*t34*t5+3*t89*t8-3*t91*t4+6*t91*t5+3*t9*t54+3*t89*t9
     #-4*t36*t5+3*t89*t1-3*t89*t4+3*t8*t54+12*t8*t34+6*t89*t5-3*t85*t4+6
     #*t85*t5
      t172 = -12*t60*t4-12*t67*t4+12*t69*t5-6*t62*t4+12*t65*t5-6*t65*t4+
     #6*t69*t8-6*t69*t4-2*t8*t36-12*t4*t34+2*t4*t36+3*t9*t85+12*t9*t34+3
     #*t1*t91+3*t1*t85+12*t1*t34
      t206 = -3*t87*t4+6*t87*t5-2*t9*t36-2*t1*t36+3*t1*t54+3*t87*t9+3*t8
     #7*t8+3*t13*t91+3*t13*t85+12*t13*t34-2*t13*t36+3*t13*t54+12*t79*t5-
     #12*t77*t4-6*t83*t4+12*t83*t5
      t234 = 6*t79*t9+3*t13*t87+6*t5*t54+3*t91*t8+12*t81*t5-6*t81*t4+t85
     #*t8+t91*t9+6*t79*t8-6*t79*t4+12*t62*t5-12*t75*t4-t54*t4+t89*t13+t8
     #7*t1+6*t83*t8
      t243 = t35-t37-t47+t55-t57-t59+t61+t63+t66+t68+t70-t72-t74+t76
      t256 = t2**2
      t258 = t78+t80+t82+t84+t86+t88+t90+t92+12*t2*t5+6*t13*t2+6*t2*t9-6
     #*t2*t4+6*t2*t8+6*t1*t2+3*t256
      t282 = t10**2
      t288 = t78+t80+t82+t84+t86+t88+t90+t92+24*t10*t5+12*t10*t9-12*t10*
     #t4+12*t10*t8+12*t282+12*t13*t10+12*t1*t10
      tncoff1 = (t1+t2-t4+t6+t8+t9+t11+t13)*NF*t18*t31/1024+NF*(t64+t95)
     #*t99*t31/3072-NF*(t103+t13+t9-t4+t6+t8+t2+t1)*t17*t99*t31/512+NF*(
     #t139+t172+t206+t234)*t21*t24*t28*t29/12288-NF*(t243+t258)*t17*t31/
     #12288+NF*t99*t18*t31/256-(t1+t103+t11+t8+t6+t9-t4+t13)*NF*t16*t99*
     #t31/256-(t243+t288)*NF*t16*t31/6144
      t1 = dlog(l4)
      t2 = 1-z
      t3 = dlog(t2)
      t5 = dlog(l3)
      t6 = 2*t5
      t7 = dlog(l2)
      t8 = -l3+2
      t9 = dlog(t8)
      t10 = 2*t9
      t12 = dlog(1-l2)
      t14 = 1-l3/2
      t15 = f2(t14,1.D0)
      t16 = dlog(t15)
      t18 = dlog(1-l4)
      t22 = -1/t8
      t23 = 1/l3
      t24 = t22*t23
      t27 = -1/t2
      t28 = s34r2(1.D0,0.D0,l2,t14,l4)
      t29 = 1/t28
      t31 = s2vr2(1.D0,0.D0,l2,t14,l4)
      t32 = t31**2
      t33 = 1/t32
      t34 = s1vr2(1.D0,0.D0,l2,t14,l4)
      t35 = t34**2
      t36 = 1/t35
      t38 = f1(t14,1.D0,0.D0)
      t39 = t33*t36*t38
      t43 = t18**2
      t45 = t1**2
      t47 = 0.3141592653589793D1**2
      t49 = t16**2
      t51 = t5**2
      t71 = 3*t43+3*t45-2*t47+3*t49+12*t51+6*t7*t12+6*t18*t1+6*t18*t12+6
     #*t18*t7+6*t1*t12+6*t1*t7+12*t5*t12+12*t18*t9-6*t1*t16
      t96 = t9**2
      t98 = t12**2
      t100 = t7**2
      t102 = 12*t1*t9-6*t18*t16-12*t5*t16+24*t9*t5+12*t12*t9-12*t16*t9-6
     #*t12*t16-6*t7*t16+12*t7*t9+12*t18*t5+12*t1*t5+12*t5*t7+12*t96+3*t9
     #8+3*t100
      t106 = t29*t33
      t113 = 1/l1
      t120 = dlog(l1)
      tncoff2 = -NF*(t1+4*t3+t6+t7+t10+t12-t16+t18)*t24*t27*t29*t39/128-
     #NF*(t71+t102)*t24*t106*t36*t38/3072+NF*t22*t106*t36*t113*t23*t27*t
     #38/128+NF*(t1+t7+t120+t10+t12-t16+t6+t18)*t22*t113*t23*t29*t39/512
      t1 = 0.3141592653589793D1**2
      t3 = f3(1.D0,1.D0)
      t6 = dlog(1-l1*t3)
      t7 = t6**2
      t9 = dlog(2.D0)
      t10 = t9**2
      t13 = dlog(1-l1)
      t14 = t13**2
      t17 = dlog(1-l2)
      t18 = t17**2
      t20 = dlog(l2)
      t21 = t20**2
      t25 = f2(1.D0,1.D0)
      t26 = dlog(t25)
      t27 = t26**2
      t31 = dlog(l1)
      t35 = dlog(1-l4)
      t36 = t35**2
      t38 = -2*t1+3*t7+12*t10+3*t14+3*t18+3*t21+12*t6*t9+3*t27-6*t26*t13
     #+6*t31*t13+3*t36
      t39 = dlog(l4)
      t40 = t39**2
      t42 = t31**2
      t64 = 3*t40+3*t42+6*t6*t13+6*t20*t13+6*t20*t6+6*t13*t39+6*t39*t6+1
     #2*t13*t9+6*t17*t6+6*t17*t13-6*t26*t6-12*t26*t9
      t88 = -6*t20*t26-6*t17*t26+12*t17*t9+6*t20*t17+12*t20*t9+12*t35*t9
     #+6*t35*t39+6*t35*t17+6*t35*t20+6*t39*t17+12*t39*t9
      t113 = 6*t39*t20-6*t39*t26-6*t35*t26+12*t31*t9+6*t35*t31+6*t31*t20
     #-6*t31*t26+6*t31*t17+6*t39*t31+6*t31*t6+6*t35*t13+6*t35*t6
      t117 = 1/l1
      t119 = s34r2(1.D0,l1,l2,1.D0,l4)
      t121 = s2vr2(1.D0,l1,l2,1.D0,l4)
      t122 = t121**2
      t125 = s1vr2(1.D0,l1,l2,1.D0,l4)
      t126 = t125**2
      t128 = f1(1.D0,1.D0,l1)
      t130 = 1/t119/t122/t126*t128
      t133 = -1+z
      t134 = 1/t133
      t137 = 1/l3*t117
      t141 = 2*t9
      t142 = dlog(l3)
      t150 = dlog(-t133)
      tncoff3 = NF*(t38+t64+t88+t113)*t117*t130/12288-NF*t134*t137*t130/
     #256-(t39+t20+t31+t13-t26+t17+t6+t141+2*t142+t35)*NF*t137*t130/1024
     #+NF*(4*t150+t35+t141-t26+t17+t13+t20+t6+t31+t39)*t117*t134*t130/51
     #2
      t1 = dlog(l4)
      t2 = dlog(l1)
      t4 = 1-l3/2
      t5 = f2(t4,1.D0)
      t6 = dlog(t5)
      t7 = f3(t4,1.D0)
      t10 = dlog(1-l1*t7)
      t12 = dlog(1-l2)
      t13 = dlog(l2)
      t15 = dlog(1-l1)
      t16 = -l3+2
      t17 = dlog(t16)
      t19 = dlog(l3)
      t22 = dlog(1-l4)
      t26 = -1/t16
      t27 = 1/l1
      t30 = 1/l3
      t31 = s34r2(1.D0,l1,l2,t4,l4)
      t32 = 1/t31
      t34 = s2vr2(1.D0,l1,l2,t4,l4)
      t35 = t34**2
      t36 = 1/t35
      t37 = s1vr2(1.D0,l1,l2,t4,l4)
      t38 = t37**2
      t39 = 1/t38
      t41 = f1(t4,1.D0,l1)
      tncoff4 = -NF*(t1+t2-t6+t10+t12+t13+t15+2*t17+2*t19+t22)*t26*t27*t
     #30*t32*t36*t39*t41/512-NF*t26*t27*t30*t36*t32*t39/(-1+z)*t41/128
      t1 = 1-z
      t2 = dlog(t1)
      t3 = 4*t2
      t5 = dlog(1-l4)
      t7 = dlog(1-l2)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t10 = dlog(l2)
      t11 = dlog(2.D0)
      t12 = 2*t11
      t13 = dlog(z)
      t14 = 2*t13
      t15 = dlog(l1)
      t16 = dlog(l4)
      t19 = 1/l1
      t21 = -1/t1
      t24 = s34r2(z,0.D0,l2,1.D0,l4)
      t26 = s2vr2(z,0.D0,l2,1.D0,l4)
      t27 = t26**2
      t30 = s1vr2(z,0.D0,l2,1.D0,l4)
      t31 = t30**2
      t33 = f1(1.D0,z,0.D0)
      t35 = 1/t24/t27/t31*t33
      t38 = dlog(l3)
      t42 = 1/l3
      t52 = t9**2
      t66 = t7**2
      t68 = t2**2
      t70 = t10**2
      t80 = t5**2
      t88 = 3*t52+6*t16*t10-48*t2*t13+24*t2*t10+24*t2*t7+48*t2*t11+6*t5*
     #t7+3*t66+48*t68+3*t70+12*t5*t11+6*t5*t16+6*t5*t10+6*t16*t7+3*t80+1
     #2*t16*t11-6*t10*t9-24*t2*t9
      t107 = t16**2
      t109 = 0.3141592653589793D1**2
      t111 = t11**2
      t113 = t13**2
      t127 = -6*t16*t9-12*t16*t13-12*t10*t13-12*t9*t11+24*t5*t2+24*t16*t
     #2+12*t7*t11+6*t10*t7+12*t10*t11+3*t107-2*t109+12*t111+12*t113-24*t
     #13*t11+12*t9*t13-12*t5*t13-6*t5*t9-6*t7*t9-12*t7*t13
      tncoff5 = -NF*(-t3-t5-t7+t9-t10-t12+t14-t15-t16)*t19*t21*t35/512-(
     #-t16-t3+t9-2*t38-t12-t10+t14-t7-t5)*NF*t42*t21*t35/256-NF*t21*t42*
     #t19*t35/256-NF*(t88+t127)*t21*t35/3072
      t1 = l3-2
      t2 = 1/t1
      t5 = 1-l3/2
      t6 = s2vr2(z,0.D0,l2,t5,l4)
      t7 = t6**2
      t8 = 1/t7
      t9 = s34r2(z,0.D0,l2,t5,l4)
      t10 = 1/t9
      t13 = s1vr2(z,0.D0,l2,t5,l4)
      t14 = t13**2
      t15 = 1/t14
      t18 = 1/l3
      t19 = -1+z
      t20 = 1/t19
      t22 = f1(t5,z,0.D0)
      t26 = dlog(l4)
      t28 = dlog(-t19)
      t30 = dlog(l2)
      t31 = dlog(l3)
      t34 = dlog(1-l2)
      t35 = f2(t5,z)
      t36 = dlog(t35)
      t38 = dlog(-t1)
      t40 = dlog(z)
      t43 = dlog(1-l4)
      tncoff6 = -NF*t2*t8*t10*t15/l1*t18*t20*t22/128-NF*(-t26-4*t28-t30-
     #2*t31-t34+t36-2*t38+2*t40-t43)*t2*t18*t20*t10*t8*t15*t22/128
      t1 = 1-z
      t2 = dlog(t1)
      t5 = dlog(1-l4)
      t7 = dlog(1-l2)
      t8 = f3(1.D0,z)
      t11 = dlog(1-l1*t8)
      t13 = dlog(1-l1)
      t14 = dlog(l1)
      t15 = dlog(z)
      t17 = dlog(2.D0)
      t19 = dlog(l2)
      t20 = f2(1.D0,z)
      t21 = dlog(t20)
      t22 = dlog(l4)
      t23 = -4*t2-t5-t7-t11-t13-t14+2*t15-2*t17-t19+t21-t22
      t25 = 1/l1
      t27 = -1/t1
      t30 = s34r2(z,l1,l2,1.D0,l4)
      t32 = s2vr2(z,l1,l2,1.D0,l4)
      t33 = t32**2
      t36 = s1vr2(z,l1,l2,1.D0,l4)
      t37 = t36**2
      t39 = f1(1.D0,z,l1)
      t41 = 1/t30/t33/t37*t39
      tncoff7 = NF*t23*t25*t27*t41/512+NF*t27/l3*t25*t41/256
      t2 = 1-l3/2
      t3 = f1(t2,z,l1)
      t5 = s2vr2(z,l1,l2,t2,l4)
      t6 = t5**2
      t8 = s34r2(z,l1,l2,t2,l4)
      t12 = s1vr2(z,l1,l2,t2,l4)
      t13 = t12**2
      tncoff8 = NF*t3/t6/t8/t13/(l3-2)/l1/l3/(-1+z)/128
      t1 = 1/l1
      t3 = -1+z
      t4 = 1/t3
      t6 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t7 = 1/t6
      t8 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t9 = t8**2
      t10 = 1/t9
      t12 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t13 = t12**2
      t14 = 1/t13
      t15 = f1(1.D0,1.D0,0.D0)
      t17 = t7*t10*t14*t15
      t21 = dlog(1-l4)
      t22 = t21**2
      t24 = dlog(l4)
      t25 = t24**2
      t27 = f2(1.D0,1.D0)
      t28 = dlog(t27)
      t29 = t28**2
      t31 = 0.3141592653589793D1**2
      t33 = dlog(2.D0)
      t34 = t33**2
      t37 = dlog(1-l2)
      t38 = t37**2
      t40 = dlog(l2)
      t41 = t40**2
      t51 = -3*t22-3*t25-3*t29+2*t31-12*t34-3*t38-3*t41-12*t37*t33+6*t37
     #*t28-12*t40*t33-6*t24*t40
      t74 = -6*t24*t37+6*t24*t28-12*t24*t33-6*t21*t40-6*t21*t37+6*t21*t2
     #8-12*t21*t33-6*t21*t24-6*t40*t37+6*t40*t28+12*t28*t33
      t82 = 1/l3
      t83 = NF*t82
      t87 = 6*t24
      t88 = 6*t21
      t89 = 6*t40
      t90 = 6*t37
      t91 = 12*t33
      t92 = 6*t28
      t94 = dlog(-t3)
      t104 = dlog(l1)
      t111 = dlog(l3)
      tncoff9 = NF*t1*t4*t17/512+NF*(t51+t74)*t7*t10*t14*t15/12288-t83*t
     #1*t17/1024+NF*(-t87-t88-t89-t90-t91+t92-24*t94)*t4*t17/3072+t83*t4
     #*t17/256-NF*(-t87-t88-t89-t90-6*t104-t91+t92)*t1*t17/12288-(-t87-t
     #88-t89-t90-t91+t92-12*t111)*NF*t82*t17/6144
      t1 = l3-2
      t2 = 1/t1
      t3 = NF*t2
      t4 = 1/l3
      t10 = 1-l3/2
      t11 = s34r2(1.D0,0.D0,l2,t10,l4)
      t13 = s2vr2(1.D0,0.D0,l2,t10,l4)
      t14 = t13**2
      t17 = s1vr2(1.D0,0.D0,l2,t10,l4)
      t18 = t17**2
      t20 = f1(t10,1.D0,0.D0)
      t22 = 1/t11/t14/t18*t20
      t26 = dlog(1-l4)
      t28 = dlog(l2)
      t31 = dlog(-t1)
      t33 = f2(t10,1.D0)
      t34 = dlog(t33)
      t37 = dlog(1-l2)
      t39 = dlog(l3)
      t41 = dlog(l4)
      tncoff10 = t3*t4/(-1+z)*t22/128-NF*(-6*t26-6*t28-12*t31+6*t34-6*t3
     #7-12*t39-6*t41)*t2*t4*t22/3072-t3/l1*t4*t22/512
      t1 = 1/l1
      t6 = s34r2(1.D0,l1,l2,1.D0,l4)
      t8 = s2vr2(1.D0,l1,l2,1.D0,l4)
      t9 = t8**2
      t12 = s1vr2(1.D0,l1,l2,1.D0,l4)
      t13 = t12**2
      t15 = f1(1.D0,1.D0,l1)
      t17 = 1/t6/t9/t13*t15
      t20 = dlog(l4)
      t23 = dlog(1-l4)
      t25 = dlog(l1)
      t27 = dlog(l2)
      t30 = dlog(1-l2)
      t32 = f3(1.D0,1.D0)
      t35 = dlog(1-l1*t32)
      t37 = f2(1.D0,1.D0)
      t38 = dlog(t37)
      t41 = dlog(1-l1)
      t43 = dlog(2.D0)
      tncoff11 = -NF*t1/(-1+z)*t17/512+NF*(-6*t20-6*t23-6*t25-6*t27-6*t3
     #0-6*t35+6*t38-6*t41-12*t43)*t1*t17/12288+NF/l3*t1*t17/1024
      t9 = 1-l3/2
      t10 = s34r2(1.D0,l1,l2,t9,l4)
      t12 = s2vr2(1.D0,l1,l2,t9,l4)
      t13 = t12**2
      t16 = s1vr2(1.D0,l1,l2,t9,l4)
      t17 = t16**2
      t19 = f1(t9,1.D0,l1)
      tncoff12 = NF/(l3-2)/l1/l3/t10/t13/t17*t19/512
      t1 = dlog(l4)
      t3 = dlog(l2)
      t6 = dlog(1-l2)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t11 = dlog(z)
      t13 = dlog(2.D0)
      t15 = 1-z
      t16 = dlog(t15)
      t19 = dlog(1-l4)
      t24 = -1/t15
      t26 = s34r2(z,0.D0,l2,1.D0,l4)
      t28 = s2vr2(z,0.D0,l2,1.D0,l4)
      t29 = t28**2
      t32 = s1vr2(z,0.D0,l2,1.D0,l4)
      t33 = t32**2
      t35 = f1(1.D0,z,0.D0)
      t37 = 1/t26/t29/t33*t35
      tncoff13 = -NF*(-6*t1-6*t3-6*t6+6*t9+12*t11-12*t13-24*t16-6*t19)*t
     #24*t37/3072-NF/l3*t24*t37/256-NF/l1*t24*t37/512
      t10 = 1-l3/2
      t11 = s34r2(z,0.D0,l2,t10,l4)
      t13 = s2vr2(z,0.D0,l2,t10,l4)
      t14 = t13**2
      t17 = s1vr2(z,0.D0,l2,t10,l4)
      t18 = t17**2
      t20 = f1(t10,z,0.D0)
      tncoff14 = -NF/(l3-2)/l3/(-1+z)/t11/t14/t18*t20/128
      t6 = s34r2(z,l1,l2,1.D0,l4)
      t8 = s2vr2(z,l1,l2,1.D0,l4)
      t9 = t8**2
      t12 = s1vr2(z,l1,l2,1.D0,l4)
      t13 = t12**2
      t15 = f1(1.D0,z,l1)
      tncoff15 = NF/l1/(-1+z)/t6/t9/t13*t15/512
      t1 = dlog(2.D0)
      t4 = dlog(1-l4)
      t6 = dlog(l4)
      t9 = dlog(1-l2)
      t11 = dlog(l2)
      t13 = f2(1.D0,1.D0)
      t14 = dlog(t13)
      t18 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t19 = 1/t18
      t21 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t22 = t21**2
      t24 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t25 = t24**2
      t28 = f1(1.D0,1.D0,0.D0)
      t29 = 1/t22/t25*t28
      tncoff16 = NF*(12*t1+6*t4+6*t6+6*t9+6*t11-6*t14)*t19*t29/12288+NF/
     #(-1+z)*t19*t29/512-NF/l3*t19*t29/1024-NF/l1*t19*t29/2048
      t7 = 1-l3/2
      t8 = s34r2(1.D0,0.D0,l2,t7,l4)
      t10 = s2vr2(1.D0,0.D0,l2,t7,l4)
      t11 = t10**2
      t14 = s1vr2(1.D0,0.D0,l2,t7,l4)
      t15 = t14**2
      t17 = f1(t7,1.D0,0.D0)
      tncoff17 = -NF/(l3-2)/l3/t8/t11/t15*t17/512
      t3 = s34r2(1.D0,l1,l2,1.D0,l4)
      t6 = s2vr2(1.D0,l1,l2,1.D0,l4)
      t7 = t6**2
      t9 = s1vr2(1.D0,l1,l2,1.D0,l4)
      t10 = t9**2
      t13 = f1(1.D0,1.D0,l1)
      tncoff18 = NF/l1/t3/t7/t10*t13/2048
      t4 = s34r2(z,0.D0,l2,1.D0,l4)
      t7 = s2vr2(z,0.D0,l2,1.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,0.D0,l2,1.D0,l4)
      t11 = t10**2
      t14 = f1(1.D0,z,0.D0)
      tncoff19 = -NF/(-1+z)/t4/t8/t11*t14/512
      t1 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t4 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t5 = t4**2
      t7 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t8 = t7**2
      t11 = f1(1.D0,1.D0,0.D0)
      tncoff20 = -NF/t1/t5/t8*t11/2048

CCC Numerator functions

      tcnum1 = A2num1e0(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = A2num1e0(var,ndim,1.D0,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = A2num1e0(var,ndim,1.D0,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = A2num1e0(var,ndim,1.D0,ux,l1,l2,1-l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = A2num1e0(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = A2num1e0(var,ndim,z,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = A2num1e0(var,ndim,z,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = A2num1e0(var,ndim,z,ux,l1,l2,1-l3/2,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = A2num1e1(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = A2num1e1(var,ndim,1.D0,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = A2num1e1(var,ndim,1.D0,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = A2num1e1(var,ndim,1.D0,ux,l1,l2,1-l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = A2num1e1(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = A2num1e1(var,ndim,z,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = A2num1e1(var,ndim,z,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = A2num1e2(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = A2num1e2(var,ndim,1.D0,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = A2num1e2(var,ndim,1.D0,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = A2num1e2(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = A2num1e3(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      A22e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16
     #+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncoff
     #20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function A22e1int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
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
      double precision A2num1e0
      external A2num1e0
      double precision A2num1e1
      external A2num1e1
      double precision A2num1e2
      external A2num1e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l5 = var(7)
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = var(9)
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
      t4 = f2(1.D0,1.D0)
      t5 = dlog(t4)
      t6 = t5**2
      t8 = 0.3141592653589793D1**2
      t11 = dlog(1-l2)
      t12 = dlog(2.D0)
      t17 = dlog(l2)
      t29 = dlog(1-l4)
      t32 = 3*t2+3*t6-2*t8+12*t11*t12-12*t5*t12+6*t1*t17+12*t17*t12-6*t1
     #1*t5+12*t1*t12-6*t1*t5+6*t29*t1
      t43 = t11**2
      t45 = t17**2
      t47 = t29**2
      t53 = t12**2
      t55 = 12*t29*t12+6*t1*t11+6*t29*t17+6*t29*t11-6*t29*t5+3*t43+3*t45
     #+3*t47+6*t17*t11-6*t17*t5+12*t53
      t58 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t59 = 1/t58
      t61 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t62 = t61**2
      t63 = 1/t62
      t64 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t65 = t64**2
      t66 = 1/t65
      t68 = f1(1.D0,1.D0,0.D0)
      t72 = 2*t12
      t73 = 1-z
      t74 = dlog(t73)
      t79 = -1/t73
      t83 = t59*t63*t66*t68
      t86 = dlog(l3)
      t90 = 1/l3
      t94 = NF*t90
      t98 = dlog(l1)
      t101 = 1/l1
      tncoff1 = -NF*(t32+t55)*t59*t63*t66*t68/12288-NF*(t17-t5+t72+t11+t
     #29+4*t74+t1)*t79*t83/512+(t29+2*t86+t11+t72+t17-t5+t1)*NF*t90*t83/
     #1024+t94*t79*t83/256+NF*(t98+t72-t5+t11+t17+t1+t29)*t101*t83/2048+
     #NF*t101*t79*t83/512-t94*t101*t83/1024
      t2 = dlog(1-l4)
      t3 = dlog(l3)
      t5 = dlog(l2)
      t6 = -l3+2
      t7 = dlog(t6)
      t10 = dlog(1-l2)
      t12 = 1-l3/2
      t13 = f2(t12,1.D0)
      t14 = dlog(t13)
      t15 = dlog(l4)
      t19 = -1/t6
      t20 = 1/l3
      t23 = s34r2(1.D0,0.D0,l2,t12,l4)
      t25 = s2vr2(1.D0,0.D0,l2,t12,l4)
      t26 = t25**2
      t29 = s1vr2(1.D0,0.D0,l2,t12,l4)
      t30 = t29**2
      t32 = f1(t12,1.D0,0.D0)
      t34 = 1/t23/t26/t30*t32
      t37 = NF*t19
      tncoff2 = NF*(t2+2*t3+t5+2*t7+t10-t14+t15)*t19*t20*t34/512+t37*t20
     #/(-1+z)*t34/128-t37/l1*t20*t34/512
      t1 = dlog(l1)
      t3 = dlog(1-l2)
      t4 = dlog(2.D0)
      t6 = f2(1.D0,1.D0)
      t7 = dlog(t6)
      t8 = f3(1.D0,1.D0)
      t11 = dlog(1-l1*t8)
      t13 = dlog(1-l1)
      t14 = dlog(l2)
      t15 = dlog(l4)
      t17 = dlog(1-l4)
      t20 = 1/l1
      t22 = s34r2(1.D0,l1,l2,1.D0,l4)
      t24 = s2vr2(1.D0,l1,l2,1.D0,l4)
      t25 = t24**2
      t28 = s1vr2(1.D0,l1,l2,1.D0,l4)
      t29 = t28**2
      t31 = f1(1.D0,1.D0,l1)
      t33 = 1/t22/t25/t29*t31
      tncoff3 = -NF*(t1+t3+2*t4-t7+t11+t13+t14+t15+t17)*t20*t33/2048-NF*
     #t20/(-1+z)*t33/512+NF/l3*t20*t33/1024
      t9 = 1-l3/2
      t10 = s34r2(1.D0,l1,l2,t9,l4)
      t12 = s2vr2(1.D0,l1,l2,t9,l4)
      t13 = t12**2
      t16 = s1vr2(1.D0,l1,l2,t9,l4)
      t17 = t16**2
      t19 = f1(t9,1.D0,l1)
      tncoff4 = NF/(l3-2)/l1/l3/t10/t13/t17*t19/512
      t1 = dlog(2.D0)
      t4 = dlog(1-l2)
      t5 = dlog(l2)
      t6 = dlog(z)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t11 = dlog(1-l4)
      t12 = 1-z
      t13 = dlog(t12)
      t15 = dlog(l4)
      t19 = -1/t12
      t21 = s34r2(z,0.D0,l2,1.D0,l4)
      t23 = s2vr2(z,0.D0,l2,1.D0,l4)
      t24 = t23**2
      t27 = s1vr2(z,0.D0,l2,1.D0,l4)
      t28 = t27**2
      t30 = f1(1.D0,z,0.D0)
      t32 = 1/t21/t24/t28*t30
      tncoff5 = -NF*(-2*t1-t4-t5+2*t6+t9-t11-4*t13-t15)*t19*t32/512-NF/l
     #3*t19*t32/256-NF/l1*t19*t32/512
      t10 = 1-l3/2
      t11 = s34r2(z,0.D0,l2,t10,l4)
      t13 = s2vr2(z,0.D0,l2,t10,l4)
      t14 = t13**2
      t17 = s1vr2(z,0.D0,l2,t10,l4)
      t18 = t17**2
      t20 = f1(t10,z,0.D0)
      tncoff6 = -NF/(l3-2)/l3/(-1+z)/t11/t14/t18*t20/128
      t6 = s34r2(z,l1,l2,1.D0,l4)
      t8 = s2vr2(z,l1,l2,1.D0,l4)
      t9 = t8**2
      t12 = s1vr2(z,l1,l2,1.D0,l4)
      t13 = t12**2
      t15 = f1(1.D0,z,l1)
      tncoff7 = NF/l1/(-1+z)/t6/t9/t13*t15/512
      t1 = dlog(l4)
      t4 = dlog(1-l4)
      t6 = f2(1.D0,1.D0)
      t7 = dlog(t6)
      t9 = dlog(2.D0)
      t12 = dlog(1-l2)
      t14 = dlog(l2)
      t18 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t19 = 1/t18
      t21 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t22 = t21**2
      t24 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t25 = t24**2
      t28 = f1(1.D0,1.D0,0.D0)
      t29 = 1/t22/t25*t28
      tncoff8 = -NF*(-6*t1-6*t4+6*t7-12*t9-6*t12-6*t14)*t19*t29/12288+NF
     #/(-1+z)*t19*t29/512-NF/l3*t19*t29/1024-NF/l1*t19*t29/2048
      t7 = 1-l3/2
      t8 = s34r2(1.D0,0.D0,l2,t7,l4)
      t10 = s2vr2(1.D0,0.D0,l2,t7,l4)
      t11 = t10**2
      t14 = s1vr2(1.D0,0.D0,l2,t7,l4)
      t15 = t14**2
      t17 = f1(t7,1.D0,0.D0)
      tncoff9 = -NF/(l3-2)/l3/t8/t11/t15*t17/512
      t3 = s34r2(1.D0,l1,l2,1.D0,l4)
      t6 = s2vr2(1.D0,l1,l2,1.D0,l4)
      t7 = t6**2
      t9 = s1vr2(1.D0,l1,l2,1.D0,l4)
      t10 = t9**2
      t13 = f1(1.D0,1.D0,l1)
      tncoff10 = NF/l1/t3/t7/t10*t13/2048
      t4 = s34r2(z,0.D0,l2,1.D0,l4)
      t7 = s2vr2(z,0.D0,l2,1.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,0.D0,l2,1.D0,l4)
      t11 = t10**2
      t14 = f1(1.D0,z,0.D0)
      tncoff11 = -NF/(-1+z)/t4/t8/t11*t14/512
      t1 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t4 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t5 = t4**2
      t7 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t8 = t7**2
      t11 = f1(1.D0,1.D0,0.D0)
      tncoff12 = -NF/t1/t5/t8*t11/2048

CCC Numerator functions

      tcnum1 = A2num1e0(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = A2num1e0(var,ndim,1.D0,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = A2num1e0(var,ndim,1.D0,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = A2num1e0(var,ndim,1.D0,ux,l1,l2,1-l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = A2num1e0(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = A2num1e0(var,ndim,z,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = A2num1e0(var,ndim,z,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = A2num1e1(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = A2num1e1(var,ndim,1.D0,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = A2num1e1(var,ndim,1.D0,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = A2num1e1(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = A2num1e2(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      A22e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function A22e2int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
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
      double precision A2num1e0
      external A2num1e0
      double precision A2num1e1
      external A2num1e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l5 = var(7)
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = var(9)
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l2)
      t5 = f2(1.D0,1.D0)
      t6 = dlog(t5)
      t7 = dlog(2.D0)
      t10 = dlog(1-l2)
      t13 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t14 = 1/t13
      t16 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t17 = t16**2
      t19 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t20 = t19**2
      t23 = f1(1.D0,1.D0,0.D0)
      t24 = 1/t17/t20*t23
      tncoff1 = NF*(t2+t3+t4-t6+2*t7+t10)*t14*t24/2048+NF/(-1+z)*t14*t24
     #/512-NF/l3*t14*t24/1024-NF/l1*t14*t24/2048
      t7 = 1-l3/2
      t8 = s34r2(1.D0,0.D0,l2,t7,l4)
      t10 = s2vr2(1.D0,0.D0,l2,t7,l4)
      t11 = t10**2
      t14 = s1vr2(1.D0,0.D0,l2,t7,l4)
      t15 = t14**2
      t17 = f1(t7,1.D0,0.D0)
      tncoff2 = -NF/(l3-2)/l3/t8/t11/t15*t17/512
      t3 = s34r2(1.D0,l1,l2,1.D0,l4)
      t6 = s2vr2(1.D0,l1,l2,1.D0,l4)
      t7 = t6**2
      t9 = s1vr2(1.D0,l1,l2,1.D0,l4)
      t10 = t9**2
      t13 = f1(1.D0,1.D0,l1)
      tncoff3 = NF/l1/t3/t7/t10*t13/2048
      t4 = s34r2(z,0.D0,l2,1.D0,l4)
      t7 = s2vr2(z,0.D0,l2,1.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,0.D0,l2,1.D0,l4)
      t11 = t10**2
      t14 = f1(1.D0,z,0.D0)
      tncoff4 = -NF/(-1+z)/t4/t8/t11*t14/512
      t1 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t4 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t5 = t4**2
      t7 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t8 = t7**2
      t11 = f1(1.D0,1.D0,0.D0)
      tncoff5 = -NF/t1/t5/t8*t11/2048

CCC Numerator functions

      tcnum1 = A2num1e0(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = A2num1e0(var,ndim,1.D0,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = A2num1e0(var,ndim,1.D0,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = A2num1e0(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = A2num1e1(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      A22e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

