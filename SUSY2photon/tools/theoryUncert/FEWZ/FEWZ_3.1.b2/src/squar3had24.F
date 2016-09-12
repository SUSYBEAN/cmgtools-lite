CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function sq34e0int(var,ndim,Tcoeff)

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
      double precision sq3num2e0
      external sq3num2e0
      double precision sq3num2e1
      external sq3num2e1
      double precision sq3num2e2
      external sq3num2e2
      double precision sq3num2e3
      external sq3num2e3
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

      t1 = f3(1.D0,z)
      t3 = dlog(1-t1)
      t4 = t3**2
      t5 = dlog(z)
      t9 = dlog(1-z)
      t10 = t9**2
      t14 = dlog(1-l4)
      t15 = dlog(l4)
      t16 = t15**2
      t23 = t5**2
      t26 = 0.3141592653589793D1**2
      t29 = dlog(2.D0)
      t30 = t29**2
      t35 = f2(1.D0,z)
      t36 = dlog(t35)
      t37 = t36**2
      t42 = -t4*t5/2048+t3*t10/256+t14*t16/4096+t14*t10/256+t14*t4/4096+
     #t3*t23/1024-t3*t26/6144+t3*t30/256-t26*t29/1536+t14*t37/4096+t30*t
     #29/192
      t43 = t14**2
      t68 = t43*t15/4096+t43*t3/4096+t43*t9/1024+t16*t9/1024-t16*t36/409
     #6+t16*t29/1024-t16*t5/2048+t15*t4/4096+t15*t37/4096+t15*t23/1024+t
     #15*t10/256+t37*t29/1024
      t91 = -t15*t26/6144+t15*t30/256+t14*t23/1024-t14*t26/6144-t10*t36/
     #256-t43*t36/4096+t43*t29/1024-t43*t5/2048+t16*t3/4096+t14*t30/256+
     #z3/768
      t92 = t14*t15
      t95 = t14*t9
      t98 = t3*t9
      t119 = t92*t29/512-t95*t5/256-t98*t36/512-t92*t5/1024+t92*t9/512+t
     #16*t15/12288+t10*t29/64+t36*t26/6144-t36*t30/256-t37*t5/2048-t36*t
     #23/1024-t9*t26/1536
      t126 = t15*t5
      t129 = t3*t36
      t132 = t9*t5
      t135 = t9*t36
      t144 = t3*t5
      t147 = t15*t9
      t150 = t9*t30/64+t98*t29/128-t126*t29/256-t129*t29/512-t132*t29/64
     #-t135*t29/128+t10*t9/192+t135*t5/256+t129*t5/1024-t144*t29/256+t14
     #7*t29/128
      t175 = -t147*t5/256+t4*t9/1024-t4*t36/4096+t4*t29/1024-t10*t5/128+
     #t9*t37/1024+t9*t23/256+t26*t5/3072+t23*t29/256-t23*t5/1536-t37*t36
     #/12288+t3*t37/4096
      t181 = t15*t36
      t192 = t14*t36
      t199 = t15*t3
      t202 = t4*t3/12288-t92*t36/2048-t181*t29/512-t98*t5/256+t181*t5/10
     #24-t147*t36/512-t5*t30/128-t192*t29/512-t95*t36/512+t192*t5/1024+t
     #199*t29/512
      t205 = t14*t3
      t214 = t14*t5
      t225 = t36*t5
      t230 = -t199*t5/1024-t205*t36/2048+t205*t9/512+t199*t9/512-t199*t3
     #6/2048-t214*t29/256-t205*t5/1024+t205*t29/512+t95*t29/128+t92*t3/2
     #048+t225*t29/256+t43*t14/12288
      t234 = Dsquar312(z,0.D0,0.D0,1.D0,l4)
      t235 = 1/t234
      t237 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t238 = t237**2
      t239 = 1/t238
      t240 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t241 = 1/t240
      t243 = f1(1.D0,z,1.D0)
      t244 = t239*t241*t243
      t246 = dlog(l3)
      t253 = 3*t16
      t254 = t246**2
      t256 = 12*t126
      t257 = 24*t147
      t259 = 96*t9*t29
      t261 = 48*t5*t29
      t262 = 2*t26
      t263 = 24*t135
      t264 = 3*t43
      t266 = 24*t36*t29
      t267 = 3*t37
      t269 = 24*t14*t29
      t270 = 6*t92
      t271 = 24*t95
      t272 = 12*t144
      t273 = -24*t246*t5+12*t14*t246+12*t246*t3+t253+12*t254-t256+t257+t
     #259-t261-t262-t263+t264-t266+t267+t269+t270+t271-t272
      t274 = 3*t4
      t276 = 24*t3*t29
      t277 = 6*t129
      t278 = 12*t23
      t279 = 48*t10
      t280 = 24*t98
      t281 = 6*t192
      t282 = 6*t199
      t283 = 12*t225
      t284 = 48*t30
      t285 = 6*t181
      t286 = 6*t205
      t287 = 12*t214
      t289 = 24*t15*t29
      t290 = 48*t132
      t299 = t274+t276-t277+t278+t279+t280-t281+t282+t283+t284-t285+t286
     #-t287+t289-t290+48*t246*t29-12*t246*t36+48*t246*t9+12*t15*t246
      t301 = 1/l3
      t306 = t253-t256+t257+t259-t261-t262-t263+t264-t266+t267+t269+t270
     #+t271-t272+t274+t276-t277+t278
      t307 = dlog(l2)
      t322 = t307**2
      t324 = t279+t280-t281+t282+t283+t284-t285+24*t307*t9-12*t307*t5+6*
     #t14*t307+6*t15*t307+6*t307*t3+t286-t287+24*t307*t29-6*t307*t36+3*t
     #322+t289-t290
      t326 = 1/l2
      t327 = (t306+t324)*t326
      t328 = Dsquar312(z,0.D0,l2,1.D0,l4)
      t329 = 1/t328
      t335 = t241*t243
      t337 = 2*t246
      t338 = 4*t9
      t339 = 4*t29
      t340 = 2*t5
      t342 = (t14+t337+t15+t3+t338-t36+t339-t340+t307)*t326
      t345 = t301*t235
      t351 = dlog(l1)
      t356 = t351**2
      t368 = t279+t280-t281+t282+t283+t284-t285+t286+48*t351*t9+12*t351*
     #t3+12*t356-12*t351*t36+48*t351*t29-24*t351*t5+12*t15*t351-t287+12*
     #t14*t351+t289-t290
      t370 = 1/l1
      t375 = 2*t351
      t384 = (t14+t375+t3+t338-t36+t339-t340+t307+t15)*t370
      tncoff1 = (t42+t68+t91+t119+t150+t175+t202+t230)*t235*t244-(t273+t
     #299)*t301*t235*t244/6144+(t327*t329/12288-t327*t235/12288)*t239*t3
     #35+(-t342*t301*t329/1024+t342*t345/1024)*t239*t335-(t306+t368)*t37
     #0*t235*t244/6144+(t375+t3+t338-t36+t339-t340+t15+t14+t337)*t370*t3
     #01*t235*t239*t335/512+(-t384*t326*t329/1024+t384*t326*t235/1024)*t
     #239*t335+(t329*t370*t326*t301/512-t370*t326*t345/512)*t239*t335
      t1 = 0.3141592653589793D1**2
      t3 = dlog(2.D0)
      t4 = t3**2
      t6 = dlog(z)
      t7 = t6**2
      t10 = dlog(1-z)
      t11 = t10**2
      t14 = dlog(1-l4)
      t15 = t14**2
      t17 = dlog(l4)
      t18 = t17**2
      t20 = dlog(l3)
      t31 = -2*t1+48*t4+12*t7+48*t11+3*t15+3*t18+12*t14*t20-24*t20*t6+12
     #*t17*t20+48*t20*t10+48*t20*t3
      t46 = t20**2
      t48 = 1-l3
      t49 = dlog(t48)
      t50 = t49**2
      t52 = f3(t48,z)
      t54 = dlog(1-t52)
      t55 = t54**2
      t57 = f2(t48,z)
      t58 = dlog(t57)
      t59 = t58**2
      t63 = 24*t17*t3-48*t10*t6-12*t14*t6-48*t6*t3+96*t3*t10-12*t17*t6+2
     #4*t17*t10+12*t46+12*t50+3*t55+3*t59+6*t14*t54
      t87 = -12*t54*t6-6*t14*t58+24*t54*t10+12*t14*t49+12*t58*t6-24*t58*
     #t3+48*t10*t49-24*t10*t58+12*t54*t49-6*t54*t58+12*t20*t54
      t112 = 24*t20*t49-12*t49*t58+24*t54*t3-24*t49*t6+48*t49*t3-12*t20*
     #t58+6*t17*t54+12*t17*t49-6*t17*t58+24*t14*t3+6*t14*t17+24*t14*t10
      t115 = 1/l3
      t117 = Dsquar312(z,0.D0,0.D0,t48,l4)
      t118 = 1/t117
      t120 = s2vr2(z,1.D0,1.D0,t48,l4)
      t121 = t120**2
      t122 = 1/t121
      t123 = s23r2(z,1.D0,1.D0,t48,l4)
      t124 = 1/t123
      t126 = f1(t48,z,1.D0)
      t130 = 2*t20
      t131 = 4*t10
      t132 = 2*t49
      t133 = dlog(l2)
      t134 = 2*t6
      t135 = 4*t3
      t137 = 1/l2
      t138 = (t14+t130+t17+t131+t132+t133-t134-t58+t135+t54)*t137
      t139 = Dsquar312(z,0.D0,l2,t48,l4)
      t140 = 1/t139
      t143 = t115*t118
      t148 = t124*t126
      t150 = dlog(l1)
      t153 = 1/l1
      tncoff2 = (t31+t63+t87+t112)*t115*t118*t122*t124*t126/6144+(t138*t
     #115*t140/1024-t138*t143/1024)*t122*t148-(t135+t131+t132-t134+t54+2
     #*t150-t58+t17+t14+t130)*t153*t115*t118*t122*t148/512+(-t140*t153*t
     #137*t115/512+t153*t137*t143/512)*t122*t148
      t1 = dlog(l1)
      t2 = dlog(2.D0)
      t5 = dlog(l4)
      t6 = dlog(z)
      t12 = dlog(1-z)
      t16 = dlog(1-l4)
      t17 = f3(1.D0,z)
      t20 = dlog(1-t17+t17*l1)
      t23 = f2(1.D0,z)
      t24 = dlog(t23)
      t27 = 1-l1
      t28 = dlog(t27)
      t33 = t20**2
      t35 = t28**2
      t37 = t24**2
      t39 = 48*t1*t2-12*t5*t6-24*t1*t6+96*t12*t2+6*t16*t20-6*t16*t24+6*t
     #16*t28-12*t20*t6+3*t33+3*t35+3*t37
      t42 = t6**2
      t46 = t12**2
      t48 = t16**2
      t50 = t5**2
      t52 = t1**2
      t64 = 12*t24*t6+12*t42-6*t5*t24+48*t46+3*t48+3*t50+12*t52-12*t1*t2
     #4+12*t5*t1+24*t16*t2+24*t16*t12+12*t16*t1
      t88 = 24*t20*t2+12*t1*t28-48*t12*t6+24*t5*t2+6*t20*t28-12*t28*t6+2
     #4*t28*t2-6*t20*t24-6*t28*t24+12*t1*t20+24*t12*t28
      t107 = 0.3141592653589793D1**2
      t109 = t2**2
      t113 = 6*t5*t20+6*t5*t28+24*t20*t12-12*t16*t6-48*t6*t2+24*t5*t12-2
     #4*t24*t2+6*t16*t5-24*t12*t24-2*t107+48*t109+48*t1*t12
      t116 = 1/l1
      t118 = Dsquar312(z,l1,0.D0,1.D0,l4)
      t119 = 1/t118
      t121 = s2vr2(z,t27,1.D0,1.D0,l4)
      t122 = t121**2
      t123 = 1/t122
      t124 = s23r2(z,t27,1.D0,1.D0,l4)
      t125 = 1/t124
      t127 = f1(1.D0,z,t27)
      t131 = 2*t1
      t132 = 4*t12
      t133 = 4*t2
      t134 = 2*t6
      t135 = dlog(l3)
      t139 = 1/l3
      t142 = t125*t127
      t143 = t119*t123*t142
      t146 = dlog(l2)
      t149 = 1/l2
      tncoff3 = (t39+t64+t88+t113)*t116*t119*t123*t125*t127/6144-(t131+t
     #20+t132-t24+t133-t134+t28+t5+t16+2*t135)*t116*t139*t143/512-(t16+t
     #131+t20+t132-t24+t133-t134+t28+t146+t5)*t116*t149*t143/1024+t116*t
     #123*t119*t149*t139*t142/512
      t1 = 1-l3
      t2 = f3(t1,z)
      t5 = dlog(1-t2+t2*l1)
      t6 = 1-l1
      t7 = dlog(t6)
      t8 = dlog(l1)
      t10 = dlog(t1)
      t12 = dlog(z)
      t14 = f2(t1,z)
      t15 = dlog(t14)
      t17 = dlog(1-z)
      t19 = dlog(2.D0)
      t21 = dlog(l4)
      t23 = dlog(1-l4)
      t24 = dlog(l3)
      t26 = t5+t7+2*t8+2*t10-2*t12-t15+4*t17+4*t19+t21+t23+2*t24
      t27 = 1/l1
      t29 = 1/l3
      t31 = Dsquar312(z,l1,0.D0,t1,l4)
      t32 = 1/t31
      t33 = s2vr2(z,t6,1.D0,t1,l4)
      t34 = t33**2
      t35 = 1/t34
      t37 = s23r2(z,t6,1.D0,t1,l4)
      t39 = f1(t1,z,t6)
      t40 = 1/t37*t39
      tncoff4 = t26*t27*t29*t32*t35*t40/512-t27*t35*t32/l2*t29*t40/512
      t2 = dlog(1-l4)
      t3 = f2(1.D0,z)
      t4 = dlog(t3)
      t5 = 1-l1
      t6 = dlog(t5)
      t8 = dlog(1-z)
      t10 = dlog(l1)
      t13 = 1-l1*l2
      t14 = dlog(t13)
      t15 = dlog(z)
      t17 = dlog(2.D0)
      t19 = f3(1.D0,z)
      t22 = dlog(1-t19+t19*l1)
      t23 = dlog(l2)
      t24 = dlog(l4)
      t25 = t2-t4+t6+4*t8+2*t10+t14-2*t15+4*t17+t22+t23+t24
      t26 = 1/l1
      t28 = 1/l2
      t30 = Dsquar312(z,l1,l2,1.D0,l4)
      t31 = 1/t30
      t32 = s2vr2(z,t5,t13,1.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = s23r2(z,t5,t13,1.D0,l4)
      t37 = 1/t36
      t38 = f1(1.D0,z,t5)
      tncoff5 = t25*t26*t28*t31*t34*t37*t38/1024-t34*t37*t31*t26*t28/l3*
     #t38/512
      t1 = 1-l1
      t3 = 1-l1*l2
      t4 = 1-l3
      t5 = s2vr2(z,t1,t3,t4,l4)
      t6 = t5**2
      t8 = s23r2(z,t1,t3,t4,l4)
      t11 = Dsquar312(z,l1,l2,t4,l4)
      t18 = f1(t4,z,t1)
      tncoff6 = 1/t6/t8/t11/l1/l2/l3*t18/512
      t2 = dlog(1-l4)
      t3 = f3(1.D0,z)
      t5 = dlog(1-t3)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t13 = dlog(1-z)
      t18 = dlog(2.D0)
      t21 = dlog(z)
      t26 = dlog(l4)
      t41 = -t5*t2/2048+t2*t9/2048-t5*t13/512+t5*t9/2048-t5*t18/512+t5*t
     #21/1024-t9*t21/1024-t26*t5/2048+t26*t9/2048-t2*t13/512-t2*t26/2048
     #+t2*t21/1024-t2*t18/512+t13*t9/512
      t50 = t5**2
      t52 = t9**2
      t54 = t21**2
      t56 = t13**2
      t58 = t2**2
      t60 = t26**2
      t68 = 0.3141592653589793D1**2
      t70 = t18**2
      t72 = t26*t21/1024-t13*t18/128-t26*t13/512+t21*t18/256-t50/4096-t5
     #2/4096-t54/1024-t56/256-t58/4096-t60/4096+t9*t18/512+t13*t21/256-t
     #26*t18/512+t68/6144-t70/256
      t74 = Dsquar312(z,0.D0,0.D0,1.D0,l4)
      t75 = 1/t74
      t77 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t78 = t77**2
      t79 = 1/t78
      t80 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t81 = 1/t80
      t83 = f1(1.D0,z,1.D0)
      t84 = t79*t81*t83
      t86 = 6*t9
      t87 = 24*t13
      t88 = 24*t18
      t89 = 12*t21
      t90 = 6*t5
      t91 = 6*t2
      t92 = 6*t26
      t93 = dlog(l3)
      t96 = 1/l3
      t101 = dlog(l2)
      t104 = 1/l2
      t105 = (t86-t87-t88+t89-t90-t91-t92-6*t101)*t104
      t106 = Dsquar312(z,0.D0,l2,1.D0,l4)
      t107 = 1/t106
      t113 = t81*t83
      t115 = t104*t96
      t122 = dlog(l1)
      t125 = 1/l1
      t134 = t125*t104
      tncoff7 = (t41+t72)*t75*t84-(t86-t87-t88+t89-t90-t91-t92-12*t93)*t
     #96*t75*t84/6144+(t105*t107/12288-t105*t75/12288)*t79*t113+(t115*t1
     #07/1024-t115*t75/1024)*t79*t113-(t86-t87-t88+t89-t90-t91-t92-12*t1
     #22)*t125*t75*t84/6144-t125*t96*t75*t84/512+(t134*t107/1024-t134*t7
     #5/1024)*t79*t113
      t1 = 1-l3
      t2 = dlog(t1)
      t4 = dlog(2.D0)
      t6 = f3(t1,z)
      t8 = dlog(1-t6)
      t10 = dlog(z)
      t12 = dlog(l4)
      t15 = dlog(1-l4)
      t17 = dlog(l3)
      t20 = dlog(1-z)
      t22 = f2(t1,z)
      t23 = dlog(t22)
      t26 = 1/l3
      t28 = Dsquar312(z,0.D0,0.D0,t1,l4)
      t29 = 1/t28
      t31 = s2vr2(z,1.D0,1.D0,t1,l4)
      t32 = t31**2
      t33 = 1/t32
      t34 = s23r2(z,1.D0,1.D0,t1,l4)
      t35 = 1/t34
      t37 = f1(t1,z,1.D0)
      t38 = t33*t35*t37
      t42 = 1/l2*t26
      t43 = Dsquar312(z,0.D0,l2,t1,l4)
      tncoff8 = (-12*t2-24*t4-6*t8+12*t10-6*t12-6*t15-12*t17-24*t20+6*t2
     #3)*t26*t29*t38/6144+(-t42/t43/1024+t42*t29/1024)*t33*t35*t37+1/l1*
     #t26*t29*t38/512
      t1 = dlog(l1)
      t3 = dlog(l4)
      t6 = dlog(1-l4)
      t9 = dlog(1-z)
      t11 = f2(1.D0,z)
      t12 = dlog(t11)
      t14 = dlog(z)
      t16 = dlog(2.D0)
      t18 = f3(1.D0,z)
      t21 = dlog(1-t18+t18*l1)
      t23 = 1-l1
      t24 = dlog(t23)
      t27 = 1/l1
      t29 = Dsquar312(z,l1,0.D0,1.D0,l4)
      t30 = 1/t29
      t32 = s2vr2(z,t23,1.D0,1.D0,l4)
      t33 = t32**2
      t35 = s23r2(z,t23,1.D0,1.D0,l4)
      t38 = f1(1.D0,z,t23)
      t39 = 1/t33/t35*t38
      tncoff9 = (-12*t1-6*t3-6*t6-24*t9+6*t12+12*t14-24*t16-6*t21-6*t24)
     #*t27*t30*t39/6144+t27/l3*t30*t39/512+t27/l2*t30*t39/1024
      t4 = 1-l3
      t5 = Dsquar312(z,l1,0.D0,t4,l4)
      t8 = 1-l1
      t9 = s2vr2(z,t8,1.D0,t4,l4)
      t10 = t9**2
      t12 = s23r2(z,t8,1.D0,t4,l4)
      t15 = f1(t4,z,t8)
      tncoff10 = -1/l1/l3/t5/t10/t12*t15/512
      t4 = Dsquar312(z,l1,l2,1.D0,l4)
      t7 = 1-l1
      t9 = 1-l1*l2
      t10 = s2vr2(z,t7,t9,1.D0,l4)
      t11 = t10**2
      t13 = s23r2(z,t7,t9,1.D0,l4)
      t16 = f1(1.D0,z,t7)
      tncoff11 = -1/l1/l2/t4/t11/t13*t16/1024
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = f2(1.D0,z)
      t7 = dlog(t6)
      t9 = dlog(2.D0)
      t11 = dlog(z)
      t14 = dlog(1-z)
      t16 = f3(1.D0,z)
      t18 = dlog(1-t16)
      t21 = Dsquar312(z,0.D0,0.D0,1.D0,l4)
      t22 = 1/t21
      t24 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t25 = t24**2
      t26 = 1/t25
      t27 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t28 = 1/t27
      t30 = f1(1.D0,z,1.D0)
      t31 = t26*t28*t30
      t37 = 1/l2
      t38 = Dsquar312(z,0.D0,l2,1.D0,l4)
      tncoff12 = (t2/2048+t4/2048-t7/2048+t9/512-t11/1024+t14/512+t18/20
     #48)*t22*t31-1/l3*t22*t31/1024+(t37/t38/2048-t37*t22/2048)*t26*t28*
     #t30-1/l1*t22*t31/1024
      t2 = 1-l3
      t3 = Dsquar312(z,0.D0,0.D0,t2,l4)
      t6 = s2vr2(z,1.D0,1.D0,t2,l4)
      t7 = t6**2
      t9 = s23r2(z,1.D0,1.D0,t2,l4)
      t12 = f1(t2,z,1.D0)
      tncoff13 = 1/l3/t3/t7/t9*t12/1024
      t2 = Dsquar312(z,l1,0.D0,1.D0,l4)
      t5 = 1-l1
      t6 = s2vr2(z,t5,1.D0,1.D0,l4)
      t7 = t6**2
      t9 = s23r2(z,t5,1.D0,1.D0,l4)
      t12 = f1(1.D0,z,t5)
      tncoff14 = 1/l1/t2/t7/t9*t12/1024
      t1 = Dsquar312(z,0.D0,0.D0,1.D0,l4)
      t3 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t4 = t3**2
      t7 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t9 = f1(1.D0,z,1.D0)
      tncoff15 = -1/t1/t4/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num2e0(var,ndim,z,ux,1.D0,1.D0,1-l3,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num2e0(var,ndim,z,ux,1-l1,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num2e0(var,ndim,z,ux,1-l1,1.D0,1-l3,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num2e0(var,ndim,z,ux,1-l1,1-l1*l2,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num2e0(var,ndim,z,ux,1-l1,1-l1*l2,1-l3,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num2e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num2e1(var,ndim,z,ux,1.D0,1.D0,1-l3,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num2e1(var,ndim,z,ux,1-l1,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = sq3num2e1(var,ndim,z,ux,1-l1,1.D0,1-l3,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = sq3num2e1(var,ndim,z,ux,1-l1,1-l1*l2,1.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = sq3num2e2(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = sq3num2e2(var,ndim,z,ux,1.D0,1.D0,1-l3,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = sq3num2e2(var,ndim,z,ux,1-l1,1.D0,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = sq3num2e3(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff15)

CCC Sum

      sq34e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+
     #tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq34e1int(var,ndim,Tcoeff)

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
      double precision sq3num2e0
      external sq3num2e0
      double precision sq3num2e1
      external sq3num2e1
      double precision sq3num2e2
      external sq3num2e2
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

      t1 = f3(1.D0,z)
      t3 = dlog(1-t1)
      t4 = t3**2
      t6 = dlog(l4)
      t7 = dlog(z)
      t10 = f2(1.D0,z)
      t11 = dlog(t10)
      t12 = t11**2
      t14 = dlog(2.D0)
      t18 = dlog(1-z)
      t22 = dlog(1-l4)
      t33 = t7**2
      t37 = t6**2
      t39 = t18**2
      t41 = -t4/4096+t6*t7/1024-t12/4096+t11*t14/512+t18*t11/512-t22*t14
     #/512-t22*t6/2048-t22*t18/512-t3*t14/512+t3*t11/2048-t33/1024-t3*t1
     #8/512-t37/4096-t39/256
      t42 = t22**2
      t68 = 0.3141592653589793D1**2
      t70 = t14**2
      t72 = -t42/4096+t22*t11/2048+t3*t7/1024-t11*t7/1024-t6*t3/2048+t6*
     #t11/2048-t6*t18/512-t22*t3/2048-t6*t14/512+t18*t7/256+t22*t7/1024+
     #t7*t14/256-t18*t14/128+t68/6144-t70/256
      t74 = Dsquar312(z,0.D0,0.D0,1.D0,l4)
      t75 = 1/t74
      t77 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t78 = t77**2
      t79 = 1/t78
      t80 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t81 = 1/t80
      t83 = f1(1.D0,z,1.D0)
      t84 = t79*t81*t83
      t86 = 4*t18
      t87 = 4*t14
      t88 = 2*t7
      t89 = dlog(l3)
      t92 = 1/l3
      t97 = dlog(l2)
      t99 = 1/l2
      t100 = (t22+t6+t3+t86-t11+t87-t88+t97)*t99
      t101 = Dsquar312(z,0.D0,l2,1.D0,l4)
      t102 = 1/t101
      t108 = t81*t83
      t110 = t99*t92
      t117 = dlog(l1)
      t120 = 1/l1
      t129 = t120*t99
      tncoff1 = (t41+t72)*t75*t84+(t6+t22+t3+t86-t11+t87-t88+2*t89)*t92*
     #t75*t84/1024+(-t100*t102/2048+t100*t75/2048)*t79*t108+(t110*t102/1
     #024-t110*t75/1024)*t79*t108+(2*t117+t3+t86-t11+t87-t88+t22+t6)*t12
     #0*t75*t84/1024-t120*t92*t75*t84/512+(t129*t102/1024-t129*t75/1024)
     #*t79*t108
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = 1-l3
      t5 = f3(t4,z)
      t7 = dlog(1-t5)
      t9 = dlog(1-z)
      t11 = f2(t4,z)
      t12 = dlog(t11)
      t13 = dlog(2.D0)
      t15 = dlog(z)
      t17 = dlog(t4)
      t19 = dlog(l3)
      t22 = 1/l3
      t24 = Dsquar312(z,0.D0,0.D0,t4,l4)
      t25 = 1/t24
      t27 = s2vr2(z,1.D0,1.D0,t4,l4)
      t28 = t27**2
      t29 = 1/t28
      t30 = s23r2(z,1.D0,1.D0,t4,l4)
      t31 = 1/t30
      t33 = f1(t4,z,1.D0)
      t34 = t29*t31*t33
      t38 = 1/l2*t22
      t39 = Dsquar312(z,0.D0,l2,t4,l4)
      tncoff2 = -(t1+t3+t7+4*t9-t12+4*t13-2*t15+2*t17+2*t19)*t22*t25*t34
     #/1024+(-t38/t39/1024+t38*t25/1024)*t29*t31*t33+1/l1*t22*t25*t34/51
     #2
      t1 = dlog(l1)
      t3 = f2(1.D0,z)
      t4 = dlog(t3)
      t5 = f3(1.D0,z)
      t8 = dlog(1-t5+t5*l1)
      t10 = dlog(1-z)
      t12 = dlog(2.D0)
      t14 = dlog(z)
      t16 = 1-l1
      t17 = dlog(t16)
      t19 = dlog(1-l4)
      t20 = dlog(l4)
      t22 = 1/l1
      t24 = Dsquar312(z,l1,0.D0,1.D0,l4)
      t25 = 1/t24
      t27 = s2vr2(z,t16,1.D0,1.D0,l4)
      t28 = t27**2
      t30 = s23r2(z,t16,1.D0,1.D0,l4)
      t33 = f1(1.D0,z,t16)
      t34 = 1/t28/t30*t33
      tncoff3 = -(2*t1-t4+t8+4*t10+4*t12-2*t14+t17+t19+t20)*t22*t25*t34/
     #1024+t22/l3*t25*t34/512+t22/l2*t25*t34/1024
      t4 = 1-l3
      t5 = Dsquar312(z,l1,0.D0,t4,l4)
      t8 = 1-l1
      t9 = s2vr2(z,t8,1.D0,t4,l4)
      t10 = t9**2
      t12 = s23r2(z,t8,1.D0,t4,l4)
      t15 = f1(t4,z,t8)
      tncoff4 = -1/l1/l3/t5/t10/t12*t15/512
      t4 = Dsquar312(z,l1,l2,1.D0,l4)
      t7 = 1-l1
      t9 = 1-l1*l2
      t10 = s2vr2(z,t7,t9,1.D0,l4)
      t11 = t10**2
      t13 = s23r2(z,t7,t9,1.D0,l4)
      t16 = f1(1.D0,z,t7)
      tncoff5 = -1/l1/l2/t4/t11/t13*t16/1024
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = f2(1.D0,z)
      t7 = dlog(t6)
      t9 = dlog(2.D0)
      t11 = dlog(z)
      t14 = dlog(1-z)
      t16 = f3(1.D0,z)
      t18 = dlog(1-t16)
      t21 = Dsquar312(z,0.D0,0.D0,1.D0,l4)
      t22 = 1/t21
      t24 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t25 = t24**2
      t26 = 1/t25
      t27 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t28 = 1/t27
      t30 = f1(1.D0,z,1.D0)
      t31 = t26*t28*t30
      t37 = 1/l2
      t38 = Dsquar312(z,0.D0,l2,1.D0,l4)
      tncoff6 = (t2/2048+t4/2048-t7/2048+t9/512-t11/1024+t14/512+t18/204
     #8)*t22*t31-1/l3*t22*t31/1024+(t37/t38/2048-t37*t22/2048)*t26*t28*t
     #30-1/l1*t22*t31/1024
      t2 = 1-l3
      t3 = Dsquar312(z,0.D0,0.D0,t2,l4)
      t6 = s2vr2(z,1.D0,1.D0,t2,l4)
      t7 = t6**2
      t9 = s23r2(z,1.D0,1.D0,t2,l4)
      t12 = f1(t2,z,1.D0)
      tncoff7 = 1/l3/t3/t7/t9*t12/1024
      t2 = Dsquar312(z,l1,0.D0,1.D0,l4)
      t5 = 1-l1
      t6 = s2vr2(z,t5,1.D0,1.D0,l4)
      t7 = t6**2
      t9 = s23r2(z,t5,1.D0,1.D0,l4)
      t12 = f1(1.D0,z,t5)
      tncoff8 = 1/l1/t2/t7/t9*t12/1024
      t1 = Dsquar312(z,0.D0,0.D0,1.D0,l4)
      t3 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t4 = t3**2
      t7 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t9 = f1(1.D0,z,1.D0)
      tncoff9 = -1/t1/t4/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num2e0(var,ndim,z,ux,1.D0,1.D0,1-l3,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num2e0(var,ndim,z,ux,1-l1,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num2e0(var,ndim,z,ux,1-l1,1.D0,1-l3,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num2e0(var,ndim,z,ux,1-l1,1-l1*l2,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num2e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num2e1(var,ndim,z,ux,1.D0,1.D0,1-l3,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num2e1(var,ndim,z,ux,1-l1,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num2e2(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      sq34e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq34e2int(var,ndim,Tcoeff)

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
      double precision sq3num2e0
      external sq3num2e0
      double precision sq3num2e1
      external sq3num2e1
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

      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = f2(1.D0,z)
      t7 = dlog(t6)
      t9 = dlog(2.D0)
      t11 = dlog(z)
      t14 = dlog(1-z)
      t16 = f3(1.D0,z)
      t18 = dlog(1-t16)
      t21 = Dsquar312(z,0.D0,0.D0,1.D0,l4)
      t22 = 1/t21
      t24 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t25 = t24**2
      t26 = 1/t25
      t27 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t28 = 1/t27
      t30 = f1(1.D0,z,1.D0)
      t31 = t26*t28*t30
      t37 = 1/l2
      t38 = Dsquar312(z,0.D0,l2,1.D0,l4)
      tncoff1 = (t2/2048+t4/2048-t7/2048+t9/512-t11/1024+t14/512+t18/204
     #8)*t22*t31-1/l3*t22*t31/1024+(t37/t38/2048-t37*t22/2048)*t26*t28*t
     #30-1/l1*t22*t31/1024
      t2 = 1-l3
      t3 = Dsquar312(z,0.D0,0.D0,t2,l4)
      t6 = s2vr2(z,1.D0,1.D0,t2,l4)
      t7 = t6**2
      t9 = s23r2(z,1.D0,1.D0,t2,l4)
      t12 = f1(t2,z,1.D0)
      tncoff2 = 1/l3/t3/t7/t9*t12/1024
      t2 = Dsquar312(z,l1,0.D0,1.D0,l4)
      t5 = 1-l1
      t6 = s2vr2(z,t5,1.D0,1.D0,l4)
      t7 = t6**2
      t9 = s23r2(z,t5,1.D0,1.D0,l4)
      t12 = f1(1.D0,z,t5)
      tncoff3 = 1/l1/t2/t7/t9*t12/1024
      t1 = Dsquar312(z,0.D0,0.D0,1.D0,l4)
      t3 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t4 = t3**2
      t7 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t9 = f1(1.D0,z,1.D0)
      tncoff4 = -1/t1/t4/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num2e0(var,ndim,z,ux,1.D0,1.D0,1-l3,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num2e0(var,ndim,z,ux,1-l1,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num2e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      sq34e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq34e3int(var,ndim,Tcoeff)

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
      double precision sq3num2e0
      external sq3num2e0
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

      t1 = Dsquar312(z,0.D0,0.D0,1.D0,l4)
      t3 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t4 = t3**2
      t7 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t9 = f1(1.D0,z,1.D0)
      tncoff1 = -1/t1/t4/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      sq34e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq34e4int(var,ndim,Tcoeff)

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

      sq34e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

