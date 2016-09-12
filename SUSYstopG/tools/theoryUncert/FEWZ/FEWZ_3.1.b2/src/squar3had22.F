CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function sq32e0int(var,ndim,Tcoeff)

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
      double precision sq3num1e0
      external sq3num1e0
      double precision sq3num1e1
      external sq3num1e1
      double precision sq3num1e2
      external sq3num1e2
      double precision sq3num1e3
      external sq3num1e3
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
      t3 = f2(0.D0,z)
      t4 = dlog(t3)
      t7 = dlog(z)
      t8 = t1*t7
      t9 = dlog(2.D0)
      t12 = f3(0.D0,z)
      t14 = dlog(1-t12)
      t15 = t1*t14
      t17 = dlog(1-z)
      t21 = dlog(1-l4)
      t22 = t21**2
      t25 = t4**2
      t30 = t7**2
      t33 = 0.3141592653589793D1**2
      t36 = t9**2
      t39 = t17**2
      t44 = -t2*t4/4096-t8*t9/256+t15*t17/512+t22*t9/1024+t1*t25/4096+t2
     #1*t25/4096+t21*t30/1024-t21*t33/6144+t21*t36/256-t39*t4/256-t1*t33
     #/6144
      t51 = t14**2
      t66 = t1*t4
      t69 = t21*t4
      t72 = t1*t36/256+t21*t2/4096+t21*t39/256+t21*t51/4096+t22*t17/1024
     #-t22*t7/2048-t14*t33/6144+t14*t36/256+t2*t14/4096-t33*t9/1536+t66*
     #t7/1024+t69*t7/1024
      t96 = t2*t9/1024-t2*t7/2048+t1*t39/256+t1*t51/4096+t1*t30/1024+t33
     #*t7/3072+t30*t9/256-t7*t36/128+t39*t9/64-t39*t7/128+t17*t25/1024
      t97 = t1*t17
      t100 = t21*t1
      t108 = t21*t14
      t119 = t14*t7
      t122 = t14*t4
      t125 = -t97*t7/256-t100*t7/1024+t36*t9/192+z3/768+t4*t33/6144-t108
     #*t7/1024+t100*t9/512-t51*t4/4096-t15*t4/2048+t15*t9/512-t119*t9/25
     #6+t122*t7/1024
      t134 = t21*t17
      t141 = t17*t7
      t152 = -t69*t9/512-t122*t9/512+t100*t17/512+t134*t9/128+t108*t17/5
     #12-t134*t4/512-t141*t9/64+t17*t30/256-t17*t33/1536+t17*t36/64+t25*
     #t9/1024
      t155 = t21*t7
      t162 = t14*t17
      t169 = t17*t4
      t180 = t51*t9/1024-t155*t9/256-t134*t7/256-t100*t4/2048-t162*t7/25
     #6-t162*t4/512+t108*t9/512-t169*t9/128-t108*t4/2048-t25*t7/2048-t4*
     #t30/1024+t51*t17/1024
      t204 = -t51*t7/2048+t14*t39/256+t14*t25/4096+t14*t30/1024+t22*t1/4
     #096+t22*t14/4096-t4*t36/256+t2*t17/1024+t22*t21/12288-t15*t7/1024+
     #t100*t14/2048
      t209 = t4*t7
      t230 = t169*t7/256+t162*t9/128+t209*t9/256+t51*t14/12288+t2*t1/122
     #88+t39*t17/192-t25*t4/12288-t30*t7/1536-t22*t4/4096-t66*t9/512-t97
     #*t4/512+t97*t9/128
      t234 = Dsquar312(z,0.D0,0.D0,0.D0,l4)
      t235 = 1/t234
      t237 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t238 = 1/t237
      t239 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t240 = t239**2
      t241 = 1/t240
      t243 = f1(0.D0,z,1.D0)
      t244 = t238*t241*t243
      t246 = 3*t22
      t247 = 2*t33
      t248 = 48*t36
      t249 = 3*t2
      t250 = 24*t162
      t251 = dlog(l3)
      t265 = 24*t1*t9
      t266 = 6*t122
      t268 = 24*t4*t9
      t269 = 48*t141
      t272 = t251**2
      t274 = 3*t51
      t275 = t246-t247+t248+t249+t250+12*t251*t14+12*t21*t251+48*t251*t9
     #-24*t251*t7+48*t251*t17-12*t251*t4+t265-t266-t268-t269+12*t1*t251+
     #12*t272+t274
      t277 = 24*t14*t9
      t278 = 12*t209
      t279 = 12*t155
      t281 = 48*t7*t9
      t282 = 24*t134
      t283 = 6*t15
      t285 = 96*t17*t9
      t286 = 6*t66
      t287 = 12*t8
      t288 = 24*t97
      t289 = 6*t69
      t290 = 6*t108
      t291 = 12*t119
      t293 = 24*t21*t9
      t294 = 6*t100
      t295 = 3*t25
      t296 = 12*t30
      t297 = 24*t169
      t298 = 48*t39
      t299 = t277+t278-t279-t281+t282+t283+t285-t286-t287+t288-t289+t290
     #-t291+t293+t294+t295+t296-t297+t298
      t301 = 1/l3
      t306 = dlog(l2)
      t311 = t306**2
      t313 = 24*t306*t9-12*t306*t7+3*t311+t246-t247+t248+t249+t250+t265-
     #t266-t268-t269+t274+t277+t278-t279-t281+t282
      t324 = t283+t285-t286-t287+t288-t289+t290-t291+t293+t294+t295+t296
     #-t297+6*t21*t306+6*t1*t306-6*t306*t4+6*t306*t14+24*t306*t17+t298
      t326 = 1/l2
      t327 = (t313+t324)*t326
      t328 = Dsquar312(z,0.D0,l2,0.D0,l4)
      t329 = 1/t328
      t335 = t241*t243
      t337 = 2*t251
      t338 = 4*t17
      t339 = 4*t9
      t340 = 2*t7
      t342 = (-t337-t1-t21-t306-t14-t338+t4-t339+t340)*t326
      t345 = t301*t235
      t351 = dlog(l1)
      t366 = t351**2
      t368 = 12*t351*t14+48*t351*t17+12*t21*t351-24*t351*t7+12*t1*t351-1
     #2*t351*t4+48*t351*t9+t246+12*t366-t247+t248+t249+t250+t265-t266-t2
     #68-t269+t274
      t370 = 1/l1
      t375 = 2*t351
      t384 = (-t21-t1+t4-t338-t375-t339-t14+t340-t306)*t370
      tncoff1 = (t44+t72+t96+t125+t152+t180+t204+t230)*t235*t244-(t275+t
     #299)*t301*t235*t244/6144+(t327*t329/12288-t327*t235/12288)*t238*t3
     #35+(t342*t301*t329/1024-t342*t345/1024)*t238*t335-(t368+t299)*t370
     #*t235*t244/6144-(-t375-t14-t338+t4-t339+t340-t337-t21-t1)*t370*t30
     #1*t235*t238*t335/512+(t384*t326*t329/1024-t384*t326*t235/1024)*t23
     #8*t335+(t329*t370*t326*t301/512-t370*t326*t345/512)*t238*t335
      t2 = dlog(1-l4)
      t3 = f3(l3,z)
      t5 = dlog(1-t3)
      t8 = dlog(l4)
      t9 = dlog(l3)
      t14 = f2(l3,z)
      t15 = dlog(t14)
      t16 = dlog(z)
      t19 = dlog(2.D0)
      t27 = dlog(1-z)
      t36 = 6*t2*t5+12*t8*t9+12*t2*t9+12*t15*t16+24*t8*t19+48*t9*t19-24*
     #t9*t16+48*t9*t27+24*t8*t27-48*t16*t19+24*t2*t19
      t48 = dlog(1-l3)
      t63 = -48*t27*t16+24*t2*t27+6*t2*t8-12*t8*t16-12*t2*t16-12*t48*t15
     #-24*t27*t15+96*t27*t19+24*t5*t19+12*t48*t5-24*t15*t19+24*t5*t27
      t87 = 48*t48*t27-12*t5*t16-6*t5*t15+48*t48*t19-24*t48*t16-12*t9*t1
     #5+6*t8*t5+12*t9*t5+12*t2*t48-6*t2*t15+24*t9*t48
      t92 = 0.3141592653589793D1**2
      t94 = t19**2
      t96 = t9**2
      t98 = t48**2
      t100 = t15**2
      t102 = t5**2
      t104 = t16**2
      t106 = t27**2
      t108 = t2**2
      t110 = t8**2
      t112 = 12*t8*t48-6*t8*t15-2*t92+48*t94+12*t96+12*t98+3*t100+3*t102
     #+12*t104+48*t106+3*t108+3*t110
      t115 = 1/l3
      t117 = Dsquar312(z,0.D0,0.D0,l3,l4)
      t118 = 1/t117
      t120 = s23r2(z,1.D0,1.D0,l3,l4)
      t121 = 1/t120
      t122 = s1vr2(z,1.D0,1.D0,l3,l4)
      t123 = t122**2
      t124 = 1/t123
      t126 = f1(l3,z,1.D0)
      t130 = 2*t9
      t131 = dlog(l2)
      t132 = 4*t27
      t133 = 4*t19
      t134 = 2*t48
      t135 = 2*t16
      t137 = 1/l2
      t138 = (t130+t8+t2+t131+t132+t133-t15+t134-t135+t5)*t137
      t139 = Dsquar312(z,0.D0,l2,l3,l4)
      t140 = 1/t139
      t143 = t115*t118
      t148 = t124*t126
      t150 = dlog(l1)
      t153 = 1/l1
      tncoff2 = (t36+t63+t87+t112)*t115*t118*t121*t124*t126/6144+(t138*t
     #115*t140/1024-t138*t143/1024)*t121*t148-(t134+2*t150-t135-t15+t5+t
     #132+t133+t130+t2+t8)*t153*t115*t118*t121*t148/512+(-t140*t153*t137
     #*t115/512+t153*t137*t143/512)*t121*t148
      t2 = dlog(1-z)
      t3 = t2**2
      t6 = dlog(1-l4)
      t7 = t6**2
      t9 = dlog(l1)
      t10 = t9**2
      t12 = dlog(l4)
      t13 = t12**2
      t15 = dlog(z)
      t16 = dlog(2.D0)
      t25 = f3(0.D0,z)
      t28 = dlog(1-t25+l1*t25)
      t29 = 1-l1
      t30 = dlog(t29)
      t37 = 48*t3+3*t7+12*t10+3*t13-48*t15*t16+24*t12*t2+48*t2*t9+12*t12
     #*t9+6*t28*t30+48*t9*t16+24*t6*t16
      t42 = f2(0.D0,z)
      t43 = dlog(t42)
      t64 = 6*t6*t12+24*t28*t16-6*t6*t43-24*t9*t15-24*t2*t43+24*t6*t2-12
     #*t6*t15-6*t28*t43-48*t2*t15+6*t6*t30+6*t12*t28+6*t6*t28
      t72 = 0.3141592653589793D1**2
      t76 = t16**2
      t86 = t28**2
      t88 = 12*t9*t28+24*t28*t2-12*t28*t15-2*t72+24*t30*t16+48*t76+24*t2
     #*t30-24*t43*t16+6*t12*t30+12*t9*t30+3*t86
      t93 = t30**2
      t109 = t15**2
      t111 = t43**2
      t113 = -12*t12*t15+96*t2*t16+3*t93-6*t12*t43+12*t6*t9+12*t43*t15-1
     #2*t9*t43+24*t12*t16-6*t30*t43-12*t30*t15+12*t109+3*t111
      t116 = 1/l1
      t118 = Dsquar312(z,l1,0.D0,0.D0,l4)
      t119 = 1/t118
      t121 = s23r2(z,t29,1.D0,0.D0,l4)
      t122 = 1/t121
      t123 = s1vr2(z,t29,1.D0,0.D0,l4)
      t124 = t123**2
      t125 = 1/t124
      t127 = f1(0.D0,z,t29)
      t131 = 2*t9
      t132 = 4*t2
      t133 = 2*t15
      t134 = 4*t16
      t135 = dlog(l3)
      t139 = 1/l3
      t142 = t125*t127
      t143 = t119*t122*t142
      t146 = dlog(l2)
      t149 = 1/l2
      tncoff3 = (t37+t64+t88+t113)*t116*t119*t122*t125*t127/6144+(-t131-
     #t132-t28+t133+t43-t30-t134-2*t135-t6-t12)*t116*t139*t143/512+(-t6-
     #t12-t131+t133-t134+t43-t30-t28-t132-t146)*t116*t149*t143/1024+t116
     #*t122*t119*t149*t139*t142/512
      t1 = dlog(z)
      t4 = dlog(1-z)
      t6 = 1-l1
      t7 = dlog(t6)
      t8 = dlog(2.D0)
      t10 = f2(l3,z)
      t11 = dlog(t10)
      t13 = dlog(1-l3)
      t15 = f3(l3,z)
      t18 = dlog(1-t15+l1*t15)
      t19 = dlog(l1)
      t21 = dlog(l3)
      t24 = dlog(1-l4)
      t25 = dlog(l4)
      t26 = -2*t1+4*t4+t7+4*t8-t11+2*t13+t18+2*t19+2*t21+t24+t25
      t27 = 1/l1
      t29 = 1/l3
      t31 = Dsquar312(z,l1,0.D0,l3,l4)
      t32 = 1/t31
      t33 = s23r2(z,t6,1.D0,l3,l4)
      t34 = 1/t33
      t36 = s1vr2(z,t6,1.D0,l3,l4)
      t37 = t36**2
      t39 = f1(l3,z,t6)
      t40 = 1/t37*t39
      tncoff4 = t26*t27*t29*t32*t34*t40/512-t27*t34*t32/l2*t29*t40/512
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t5 = dlog(1-z)
      t7 = dlog(z)
      t9 = dlog(2.D0)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t13 = f3(0.D0,z)
      t16 = dlog(1-t13+l1*t13)
      t17 = 1-l1
      t18 = dlog(t17)
      t20 = 1-l1*l2
      t21 = dlog(t20)
      t22 = dlog(l1)
      t24 = dlog(l2)
      t25 = -t2-t3-4*t5+2*t7-4*t9+t12-t16-t18-t21-2*t22-t24
      t26 = 1/l1
      t28 = 1/l2
      t30 = Dsquar312(z,l1,l2,0.D0,l4)
      t31 = 1/t30
      t32 = s23r2(z,t17,t20,0.D0,l4)
      t33 = 1/t32
      t35 = s1vr2(z,t17,t20,0.D0,l4)
      t36 = t35**2
      t37 = 1/t36
      t38 = f1(0.D0,z,t17)
      tncoff5 = -t25*t26*t28*t31*t33*t37*t38/1024-t37*t33*t31*t26*t28/l3
     #*t38/512
      t1 = 1-l1
      t3 = 1-l1*l2
      t4 = s1vr2(z,t1,t3,l3,l4)
      t5 = t4**2
      t7 = s23r2(z,t1,t3,l3,l4)
      t10 = Dsquar312(z,l1,l2,l3,l4)
      t17 = f1(l3,z,t1)
      tncoff6 = 1/t5/t7/t10/l1/l2/l3*t17/512
      t1 = dlog(z)
      t2 = t1**2
      t5 = dlog(1-l4)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t10 = t7**2
      t12 = f3(0.D0,z)
      t14 = dlog(1-t12)
      t15 = dlog(2.D0)
      t23 = dlog(1-z)
      t26 = dlog(l4)
      t41 = -t2/1024+t5*t7/2048-t10/4096-t14*t15/512+t14*t7/2048+t7*t15/
     #512+t23*t1/256-t26*t15/512-t23*t15/128-t26*t14/2048-t14*t23/512+t1
     #4*t1/1024+t26*t7/2048-t5*t14/2048
      t58 = t23**2
      t60 = t26**2
      t62 = t14**2
      t64 = t5**2
      t68 = 0.3141592653589793D1**2
      t70 = t15**2
      t72 = t1*t15/256+t23*t7/512-t5*t23/512-t5*t26/2048-t5*t15/512-t26*
     #t23/512+t26*t1/1024+t5*t1/1024-t58/256-t60/4096-t62/4096-t64/4096-
     #t7*t1/1024+t68/6144-t70/256
      t74 = Dsquar312(z,0.D0,0.D0,0.D0,l4)
      t75 = 1/t74
      t77 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t78 = 1/t77
      t79 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t80 = t79**2
      t81 = 1/t80
      t83 = f1(0.D0,z,1.D0)
      t84 = t78*t81*t83
      t86 = dlog(l3)
      t88 = 24*t15
      t89 = 12*t1
      t90 = 6*t5
      t91 = 6*t26
      t92 = 6*t14
      t93 = 24*t23
      t94 = 6*t7
      t96 = 1/l3
      t101 = dlog(l2)
      t104 = 1/l2
      t105 = (-t88+t89-t90-t91-t92-t93+t94-6*t101)*t104
      t106 = Dsquar312(z,0.D0,l2,0.D0,l4)
      t107 = 1/t106
      t113 = t81*t83
      t115 = t104*t96
      t122 = dlog(l1)
      t125 = 1/l1
      t134 = t125*t104
      tncoff7 = (t41+t72)*t75*t84-(-12*t86-t88+t89-t90-t91-t92-t93+t94)*
     #t96*t75*t84/6144+(t105*t107/12288-t105*t75/12288)*t78*t113+(t115*t
     #107/1024-t115*t75/1024)*t78*t113-(-t88+t89-t90-t91-t92-t93+t94-12*
     #t122)*t125*t75*t84/6144-t125*t96*t75*t84/512+(t134*t107/1024-t134*
     #t75/1024)*t78*t113
      t1 = dlog(z)
      t3 = f3(l3,z)
      t5 = dlog(1-t3)
      t8 = dlog(1-z)
      t11 = dlog(1-l3)
      t13 = dlog(l4)
      t16 = dlog(1-l4)
      t18 = dlog(l3)
      t20 = dlog(2.D0)
      t22 = f2(l3,z)
      t23 = dlog(t22)
      t26 = 1/l3
      t28 = Dsquar312(z,0.D0,0.D0,l3,l4)
      t29 = 1/t28
      t31 = s23r2(z,1.D0,1.D0,l3,l4)
      t32 = 1/t31
      t33 = s1vr2(z,1.D0,1.D0,l3,l4)
      t34 = t33**2
      t35 = 1/t34
      t37 = f1(l3,z,1.D0)
      t38 = t32*t35*t37
      t42 = 1/l2*t26
      t43 = Dsquar312(z,0.D0,l2,l3,l4)
      tncoff8 = (12*t1-6*t5-24*t8-12*t11-6*t13-6*t16-12*t18-24*t20+6*t23
     #)*t26*t29*t38/6144+(-t42/t43/1024+t42*t29/1024)*t32*t35*t37+1/l1*t
     #26*t29*t38/512
      t1 = dlog(l1)
      t4 = dlog(1-l4)
      t6 = dlog(l4)
      t8 = 1-l1
      t9 = dlog(t8)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t15 = dlog(1-z)
      t17 = f3(0.D0,z)
      t20 = dlog(1-t17+l1*t17)
      t22 = dlog(2.D0)
      t24 = dlog(z)
      t27 = 1/l1
      t29 = Dsquar312(z,l1,0.D0,0.D0,l4)
      t30 = 1/t29
      t32 = s23r2(z,t8,1.D0,0.D0,l4)
      t34 = s1vr2(z,t8,1.D0,0.D0,l4)
      t35 = t34**2
      t38 = f1(0.D0,z,t8)
      t39 = 1/t32/t35*t38
      tncoff9 = (-12*t1-6*t4-6*t6-6*t9+6*t12-24*t15-6*t20-24*t22+12*t24)
     #*t27*t30*t39/6144+t27/l3*t30*t39/512+t27/l2*t30*t39/1024
      t4 = Dsquar312(z,l1,0.D0,l3,l4)
      t7 = 1-l1
      t8 = s23r2(z,t7,1.D0,l3,l4)
      t10 = s1vr2(z,t7,1.D0,l3,l4)
      t11 = t10**2
      t14 = f1(l3,z,t7)
      tncoff10 = -1/l1/l3/t4/t8/t11*t14/512
      t4 = Dsquar312(z,l1,l2,0.D0,l4)
      t7 = 1-l1
      t9 = 1-l1*l2
      t10 = s23r2(z,t7,t9,0.D0,l4)
      t12 = s1vr2(z,t7,t9,0.D0,l4)
      t13 = t12**2
      t16 = f1(0.D0,z,t7)
      tncoff11 = -1/l1/l2/t4/t10/t13*t16/1024
      t1 = dlog(l4)
      t4 = dlog(1-l4)
      t6 = f3(0.D0,z)
      t8 = dlog(1-t6)
      t11 = dlog(1-z)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t16 = dlog(2.D0)
      t18 = dlog(z)
      t21 = Dsquar312(z,0.D0,0.D0,0.D0,l4)
      t22 = 1/t21
      t24 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t25 = 1/t24
      t26 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t27 = t26**2
      t28 = 1/t27
      t30 = f1(0.D0,z,1.D0)
      t31 = t25*t28*t30
      t37 = 1/l2
      t38 = Dsquar312(z,0.D0,l2,0.D0,l4)
      tncoff12 = (t1/2048+t4/2048+t8/2048+t11/512-t14/2048+t16/512-t18/1
     #024)*t22*t31-1/l3*t22*t31/1024+(t37/t38/2048-t37*t22/2048)*t25*t28
     #*t30-1/l1*t22*t31/1024
      t2 = Dsquar312(z,0.D0,0.D0,l3,l4)
      t5 = s23r2(z,1.D0,1.D0,l3,l4)
      t7 = s1vr2(z,1.D0,1.D0,l3,l4)
      t8 = t7**2
      t11 = f1(l3,z,1.D0)
      tncoff13 = 1/l3/t2/t5/t8*t11/1024
      t2 = Dsquar312(z,l1,0.D0,0.D0,l4)
      t5 = 1-l1
      t6 = s23r2(z,t5,1.D0,0.D0,l4)
      t8 = s1vr2(z,t5,1.D0,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,z,t5)
      tncoff14 = 1/l1/t2/t6/t9*t12/1024
      t1 = Dsquar312(z,0.D0,0.D0,0.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t6 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,1.D0)
      tncoff15 = -1/t1/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num1e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num1e0(var,ndim,z,ux,1-l1,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num1e0(var,ndim,z,ux,1-l1,1.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num1e0(var,ndim,z,ux,1-l1,1-l1*l2,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num1e0(var,ndim,z,ux,1-l1,1-l1*l2,l3,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num1e1(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num1e1(var,ndim,z,ux,1-l1,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = sq3num1e1(var,ndim,z,ux,1-l1,1.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = sq3num1e1(var,ndim,z,ux,1-l1,1-l1*l2,0.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = sq3num1e2(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = sq3num1e2(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = sq3num1e2(var,ndim,z,ux,1-l1,1.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = sq3num1e3(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff15)

CCC Sum

      sq32e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+
     #tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq32e1int(var,ndim,Tcoeff)

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
      double precision sq3num1e0
      external sq3num1e0
      double precision sq3num1e1
      external sq3num1e1
      double precision sq3num1e2
      external sq3num1e2
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

      t1 = dlog(z)
      t2 = t1**2
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t6 = t5**2
      t9 = dlog(1-l4)
      t12 = f3(0.D0,z)
      t14 = dlog(1-t12)
      t19 = dlog(2.D0)
      t22 = dlog(l4)
      t26 = dlog(1-z)
      t31 = t26**2
      t33 = t14**2
      t35 = t9**2
      t41 = -t2/1024-t6/4096+t9*t5/2048-t9*t14/2048+t14*t1/1024-t9*t19/5
     #12-t9*t22/2048-t9*t26/512+t26*t5/512-t31/256-t33/4096-t35/4096-t22
     #*t19/512+t14*t5/2048
      t48 = t22**2
      t50 = 0.3141592653589793D1**2
      t52 = t19**2
      t72 = t5*t19/512+t26*t1/256-t5*t1/1024-t48/4096+t50/6144-t52/256-t
     #14*t19/512-t14*t26/512+t9*t1/1024+t1*t19/256-t22*t14/2048-t26*t19/
     #128+t22*t5/2048+t22*t1/1024-t22*t26/512
      t74 = Dsquar312(z,0.D0,0.D0,0.D0,l4)
      t75 = 1/t74
      t77 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t78 = 1/t77
      t79 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t80 = t79**2
      t81 = 1/t80
      t83 = f1(0.D0,z,1.D0)
      t84 = t78*t81*t83
      t86 = 4*t26
      t87 = 4*t19
      t88 = 2*t1
      t89 = dlog(l3)
      t92 = 1/l3
      t97 = dlog(l2)
      t99 = 1/l2
      t100 = (-t14-t86+t5-t87+t88-t97-t22-t9)*t99
      t101 = Dsquar312(z,0.D0,l2,0.D0,l4)
      t102 = 1/t101
      t108 = t81*t83
      t110 = t99*t92
      t117 = dlog(l1)
      t120 = 1/l1
      t129 = t120*t99
      tncoff1 = (t41+t72)*t75*t84-(-t14-t86+t5-t87+t88-t9-t22-2*t89)*t92
     #*t75*t84/1024+(t100*t102/2048-t100*t75/2048)*t78*t108+(t110*t102/1
     #024-t110*t75/1024)*t78*t108-(-t22-2*t117-t14-t86+t5-t87+t88-t9)*t1
     #20*t75*t84/1024-t120*t92*t75*t84/512+(t129*t102/1024-t129*t75/1024
     #)*t78*t108
      t2 = dlog(1-l3)
      t4 = f3(l3,z)
      t6 = dlog(1-t4)
      t8 = dlog(1-z)
      t10 = f2(l3,z)
      t11 = dlog(t10)
      t12 = dlog(2.D0)
      t14 = dlog(z)
      t17 = dlog(1-l4)
      t18 = dlog(l4)
      t19 = dlog(l3)
      t22 = 1/l3
      t24 = Dsquar312(z,0.D0,0.D0,l3,l4)
      t25 = 1/t24
      t27 = s23r2(z,1.D0,1.D0,l3,l4)
      t28 = 1/t27
      t29 = s1vr2(z,1.D0,1.D0,l3,l4)
      t30 = t29**2
      t31 = 1/t30
      t33 = f1(l3,z,1.D0)
      t34 = t28*t31*t33
      t38 = 1/l2*t22
      t39 = Dsquar312(z,0.D0,l2,l3,l4)
      tncoff2 = -(2*t2+t6+4*t8-t11+4*t12-2*t14+t17+t18+2*t19)*t22*t25*t3
     #4/1024+(-t38/t39/1024+t38*t25/1024)*t28*t31*t33+1/l1*t22*t25*t34/5
     #12
      t1 = dlog(l4)
      t2 = dlog(l1)
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t6 = f3(0.D0,z)
      t9 = dlog(1-t6+l1*t6)
      t11 = dlog(1-z)
      t13 = dlog(2.D0)
      t15 = dlog(z)
      t17 = 1-l1
      t18 = dlog(t17)
      t20 = dlog(1-l4)
      t22 = 1/l1
      t24 = Dsquar312(z,l1,0.D0,0.D0,l4)
      t25 = 1/t24
      t27 = s23r2(z,t17,1.D0,0.D0,l4)
      t29 = s1vr2(z,t17,1.D0,0.D0,l4)
      t30 = t29**2
      t33 = f1(0.D0,z,t17)
      t34 = 1/t27/t30*t33
      tncoff3 = (-t1-2*t2+t5-t9-4*t11-4*t13+2*t15-t18-t20)*t22*t25*t34/1
     #024+t22/l3*t25*t34/512+t22/l2*t25*t34/1024
      t4 = Dsquar312(z,l1,0.D0,l3,l4)
      t7 = 1-l1
      t8 = s23r2(z,t7,1.D0,l3,l4)
      t10 = s1vr2(z,t7,1.D0,l3,l4)
      t11 = t10**2
      t14 = f1(l3,z,t7)
      tncoff4 = -1/l1/l3/t4/t8/t11*t14/512
      t4 = Dsquar312(z,l1,l2,0.D0,l4)
      t7 = 1-l1
      t9 = 1-l1*l2
      t10 = s23r2(z,t7,t9,0.D0,l4)
      t12 = s1vr2(z,t7,t9,0.D0,l4)
      t13 = t12**2
      t16 = f1(0.D0,z,t7)
      tncoff5 = -1/l1/l2/t4/t10/t13*t16/1024
      t1 = dlog(l4)
      t4 = dlog(1-l4)
      t6 = f3(0.D0,z)
      t8 = dlog(1-t6)
      t11 = dlog(1-z)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t16 = dlog(2.D0)
      t18 = dlog(z)
      t21 = Dsquar312(z,0.D0,0.D0,0.D0,l4)
      t22 = 1/t21
      t24 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t25 = 1/t24
      t26 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t27 = t26**2
      t28 = 1/t27
      t30 = f1(0.D0,z,1.D0)
      t31 = t25*t28*t30
      t37 = 1/l2
      t38 = Dsquar312(z,0.D0,l2,0.D0,l4)
      tncoff6 = (t1/2048+t4/2048+t8/2048+t11/512-t14/2048+t16/512-t18/10
     #24)*t22*t31-1/l3*t22*t31/1024+(t37/t38/2048-t37*t22/2048)*t25*t28*
     #t30-1/l1*t22*t31/1024
      t2 = Dsquar312(z,0.D0,0.D0,l3,l4)
      t5 = s23r2(z,1.D0,1.D0,l3,l4)
      t7 = s1vr2(z,1.D0,1.D0,l3,l4)
      t8 = t7**2
      t11 = f1(l3,z,1.D0)
      tncoff7 = 1/l3/t2/t5/t8*t11/1024
      t2 = Dsquar312(z,l1,0.D0,0.D0,l4)
      t5 = 1-l1
      t6 = s23r2(z,t5,1.D0,0.D0,l4)
      t8 = s1vr2(z,t5,1.D0,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,z,t5)
      tncoff8 = 1/l1/t2/t6/t9*t12/1024
      t1 = Dsquar312(z,0.D0,0.D0,0.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t6 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,1.D0)
      tncoff9 = -1/t1/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num1e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num1e0(var,ndim,z,ux,1-l1,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num1e0(var,ndim,z,ux,1-l1,1.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num1e0(var,ndim,z,ux,1-l1,1-l1*l2,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num1e1(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num1e1(var,ndim,z,ux,1-l1,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num1e2(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      sq32e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq32e2int(var,ndim,Tcoeff)

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
      double precision sq3num1e0
      external sq3num1e0
      double precision sq3num1e1
      external sq3num1e1
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
      t4 = dlog(1-l4)
      t6 = f3(0.D0,z)
      t8 = dlog(1-t6)
      t11 = dlog(1-z)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t16 = dlog(2.D0)
      t18 = dlog(z)
      t21 = Dsquar312(z,0.D0,0.D0,0.D0,l4)
      t22 = 1/t21
      t24 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t25 = 1/t24
      t26 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t27 = t26**2
      t28 = 1/t27
      t30 = f1(0.D0,z,1.D0)
      t31 = t25*t28*t30
      t37 = 1/l2
      t38 = Dsquar312(z,0.D0,l2,0.D0,l4)
      tncoff1 = (t1/2048+t4/2048+t8/2048+t11/512-t14/2048+t16/512-t18/10
     #24)*t22*t31-1/l3*t22*t31/1024+(t37/t38/2048-t37*t22/2048)*t25*t28*
     #t30-1/l1*t22*t31/1024
      t2 = Dsquar312(z,0.D0,0.D0,l3,l4)
      t5 = s23r2(z,1.D0,1.D0,l3,l4)
      t7 = s1vr2(z,1.D0,1.D0,l3,l4)
      t8 = t7**2
      t11 = f1(l3,z,1.D0)
      tncoff2 = 1/l3/t2/t5/t8*t11/1024
      t2 = Dsquar312(z,l1,0.D0,0.D0,l4)
      t5 = 1-l1
      t6 = s23r2(z,t5,1.D0,0.D0,l4)
      t8 = s1vr2(z,t5,1.D0,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,z,t5)
      tncoff3 = 1/l1/t2/t6/t9*t12/1024
      t1 = Dsquar312(z,0.D0,0.D0,0.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t6 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,1.D0)
      tncoff4 = -1/t1/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num1e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num1e0(var,ndim,z,ux,1-l1,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      sq32e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq32e3int(var,ndim,Tcoeff)

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
      double precision sq3num1e0
      external sq3num1e0
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

      t1 = Dsquar312(z,0.D0,0.D0,0.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t6 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,1.D0)
      tncoff1 = -1/t1/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      sq32e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq32e4int(var,ndim,Tcoeff)

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

      sq32e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

